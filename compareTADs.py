## script to take in computed TAD sets on multiple cell types, compare the reproducibility scores for all pairs of replicates, and compare across resolutions, (to add: comparing across subsampled data)

import argparse
import sys
import glob
import numpy as np
import pandas as pd
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import scipy.stats
import collections
import seaborn as sns
import math

def readTADFile(filename, res):

    tadset = []
    with open(filename, 'r') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if len(line) == 3: # armatus output
                offset = 1
            else:
                offset = 0
            #print(line,len(line),offset)
            try:
                # how to know if TADs should be divided by res?
                tadset.append([int(line[offset])/res, (int(line[offset+1])+1)/res - 1])
            except:
                continue

    if len(tadset) > 0 and tadset[0][0] > tadset[-1][0]:
        tadset.reverse()
    return tadset

def readHiCFile(filename, res):

    # should be able to read matrices in both sparse and full formats eventually, just starting w/ reading sparse format
    hicidx = []
    hicdata = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            hicidx.append([int(line[0])/res, int(line[1])/res])
            hicdata.extend([float(line[2])])
    hicidx = np.array(hicidx)
    chrlength = np.max(hicidx)
    hicdata = np.array(hicdata)
    hicsparse = coo_matrix( (hicdata, (hicidx[:,0], hicidx[:,1]) ))

    # compute & return chrlength
    return hicsparse, chrlength

def readScores(filename):

    #read scores into dictionary with tad as key
    scoredict = {}
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        freader.next()
        maxscores = [-float('inf')]*4
        for line in freader:
            dictkey = (int(line[0]),int(line[1]))
            scoredict[dictkey] = [float(x) for x in line[2:]]
            for i,currmax in enumerate(maxscores):
                if float(line[2+i]) > currmax:
                    maxscores[i] = float(line[2+i])
        scoredict['maxscores'] = maxscores
        #print(maxscores)
        #sys.exit()
    return scoredict


def plotTADs(tadsets,hicdata,filename):

    # generate full matrix from sparse format for heatmap
    datamat = hicdata.todense()
    if datamat.shape[0] > datamat.shape[1]:
        # add column(s) of zeros to make square matrix
        ncols = datamat.shape[0] - datamat.shape[1]
        sqmat = np.zeros((datamat.shape[0],datamat.shape[0]))
        sqmat[:,:-1*ncols] = datamat
        datamat = sqmat
    elif datamat.shape[1] > datamat.shape[0]:
        # add row(s) of zeros to make square matrix
        nrows = datamat.shape[1] - datamat.shape[0]
        sqmat = np.zeros((datamat.shape[1],datamat.shape[1]))
        sqmat[:-1*nrows,:] = datamat
        datamat = sqmat
    datamat = datamat + np.transpose(datamat) - np.diagonal(datamat)*np.identity(len(datamat))
    logdata = np.ma.log(datamat) # log transform data to show up better in figure
    logdata = logdata.filled(0)
    labelspacing = int(math.floor(round(len(logdata),-int(math.floor(math.log10(len(logdata)))))/10))

    c = np.linspace(0,1,101)
    colors = plt.get_cmap("YlOrRd",101)(c)
    colors[0,:] = np.array([1,1,1,1])
    customcmap = matplotlib.colors.ListedColormap(colors)

    plt.figure(figsize=(10,10))
    ax = sns.heatmap(logdata,cbar=False,xticklabels=labelspacing,yticklabels=labelspacing,cmap=customcmap,square=True)

    #plot TADs on top of heatmap
    colors = ['c','m','g','k']
    for i,tset in enumerate(tadsets):
        for interval in tset: # plot outline of each domain
            plt.plot((interval[0]-1,interval[1]),(interval[0]+1,interval[0]+1),colors[i])
            plt.plot((interval[1],interval[1]),(interval[0]+1,interval[1]),colors[i])

    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    print('Saved plot of Hi-C with TADs to',filename)

def calcJI(tadset1,tadset2):
    # have to return number of overlaps and total number by chromosome, then later add them up for genome-wide value
    tadbdys1 = tadset1.flatten()
    tadbdys2 = tadset2.flatten()
    intersectsize = len(set(tadbdys1) & set(tadbdys2))
    unionsize = len(set(tadbdys2) | set(tadbdys2))
    return intersectsize,unionsize

def calcApproxJI(tadset1,tadset2,window):
    tadbdys1 = tadset1.flatten()
    tadbdys2 = tadset2.flatten()
    intersectsize = 0
    unionlist = list(tadbdys1)+list(tadbdys2)
    for bdy1 in tadbdys1:
        bdy1matched = False
        for bdy2 in tadbdys2:
            if not bdy1matched:
                if abs(bdy1-bdy2) <= window:
                    intersectsize += 1
                    unionlist.remove(bdy1)
                    bdy1matched = True
    unionsize = len(unionlist)
    return intersectsize,unionsize

def computeReproducibility(celltypelist,filepath):

    reproducibility = []
    for celltype in celltypelist:
        if celltype == 'NHEK': continue
        celltyperepro = []
        for chrnum in range(1,23):
            chrstr = 'chr'+str(chrnum)
            repfiles = glob.glob(filepath+'*'+celltype+'*rep*'+chrstr+'.txt')
            if len(repfiles) == 0:
                repfiles = glob.glob(filepath+'*'+celltype+'*rep*'+chrstr+'_clean.txt')
            if len(repfiles) == 0:
                repfiles = glob.glob(filepath+'*'+celltype+'*rep*.txt')
                repfiles = [x for x in repfiles if 'hr'+str(chrnum)+'_' in x]
            if len(repfiles) == 0 and celltype == 'HFFc6':
                repfiles = glob.glob(filepath+'4DNChr'+str(chrnum)+'_HFF-c6_rep*')
            if len(repfiles) == 0 and celltype == 'SkeletalMuscle':
                repfiles = glob.glob(filepath+'*hr'+str(chrnum)+'_*SkelMuscle*')
            if len(repfiles) == 0:
                print(filepath)
                print(filepath+'*'+celltype+'*rep*'+chrstr+'.txt')
                sys.exit()
            #if 'armatus' in repfiles[0]:
            #    print(repfiles)#,filepath+'*'+celltype+'*rep*.txt')
                #sys.exit()
            chrrepdata = []
            for rf in repfiles:
                # figure out resolution from file name
                if '40k' not in rf: continue
                res = 40000
                chrrepdata.append(readTADFile(rf,res))
            for idx,rep1 in enumerate(chrrepdata):
                for _,rep2 in enumerate(chrrepdata[idx:]):
                    jiint,jiunion = calcJI(np.array(rep1), np.array(rep2))
                    celltyperepro.append([jiint,jiunion])
        celltyperepro = np.array(celltyperepro)
        jival = np.sum(celltyperepro[:,0])/float(np.sum(celltyperepro[:,1]))
        reproducibility.extend([jival])
    return reproducibility

def computeSubsampling(fullfilepath,subsampfilepath):

    #print(fullfilepath)
    subsampji = {}
    for chrnum in range(1,23):
        chrstr = 'chr'+str(chrnum)
        fullfile = glob.glob(fullfilepath+'*'+chrstr+'.txt')
        #print(fullfile)
        #print(fullfilepath+'*'+chrstr+'.txt')
        if len(fullfile) == 0:
            fullfile = glob.glob(fullfilepath+'*'+chrstr+'_*.txt')
        if len(fullfile) == 0:
            fullfile = glob.glob(fullfilepath)
            fullfile = [x for x in fullfile if 'hr'+str(chrnum)+'_' in x]
            #print(fullfile, fullfilepath)
        if len(fullfile) != 1:
            print('couldnt find unique file for subsample comparison')
            print(fullfile)
            sys.exit()
        #print('full file',fullfile[0])
        #print(fullfilepath+'*'+chrstr+'.txt')
        fulldata = readTADFile(fullfile[0],40000)
        subfiles = glob.glob(subsampfilepath+'*'+chrstr+'.txt')
        #if 'armatus' in fullfile[0]:
        #    print('firstpass')
        #    print(subfiles)
        if len(subfiles) == 0:
            subfiles = glob.glob(subsampfilepath+chrstr+'_*')
            #if 'armatus' in fullfile[0]:
            #    print('secondpass')
            #    print(subfiles)
        if len(subfiles) == 0:
            subfiles = glob.glob(subsampfilepath+'hr'+str(chrnum)+'_*')
            #if 'armatus' in fullfile[0]:
            #    print('lastpass')
            #    print(subfiles)
            #sys.exit()
        if len(subfiles) != 3:
            print('couldnt find correct number of subsampled files')
            print(subfiles)
            sys.exit()
        #print('full', len(fulldata))
        subdata = {}
        #print(subfiles)
        #print(fullfile)
        for sf in subfiles:
            sfsplit = sf.split('_')
            sublevel = [x[7:] for x in sfsplit if 'subsamp' in x][1]
            if '40k' not in sf: continue
            res = 40000
            #print(sublevel,sf)
            subdata[sublevel] = readTADFile(sf,res)
            #print(sublevel,len(subdata[sublevel]))
            #print(sublevel,subdata[sublevel])
        #if 'hicseg' in fullfile[0]:
        #    print(subdata)
        for sublevel,subsamp in subdata.items():
            jiint,jiunion = calcJI(np.array(fulldata), np.array(subsamp))
            if sublevel in subsampji:
                subsampji[sublevel].append([jiint,jiunion])
            else:
                subsampji[sublevel] = [[jiint,jiunion]]
            #if 'hicseg' in fullfile[0]:
            #    print(chrstr, sublevel, jiint,jiunion)
    #if 'hicseg' in fullfile[0]:
    #    sys.exit()
    #print(subsampji)
    #if 'armatus' in fullfile[0]:
    #    sys.exit()
    for key,subsampdata in subsampji.items():
        ssji = np.array(subsampdata)
        jival = np.sum(ssji[:,0])/float(np.sum(ssji[:,1]))
        subsampji[key] = jival
        #print(key,jival)
    return subsampji

def compareResolutions(celltypelist,filepath,window):

    rescomp = []
    for celltype in celltypelist:
        celltyperescomp = []
        for chrnum in range(1,23):
            chrstr = 'chr'+str(chrnum)
            chrresdata = []
            resfiles = glob.glob(filepath+'*'+celltype+'*defaultparams*'+chrstr+'.txt')
            for rf in resfiles:
                if 'rep' in rf: continue
                chrresdata.append(readTADFile(rf,1))
            for idx,rep1 in enumerate(chrresdata):
                for _,rep2 in enumerate(chrresdata[idx:]):
                    jiint,jiunion = calcApproxJI(np.array(rep1), np.array(rep2), window)
                    celltyperescomp.append([jiint,jiunion])
        celltyperescomp = np.array(celltyperescomp)
        jival = np.sum(celltyperescomp[:,0])/float(np.sum(celltyperescomp[:,1]))
        rescomp.extend([jival])
        #print(celltype,jival)
    return rescomp
    #print(calcApproxJI(np.array([[1,3],[6,10],[14,30]]), np.array([[1,2],[8,11],[14,20]]), 1))
    

def compareToSingleFcns(celltypelist, filepath, hicfile, res,outfile):

    for celltype in celltypelist:
        tadsfoundbyall = 0
        tadsonlyfrankentad = 0
        totaltads = 0
        #datadf = pd.DataFrame(columns=['TAD','chr','GDscore','BIscore','BDscore','fTscore','unique'])
        datalist = []
        for chrnum in range(1,23):
            chrstr = 'chr'+str(chrnum)
            frankentadfile = glob.glob(filepath+'*'+celltype+'_40k_dp_defaultparams_'+chrstr+'.txt')
            if len(frankentadfile) != 1:
                frankentadfile = glob.glob(filepath+'*'+celltype+'_40kb_dp_defaultparams_'+chrstr+'.txt')
                #print('couldnt find unique file for frankentads')
                #print(celltype,chrnum, filepath+'*'+celltype+'_40k_dp_defaultparams_'+chrstr+'.txt')
                #sys.exit()
            frankentads = readTADFile(frankentadfile[0],40000)
            frankentadset = [tuple(tad) for tad in frankentads]
            gdonlyfile = glob.glob(filepath+'*'+celltype+'*gdonly*'+chrstr+'.txt')
            gdtads = readTADFile(gdonlyfile[0],40000)
            gdtadset = [tuple(tad) for tad in gdtads]
            bdonlyfile = glob.glob(filepath+'*'+celltype+'*bdonly*'+chrstr+'.txt')
            if len(bdonlyfile) < 1:
                print('couldnt find unique file for bias direction')
                print(celltype,chrnum, filepath+'*'+celltype+'_40k_dp_defaultparams_'+chrstr+'.txt')
                sys.exit()
            bdtads = readTADFile(bdonlyfile[0],40000)
            bdtadset = [tuple(tad) for tad in gdtads]
            bionlyfile = glob.glob(filepath+'*'+celltype+'*bionly*'+chrstr+'.txt')
            #if len(bionlyfile) > 0:
            bitads = readTADFile(bionlyfile[0],40000)
            bitadset = [tuple(tad) for tad in gdtads]
            #else:
            #bitads = bdtads

            #hicfilename = glob.glob(hicfile+chrstr+'.matrix')
            #hicdata,_ = readHiCFile(hicfilename[0],res)

            #plotTADs([gdtads,bdtads,bitads,frankentads],hicdata,outfile+celltype+'_'+chrstr+'_alltads.png')

            #tadsfoundbyall += len(set(frankentadset) & set(gdtadset) & set(bdtadset) & set(bitadset))
            #singlefcntads = gdtadset+bdtadset+bitadset
            #ftads = np.array(frankentadset)
            #ftadrows = ftads.view([('',ftads.dtype)]*ftads.shape[1])
            #nonftads = np.array(singlefcntads)
            #nonftadrows = nonftads.view([('',nonftads.dtype)]*nonftads.shape[1])
            #frankentadunique = np.setdiff1d(ftadrows,nonftadrows).view(ftads.dtype).reshape(-1,ftads.shape[1])
            
            frankentadunique = []
            tadsfoundbyall = []
            for ftad in frankentadset:
                if ftad in gdtadset and ftad in bdtadset and ftad in bitadset:
                    tadsfoundbyall.append(ftad)
                if ftad not in gdtadset and ftad not in bdtadset and ftad not in bitadset:
                    frankentadunique.append(ftad)
            #print(frankentadset)
            #print(gdtads)
            #print(bdtads)
            #print(bitads)
            #sys.exit()

            #plotTADs([frankentadunique],hicdata,outfile+celltype+'_'+chrstr+'_uniquetads.png')

            scorefile = glob.glob(filepath+'*'+celltype+'_40k_dp_defaultparams_'+chrstr+'_tadscores.txt')
            scores = readScores(scorefile[0])
            #print(scores)

            #print(datalist)
            chrdata =  updateDataFrame(frankentadset,frankentadunique,scores,chrnum)
            datalist.extend(chrdata)
            
            datadf = pd.DataFrame(chrdata)
            scoresunique = datadf.loc[(datadf['label'] =='frankenTAD') & (datadf['unique'] == True )]
            scoresnotunique = datadf.loc[(datadf['label'] =='frankenTAD') & (datadf['unique'] == False )]
            scoresunique = scoresunique['score'].tolist()
            scoresnotunique = scoresnotunique['score'].tolist()
            print(chrnum,scipy.stats.ks_2samp(scoresunique,scoresnotunique))

            scoresunique = datadf.loc[(datadf['label'] =='GD') & (datadf['unique'] == True )]
            scoresnotunique = datadf.loc[(datadf['label'] =='GD') & (datadf['unique'] == False )]
            scoresunique = scoresunique['score'].tolist()
            scoresnotunique = scoresnotunique['score'].tolist()
            print('GD',chrnum,scipy.stats.ks_2samp(scoresunique,scoresnotunique))

            scoresunique = datadf.loc[(datadf['label'] =='BD') & (datadf['unique'] == True )]
            scoresnotunique = datadf.loc[(datadf['label'] =='BD') & (datadf['unique'] == False )]
            scoresunique = scoresunique['score'].tolist()
            scoresnotunique = scoresnotunique['score'].tolist()
            print('BD',chrnum,scipy.stats.ks_2samp(scoresunique,scoresnotunique))

            scoresunique = datadf.loc[(datadf['label'] =='BI') & (datadf['unique'] == True )]
            scoresnotunique = datadf.loc[(datadf['label'] =='BI') & (datadf['unique'] == False )]
            scoresunique = scoresunique['score'].tolist()
            scoresnotunique = scoresnotunique['score'].tolist()
            print('BI',chrnum,scipy.stats.ks_2samp(scoresunique,scoresnotunique))


            #print("got through once")
            #print(celltype,chrnum,datalist)
            #print(celltype, chrnum, frankentadunique)
            tadsonlyfrankentad += len(frankentadunique)
            totaltads += len(frankentads)
        #if celltype == 'GM12878':
        #    sys.exit()
        datadf = pd.DataFrame(datalist)
        scoresunique = datadf.loc[(datadf['label'] =='frankenTAD') & (datadf['unique'] == True )]
        scoresnotunique = datadf.loc[(datadf['label'] =='frankenTAD') & (datadf['unique'] == False )]
        scoresunique = scoresunique['score'].tolist()
        scoresnotunique = scoresnotunique['score'].tolist()
        print(scipy.stats.ks_2samp(scoresunique,scoresnotunique))
        sys.exit()
        plotScores(datadf,outfile+'_'+celltype+'_tadscoresviolinplot.png')
        sys.exit()
        print(celltype)
        print('common to all fcns:',len(tadsfoundbyall))
        print('unique to frankentad:',tadsonlyfrankentad)
        print('nonunique TADs:',totaltads-tadsonlyfrankentad)
        print('total number of tads from frankentad:',totaltads)

def updateDataFrame(ftads,ftadunique,scores,chrnum):

    #['TAD','chr','GDscore','BIscore','BDscore','fTscore','unique']
    datatoadd = []
    for tad in ftads:
        if tad in ftadunique:
            isunique = True
        else:
            isunique = False
        if tad in scores:
            for i,score in enumerate(scores[tad]):
                if i == 0: 
                    label = 'GD'
                elif i == 1: 
                    label = 'BI'
                elif i == 2: 
                    label = 'BD'
                elif i == 3: 
                    label = 'frankenTAD'
                datatoadd.append({'TAD':tad, 'chr': chrnum, 'score': score/scores['maxscores'][i], 'label': label, 'unique': isunique})
        else:
            continue
            #print('couldnt find scores for tad',tad)
    #print(datalist,chrnum)
    #datalist.extend(datatoadd)
    return datatoadd
    

def plotScores(datadf,filename):

    # want to plot scores of TADs - all vs unique to frankenTAD?
    #ax = sns.violinplot(x='label',y='score',hue='unique',split=True,data=datadf)
    ax = sns.violinplot(x='unique',y='score',data=datadf.loc[(datadf['label'] == 'frankenTAD')])
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    print('Saved plot of TAD scores to',filename)


def main(celltypelist,filepath,fulldatapath,subsampfilepath,hicfiles,outpath):

    computeSubsampling(fulldatapath,subsampfilepath)
    sys.exit()

    reproducibility = computeReproducibility(celltypelist,filepath)
    print("reproducibility scores: ",reproducibility)
    sys.exit()
    window = 30000
    #resolutions = compareResolutions(celltypelist,filepath,window)
    #print("consistency across resolutions",resolutions)

    compareToSingleFcns(celltypelist,filepath,hicfiles,40000,outpath)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i',type=str, nargs='+', help='list of cell type names to search filepath for') # list of cell type names
    parser.add_argument('-p',type=str, help='Path to search for TAD files in') # file path to search for TAD files
    parser.add_argument('-o', type=str, help='File string for outputs') # file name/path for outputs
    parser.add_argument('-hic',type=str, help='HiC data')
    parser.add_argument('-ss',type=str,help='Path with subsampled data results')
    parser.add_argument('-fp',type=str,help='Path to results corresponding to full data')

    args=parser.parse_args()

    main(args.i, args.p, args.fp, args.ss, args.hic, args.o)
