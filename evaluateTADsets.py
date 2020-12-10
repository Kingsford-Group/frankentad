import argparse
import math
import sys
import glob
import numpy as np
import pandas as pd
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import seaborn as sns
import scipy.stats
import collections

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

def readCTCFfile(res,chrnum,filename = '/mnt/disk70/user/nsauerwa/frankenTAD/go/inputdata/CTCFBSDB_all_exp_sites_Sept12_2012_hg19.txt'):
    
    # for ctcf and histone markers
    if chrnum == 'all':
        ctcfpeaks = {}
    else:
        ctcfpeaks = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if line[4] == 'chr'+chrnum:
                ctcfpeaks.append(int((int(line[5])+int(line[6]))/2/res))
            elif chrnum == 'all':
                try:
                    chrkey = int(line[4][3:])
                except:
                   chrkey = line[4][3:] 
                if chrkey not in ctcfpeaks:
                    ctcfpeaks[chrkey] = []
                ctcfpeaks[chrkey].append(int((int(line[5])+int(line[6]))/2/res))
    if chrnum == 'all':
        for key,ctcfdata in ctcfpeaks.iteritems():
            ctcfdata.sort()
            ctcfpeaks[key] = collections.Counter(ctcfdata)
    else:
        ctcfpeaks.sort()
        ctcfpeaks = collections.Counter(ctcfpeaks)
    return ctcfpeaks

def readBEDfile(filename,res,chrnum):

    if chrnum == 'all':
        peaks = {}
    else:
        peaks = []
    with open(filename,'r') as f:
        freader = csv.reader(f,delimiter='\t')
        row1 = next(freader)
        if len(row1) == 1: freader = csv.reader(f,delimiter=' ')
        for line in freader:
            if line[0] == 'chr'+chrnum:
                peaks.append(int((int(line[1])+int(line[2]))/2/res))
            elif chrnum == 'all':
                try:
                    chrkey = int(line[0][3:])
                except:
                    chrkey = line[0][3:]
                if chrkey not in peaks:
                    peaks[chrkey] = []
                #print(line)
                peaks[chrkey].append(int((int(line[1])+int(line[2]))/2/res))
    if chrnum == 'all':
        for key,data in peaks.iteritems():
            data.sort()
            peaks[key] = collections.Counter(data)
    else:
        peaks.sort()
        peaks = collections.Counter(peaks)
    return peaks

def calcJI(tadset1,tadset2):
    # have to return number of overlaps and total number by chromosome, then later add them up for genome-wide value
    tadbdys1 = tadset1.flatten()
    tadbdys2 = tadset2.flatten()
    intersectsize = len(set(tadbdys1) & set(tadbdys2))
    unionsize = len(set(tadbdys2) | set(tadbdys2))
    return intersectsize,unionsize

def plotTADs(tadset,hicdata,filename):
    
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
    for interval in tadset: # plot outline of each domain
        plt.plot((interval[0]-1,interval[1]),(interval[0]+1,interval[0]+1),'k')
        plt.plot((interval[1],interval[1]),(interval[0]+1,interval[1]),'k')
        
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    print('Saved plot of Hi-C with TADs to',filename)
    
def preprocessEpigenetics(tadset, markerlocs, res, chrlength):

    # need dataframe with 11*2 rows for each TAD (-500kb, 400kb, -300kb, ... 0kb, 100kb, 200kb, ... 500kb) * 2 (one centered at bdy and one centered at midpt)
    # create list of dicts and convert to dataframe
    binrange = 500000/res
    data = []
    for tad in tadset:
        tadmidpt = (tad[0]+tad[1])/2
        for binloc in range(-binrange,binrange+1): #range(max([0,tadmidpt-binrange]), min([tadmidpt+binrange+1,chrlength])):
            try:
                centerval = markerlocs[tadmidpt+binloc]
            except:
                centerval = 0
            #data.extend([{'dist':binloc-tadmidpt, 'markerval':markerval, 'type':'TAD center'}])
            bdylocs = [tad[0],tad[1]]
            for bdy in bdylocs:
            #for binloc in range(max([0,bdy-binrange]), min([bdy+binrange+1,chrlength])):
                try:
                    bdyval = markerlocs[bdy+binloc]
                except:
                    bdyval = 0
                if centerval > 0:
                    fc = bdyval/float(centerval)
                else:
                    fc = float('NaN')
                data.extend([{'dist': binloc, 'foldchange':fc, 'boundaryval':bdyval }])
    datadf = pd.DataFrame(data)
    return datadf

def plotEpigenetics(datadf,label,filename):
    
    # plot lines of average # of ctcf sites from -500kb to 500kb of (1) TAD boundaries and (2) TAD centers, with shaded error regions
    #print(datadf.columns)
    ax = sns.lineplot(x="dist", y="boundaryval", data=datadf)
    plt.ylabel(label)
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    print('Saved plot of epigenetics at TAD boundaries vs centers to',filename)
    
    # compute difference between average CTCF value at boundary and average CTCF value at midpt
    bdyvals = datadf.loc[(datadf['dist'] == 0)]['boundaryval']
    avgbdyval = np.mean(bdyvals)#['markerval'].tolist())
    midptvals = datadf.loc[(datadf['dist'] == 0) & (datadf['type'] == 'TAD center')]
    avgmidptval = np.mean(midptvals['markerval'].tolist())
    return avgbdyval - avgmidptval

def computeInterIntraTADcts(tadset,hicdata):
    
    hicdata = hicdata.tocsr()

    intertadcts = []
    intratadcts = []
    for i,tad in enumerate(tadset):
        for cidx in tad:
            for ridx in range(cidx+1,tad[1]):
                try:
                    intratadcts.extend([math.log(hicdata[cidx,ridx],2)])
                except:
                    continue
            if i > 0:
                for ridx in range(tadset[i-1][0],tadset[i-1][1]):
                    try: intertadcts.extend([math.log(hicdata[ridx,cidx],2)])
                    except: continue
                    #print(ridx,cidx)
    return intertadcts, intratadcts

def plotInterIntraTADcts(intertadcts, intratadcts, filename):

    muval,pval = scipy.stats.mannwhitneyu(intratadcts, intertadcts, alternative='greater')
    print('Mann Whitney test results for intra vs inter TAD counts: mu =',muval,'p =',pval)

    fig = plt.figure()
    ax = plt.gca()

    alldata = intertadcts + intratadcts
    alllabels = ['inter-TAD values']*len(intertadcts) + ['intra-TAD values']*len(intratadcts)
    datadict = {'labels':alllabels, 'values': alldata}
    alldatadf = pd.DataFrame(datadict)
    sns.boxplot(x=alldatadf['labels'],y=alldatadf['values'])

    plt.ylabel('Log2 of interaction frequency')
    plt.xlabel('')
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('Saved plot of inter vs intra TAD contact values to',filename)

    return np.mean(intratadcts) - np.mean(intertadcts)


def writeTADQualReport(inputtadfile, inputhicfile, res, tadset, genomeji, ctcfdiff, histdiff, countdiff, filename):
    
    if isinstance(tadset,dict):
        numtads = np.mean([len(t) for _,t in tadset.iteritems()])
        tadset = np.array([x for _,t in tadset.iteritems() for x in t])
        avgtadlen = res*np.mean(tadset[:,1] - tadset[:,0])
    else:
        numtads = len(tadset)
        tadset = np.array(tadset)
        avgtadlen = res*np.mean(tadset[:,1] - tadset[:,0])
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        fwriter.writerow(['inputs:', inputtadfile, inputhicfile])
        fwriter.writerow(['number of TADs (average if full genome)', numtads])
        fwriter.writerow(['average TAD length (bp)', avgtadlen])
        if genomeji > -1:
            fwriter.writerow(['JI between two input TAD sets',genomeji])
        fwriter.writerow(['difference in avg CTCF at TAD bdys vs centers',ctcfdiff])
        if len(histdiff) > 0:
            for hlabel,hdiff in histdiff.iteritems():
                fwriter.writerow(['difference in avg '+hlabel+' at TAD bdys vs centers',hdiff])
        fwriter.writerow(['difference in avg intraTAD vs interTAD contacts', countdiff])
    print('Wrote evaluation report to',filename)

def computeGenomeJI(comptadfiles,tadsets,res):

    jivals = []
    comptadlist = glob.glob(comptadfiles)
    for ctadf in comptadlist:
        try:
            chridx = ctadf.index('chr')
        except:
            chridx = ctadf.index('Chr')
        try:
            chrval = int(ctadf[chridx+3:chridx+5])
        except:
            chrval = int(ctadf[chridx+3])
        if len(res) == 1: 
            chrtads = readTADFile(ctadf,res[0])
        else: 
            chrtads = readTADFile(ctadf,res[1])
        jiint,jiunion = calcJI(np.array(chrtads), np.array(tadsets[chrval]))
        jivals.append([jiint,jiunion])
    jivals = np.array(jivals)
    return np.sum(jivals[:,0])/float(np.sum(jivals[:,1]))
        

def combineAndWriteData(histdf, ctcfdf, intertadcts, intratadcts, filename):
    
    # need to add a markertype column to dataframe, put histones and ctcf all together
    masterdf = ctcfdf.assign(markertype = ["CTCF"]*len(ctcfdf))
    for key,df in histdf.iteritems():
        label = [key]*len(df)
        df['markertype'] = label 
        masterdf = pd.concat([masterdf, df], ignore_index=True)

    masterdf.to_csv(filename+'_epigendata.txt',sep='\t')
    print('wrote epigenetic data to file',filename+'_epigendata.txt')

    #write intra/intertad counts to separate file?
    with open(filename+'_interintravals.txt','w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        fwriter.writerow(['interTADvals']+intertadcts)
        fwriter.writerow(['intraTADvals']+intratadcts)
    print('wrote inter/intra tad counts to',filename+'_interintravals.txt')


def main(tadfile, res, hicfile, ctcffile, histonefiles, histlabels, chrnum, outfile):

    # read in TAD files and hic files
    # add possibility for full genome (if chrnum == 'all')
    if chrnum == 'all':
        tadsets = {}
        hicdata = {}
        chrlengths = {}
        tadfilelist = glob.glob(tadfile[0])
        #print(len(tadfilelist))
        #sys.exit()
        for tadf in tadfilelist:
            if 'tadscores' in tadf: continue
            #print(tadf)
            try:
                chridx = tadf.index('chr')
            except:
                chridx = tadf.index('Chr')
            try:
                chrval = int(tadf[chridx+3:chridx+5])
            except:
                chrval = int(tadf[chridx+3])
            tadsets[chrval] = readTADFile(tadf,res[0])
        print('finished reading TAD sets')
        #sys.exit()
        if len(tadfile) > 1:
            print('computing JI between input sets')
            genomeji = computeGenomeJI(tadfile[1],tadsets,res)
        else:
            genomeji = -1

        hicfilelist = glob.glob(hicfile)
        for hicf in hicfilelist:
            if "full" in hicf: continue
            chridx = hicf.index('chr')
            try:
                chrval = int(hicf[chridx+3:chridx+5])
            except:
                try:
                    chrval = int(hicf[chridx+3])
                except:
                    chrval = hicf[chridx+3]
            hicdata[chrval],chrlengths[chrval] = readHiCFile(hicf,res[0])
        print('finished reading HiC data')
        if len(ctcffile) == 0:
            ctcfcounts = readCTCFfile(res[0],chrnum)
        elif ctcffile[-4:] == "Peak" or ctcffile[-3:] == "bed":
            ctcfcounts = readBEDfile(ctcffile,res[0],chrnum)
        else:
            ctcfcounts = readCTCFfile(res[0],chrnum,ctcffile)
        print('finished reading CTCF file')

        histdata = {}
        if len(histonefiles) != len(histlabels):
            print('Please provide labels for each histone marker bed file provided')
            sys.exit()
        elif len(histonefiles) > 0:
            histdf = {}
            histdiff = {}
            for idx,hf in enumerate(histonefiles):
                #print(hf)
                histdata[histlabels[idx]] = readBEDfile(hf,res[0],chrnum)
            print('read',histlabels[idx],'data')
        else:
            histdf = {}
            histdiff = []

        intertadcts = []
        intratadcts = []
        for chrval in range(1,23):
            #plotTADs(tadsets[chrval],hicdata[chrval],outfile+'_chr'+str(chrval)+'_hictadplot.png')
            intertadctschr, intratadctschr = computeInterIntraTADcts(tadsets[chrval],hicdata[chrval])
            #print(intertadctschr,intratadctschr)
            intertadcts.extend(intertadctschr)
            intratadcts.extend(intratadctschr)

            chrdf = preprocessEpigenetics(tadsets[chrval], ctcfcounts[chrval], res[0], chrlengths[chrval]) # concatenate dfs from each chr
            if chrval > 1:
                ctcfdf = pd.concat([ctcfdf, chrdf], ignore_index=True)
            else:
                ctcfdf = chrdf
                
            if len(histonefiles) > 0:
                for hlabel,hdata in histdata.iteritems():
                    histdfchr = preprocessEpigenetics(tadsets[chrval],hdata[chrval],res[0],chrlengths[chrval])
                    if chrval > 1:
                        histdf[hlabel] = pd.concat([histdf[hlabel], histdfchr], ignore_index=True)
                    else:
                        histdf[hlabel] = histdfchr
            print(chrval,"done")
        
        #for hlabel,hdf in histdf.iteritems():
        #    histdiff[histlabels[idx]] = plotEpigenetics(hdf,hlabel,outfile+'_fullgenome_'+hlabel+'atbdys.png')

        #tadctdiff = plotInterIntraTADcts(intertadcts, intratadcts, outfile+'_fullgenome_interVSintraTADcounts.png')
        #ctcfdiff = plotEpigenetics(ctcfdf,'CTCF', outfile+'_fullgenome_ctcfatbdys.png')

        combineAndWriteData(histdf, ctcfdf, intertadcts, intratadcts, outfile)
        sys.exit()

        writeTADQualReport(tadfile,hicfile,res[0],tadsets,genomeji,ctcfdiff,histdiff,tadctdiff,outfile+'_fullgenome_evaluationreport.txt')
    
    else:
        tadset = readTADFile(tadfile[0],res[0])
        #hicdata,chrlength = readHiCFile(hicfile,res[0])

        if len(tadfile) > 1:
            comptadset = readTADFile(tadfile[1],res[-1])
            jiint,jiunion = calcJI(np.array(tadset),np.array(comptadset))
            jival = jiint/float(jiunion)
        else:
            jival = -1

        print(jival)
        sys.exit()
        plotTADs(tadset,hicdata,outfile+'_hictadplot.png')
        sys.exit()

        intertadcts, intratadcts = computeInterIntraTADcts(tadset,hicdata)
        tadctdiff = plotInterIntraTADcts(intertadcts, intratadcts, outfile+'_interVSintraTADcounts.png')

        if len(histonefiles) != len(histlabels):
            print('Please provide labels for each histone marker bed file provided')
            sys.exit()
        elif len(histonefiles) > 0:
            histdiff = {}
            for idx,hf in enumerate(histonefiles):
                histdata = readBEDfile(hf,res[0],chrnum)
                histdf = preprocessEpigenetics(tadset,histdata,res[0],chrlength)
                histdiff[histlabels[idx]] = plotEpigenetics(histdf,histlabels[idx],outfile+'_'+histlabels[idx]+'atbdys.png')
        else:
            histdiff = []
        
        ctcfcounts = readCTCFfile(res[0],chrnum)
        ctcfdf = preprocessEpigenetics(tadset, ctcfcounts, res[0], chrlength)
        ctcfdiff = plotEpigenetics(ctcfdf, 'CTCF', outfile+'_ctcfatbdys.png')

        writeTADQualReport(tadfile,hicfile,res[0],tadset,jival,ctcfdiff,histdiff,tadctdiff,outfile+'_evaluationreport.txt') # inputtadfile, inputhicfile, res, tadset,ctcfdiff,countdiff,filename)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-t',type=str,nargs='+',help='TAD set file (can include multiple to compare replicates or resolutions if desired)')
    parser.add_argument('-r', type=int,nargs='+', help='Hi-C and/or TAD resolution. If TAD sets of different resolutions are input, give all resolution values here in same order')
    parser.add_argument('-hic', type=str, help='Hi-C data file')
    parser.add_argument('-ctcf',type=str, default='',help='Peak file with CTCF binding sites')
    parser.add_argument('-histone',type=str,default='',nargs='+',help='ChIPseq files for histone markers (optional)')
    parser.add_argument('-hl',type=str,default='',nargs='+',help='Labels of histone marker data used')
    parser.add_argument('-c', type=str, help='Chromosome number')
    parser.add_argument('-o', type=str, help='File string for outputs')
    

    args=parser.parse_args()

    main(args.t, args.r, args.hic,args.ctcf,args.histone, args.hl, args.c, args.o)
