package main

import (
	"flag"
	"hicutil"
	"ftutil"
//	"math/rand"
	"fmt"
	"math"
//	"os"
	"path/filepath"
	"strconv"
	"log"
	"strings"
//	"onlinestats"
//	"runtime"
//	"bufio"
	"time"
//	"github.com/d4l3k/go-bayesopt"
	"github.com/c-bata/goptuna"
	"github.com/c-bata/goptuna/tpe"
//	"context"
//	"golang.org/x/sync/errgroup"
//	"gonum.org/v1/gonum/stat"
	//"debugtools"
//	"sync"
)

// main file for frakenTAD

// define a pair (tuple) type for keys to Hi-C maps
type Pair struct {
     A,B interface{}
}


func main() {

	//runtime.GOMAXPROCS(1)

	//rand.Seed(5)

	hicfile := flag.String("hic","","Hi-C file(s)") // if in train mode, this should be path to Hi-C files of full genome
	res := flag.Int("res",1,"resolution of Hi-C data")

	trainflag := flag.Bool("train",false,"training mode to learn parameters")
	alpha_is_in := flag.Float64("alpha_is",-1.0,"TopDom/InsScore internal parameter")
	gamma_armatus_in := flag.Float64("gamma_armatus",-1.0,"Armatus internal parameter")
	//gamma_tt_in := flag.Float64("gamma_tt",-1.0,"TADTree internal parameter")
	lambda_armatus_in := flag.Float64("lambda_arm",-1.0,"weight for Armatus function")
	lambda_is_in := flag.Float64("lambda_is",-1.0,"weight for TopDom/InsScore function")
	//lambda_tt_in := flag.Float64("lambda_tt",-1.0,"weight for HiCSeg function")
	//lambda_hs_in := flag.Float64("lambda_hs",-1.0,"weight for HiCSeg function")
	lambda_di_in := flag.Float64("lambda_di",-1.0,"weight for DI function")

	numtrainsteps := flag.Int("bosteps",100,"number of Bayesian optimization steps")

	ctcffile := flag.String("ctcf","inputdata/CTCFBSDB_all_exp_sites_Sept12_2012_hg19.txt","file with locations of CTCF binding sites")
	outfilename := flag.String("o","","Fileseed or filename for output TAD file") 

	flag.Parse()

	var orighicmaps = map[string]map[hicutil.Pair]float64{}
	var origchrlengths = map[string]int{} 
	var newhicmaps = map[string]map[hicutil.Pair]float64{}
	var newchrlengths = map[string]int{}
	var centromerelocs = map[string][]int{}
	var allgdscores = map[string]map[hicutil.Pair]float64{}
	var meangdscores = map[string][][]float64{}
	var allisscores = map[string][]float64{}
	var mindists = map[string][]int{}
	var chrhicfile string
	var newtads = map[string][][]int{}
	
	var maxsteps int


	allhicfiles,err := filepath.Glob(*hicfile)
	if err != nil{
		fmt.Println("couldnt find files matching filepath")
		log.Fatal(err)
	}
	
	ctcflocs := ftutil.ReadCTCFlocs(*ctcffile, *res)
	ctcfweights := ftutil.CountCTCFweights(ctcflocs)
	//fmt.Println(ctcflocs)
	
	windowsize := 5
	p := 3
	q := 12

	//var diridxtt = map[string][]float64{}
	var avgfreq = map[string][]float64{}
	var exptadfreq = map[string]map[hicutil.Pair]float64{}
	
	var tadmeans = map[string]map[hicutil.Pair]float64{}
	var nontadmean = map[string]float64{}
	//var nontadmeans = map[string]map[hicutil.Pair]float64{}

	var discores = map[string][]float64{}
	var didists = map[string][]int{}

	starttime := time.Now()

	for chrnum := 1; chrnum < 23; chrnum += 1 {
		chrkey := "chr"+strconv.Itoa(chrnum)
		fmt.Println(chrkey)
		for _,f := range allhicfiles {
			if strings.Contains(f,chrkey+".") {
				chrhicfile = f
				break
			}
		}			
		orighicmaps[chrkey],origchrlengths[chrkey],centromerelocs[chrkey] = hicutil.ReadHiCFile(chrhicfile,*res, []float64{})
		newhicmaps[chrkey],newchrlengths[chrkey] = hicutil.RemoveCentromere(orighicmaps[chrkey],centromerelocs[chrkey],origchrlengths[chrkey])
		//fmt.Println(centromerelocs[chrkey])

		//allgdscores[chrkey], meangdscores[chrkey], allisscores[chrkey], mindists[chrkey], diridxtt[chrkey], avgfreq[chrkey], exptadfreq[chrkey], tadmeans[chrkey], nontadmean[chrkey], discores[chrkey], didists[chrkey] = ftutil.PrecomputeAll_nott(newhicmaps[chrkey], newchrlengths[chrkey], windowsize, p, q, *res)
		//fmt.Println(discores[chrkey])

		allgdscores[chrkey], meangdscores[chrkey], allisscores[chrkey], mindists[chrkey], _, avgfreq[chrkey], exptadfreq[chrkey], tadmeans[chrkey], nontadmean[chrkey], discores[chrkey], didists[chrkey] = ftutil.PrecomputeAll_nott(newhicmaps[chrkey], newchrlengths[chrkey], windowsize, p, q, *res)


		newtads[chrkey] = ftutil.CreateRandomTADs(newchrlengths[chrkey])//, centromerelocs[chrkey])
		//fmt.Println(chrkey,centromerelocs[chrkey])
		//if chrnum > 25 {fmt.Println(p,q,diridxtt, avgfreq, exptadfreq) }
	}
	fmt.Println("preprocessing took", time.Since(starttime))
	//os.Exit(1)
	//fmt.Println(p,q,windowsize)

	var optgamma_a, optalpha_is float64 //optgamma_tt float64
	var optlambda_a, optlambda_is, optlambda_di float64 //optlambda_tt, optlambda_hs, optlambda_di float64
	var gamma_a, alpha_is float64 //gamma_tt float64
	var lambda_a, lambda_is, lambda_di float64 //lambda_tt, lambda_hs, lambda_di float64

	if *trainflag {
		fmt.Println("training parameters")
		trialchan := make(chan goptuna.FrozenTrial, 8)
		study, err := goptuna.CreateStudy(
			"goptuna-testing",
			goptuna.StudyOptionSampler(tpe.NewSampler()),
			//goptuna.StudyOptionIgnoreObjectiveErr(true),
			goptuna.StudyOptionSetTrialNotifyChannel(trialchan),
		)
		if err != nil { fmt.Println("error initializing study") }

		//define objective function
		objective := func(trial goptuna.Trial) (float64,error) {

			// NEED TO ADD CONDITION IF TRAINING ALL PARAMS, ONE MUST = 1

			if *gamma_armatus_in < 0 {
				gamma_a,_ = trial.SuggestUniform("gamma_a",0,3)
			} else { 
				gamma_a = *gamma_armatus_in
				optgamma_a = *gamma_armatus_in
			}
			if *alpha_is_in < 0 {
				alpha_is,_ = trial.SuggestUniform("alpha_is",0,200)
			} else { 
				alpha_is = *alpha_is_in
				optalpha_is = *alpha_is_in 
			}
			/*if *gamma_tt_in < 0 {
				gamma_tt,_ = trial.SuggestUniform("gamma_tt",0,1)
			} else { 
				gamma_tt = *gamma_tt_in
				optgamma_tt = *gamma_tt_in 
			}*/
			
			if *lambda_armatus_in < 0 {
				lambda_a = 50.0 //- lambda_is - lambda_tt - lambda_hs - lambda_di
				optlambda_a = 50.0 //- optlambda_is - optlambda_tt - optlambda_hs - optlambda_di
				//lambda_a,_ = trial.SuggestUniform("lambda_a",0,100)
			} else {
				lambda_a = *lambda_armatus_in
				optlambda_a = *lambda_armatus_in 
			}
			if *lambda_is_in < 0 {
				//lambda_is = 1.0
				//optlambda_is = 1.0
				lambda_is,_ = trial.SuggestUniform("lambda_is",0,100)
			} else {
				lambda_is = *lambda_is_in
				optlambda_is = *lambda_is_in 
			}
			/*if *lambda_tt_in < 0 {
				lambda_tt,_ = trial.SuggestUniform("lambda_tt",0,100)
			} else {
				lambda_tt = *lambda_tt_in
				optlambda_tt = *lambda_tt_in
			}
			if *lambda_hs_in < 0 {
				lambda_hs,_ = trial.SuggestUniform("lambda_hs",0,100)
			} else {
				lambda_hs = *lambda_hs_in
				optlambda_hs = *lambda_hs_in
			}*/
			if *lambda_di_in < 0 {
				lambda_di,_ = trial.SuggestUniform("lambda_di",0,100)
			} else { 
				lambda_di = *lambda_di_in
				optlambda_di = *lambda_di_in
			}
		
			avgctcf := 0.0
			numchr := 0
			for chrnum := 1; chrnum < 23; chrnum += 2 { //should be able to make this a go routine to speed things up
				chrkey := "chr"+strconv.Itoa(chrnum)
				fmt.Println(chrkey)
			
				maxsteps = newchrlengths[chrkey]/2
				//newtads[chrkey] = runMCMC(newtads[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], gamma_a, allisscores[chrkey], mindists[chrkey], alpha_is, lambda_a, lambda_is, lambda_tt, lambda_hs, lambda_di, maxsteps, gamma_tt, diridxtt[chrkey], avgfreq[chrkey], exptadfreq[chrkey], nontadmean[chrkey], tadmeans[chrkey], discores[chrkey], didists[chrkey], newhicmaps[chrkey]) 
				//newtads[chrkey],_ = ftutil.RunDynamicProgram(newhicmaps[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], gamma_a, allisscores[chrkey], mindists[chrkey], alpha_is, lambda_a, lambda_is, lambda_tt, lambda_hs, lambda_di, maxsteps, gamma_tt, diridxtt[chrkey], avgfreq[chrkey], exptadfreq[chrkey], nontadmean[chrkey], tadmeans[chrkey], discores[chrkey], didists[chrkey])

				newtads[chrkey] = ftutil.RunDynamicProgram(newhicmaps[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], gamma_a, allisscores[chrkey], mindists[chrkey], alpha_is, lambda_a, lambda_is, lambda_di, discores[chrkey], didists[chrkey])

				avgctcf += ftutil.AvgCTCFload(newtads[chrkey], ctcfweights[chrkey], centromerelocs[chrkey])
				numchr +=1
			}

			//compute avg ctcf load
			return -avgctcf/float64(numchr), nil
		}

		// to save data from optimization:
		/*var wg sync.WaitGroup
		var trialdetails = make([]goptuna.FrozenTrial,*numtrainsteps)
		tnum := 0
		wg.Add(2)
		go func() {
			defer wg.Done()
			err = study.Optimize(objective, *numtrainsteps)
			close(trialchan)
		}()
		go func() {
			defer wg.Done()
			for t := range trialchan {
				trialdetails[tnum] = t
				tnum +=1
				//log.Println("trial", t)
			}
		}()
		wg.Wait()*/
		
		// Run an objective function n times to find a global minimum. - coment these lines to write BO details to file
		err = study.Optimize(objective, *numtrainsteps)
		if err != nil { fmt.Println("error optimizing") }
		

		// Print the best evaluation value and the parameters.
		v, _ := study.GetBestValue()
		optparams, _ := study.GetBestParams()
		fmt.Println(v,optparams)
	
		if *gamma_armatus_in < 0 { optgamma_a = optparams["gamma_a"].(float64) }
		if *alpha_is_in < 0 { optalpha_is = optparams["alpha_is"].(float64) }
		//if *gamma_tt_in < 0 { optgamma_tt = optparams["gamma_tt"].(float64) }
		//if *lambda_armatus_in < 0 { optlambda_a = optparams["lambda_a"].(float64) }
		if *lambda_is_in < 0 { optlambda_is = optparams["lambda_is"].(float64) }
		//if *lambda_tt_in < 0 { optlambda_tt = optparams["lambda_tt"].(float64) }
		//if *lambda_hs_in < 0 { optlambda_hs = optparams["lambda_hs"].(float64) }
		if *lambda_di_in < 0 { optlambda_di = optparams["lambda_di"].(float64) }

		//optlambda_a = 1.0 - optlambda_is - optlambda_tt - optlambda_hs - optlambda_di
		
		log.Printf("Best evaluation value=%f (alpha_is=%f, gamma_a=%f, lambda_a=%f, lambda_is=%f, lambda_di=%f)", v, optalpha_is, optgamma_a, optlambda_a, optlambda_is, optlambda_di )
		//fmt.Println(trialdetails)
		//filename := *outfilename+strconv.Itoa(*numtrainsteps)+"steps_BOdetails.txt"
		//ftutil.WriteBOInfoToFile(trialdetails,filename)
		//os.Exit(1)
	
		// need to make sure this runs all functions - uncomment out HiCSeg and TADTree*/
	} else {
		fmt.Println("using given parameters")

		if *gamma_armatus_in > -1 {
			optgamma_a = *gamma_armatus_in // armatus param
		} else { optgamma_a = 1.2 } 
		if *alpha_is_in > -1 {
			optalpha_is = *alpha_is_in // IS score internal param
		} else { optalpha_is = 2.5 }
		/*if *gamma_tt_in > -1 {
			optgamma_tt = *gamma_tt_in // TADTree internal weight
		} else { optgamma_tt = 0.0 }*/
		if *lambda_armatus_in > -1 {
			optlambda_a = *lambda_armatus_in // armatus score weight
		} else { optlambda_a = 50.0 }
		if *lambda_is_in > -1 {
			optlambda_is = *lambda_is_in // IS score weight
		} else { optlambda_is = 36.0 } // old default: 15
		/*if *lambda_hs_in > -1 {
			optlambda_hs = *lambda_hs_in // HiCSeg score weight
		} else { optlambda_hs = 0.0 }
		if *lambda_tt_in > -1 {
			optlambda_tt = *lambda_tt_in // TTscore weight
		} else { optlambda_tt = 0.0 }*/
		if *lambda_di_in > -1 {
			optlambda_di = *lambda_di_in // DI score weight
		} else { optlambda_di = 0.62 } // old default: 2

	}

	//compute TADs with opt params and write them to file
	//for randseed := 1; randseed < 11; randseed++ {
	//	rand.Seed(int64(randseed))
	for chrnum := 1; chrnum < 23; chrnum += 1 {
		chrkey := "chr"+strconv.Itoa(chrnum)
		maxsteps = 3*newchrlengths[chrkey]
		//filename := *outfilename+chrkey+"_seed"+strconv.Itoa(randseed)+".txt"
		filename := *outfilename+chrkey+".txt"

		fmt.Println(chrkey)

		//opttads := runMCMC(newtads[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], optgamma_a, allisscores[chrkey], mindists[chrkey], optalpha_is, optlambda_a, optlambda_is, optlambda_tt, optlambda_hs, optlambda_di, maxsteps, optgamma_tt, diridxtt[chrkey], avgfreq[chrkey], exptadfreq[chrkey], nontadmean[chrkey], tadmeans[chrkey], discores[chrkey], didists[chrkey], newhicmaps[chrkey])

		opttads := ftutil.RunDynamicProgram(newhicmaps[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], optgamma_a, allisscores[chrkey], mindists[chrkey], optalpha_is, optlambda_a, optlambda_is, optlambda_di, discores[chrkey], didists[chrkey])

		//fmt.Println(opttads)
		//ftutil.WriteScoresToFile(alltadscores,*outfilename+chrkey+"_tadscores.txt")
		//os.Exit(1)

		cleantads := hicutil.PostprocessTADs(opttads, allgdscores[chrkey], centromerelocs[chrkey])

		fmt.Println(cleantads)

		// also write param values to file? or would this complicate reading tad files for other programs?
		ftutil.WriteTADsToFile(cleantads,*res,filename)
		log.Printf("Parameter values used (alpha_is=%f, gamma_armatus=%f, lambda_armatus=%f, lambda_IS=%f, lambda_DI=%f)", optalpha_is, optgamma_a, optlambda_a, optlambda_is, optlambda_di )

	}
}


func arange(start, stop, step float64) []float64 {
    N := int(math.Ceil((stop - start+1) / step));
    rnge := make([]float64, N, N)
    i := 0
    for x := start; x < stop; x += step {
        rnge[i] = x;
        i += 1
    }
    return rnge
}

/*func subsampHiC(hicmap map[hicutil.Pair]float64, frac float64) (map[hicutil.Pair]float64) {
	// frac is a number between 0 and 1 representing how much of the matrix to keep

	//runtime.GOMAXPROCS(1)

	totalnonzero := len(hicmap)
	hicmap_subsamp := hicmap

	keys := make([]hicutil.Pair, totalnonzero)
	i := 0
	for k := range hicmap {
		keys[i] = k
		i++
	}
	numtoremove := int((1.0-frac)*float64(totalnonzero))
	// randomly shuffle keys in case iteration order isn't random
	rand.Shuffle(totalnonzero, func(i,j int) {keys[i],keys[j] = keys[j], keys[i]})
	for i := 0; i < numtoremove; i++ {
		delete(hicmap_subsamp, keys[i])
	}

	return hicmap_subsamp
}*/

/*func precomputeAll(hicmap map[hicutil.Pair]float64, chrlength int, iswindowsize int, p int, q int,res int) (map[hicutil.Pair]float64, [][]float64, []float64, []int, []float64, []float64, map[hicutil.Pair]float64, map[hicutil.Pair]float64, float64, []float64) {

	//res int, diwindow int)
	
	//2nd to last return: map[hicutil.Pair]

	// loop through Hi-C matrix once to precompute everything needed
	//for Armatus scoring
	gdscores := make(map[hicutil.Pair]float64) // contain scores for each (a,b) interval on HiC matrix NOT NORMALIZED BY LENGTH^GAMMA 
	meangdscores := make([][]float64, chrlength) // sum of scores for all tads of given length (length = array index) and number of intervals of that length (to compute mean)

	//for TopDom/Ins Score scoring
	allisscores := make([]float64, chrlength)
	sumisscores := 0.0
	minlocs := make([]int, chrlength-2*iswindowsize)

	//for TADTree scoring
	diridxtt := make([]float64, chrlength)
	avgfreq := make([]float64, chrlength)
	exptadfreq := make(map[hicutil.Pair]float64)
	linregvals :=make(map[hicutil.Pair]LinRegIntermed)
	
	//for HiCSeg scoring
	nontadmean := 0.0
	tadmeans := make(map[hicutil.Pair]float64)
	nontadmeans := make(map[hicutil.Pair]float64)
	//nontadmeans2 := make(map[hicutil.Pair]float64)
	colsums := make([][]float64,chrlength)
	
	//for DI scoring
	diridx := make([]float64,chrlength)
	diridxscores := make([]float64,chrlength)
	dilen := 2000000/res
	//diwindow := 5

	cornercount := 0
	for a := 0; a < chrlength/4; a++ {
		for b:= 3*chrlength/4; b < chrlength; b++ {
			nontadmean += hicmap[hicutil.Pair{a,b}]
			cornercount++
		}
	}
	nontadmean = nontadmean/float64(cornercount)
	//fmt.Println(nontadmean)

	meangdscores[0] = make([]float64,2)
	colsums[0] = make([]float64,chrlength)
	for a := 0; a < chrlength; a++ {
		// initialize armatus scoring
		gdscores[hicutil.Pair{a,a}] = hicmap[hicutil.Pair{a,a}]
		meangdscores[0][0] += gdscores[hicutil.Pair{a,a}]
		meangdscores[0][1] += 1.0

		//tad means for hicseg scoring
		tadmeans[hicutil.Pair{a,a}] = hicmap[hicutil.Pair{a,a}]

		colsums[0][a] = hicmap[hicutil.Pair{0,a}]

		//compute directionality index for TADtree objective 
		diridxtt[a] = computeDirIdxTT(hicmap,a,p,q,chrlength)
		avgfreq[0] += hicmap[hicutil.Pair{a,a}] * (1.0/float64(chrlength))
		//if hicmap[hicutil.Pair{a,a}] > 0 {fmt.Println(hicmap[hicutil.Pair{a,a}])}
	}
	//fmt.Println("done initializing")
	colsum := 0.0
	minidx := 0
	for b := 1; b < chrlength; b++ {
		colsums[b] = make([]float64,chrlength)
		//directionality index
		if b > dilen && b < chrlength-dilen {
			diridx[b] = computeDI(b,dilen,hicmap)
		} else {
			diridx[b] = -1
		}
		//fmt.Println(b,diridxscores[b])
		//fmt.Println(b)
		colsum = hicmap[hicutil.Pair{b,b}]
		if b >= iswindowsize/2 && b < chrlength-iswindowsize/2+1 {
			for c := b; c < b+iswindowsize/2+1; c++ {
				allisscores[b] += hicmap[hicutil.Pair{b,c}]
			}
		} 
		colsums[b][chrlength-1] = colsums[b-1][chrlength-1]+hicmap[hicutil.Pair{b,chrlength-1}]
		for a := b-1; a > -1; a-- {
			//nontadmeans[hicutil.Pair{a,b}] += rowsum/float64((a-1)*b)
			x := int(math.Abs(float64(a-b+1)))

			//hicseg nontadmeans computed per rectangle instead of 1 overall
			if x > 0 {colsums[x][b] = colsums[x-1][b] + hicmap[hicutil.Pair{x,b}]}
			if b > 2 && a > 1 {nontadmeans[hicutil.Pair{a,b-1}] = (nontadmeans[hicutil.Pair{a,b-2}]*float64((a-1)*(b-2)) + colsums[a-1][b-1])/float64((a-1)*(b-1))}
			if b == chrlength-1 && a > 0 { nontadmeans[hicutil.Pair{a,b}] = (nontadmeans[hicutil.Pair{a,b-1}]*float64((a-1)*(b-1)) + colsums[a-1][b])/float64((a-1)*b) }

			/*for row := 0; row < a; row++ {
				for col := a; col <= b; col++ {
					nontadmeans[hicutil.Pair{a,b}] += hicmap[hicutil.Pair{row,col}]//float64((a-1)*b)
				}
			}
			if b > 1 && b < 100 {
				fmt.Println(a, b-1, nontadmeans[hicutil.Pair{a,b-1}],nontadmeans2[hicutil.Pair{a,b-1}])
			} else if b > 100 {
				os.Exit(1)
			}*/


			//if b-a > chrlength-100 { continue }
/*			avgfreq[a] += hicmap[hicutil.Pair{b,b+a}] * (1.0/float64(chrlength-a))

			// Armatus scoring
			colsum += hicmap[hicutil.Pair{a,b}]
			gdscores[hicutil.Pair{a,b}] = gdscores[hicutil.Pair{a,b-1}] + colsum
			
			//fmt.Println("a=",a,"b=",b,"score=",allscores[hicutil.Pair{a,b}])
			if len(meangdscores[b-a]) == 0 {
				meangdscores[b-a] = make([]float64,2)
			}
			meangdscores[b-a][0] += gdscores[hicutil.Pair{a,b}]
			meangdscores[b-a][1] += 1.0
		
			// hicseg scoring
			tadmeans[hicutil.Pair{a,b}] = gdscores[hicutil.Pair{a,b}]/triangleNum(b-a+1)
			//fmt.Println(a,b,tadmeans[hicutil.Pair{a,b}],triangleNum(b-a+1))
			//fmt.Println(a,b,"hicval =",hicmap[hicutil.Pair{a,b}])

			// insulation score/ topdom scoring
			if b >= iswindowsize/2 && b < chrlength-iswindowsize/2+1 {
				if a >= b-iswindowsize/2 {
					for c := b; c < b+iswindowsize/2+1; c++{
						allisscores[b] += hicmap[hicutil.Pair{a,c}]
						//if b == 20 {fmt.Println(a,c,hicmap[hicutil.Pair{a,c}])}
					} // maybe missing the row where a = b because a starts at b-1
					allisscores[b] += hicmap[hicutil.Pair{b,a+iswindowsize/2+1}]
					//if b == 20 {fmt.Println(b,a+iswindowsize/2)}
				}
				if allisscores[b] < 0 {
					fmt.Println("got negative ins score")
					fmt.Println(allisscores[b],b)
					os.Exit(1)
				}
				sumisscores += allisscores[b]
			}
		}
		if b > 1 &&  allisscores[b] > allisscores[b-1] && allisscores[b-2] > allisscores[b-1] {
			minlocs[minidx] = b-1
			minidx += 1
			//fmt.Println("local min found at",b)
		}
	}
	//for i:= 0; i < len(avgfreq); i++ {
	//	fmt.Println(i,avgfreq[i])
	//}
	minlocs[minidx] = chrlength-1
	minlocs = minlocs[:minidx+1] // note this just finds all things that are less than the two before and 1 after -- is this the best way to identify minima? maybe consider changing or at least testing to see if it works decently
	//fmt.Println("done with primary calculations")
	// normalize all scores as log_2 (binscore/meanscores)
	normfactor := sumisscores/float64(chrlength - 2*iswindowsize)
	mindists  := make([]int, chrlength)
	lastmin := 0
	nextmin := minlocs[0]
	minidx = 0
	//var xvals []float64
	//var yvals []float64
	//disum := 0.0
	
	for i := 0; i < chrlength; i++ {
		// compute differences in DI for DI scoring
		if i > 0 {
			diridxscores[i] = diridx[i] - diridx[i-1]
			//for j := 0; j < diwindow; j++{
				//diridxscores[i] +=  diridx[i+1+j] - diridx[i-j]
			//}
			//disum += diridxscores[i]
		} else { diridxscores[i] = 0 }
		//fmt.Println(i,diridx[i],diridxscores[i])
		//fmt.Println(i)
		if i >= iswindowsize/2 && i < chrlength-iswindowsize/2+1 {
			if allisscores[i] != 0 {
				allisscores[i] = math.Log2(allisscores[i]/normfactor)
				if math.IsNaN(allisscores[i]) {
					fmt.Println("uh oh, an insulation score is NaN")
					os.Exit(1)
				}
			} else if i+1 < chrlength && allisscores[i+1] != 0 {
				allisscores[i] = math.Log2(allisscores[i+1]/normfactor)-1
			} else {
				allisscores[i] = -10000
			}
			if math.IsInf(allisscores[i],0) {
				fmt.Println("IS score of infinity")
				fmt.Println(i, normfactor)
				os.Exit(1)
			}
		}
		
		/*xvals = []float64{0.0}
		if avgfreq[0] > 0 {
			yvals = []float64{hicmap[hicutil.Pair{i,i}] / avgfreq[0]}
		} else {
			yvals = []float64{0.0}
		}*/
		//fmt.Println("done with ins score")
/*		mindists,lastmin,nextmin,minidx = updateMinDists(minlocs,mindists,i,lastmin,nextmin,minidx,chrlength)
		//fmt.Println("done updating min")
		// compute Atilde from TADtree - linear regression for every tad is crazy slow, how to make this better??
		//exptadfreq[i],prevxvals,prevyvals = calcExpectedContactFreq(hicmap,prevxvals,prevyvals,i,chrlength,avgfreq)
		for j := i+1; j < chrlength; j++ {
			//fmt.Println(j)

			//if j-i > chrlength-100 { continue }
			//exptadfreq[hicutil.Pair{i,j}] = calcExpectedContactFreq(hicmap,i,j,avgfreq)
			// update linregvals for all new points
			linregvals[hicutil.Pair{i,j}] = linregvals[hicutil.Pair{i,j-1}]
			newvals := linregvals[hicutil.Pair{i,j}]
			//var avgfreqval float64
			for k := 1; k <= j - i; k++ {
				var newy float64
				//xvals = append(xvals,float64(k))
				if avgfreq[k] == 0 {
					newy = 0
				} else {
					newy = hicmap[hicutil.Pair{j-k,j}] / avgfreq[k]
				}
				if math.IsNaN(newy) { 
					fmt.Println(i,j,k,hicmap[hicutil.Pair{j-k,j}],avgfreq[k])
					os.Exit(1)
				}
				//avgfreqval = avgfreq[k]
				//yvals = append(yvals,newy) // change to be 0 if avgfreq == 0?
				//fmt.Println(k,avgfreq[k],hicmap[hicutil.Pair{j-k,j}])
				//fmt.Println(hicmap[hicutil.Pair{j-k,j}] / avgfreq[k])

				//newy := hicmap[hicutil.Pair{j-k,j}]/avgfreq[k]
				newvals.n += 1
				newvals.dx = float64(k) - newvals.meanX
				newvals.dy = newy - newvals.meanY
				newvals.varX += (((newvals.n-1)/newvals.n)*newvals.dx*newvals.dx - newvals.varX)/newvals.n
				newvals.covXY += (((newvals.n-1)/newvals.n)*newvals.dx*newvals.dy - newvals.covXY)/newvals.n
				newvals.meanX += newvals.dx/newvals.n
				newvals.meanY += newvals.dy/newvals.n
			}
			
			//fmt.Println(xvals,yvals)
			//if j > 5 {os.Exit(1)}
			// compute lin reg params
			delta := newvals.getA()
			//if math.IsNaN(delta) {
			//	delta = 0
			//	newvals
			//}
			beta := newvals.getB()
			if math.IsNaN(delta) || math.IsNaN(beta) {
				fmt.Println(delta,beta)
				fmt.Println(newvals.n)
				fmt.Println(newvals.dx)
				fmt.Println(newvals.dy)
				fmt.Println(newvals.varX)
				fmt.Println(newvals.covXY)
				fmt.Println(newvals.meanX)
				fmt.Println(newvals.meanY)
			}
			// need to compare to delta, beta from normal linear regression
			/*beta2, delta2 := stat.LinearRegression(xvals, yvals, nil, false)
			if math.IsNaN(delta2) || math.IsNaN(beta2) {
				fmt.Println("linear regression returned nan")
				fmt.Println(xvals,yvals,avgfreqval,hicmap[hicutil.Pair{i,j}])
				os.Exit(1)
			}
			if delta - delta2 > 0.1 || beta - beta2 > 0.1 {
				fmt.Println("linear regression results dont match up")
				fmt.Println(i,j)
				fmt.Println(delta,delta2)
				fmt.Println(beta,beta2)
			}*/
/*			linregvals[hicutil.Pair{i,j}] = newvals
			// calc actual exptadfreq value
			exptadfreq[hicutil.Pair{i,j}] = (float64(j-i)*delta + beta)*avgfreq[j-i] // should this be i-j or i-j+1???
			if math.IsNaN(exptadfreq[hicutil.Pair{i,j}]) {
				fmt.Println("exp tad freq NaN")
				fmt.Println("HiC coords",i,j)
				fmt.Println(delta,beta)
				fmt.Println(avgfreq[j-i])
				os.Exit(1)
			}
		}
		// do we need ANOTHER loop through everything to compute full TAD scores from linear regression values? the TADTree scoring function is the absolute worst. 3 quadratic loops -_-
		
		//if i > 1 {os.Exit(1)}
	}
	return gdscores, meangdscores, allisscores, mindists, diridxtt, avgfreq, exptadfreq, tadmeans, nontadmean, diridx//scores
}*/



/*func calcExpectedContactFreq(hicmap map[hicutil.Pair]float64, xvals []float64, yvals []float64, tadlen int, chrlength int, avgfreq []float64) (float64, []float64, []float64) {

	//xvals := make([]float64, int(math.Pow(float64(chrlength),2)))
	//yvals := make([]float64, int(math.Pow(float64(chrlength),2)))
	//a := 0
	for idx := 1; idx < chrlength-tadlen; idx++ {
		xvals = append(xvals, float64(idx+tadlen))
		yvals = append(yvals, hicmap[hicutil.Pair{idx,idx+tadlen}] / avgfreq[tadlen])
		//for idx2 := idx1; idx2 < idx1+tadlen+1; idx2++ {
		//xvals[a] = float64(idx2 - idx1)
		//yvals[a] = hicmap[hicutil.Pair{idx1,idx2}] / avgfreq[idx2-idx1]
		//	a++
		//}
	}
	//xvals = xvals[:a]
	//yvals = yvals[:a]
	// linear regression to estimate parameters
	delta, beta := stat.LinearRegression(xvals, yvals, nil, false)
	//put together (eqn 2 from TADtree paper)
	expval := (float64(tadlen)*delta + beta)*avgfreq[tadlen]
	return expval,xvals,yvals
}*/

/*func calcExpectedContactFreq(hicmap map[hicutil.Pair]float64, i int, j int, avgfreq []float64) float64 {

	tadlen := j-i+1
	xvals := make([]float64, tadlen*(tadlen+1)/2)
	yvals := make([]float64, tadlen*(tadlen+1)/2)
	a := 0
	for idx1 := i; idx1 < j; idx1++ {
		for idx2 := idx1; idx2 < j+1; idx2++ {
			xvals[a] = float64(idx2 - idx1)
			yvals[a] = hicmap[hicutil.Pair{idx1,idx2}] / avgfreq[idx2-idx1]
			a++
		}
	}
	// linear regression to estimate parameters
	delta, beta := stat.LinearRegression(xvals, yvals, nil, false)
	//put together (eqn 2 from TADtree paper)
	expval := (float64(tadlen-1)*delta + beta)*avgfreq[tadlen-1]
	return expval
}*/


/*func computeAllGDScores(hicmap map[hicutil.Pair]float64, chrlength int) (map[hicutil.Pair]float64, [][]float64) {
	// precompute values needed for graph density (Armatus) score

	//runtime.GOMAXPROCS(1)

	allscores := make(map[hicutil.Pair]float64) // contain scores for each (a,b) interval on HiC matrix NOT NORMALIZED BY LENGTH^GAMMA 
	meanscores := make([][]float64, chrlength-99) // sum of scores for all tads of given length (length = array index) and number of intervals of that length (to compute mean)

	meanscores[0] = make([]float64,2)
	for a := 0; a < chrlength; a++ {
		allscores[hicutil.Pair{a,a}] = hicmap[hicutil.Pair{a,a}]
		meanscores[0][0] += allscores[hicutil.Pair{a,a}]
		meanscores[0][1] += 1.0
	}
	
	colsum := 0.0
	for b := 1; b < chrlength; b++ {
		colsum = hicmap[hicutil.Pair{b,b}]
		for a := b-1; a > -1; a-- {
			//if b-a > chrlength-100 { continue }
			colsum += hicmap[hicutil.Pair{a,b}]
			allscores[hicutil.Pair{a,b}] = allscores[hicutil.Pair{a,b-1}] + colsum
			//fmt.Println("a=",a,"b=",b,"score=",allscores[hicutil.Pair{a,b}])
			if len(meanscores[b-a]) == 0 {
				meanscores[b-a] = make([]float64,2)
			}
			meanscores[b-a][0] += allscores[hicutil.Pair{a,b}]
			meanscores[b-a][1] += 1.0
		}
	}
	
	return allscores, meanscores
}*/

/*func computeParticularGDScore(allgdscores map[hicutil.Pair]float64, meanscores *[][]float64, tadset *[][]int, gamma float64) (float64,time.Duration) {
	// compute graph density (Armatus) score for particular tad set

	//runtime.GOMAXPROCS(1)
	timestart := time.Now()
	var tadscore,normfactor float64

	totalscore := 0.0

	for _, tad := range *tadset {
		//if tad[0] == tad[1] {continue} //is this right?

		tadscore = allgdscores[hicutil.Pair{tad[0],tad[1]}]
		//fmt.Println("tadscore =",tadscore)
		//if tad[0] == tad[1] {
		//	normfactor = 1.0
		//} else {
		normfactor = math.Pow(float64(tad[1]-tad[0]+1), gamma) // armatus code has b-a+1 
		//}
		//fmt.Println("normfactor =",normfactor)
		if tad[1]-tad[0] < 0 {
			fmt.Println("mistake in TAD set detected")
			fmt.Println(tadset, tad)
			os.Exit(1)
		}
		if tad[1]-tad[0] < len(*meanscores) {
			totalscore += tadscore/normfactor - ( (*meanscores)[tad[1]-tad[0]][0] / ( (*meanscores)[tad[1]-tad[0]][1] * normfactor ) )
		} // tads longer than chrlength-99 get scores of 0 this way - should it be negative?
		//fmt.Println("totalscore =",totalscore)
	}
	return totalscore,time.Since(timestart)
}*/

/*func computeParticularGDScoreByBdy(allgdscores map[hicutil.Pair]float64, meanscores [][]float64, tadset [][]int, gamma float64) ([]float64) {
	// compute graph density (Armatus) score for particular tad set, returning score per boundary for comparison w/ ctcf sites

	//runtime.GOMAXPROCS(1)

	var tadscore,normfactor,prevtadscore float64

	bdyscores := make([]float64,len(tadset)-2)
	//what to do about centromeres?

	for tadidx, tad := range tadset {
		tadscore = allgdscores[hicutil.Pair{tad[0],tad[1]}]
		normfactor = math.Pow(float64(tad[1]-tad[0]+1), gamma) // armatus code has b-a+1 
		tadscore = tadscore/normfactor - ( meanscores[tad[1]-tad[0]][0] / ( meanscores[tad[1]-tad[0]][1] * normfactor ) )
		
		if tadidx > 0 && tadidx < len(bdyscores) {
			bdysc := (tadscore + prevtadscore) / 2
			bdyscores[tadidx-1] = bdysc
			prevtadscore = tadscore
		}
	}
	return bdyscores
}*/

/*func computeAllInsulationVals(hicmap map[hicutil.Pair]float64, chrlength int, windowsize int) ([]float64, []int) {

	//runtime.GOMAXPROCS(1)

	// precompute values needed for insulation score

	//might need to rethink this formulation since alpha --> 0

	allisscores := make([]float64, chrlength)
	sumisscores := 0.0
	minlocs := make([]int, chrlength-2*windowsize)

	/*for a := 0; a < windowsize; a++ {
		for b := 0; b < windowsize; b++ {
			allisscores[windowsize-1] = hicmap[hicutil.Pair{a,b}]
			sumisscores += allisscores[windowsize-1]
		}
	}
	minidx := 0
	for i := windowsize; i < chrlength-windowsize; i++ {
		toprow := 0.0
		for b := i-windowsize; b < i; b++ {
			toprow += hicmap[hicutil.Pair{i-windowsize,b}]
		}
		newcol := 0.0
		for a := i-windowsize+1; a < i+1; a ++ {
			newcol += hicmap[hicutil.Pair{a,i}]
		}
		allisscores[i] = allisscores[i-1] - toprow + newcol
		if allisscores[i] < 0 {
			fmt.Println(allisscores[i-1],i,toprow,newcol)
		}
		sumisscores += allisscores[i]
		if allisscores[i] > allisscores[i-1] && allisscores[i-2] > allisscores[i-1] {
			minlocs[minidx] = i-1
			minidx += 1
		}
	}*/
	/*minidx := 0
	for i := windowsize/2; i < chrlength-windowsize/2+1; i++ {
		for a := i-windowsize/2; a < i+windowsize/2+1; a++ {
			for b := i-windowsize/2; b < i+windowsize/2+1; b++{
				allisscores[i] += hicmap[hicutil.Pair{a,b}]
			}
		} // this could be made more efficient by using idea from above
		if allisscores[i] < 0 {
			fmt.Println("got negative ins score here")
			fmt.Println(allisscores[i],i)
			for a := i-windowsize/2; a < i+windowsize/2+1; a++ {
				for b := i-windowsize/2; b < i+windowsize/2+1; b++{
					fmt.Println(a,b,hicmap[hicutil.Pair{a,b}])
				}
			}
			os.Exit(1)
		}
		sumisscores += allisscores[i]
		if allisscores[i] > allisscores[i-1] && allisscores[i-2] > allisscores[i-1] {
			minlocs[minidx] = i-1
			minidx += 1
		}
	}
	minlocs = minlocs[:minidx] // note this just finds all things that are less than the two before and 1 after -- is this the best way to identify minima? maybe consider changing or at least testing to see if it works decently

	// normalize all scores as log_2 (binscore/meanscores)
	normfactor := sumisscores/float64(chrlength - 2*windowsize)
	for i := windowsize/2; i < chrlength-windowsize/2+1; i++ {
		if allisscores[i] != 0 {
			allisscores[i] = math.Log2(allisscores[i]/normfactor)
			if math.IsNaN(allisscores[i]) {
				fmt.Println("uh oh, an insulation score is NaN")
				os.Exit(1)
			}
		}
	}

	mindists  := make([]int, chrlength)
	lastmin := 0
	nextmin := minlocs[0]
	minidx = 0
	// vector of distance to nearest minima
	for i := 0; i < chrlength; i++ {
		if lastmin == 0 {
			mindists[i] = nextmin - i
		}
		if i == nextmin {
			mindists[i] = 0
			lastmin = i
			if len(minlocs) > minidx+1 {
				nextmin = minlocs[minidx+1]
			} else {
				nextmin = chrlength*2
			}
			minidx += 1
		} else if (i - lastmin >= nextmin - i) && lastmin != 0 {
			mindists[i] = nextmin - i
		} else if (nextmin - i > i - lastmin) && lastmin != 0 {
			mindists[i] = i - lastmin
		}
		//fmt.Println(lastmin,nextmin,i,mindists[i])
	}

	// aren't we trying to minimize this score? do we need to reverse somehow for maximizing?? ********************************************
	return allisscores, mindists
}*/

/*func computeParticularInsScore(allisscores *[]float64, mindists *[]int, tadset *[][]int, alpha float64) (float64,time.Duration) {

	//runtime.GOMAXPROCS(1)
	starttime := time.Now()

	// compute insulation score for a particular tad set
	totalscore := 0.0
	steepscore := 0.0
	chrlength := len(*allisscores)

	closestbdy := 0
	tadidx := 0
	//weight := 0.0
	for i := 0; i < chrlength; i++ {
		score := math.Abs(math.Abs(float64(i-closestbdy)) - float64((*mindists)[i]))
		// add weight for slope around i - worse to be off on strong boundaries than weak ones
		// this doesn't make sense for all points - only want steepness to matter around actual boundary points
		if (i == closestbdy) && i > 9 && i < chrlength-11 {
			centerval := (*allisscores)[i]
			for j := i-10; j < i+11; j++{
				if j < i {
					steepscore += ((*allisscores)[j] - centerval)/10
					//fmt.Println(i,j,((*allisscores)[j] - centerval)/10)
				} else if j > i {
					steepscore += (centerval - (*allisscores)[j])/10
					//fmt.Println(i,j,(centerval-(*allisscores)[j])/10)
				}
			}
			//fmt.Println(i,steepscore)
		}
		if math.IsNaN(score) {
			fmt.Println("insulation score of NaN from dists")
			fmt.Println(i, closestbdy, (*mindists)[i])
			os.Exit(1)
		}
		if math.IsNaN(steepscore){
			fmt.Println("insulation score of NaN from slope")
			fmt.Println(i, (*allisscores)[i-10:i+11])
			os.Exit(1)
		}
		totalscore -= score
		//fmt.Println(i,closestbdy,math.Abs(float64(i-closestbdy)),(*mindists)[i])
		if closestbdy == (*tadset)[tadidx][0] && i+1 >= ((*tadset)[tadidx][0]+(*tadset)[tadidx][1])/2 {
			// switching from closest boundary being prev bdy to next bdy
			closestbdy = (*tadset)[tadidx][1]
		} else if i == (*tadset)[tadidx][1] && tadidx+1 < len(*tadset) {
			// switching to next tad
			closestbdy = (*tadset)[tadidx+1][0]
			tadidx +=1
		}
	}
	/*for _, tad := range *tadset {
		distscore += (*mindists)[tad[0]] + (*mindists)[tad[1]]
		
		//startis := (*allisscores)[tad[0]]
		slopeavg := 0.0
		startbin := tad[0] - 5
		if startbin < 0 {startbin = 0}
		endbin := tad[0] + 6
		if endbin > chrlength - 1 {endbin = chrlength-1}
		numbins := 0
		for i := startbin; i < endbin; i++ {
			//fmt.Println(i, (*allisscores)[i])
			//slopeavg += (*allisscores)[i] - startis
			numbins += 1
			if i < tad[0] {
				slopeavg += (*allisscores)[i] - (*allisscores)[i+1]
			} else {
				slopeavg += (*allisscores)[i+1] - (*allisscores)[i]
			}
		}
		steepscore += slopeavg/float64(numbins)
		//fmt.Println("avg slope around tad bdy:",tad[0],slopeavg/float64(numbins))
		numbins = 0
		slopeavg = 0.0
		//endis := (*allisscores)[tad[1]]
		startbin = tad[1] - 5
		if startbin < 0 {startbin = 0}
		endbin = tad[1] + 6
		if endbin > chrlength - 1 {endbin = chrlength-1}
		for i := startbin; i < endbin; i++ {
			//fmt.Println(i, (*allisscores)[i])
			//slopeavg += (*allisscores)[i] - endis
			numbins += 1
			if i < tad[1] {
				slopeavg += (*allisscores)[i] - (*allisscores)[i+1]
			} else {
				slopeavg += (*allisscores)[i+1] - (*allisscores)[i]
			}
		}
		steepscore += slopeavg/float64(numbins)
		//fmt.Println("avg slope around tad bdy:",tad[1],slopeavg/float64(numbins))
		//fmt.Println(tad)
		//fmt.Println((*mindists)[tad[0]] + (*mindists)[tad[1]], slopeavg/float64(numbins))
		//if tad[0] > 100 { os.Exit(1)}
	}
	totalscore := 1.0/float64(distscore) + alpha*steepscore*/
/*	totalscore += alpha*steepscore
	//fmt.Println(totalscore)
	//os.Exit(1)
	return totalscore,time.Since(starttime)
}*/


/*func computeParticularInsScoreByBdy(allisscores []float64, mindists []int, tadset [][]int, alpha float64) ([]float64) {

	//runtime.GOMAXPROCS(1)

	// compute insulation score for a particular tad set
	//steepscore := 0.0
	chrlength := len(allisscores)

	bdyscores := make([]float64, len(tadset)-2)

	for tadidx, tad := range tadset {
		if tadidx > 0 && tadidx < len(bdyscores) {
			distscore := mindists[tad[0]] + mindists[tad[1]]
		
			startis := allisscores[tad[0]]
			slopeavg := 0.0
			startbin := tad[0] - 5
			if startbin < 0 {startbin = 0}
			endbin := tad[0] + 6
			if endbin > chrlength - 1 {endbin = chrlength-1}
			numbins := 0
			for i := startbin; i < endbin; i++ {
				slopeavg += allisscores[i] - startis
				numbins += 1
			}
			endis := allisscores[tad[1]]
			startbin = tad[1] - 5
			if startbin < 0 {startbin = 0}
			endbin = tad[1] + 6
			if endbin > chrlength - 1 {endbin = chrlength-1}
			for i := startbin; i < endbin; i++ {
				slopeavg += allisscores[i] - endis
				numbins += 1
			}
			steepscore := slopeavg/float64(numbins)
	
			bdysc := 1.0/float64(distscore) + alpha*steepscore
			bdyscores[tadidx-1] = bdysc
		}
	}
	return bdyscores
}*/

/*func computeTTScore(tadset *[][]int, chrlength int, gamma float64, diridxtt []float64, avgfreq []float64, atilde map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, tadtreescores map[hicutil.Pair]float64) (float64, map[hicutil.Pair]float64, time.Duration) {

	starttime := time.Now()
	score := 0.0

	for _, tad := range *tadset {
		if tadscore, ok := tadtreescores[hicutil.Pair{tad[0],tad[1]}]; ok {
			score += tadscore
		} else {
			tadscore := 0.0
			//tad[0] = 547
			//tad[1] = 600
			tadscore += gamma * diridxtt[tad[0]] + gamma * diridxtt[tad[1]]
			//fmt.Println(score,diridxtt[tad[0]],diridxtt[tad[1]],gamma)
			//fmt.Println(tad)
			for a := tad[0]; a <= tad[1]; a++ {
				for b := a; b < chrlength; b++{
					//fmt.Println(a,b)
					if b <= tad[1] {
						//fmt.Println("b <= tad[1]")
						//fmt.Println(math.Pow((atilde[hicutil.Pair{a,b}] - hicmap[hicutil.Pair{a,b}]),2))
						//fmt.Println(atilde[hicutil.Pair{a,b}],hicmap[hicutil.Pair{a,b}])
						tadscore -= (1.0-gamma)*math.Pow((atilde[hicutil.Pair{a,b}] - hicmap[hicutil.Pair{a,b}]),2)
					} else {
						//fmt.Println("else")
						//fmt.Println(math.Pow(avgfreq[b-a] - hicmap[hicutil.Pair{a,b}],2))
						//fmt.Println(avgfreq[b-a], hicmap[hicutil.Pair{a,b}])
						//fmt.Println(score, 1.0-gamma, (1.0-gamma)*math.Pow(avgfreq[b-a] - hicmap[hicutil.Pair{a,b}],2))
						tadscore -= (1.0-gamma)*math.Pow(avgfreq[b-a] - hicmap[hicutil.Pair{a,b}],2)
					}
					if math.IsNaN(tadscore) {
						fmt.Println(tad,a,b)
						os.Exit(1)
					}
				}
			}
			score += tadscore
			tadtreescores[hicutil.Pair{tad[0],tad[1]}] = tadscore
		}
		//fmt.Println(tad,score)
		//os.Exit(1)
	}
	//fmt.Println(score)
	//os.Exit(1)
	return score,tadtreescores,time.Since(starttime)
}*/


/*func computeHiCSegScore(tadset *[][]int, chrlength int, nontadmean float64, tadmeans map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, hicsegscores map[hicutil.Pair]float64) (float64,map[hicutil.Pair]float64,time.Duration) {
//s map[hicutil.Pair]
	starttime := time.Now()
	score := 0.0

	for _, tad := range *tadset {
		//fmt.Println(tad)
		if tadscore, ok := hicsegscores[hicutil.Pair{tad[0],tad[1]}]; ok {
			score += tadscore
		} else {
			tadscore := 0.0
			//fmt.Println(tad)
			for a := 0; a <= tad[1]; a++ { 
				
				//for a := tad[0]; a <= tad[1]; a++ {
				bstart := math.Max(float64(a),float64(tad[0]))
				//fmt.Println("row =",a,"mincol =",bstart)
				for b := int(bstart); b <= tad[1]; b++ { 
				//for b := a; b < chrlength; b++ {
					if a >= tad[0] { 
					//if b <= tad[1] {
						
						//fmt.Println("in tad")
						//fmt.Println(a,b,hicmap[hicutil.Pair{a,b}],tadmeans[hicutil.Pair{tad[0],tad[1]}],math.Pow(hicmap[hicutil.Pair{a,b}] - tadmeans[hicutil.Pair{tad[0],tad[1]}],2))
						//scorediff := math.Pow(hicmap[hicutil.Pair{a,b}] - tadmeans[hicutil.Pair{tad[0],tad[1]}],2)
						//if scorediff > 0 {fmt.Println("in tad",a,b,scorediff,hicmap[hicutil.Pair{a,b}],tadmeans[hicutil.Pair{tad[0],tad[1]}])}
						tadscore -= math.Pow(hicmap[hicutil.Pair{a,b}] - tadmeans[hicutil.Pair{tad[0],tad[1]}],2)
					} else {
						//fmt.Println("outside of tad",a,b,hicmap[hicutil.Pair{a,b}],nontadmean,math.Pow(hicmap[hicutil.Pair{a,b}] - nontadmean,2))
						//scorediff := math.Pow(hicmap[hicutil.Pair{a,b}] - nontadmean,2)
						//if scorediff > 0 {fmt.Println("out of tad",a,b,scorediff)}
						tadscore -= math.Pow(hicmap[hicutil.Pair{a,b}] - nontadmean,2)
						//tadscore -= math.Pow(hicmap[hicutil.Pair{a,b}] - nontadmeans[hicutil.Pair{tad[0],tad[1]}],2)
					}
				}
			}
			hicsegscores[hicutil.Pair{tad[0],tad[1]}] = tadscore
			score += tadscore
			//fmt.Println(tad,tadscore)
		}
		//os.Exit(1)
	}
	//fmt.Println(score)
	//os.Exit(1)
	return score,hicsegscores,time.Since(starttime)
}*/

/*func computeDIScore(discores []float64, tadset *[][]int) (float64,time.Duration) {

	starttime := time.Now()
	score := 0.0
	for _,tad := range *tadset {
		tadmidpt := (tad[0]+tad[1])/2
		for _,i := range tad {
			if i > tadmidpt {
				score -= discores[i]
			} else {
				score += discores[i]
			}
		}
		//fmt.Println(tad,score)
		//for _,bdy := range tad {
		//	score += discores[bdy]
		//}
	}
	return score,time.Since(starttime)
}*/


/*func calcSumSqErrs(tadscores []float64) (float64) {

	//runtime.GOMAXPROCS(1)

	// first convert tad scores to be within [-1,1]
	maxscore := -1.0
	for _,s := range tadscores {
		if math.Abs(s) > maxscore { maxscore = math.Abs(s)}
	}

	for i,s := range tadscores {
		tadscores[i] = s/maxscore
	}

	sumsqerr := 0.0
	for i:=0; i < len(tadscores)-1; i++ {
		for j:=i; j < len(tadscores); j++ {
			//sumsqerr += math.Pow(math.Abs(tadscores[i]-tadscores[j]),2)
			sumsqerr += math.Abs(tadscores[i]-tadscores[j])
		}
	}
	return sumsqerr
}*/

/*func runMCMC(tadset [][]int, chrlength int, allscores map[hicutil.Pair]float64, meanscores [][]float64, gamma_a float64, allisscores []float64, mindists []int, alpha_is float64, lambda_a float64, lambda_is float64, lambda_tt float64, lambda_hs float64, lambda_di float64, maxsteps int, gamma_tt float64, diridxtt []float64, avgfreq []float64, atilde map[hicutil.Pair]float64, nontadmean float64, tadmeans map[hicutil.Pair]float64, discores []float64, didists []int, hicmap map[hicutil.Pair]float64) ([][]int) {
//map[hicutil.Pair]
	// note this is simulated annealing not MCMC, (sim annealing is an adaptation of metropolis-hastings, which is a monte carlo method so MCMC isn't totally innacurate)

	currtadset := tadset
	stepssinceimprov := 0
	var newloc int
	var newtadset [][]int
	//var idx int
	var tad []int

	if checkForDuplicateTADs(&currtadset) {
			fmt.Println(currtadset)
			fmt.Println("created duplicate in initial TAD set")
			os.Exit(1)
	}

	diwindowsize := 3

	itercount := 0
	temp := 1000000.0

	gdscore_start,gdtime := ftutil.ComputeGraphDensScore(allscores, &meanscores, &tadset, gamma_a)
	isscore_start,istime := ftutil.ComputeBdyInsScore(&allisscores, &mindists, &tadset, alpha_is)
	ttscore_start := 0.0
	//hicsegscore_start := 0.0
	//ttscore_start,tadtreescores,tttime := ftutil.ComputeLinEnrichScore(&tadset, chrlength, gamma_tt, diridxtt, avgfreq, atilde, hicmap, map[hicutil.Pair]float64{})
	hicsegscore_start,hicsegscores,hstime := ftutil.ComputeBlockConstScore(&tadset, chrlength, nontadmean, tadmeans, hicmap, map[hicutil.Pair]float64{}) // start w/ empty map for hicseg scores
	discore_start,ditime := ftutil.ComputeBiasDirScore(discores, didists, &tadset, diwindowsize)

	//fmt.Println("armatus score on random TADs:",gdscore_start)
	//fmt.Println("is score on random TADs:",isscore_start)
	//fmt.Println("tadtree score on random TADs:",ttscore_start)
	//fmt.Println("hicseg score on random TADs:",hicsegscore_start)
	//fmt.Println("DI score on random TADs:",discore_start)
	//os.Exit(1)

	//currscore := (1-lambda1-lambda2-lambda3-lambda4)*gdscore_start + lambda1*isscore_start + lambda2*hicsegscore_start+lambda3*ttscore_start+lambda4*discore_start
	currscore := lambda_a*gdscore_start + lambda_is*isscore_start + lambda_hs*hicsegscore_start+lambda_tt*ttscore_start+lambda_di*discore_start
	fmt.Println("full score on random TADs:",currscore)

	
	currbestscore := currscore
	currbestset := tadset

	//var hicsegscore,ttscore float64
	var hslen time.Duration //,ttlen time.Duration

	// make a for loop w/ stopping condition (create a "steps since improvement" variable?)
	for stepssinceimprov < maxsteps {
		var p int

		p = rand.Intn(2)
		if len(currtadset) == 1 { p = 0 }
		// also need to check if all tads are length 1 (or 0)
		tadstoosmall := 0
		for _,tad = range currtadset {
			if tad[1]-tad[0] < 4 {
				tadstoosmall += 1
			}
		}
		//if checkForDuplicateTADs(&currtadset) {
		//	fmt.Println(currtadset)
		//	fmt.Println("created duplicate TAD in MCMC")
		//	os.Exit(1)
		//}

		if tadstoosmall >= len(currtadset)-2 && len(currtadset) > 1 && tadstoosmall > 0 { p = 1 }
		if p == 0 {
			newtadset = ftutil.SplitRandomTAD(currtadset,chrlength)
		} else {
			newtadset = ftutil.MergeRandomTAD(currtadset)
		}
		if checkForSingleBinTADs(&newtadset) {
			fmt.Println("found single bin TAD")
			fmt.Println(newtadset)
			if p == 0 {
				fmt.Println("created new TAD")
				fmt.Println(newloc,tad)
			} else {
				fmt.Println("merged 2 TADs")
			}
			os.Exit(1)
		}
		//if newtadset[0] == []int{0,-1} {
		//	os.Exit(1)
		//}
		//fmt.Println("newtadset",newtadset)
		// rescore new tad set
		gdscore,gdlen := ftutil.ComputeGraphDensScore(allscores, &meanscores, &newtadset, gamma_a)
		isscore,islen := ftutil.ComputeBdyInsScore(&allisscores, &mindists, &newtadset, alpha_is)
		hicsegscore := 0.0
		ttscore := 0.0
		hicsegscore,hicsegscores,hslen = ftutil.ComputeBlockConstScore(&newtadset, chrlength, nontadmean, tadmeans, hicmap,hicsegscores)
		//ttscore,tadtreescores,ttlen = ftutil.ComputeLinEnrichScore(&newtadset, chrlength, gamma_tt, diridxtt, avgfreq, atilde, hicmap,tadtreescores)
		discore,dilen := ftutil.ComputeBiasDirScore(discores, didists, &newtadset, diwindowsize)
		gdtime += gdlen
		istime += islen
		hstime += hslen
		//tttime += ttlen
		ditime += dilen
		//newscore := (1-lambda1-lambda2-lambda3-lambda4)*gdscore + lambda1*isscore + lambda2*hicsegscore + lambda3*ttscore+lambda4*discore
		newscore := lambda_a*gdscore + lambda_is*isscore + lambda_hs*hicsegscore + lambda_tt*ttscore+lambda_di*discore

		var acceptprob float64
		temp = 0.95*temp
		if newscore > currscore {
			acceptprob = 1.0
		} else {
			acceptprob = math.Exp( -math.Abs(currscore-newscore) / temp )
		}

		// generate random number between [0,1] and accept/reject based on probability
		randnum := rand.Float64()
		//fmt.Println(randnum)
		if randnum <= acceptprob {
			//fmt.Println("accept")
			//fmt.Println(newscore,currscore)
			if newscore > currscore {
				stepssinceimprov = 0
			} else {
				stepssinceimprov += 1 
			}
			currtadset = newtadset
			currscore = newscore
		} else {
			//fmt.Println("reject")
			// keep current tad set, continue on
			stepssinceimprov += 1
		}
		if currscore > currbestscore {
			currbestscore = currscore
			currbestset = currtadset
		}
		//fmt.Println(currscore)
		itercount += 1
		if itercount % 1000 == 0 { 
			fmt.Println("MCMC iteration:",itercount)
			fmt.Println("current best score:", currbestscore)
			//fmt.Println(currbestset)
			//os.Exit(1)
		}
	}
	fmt.Println("number of simulated annealing steps", itercount) //,"final TAD set score =",currscore)
	fmt.Println("full score on final optimized TADs:",currbestscore)
	
	fmt.Println("total time spent in Armatus scoring =", gdtime)
	fmt.Println("total time spent in IS scoring =", istime)
	//fmt.Println("total time spent in TADtree scoring =", tttime)
	fmt.Println("total time spent in HiCSeg scoring =", hstime)
	fmt.Println("total time spent in DI scoring =", ditime)
	
	return currbestset
}*/


/*func paramOptObj(chrlengths map[string]int, centromerelocs map[string][]int, allgdscores map[string]map[hicutil.Pair]float64, meangdscores map[string][][]float64, allisscores map[string][]float64, mindists map[string][]int, alpha float64, gamma float64, lambda1 float64, lambda2 float64, lambda3 float64, lambda4 float64, maxsteps int, ctcfweights map[string][]int, gammaTT float64, diridxtt map[string][]float64, avgfreq map[string][]float64, atilde map[string]map[hicutil.Pair]float64, nontadmean map[string]float64, tadmeans map[string]map[hicutil.Pair]float64, discores map[string][]float64, hicmap map[string]map[hicutil.Pair]float64) float64 {
//map[hicutil.Pair]
	//change inputs to be pointers wherever possible/reasonable

	avgctcf := 0.0
	numchr := 0
	for chrnum := 1; chrnum < 23; chrnum += 2 {
		chrkey := "chr"+strconv.Itoa(chrnum)
		fmt.Println(chrkey)
		newtads := ftutil.CreateRandomTADs(chrlengths[chrkey])

		tads := runMCMC(newtads, chrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], gamma, allisscores[chrkey], mindists[chrkey], alpha, lambda1, lambda2, lambda3, lambda4,  maxsteps, gammaTT, diridxtt[chrkey], avgfreq[chrkey], atilde[chrkey], nontadmean[chrkey], tadmeans[chrkey], discores[chrkey], hicmap[chrkey])
		avgctcf += ftutil.AvgCTCFload(tads, ctcfweights[chrkey], centromerelocs[chrkey])
		numchr +=1
	}

	//compute avg ctcf load
	return avgctcf/float64(numchr)

}*/

/*func paramOptObj(chrlengths map[string]int, centromerelocs map[string][]int, allgdscores map[string]map[hicutil.Pair]float64, meangdscores map[string][][]float64, allisscores map[string][]float64, mindists map[string][]int, params map[bayesopt.Param]float64, alpha bayesopt.UniformParam, gamma bayesopt.UniformParam, lambda bayesopt.UniformParam, maxsteps int, ctcfweights map[string][]int) float64 {

	//change inputs to be pointers wherever possible/reasonable

	avgctcf := 0.0
	numchr := 0
	for chrnum := 1; chrnum < 23; chrnum += 2 {
		chrkey := "chr"+strconv.Itoa(chrnum)
		fmt.Println(chrkey)
		newtads := ftutil.CreateRandomTADs(chrlengths[chrkey], centromerelocs[chrkey])

		tads := runMCMC(newtads, chrlengths[chrkey], centromerelocs[chrkey], allgdscores[chrkey], meangdscores[chrkey], params[gamma], allisscores[chrkey], mindists[chrkey], params[alpha], params[lambda], maxsteps)
		avgctcf += ftutil.AvgCTCFload(tads, ctcfweights[chrkey])
		numchr +=1
	}

	//compute avg ctcf load
	return avgctcf/float64(numchr)

}*/


/*func removeZeros(x []int) []int {

	var xclean []int
	for _,val := range x {
		if val != 0 {
			xclean = append(xclean,val)
		}
	}
	return xclean
}*/

/*func checkTADinCentromere(tadset [][]int, centromereloc []int) bool {

	for _,tad := range tadset {
		if (tad[0] > centromereloc[0] && tad[0] < centromereloc[1]) || (tad[1] < centromereloc[1] && tad[1] > centromereloc[0]) {
			return true
		}
	}
	return false

}*/

/*func checkForDuplicateTADs(tadset *[][]int) bool {

	prevtad := []int{-1,-1}
	for _,tad := range *tadset {
		if tad[0] == prevtad[0] && tad[1] == prevtad[1] {
			return true
		}
		prevtad = tad
	}
	return false
}*/

/*func checkForSingleBinTADs(tadset *[][]int) bool {

	for _,tad := range *tadset {
		if tad[0] == tad[1] {
			return true
		}
	}
	return false
}*/


/*func writeSliceToFile(x []float64, filename string) {

	f,err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	w := bufio.NewWriter(f)
	for _,y := range x {
		fmt.Fprintf(w,fmt.Sprintf("%f",y)+"\n")
	}
	w.Flush()
	f.Close()
	fmt.Println("wrote slices values to",filename)

}*/




/*// need to come up with stopping criteria to put this in for/while loop
				//newtad set needs to be in a dictionary like chrlengths
				newtadset[chrkey] = runMCMC(newtadset[chrkey], chrlengths[chrkey], centromerelocs[chrkey], allgdscores[chrkey], meangdscores[chrkey], optgamma, allisscores[chrkey], mindists[chrkey], optalpha, lambda, maxsteps)
				//fmt.Println(newtadset)
				//fmt.Println("done computing new tad set")
				//avgctcf := ftutil.AvgCTCFload(newtadset, ctcfweights[chrkey])
v				avgctcfbdy = append(avgctcfbdy, ftutil.AvgCTCFloadByBdy(newtadset[chrkey],ctcfweights[chrkey])...)
				//fmt.Println("avgctcf for",chrkey,avgctcfbdy)
				//os.Exit(1)

				// store scores per gamma/alpha value w/ ctcf value
				for gamma := gmin; gamma < gmax+gstep; gamma += gstep {
					//score := computeParticularGDScore(allgdscores[chrkey], meangdscores[chrkey], newtadset, gamma)
					scores := computeParticularGDScoreByBdy(allgdscores[chrkey], meangdscores[chrkey], newtadset[chrkey], gamma)
					//rval,_ := onlinestats.Spearman(avgctcfbdy,scores)
					if vals,ok := gammamaps[gamma]; ok {
						gammamaps[gamma] = append(vals,scores...)
						//gammamaps[gamma] = append(vals,rval)
						//gammamaps[gamma][0] = append(val[0],avgctcf)
						//gammamaps[gamma][1] = append(val[1],score)
					} else {
						gammamaps[gamma] = scores
						//gammamaps[gamma] = make([]float64,1)
						//gammamaps[gamma] = []float64{rval}
						//gammamaps[gamma] = make([][]float64,2)
						//gammamaps[gamma][0] = []float64{avgctcf}
						//gammamaps[gamma][1] = []float64{score}
					}
				}
				//fmt.Println(gammamaps)
				
				for alpha := amin; alpha < amax+astep; alpha += astep {
					//score := computeParticularInsScore(allisscores[chrkey], mindists[chrkey], newtadset, alpha)
					scores := computeParticularInsScoreByBdy(allisscores[chrkey], mindists[chrkey], newtadset[chrkey], alpha)
					//rval,_ := onlinestats.Spearman(avgctcfbdy,scores)
					if vals,ok := alphamaps[alpha]; ok {
v						alphamaps[alpha] = append(vals,scores...)
						//alphamaps[alpha] = append(vals,rval)
						//alphamaps[alpha][0] = append(val[0],avgctcf)
						//alphamaps[alpha][1] = append(val[1],score)
					} else {
						alphamaps[alpha] = scores
						//alphamaps[alpha] = []float64{rval}
						//alphamaps[alpha] = make([][]float64,2)
						//alphamaps[alpha][0] = []float64{avgctcf}
						//alphamaps[alpha][1] = []float64{score}
					}
					if len(alphamaps[alpha]) == 0 {
						fmt.Println(alpha,scores)
						os.Exit(1)
					}
				}
				//fmt.Println(alphamaps)
			}
			//fmt.Println(gammamaps)
			//fmt.Println(alphamaps)
			//os.Exit(1)

			//write gamma maps and ctcf boundaries to files to plot
			//writeSliceToFile(avgctcfbdy,"outfiles/ctcfbybdy.txt")

			optrval := -100.0
			for gamma := gmin; gamma < gmax+gstep; gamma += gstep {
				//writeSliceToFile(gammamaps[gamma],"outfiles/armatusscores_gamma"+fmt.Sprintf("%f",gamma)+".txt")
				rval,_ := onlinestats.Spearman(avgctcfbdy,gammamaps[gamma])
				/*avgrval := 0.0
				numrvals := len(gammamaps[gamma])
				for _,chrrval := range gammamaps[gamma] {
					//fmt.Println(chrrval)
					avgrval += float64(chrrval)
				}
				avgrval = avgrval / float64(numrvals)*/
				/*fmt.Println(gamma, rval)
				//data := gammamaps[gamma]
				//rval,_ := onlinestats.Spearman(data[0],data[1])
				if rval > optrval {
					optgamma = gamma
					optrval = rval
				}
			}
			fmt.Println("optgamma =",optgamma,"optrval =",optrval)
			//fmt.Println(optgamma,optrval)
			//os.Exit(1)
			//alpha := amin
			optrval = -100.0
			//prevrval := -100.0
			//numsamerval := 0
			//for numsamerval < 10  {
			for alpha := amin; alpha < amax+astep; alpha += astep {
				rval,_ := onlinestats.Spearman(avgctcfbdy,alphamaps[alpha])
				if rval > 1 {
					fmt.Println(avgctcfbdy)
					fmt.Println(alphamaps[alpha])
					fmt.Println(len(avgctcfbdy),len(alphamaps[alpha]))
					os.Exit(1)
				}
				/*avgrval := 0.0
				numrvals := len(alphamaps[alpha])
				for _,chrrval := range alphamaps[alpha] {
					//fmt.Println("chrrval",chrrval)
					avgrval += float64(chrrval)
				}
				//fmt.Println(numrvals)
				avgrval = avgrval / float64(numrvals)*/
				//fmt.Println("avgrval",avgrval)
				//os.Exit(1)
				//rval := alphamaps[alpha][0]
				//data := alphamaps[alpha]
				//rval,_ := onlinestats.Spearman(data[0],data[1])
				/*if rval > optrval {
					optalpha = alpha
					optrval = rval
				}
				//alpha += astep
				//if prevrval == avgrval { numsamerval += 1 }
				//prevrval = avgrval
			}
			fmt.Println("optalpha =",optalpha,"optrval =",optrval)
			if prevoptalpha == optalpha && prevoptgamma == optgamma {
				sameparams += 1
				fmt.Println("optimal gammas and alphas unchanged for 2 iterations")
			} else {
				fmt.Println("previous opt alpha:", prevoptalpha, "current opt alpha:", optalpha)
				fmt.Println("previous opt gamma:", prevoptgamma, "current opt gamma:", optgamma)
				sameparams = 0
			}
			if sameparams == 2 {
				keeptraining = false
			}
			//keeptraining = false
			prevoptalpha = optalpha
			prevoptgamma = optgamma
			//fmt.Println("max alpha checked =",alpha)
			iter++
*/


/*optgamma = 0.0
	optalpha = 0.0
	optlambda = 0.0
	maxavgctcf := 0.0
	var avgctcf float64
	for gamma := gmin; gamma < gmax+gstep; gamma += gstep {
		avgctcf := 0.0
		for _,ctcfval := range gammamaps[gamma] {
			avgctcf += ctcfval
		}
		avgctcf = avgctcf / float64(len(gammamaps[gamma]))
		fmt.Println("gamma =",gamma,"avg ctcf =",avgctcf)
		if avgctcf > maxavgctcf {
			maxavgctcf = avgctcf
			optgamma = gamma
		}
	}
	maxavgctcf = 0.0
	for alpha := amin; alpha < amax+astep; alpha += astep {
		avgctcf = 0.0
		for _,ctcfval := range alphamaps[alpha] {
			avgctcf += ctcfval
		}
		avgctcf = avgctcf / float64(len(alphamaps[alpha]))
		fmt.Println("alpha =",alpha,"avg ctcf =",avgctcf)
		if avgctcf > maxavgctcf {
			maxavgctcf = avgctcf
			optalpha = alpha
		}
	}
	
	//optparam = param value w/ highest average CTCF load across chr
	//fmt.Println(newtadset)
	fmt.Println("optalplha =", optalpha, "optgamma =",optgamma)
	// train lambda parameter
	maxavgctcf = 0.0
	for lambda := lmin; lambda < lmax+lstep; lambda += lstep {
		avgctcf = 0.0
		for _,ctcfval := range lambdamaps[lambda] {
			avgctcf += ctcfval
		}
		avgctcf = avgctcf / float64(len(lambdamaps[lambda]))
		fmt.Println("lambda =",lambda,"avg ctcf =",avgctcf)
		if avgctcf > maxavgctcf {
			maxavgctcf = avgctcf
			optlambda = lambda
		}
	}

	fmt.Println("optalplha =", optalpha, "optgamma =",optgamma, "optlambda =", optlambda)
	os.Exit(1)*/




//	tadset := [][]int{[]int{0,28},[]int{29,33},[]int{34,83},[]int{84,111},[]int{112,151},[]int{152,186},[]int{187,199},[]int{200,245},[]int{246,273},[]int{274,310},[]int{311,336},[]int{337,355},[]int{356,376},[]int{377,420},[]int{421,439},[]int{440,443},[]int{444,477},[]int{478,513},[]int{514,529},[]int{530,566},[]int{567,609},[]int{610,652},[]int{653,695},[]int{696,711},[]int{712,726},[]int{727,739},[]int{740,768},[]int{769,801},[]int{802,816},[]int{817,845},[]int{846,893},[]int{894,901},[]int{902,913},[]int{914,917},[]int{918,921},[]int{922,950},[]int{951,979},[]int{980,1027}}

//	tadset := [][]int{[]int{0,17},[]int{18,64},[]int{65,67},[]int{68,98},[]int{99,139},[]int{140,179},[]int{180,229},[]int{230,245},[]int{246,270},[]int{271,317},[]int{318,364},[]int{365,372},[]int{373,385},[]int{386,416},[]int{417,424},[]int{425,438},[]int{439,486},[]int{487,522},[]int{523,553},[]int{554,580},[]int{581,626},[]int{627,632},[]int{633,649},[]int{650,684},[]int{685,711},[]int{712,719},[]int{720,743},[]int{744,751},[]int{752,754},[]int{755,763},[]int{764,767},[]int{768,814},[]int{815,822},[]int{823,845},[]int{846,879},[]int{880,910},[]int{911,947},[]int{948,959},[]int{960,990},[]int{991,998}}

 //	tadset := [][]int{[]int{7,13},[]int{14,17},[]int{18,20},[]int{21,24},[]int{25,39},[]int{40,62},[]int{63,67},[]int{68,80},[]int{81,83},[]int{84,88},[]int{89,96},[]int{97,103},[]int{104,110},[]int{111,119},[]int{120,123},[]int{124,126},[]int{127,130},[]int{131,133},[]int{127,154},[]int{155,158},[]int{159,164},[]int{165,167},[]int{168,173},[]int{174,179},[]int{180,191},[]int{192,196},[]int{197,200},[]int{201,208},[]int{209,216},[]int{217,219},[]int{220,224},[]int{225,230},[]int{231,232},[]int{233,239},[]int{240,243},[]int{244,251},[]int{252,255},[]int{256,264},[]int{265,267},[]int{268,272},[]int{273,281},[]int{282,283},[]int{284,289},[]int{290,296},[]int{297,299},[]int{300,311},[]int{312,314},[]int{315,320},[]int{321,325},[]int{326,330},[]int{331,336},[]int{337,341},[]int{342,352},[]int{353,364},[]int{365,370},[]int{371,380},[]int{381,385},[]int{386,393},[]int{394,402},[]int{403,405},[]int{406,433},[]int{434,438},[]int{430,446},[]int{447,452},[]int{453,461},[]int{462,468},[]int{469,470},[]int{471,492},[]int{493,508},[]int{509,518},[]int{519,526},[]int{527,536},[]int{537,543},[]int{544,554},[]int{555,572},[]int{573,592},[]int{593,603,},[]int{604,628},[]int{629,632},[]int{633,639},[]int{640,647},[]int{648,653},[]int{654,671},[]int{672,680},[]int{681,690},[]int{691,712},[]int{713,760},[]int{761,777},[]int{778,783},[]int{784,791},[]int{792,844},[]int{845,854},[]int{855,868},[]int{869,893},[]int{894,918},[]int{919,931},[]int{932,935},[]int{936,943},[]int{944,950},[]int{951,955},[]int{956,1004},[]int{1005,1013},[]int{1014,1017},[]int{1018,1037},[]int{1040,1042},[]int{1043,1076},[]int{1077,1091},[]int{1092,1097},[]int{1098,1100},[]int{1101,1107},[]int{1108,1111},[]int{1112,1129},[]int{1130,1134},[]int{1135,1143},[]int{1144,1152},[]int{1153,1163},[]int{1164,1172},[]int{1173,1175},[]int{1176,1184},[]int{1185,1200},[]int{1201,1207},[]int{1208,1209},[]int{1212,1214},[]int{1430,1432},[]int{1439,1441},[]int{1445,1447},[]int{1448,1454},[]int{1455,1459},[]int{1464,1474},[]int{1478,1481},[]int{1484,1487},[]int{1488,1492},[]int{1493,1495},[]int{1496,1497},[]int{1498,1501},[]int{1502,1509},[]int{1510,1512},[]int{1513,1519},[]int{1520,1533},[]int{1534,1536},[]int{1537,1540},[]int{1541,1544},[]int{1545,1551},[]int{1552,1555},[]int{1556,1561},[]int{1562,1566},[]int{1567,1569},[]int{1570,1610},[]int{1611,1613},[]int{1614,1618},[]int{1619,1623},[]int{1624,1628},[]int{1629,1633},[]int{1634,1651},[]int{1652,1676},[]int{1677,1680},[]int{1681,1696},[]int{1697,1715},[]int{1716,1751},[]int{1752,1767},[]int{1768,1779},[]int{1780,1792},[]int{1793,1796},[]int{1797,1804},[]int{1805,1810},[]int{1811,1828},[]int{1829,1835},[]int{1836,1846},[]int{1847,1862},[]int{1863,1929},[]int{1930,1966},[]int{1967,1975},[]int{1976,1985},[]int{1986,2005},[]int{2006,2008},[]int{2009,2012},[]int{2013,2018},[]int{2019,2022},[]int{2023,2029},[]int{2030,2037},[]int{2038,2043},[]int{2044,2048},[]int{2049,2054},[]int{2055,2065},[]int{2066,2069},[]int{2070,2079},[]int{2080,2100},[]int{2101,2116},[]int{2117,2127},[]int{2128,2132},[]int{2133,2144},[]int{2145,2159},[]int{2160,2184},[]int{2185,2200},[]int{2201,2203},[]int{2204,2208},[]int{2209,2226},[]int{2227,2239},[]int{2240,2246},[]int{2247,2251},[]int{2252,2255},[]int{2256,2263},[]int{2264,2266},[]int{2267,2277},[]int{2278,2282},[]int{2283,2285},[]int{2286,2297},[]int{2298,2304},[]int{2305,2314},[]int{2315,2334},[]int{2335,2344},[]int{2345,2351},[]int{2352,2363},[]int{2364,2378},[]int{2379,2401},[]int{2402,2409},[]int{2410,2432},[]int{2433,2448},[]int{2449,2470},[]int{2471,2489},[]int{2490,2491}}



/*		meangdchr := meangdscores[chrkey]
		isscoreschr := allisscores[chrkey]
		mindistchr := mindists[chrkey]
	gdscore,_ := computeParticularGDScore(allgdscores[chrkey], &meangdchr, &tadset, optgamma)
	isscore,_ := computeParticularInsScore(&isscoreschr, &mindistchr, &tadset, optalpha)
	//ttscore,_ := computeTTScore(&tadset, chrlengths[chrkey], optgammaTT, diridxtt[chrkey], avgfreq[chrkey], exptadfreq[chrkey], allhicmaps[chrkey])
	hicsegscore,_,_ := computeHiCSegScore(&tadset, chrlengths[chrkey], nontadmean[chrkey], tadmeans[chrkey], allhicmaps[chrkey],map[hicutil.Pair]float64{})
		discore,_ := computeDIScore(discores[chrkey],&tadset)

	fmt.Println("armatus score on true TADs:",gdscore)
	fmt.Println("is score on true TADs:",isscore)
	//fmt.Println("tadtree score on true TADs:",ttscore)
	fmt.Println("hicseg score on true TADs:",hicsegscore)
		fmt.Println("DI score on true TADs:",discore)
	//os.Exit(1)
*/
