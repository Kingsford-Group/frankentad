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
	"github.com/c-bata/goptuna"
	"github.com/c-bata/goptuna/tpe"
//	"context"
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
	//rep1file := flag.String("rep1","","Hi-C file(s) of replicate 1")
	//rep2file := flag.String("rep2","","Hi-C file(s) of replicate 2")
	res := flag.Int("res",1,"resolution of Hi-C data")

	trainflag := flag.Bool("train",false,"training mode to learn parameters")
	alpha_is_in := flag.Float64("alpha_is",-1.0,"TopDom/InsScore internal parameter")
	gamma_armatus_in := flag.Float64("gamma_armatus",-1.0,"Armatus internal parameter")
	lambda_armatus_in := flag.Float64("lambda_arm",-1.0,"weight for Armatus function")
	lambda_is_in := flag.Float64("lambda_is",-1.0,"weight for TopDom/InsScore function")
	lambda_di_in := flag.Float64("lambda_di",-1.0,"weight for DI function")

	numtrainsteps := flag.Int("bosteps",100,"number of Bayesian optimization steps")

	ctcffile := flag.String("ctcf","inputdata/CTCFBSDB_all_exp_sites_Sept12_2012_hg19.txt","file with locations of CTCF binding sites")
	//bedfile := flag.String("bed","","BED file of ChIPseq data to optimize parameters for")
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
	
	//var maxsteps int

	/*var orighicmapsrep1 = map[string]map[hicutil.Pair]float64{}
	//var origchrlengths = map[string]int{} 
	var newhicmapsrep1 = map[string]map[hicutil.Pair]float64{}
	//var newchrlengths = map[string]int{}
	//var centromerelocs = map[string][]int{}
	var allgdscoresrep1 = map[string]map[hicutil.Pair]float64{}
	var meangdscoresrep1 = map[string][][]float64{}
	var allisscoresrep1 = map[string][]float64{}
	var mindistsrep1 = map[string][]int{}
	var chrrep1file string
	var newtadsrep1 = map[string][][]int{}

	var orighicmapsrep2 = map[string]map[hicutil.Pair]float64{}
	//var origchrlengths = map[string]int{} 
	var newhicmapsrep2 = map[string]map[hicutil.Pair]float64{}
	//var newchrlengths = map[string]int{}
	//var centromerelocs = map[string][]int{}
	var allgdscoresrep2 = map[string]map[hicutil.Pair]float64{}
	var meangdscoresrep2 = map[string][][]float64{}
	var allisscoresrep2 = map[string][]float64{}
	var mindistsrep2 = map[string][]int{}
	var chrrep2file string
	var newtadsrep2 = map[string][][]int{}*/


	allhicfiles,err := filepath.Glob(*hicfile)
	//allrep1files,err := filepath.Glob(*rep1file)
	//allrep2files,err := filepath.Glob(*rep2file)
	if err != nil{
		fmt.Println("couldnt find files matching filepath")
		log.Fatal(err)
	}
	
	ctcflocs := ftutil.ReadCTCFlocs(*ctcffile, *res)
	//chiplocs := ftutil.ReadBEDfile(*bedfile,*res)
	ctcfweights := ftutil.CountCTCFweights(ctcflocs)
	//chipweights := ftutil.CountCTCFweights(chiplocs)
	//fmt.Println(ctcflocs)
	
	windowsize := 5
	p := 3
	q := 12

	var avgfreq = map[string][]float64{}
	var exptadfreq = map[string]map[hicutil.Pair]float64{}
	
	var tadmeans = map[string]map[hicutil.Pair]float64{}
	var nontadmean = map[string]float64{}

	var discores = map[string][]float64{}
	var didists = map[string][]int{}


	/*var avgfreqrep1 = map[string][]float64{}
	var exptadfreqrep1 = map[string]map[hicutil.Pair]float64{}
	
	var tadmeansrep1 = map[string]map[hicutil.Pair]float64{}
	var nontadmeanrep1 = map[string]float64{}

	var discoresrep1 = map[string][]float64{}
	var didistsrep1 = map[string][]int{}

	var avgfreqrep2 = map[string][]float64{}
	var exptadfreqrep2 = map[string]map[hicutil.Pair]float64{}
	
	var tadmeansrep2 = map[string]map[hicutil.Pair]float64{}
	var nontadmeanrep2 = map[string]float64{}

	var discoresrep2 = map[string][]float64{}
	var didistsrep2 = map[string][]int{}*/


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

		/*for _,f := range allrep1files {
			if strings.Contains(f,chrkey+".") {
				chrrep1file = f
				break
			}
		}
		orighicmapsrep1[chrkey],origchrlengths[chrkey],centromerelocs[chrkey] = hicutil.ReadHiCFile(chrrep1file,*res, []float64{})
		newhicmapsrep1[chrkey],newchrlengths[chrkey] = hicutil.RemoveCentromere(orighicmapsrep1[chrkey],centromerelocs[chrkey],origchrlengths[chrkey])
		for _,f := range allrep2files {
			if strings.Contains(f,chrkey+".") {
				chrrep2file = f
				break
			}
		}
		orighicmapsrep2[chrkey],_,_ = hicutil.ReadHiCFile(chrrep2file,*res, []float64{})
		newhicmapsrep2[chrkey],_ = hicutil.RemoveCentromere(orighicmapsrep2[chrkey],centromerelocs[chrkey],origchrlengths[chrkey])*/


		allgdscores[chrkey], meangdscores[chrkey], allisscores[chrkey], mindists[chrkey], _, avgfreq[chrkey], exptadfreq[chrkey], tadmeans[chrkey], nontadmean[chrkey], discores[chrkey], didists[chrkey] = ftutil.PrecomputeAll_nott(newhicmaps[chrkey], newchrlengths[chrkey], windowsize, p, q, *res)


		/*allgdscoresrep1[chrkey], meangdscoresrep1[chrkey], allisscoresrep1[chrkey], mindistsrep1[chrkey], _, avgfreqrep1[chrkey], exptadfreqrep1[chrkey], tadmeansrep1[chrkey], nontadmeanrep1[chrkey], discoresrep1[chrkey], didistsrep1[chrkey] = ftutil.PrecomputeAll_nott(newhicmapsrep1[chrkey], newchrlengths[chrkey], windowsize, p, q, *res)

		allgdscoresrep2[chrkey], meangdscoresrep2[chrkey], allisscoresrep2[chrkey], mindistsrep2[chrkey], _, avgfreqrep2[chrkey], exptadfreqrep2[chrkey], tadmeansrep2[chrkey], nontadmeanrep2[chrkey], discoresrep2[chrkey], didistsrep2[chrkey] = ftutil.PrecomputeAll_nott(newhicmapsrep2[chrkey], newchrlengths[chrkey], windowsize, p, q, *res)*/

	}
	fmt.Println("preprocessing took", time.Since(starttime))

	var optgamma_a, optalpha_is float64 //optgamma_tt float64
	var optlambda_a, optlambda_is, optlambda_di float64 //optlambda_tt, optlambda_hs, optlambda_di float64
	var gamma_a, alpha_is float64 //gamma_tt float64
	var lambda_a, lambda_is, lambda_di float64 //lambda_tt, lambda_hs, lambda_di float64


	gamma_a = *gamma_armatus_in

	if *trainflag {
		fmt.Println("training parameters")
		//trialchan := make(chan goptuna.FrozenTrial, 8)
		//relativeSampler := bayesopt.NewSampler()
		study, err := goptuna.CreateStudy(
			"goptuna-testing",
			//goptuna.StudyOptionRelativeSampler(relativeSampler),
			goptuna.StudyOptionSampler(tpe.NewSampler()),
			//goptuna.StudyOptionIgnoreObjectiveErr(true),
			//goptuna.StudyOptionSetTrialNotifyChannel(trialchan),
		)
		if err != nil { fmt.Println("error initializing study") }

		//define objective function
		objective := func(trial goptuna.Trial) (float64,error) {

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
			if *lambda_di_in < 0 {
				lambda_di,_ = trial.SuggestUniform("lambda_di",0,100)
			} else { 
				lambda_di = *lambda_di_in
				optlambda_di = *lambda_di_in
			}
		
			avgctcf := 0.0
			//avgdiff := 0.0
			//avgchip := 0.0
			//avgji := 0.0
			numchr := 0
			for chrnum := 1; chrnum < 23; chrnum += 2 { //should be able to make this a go routine to speed things up
				chrkey := "chr"+strconv.Itoa(chrnum)
				fmt.Println(chrkey)
			

				newtads[chrkey] = ftutil.RunDynamicProgram(newhicmaps[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], gamma_a, allisscores[chrkey], mindists[chrkey], alpha_is, lambda_a, lambda_is, lambda_di, discores[chrkey], didists[chrkey])

				/*newtadsrep1[chrkey] = ftutil.RunDynamicProgram(newhicmapsrep1[chrkey], newchrlengths[chrkey], allgdscoresrep1[chrkey], meangdscoresrep1[chrkey], gamma_a, allisscoresrep1[chrkey], mindistsrep1[chrkey], alpha_is, lambda_a, lambda_is, lambda_di, discoresrep1[chrkey], didistsrep1[chrkey])
				
				newtadsrep2[chrkey] = ftutil.RunDynamicProgram(newhicmapsrep2[chrkey], newchrlengths[chrkey], allgdscoresrep2[chrkey], meangdscoresrep2[chrkey], gamma_a, allisscoresrep2[chrkey], mindistsrep2[chrkey], alpha_is, lambda_a, lambda_is, lambda_di, discoresrep2[chrkey], didistsrep2[chrkey])*/

				/*if len(newtads[chrkey]) == 0 {
					fmt.Println("DP returned no TADs")
					os.Exit(1)
				}*/

				/*tadswithcmrep1 := hicutil.PostprocessTADs(newtadsrep1[chrkey], allgdscoresrep1[chrkey], centromerelocs[chrkey])
				
				tadswithcmrep2 := hicutil.PostprocessTADs(newtadsrep2[chrkey], allgdscoresrep2[chrkey], centromerelocs[chrkey])*/

				/*if len(tadswithcm) == 0 {
					fmt.Println("postprocessing removed all TADs")
					os.Exit(1)
				}*/
				
				tadswithcm := hicutil.PostprocessTADs(newtads[chrkey], allgdscores[chrkey], centromerelocs[chrkey])

				//avgji += hicutil.ComputeJaccardIndex(tadswithcmrep1,tadswithcmrep2)
				//avgdiff += ftutil.InterIntraDiff(newtads[chrkey], newhicmaps[chrkey])
				//avgchip += ftutil.AvgCTCFload(tadswithcm,chipweights[chrkey])
				avgctcf += ftutil.AvgCTCFload(tadswithcm, ctcfweights[chrkey])
				//avgctcf += ftutil.AvgCTCFload(newtads[chrkey],ctcfweights[chrkey],centromerelocs[chrkey])
				numchr +=1
			}

			// print interation number, param values, tad set, ctcf value to file for debugging
			//fmt.Println("avg inter/intra freq diff:",avgdiff)
			//compute avg ctcf load
			return -avgctcf/float64(numchr), nil
			//return -avgchip/float64(numchr), nil
			//return -avgji/float64(numchr), nil
			//return avgdiff, nil ///float64(numchr),nil
		}

		// Run an objective function n times to find a global minimum. - coment these lines to write BO details to file
		err = study.Optimize(objective, *numtrainsteps)
		if err != nil { fmt.Println("error optimizing") }
		

		// Print the best evaluation value and the parameters.
		v, _ := study.GetBestValue()
		optparams, _ := study.GetBestParams()
		fmt.Println(v,optparams)
	
		if *gamma_armatus_in < 0 { optgamma_a = optparams["gamma_a"].(float64) }
		if *alpha_is_in < 0 { optalpha_is = optparams["alpha_is"].(float64) }
		if *lambda_is_in < 0 { optlambda_is = optparams["lambda_is"].(float64) }
		if *lambda_di_in < 0 { optlambda_di = optparams["lambda_di"].(float64) }

		
		log.Printf("Best evaluation value=%f (alpha_is=%f, gamma_a=%f, lambda_a=%f, lambda_is=%f, lambda_di=%f)", v, optalpha_is, optgamma_a, optlambda_a, optlambda_is, optlambda_di )
	} else {
		fmt.Println("using given parameters")

		if *gamma_armatus_in > -1 {
			optgamma_a = *gamma_armatus_in // armatus param
		} else { optgamma_a = 0.53 }  
		if *alpha_is_in > -1 {
			optalpha_is = *alpha_is_in // IS score internal param
		} else { optalpha_is = 56.54 }
		if *lambda_armatus_in > -1 {
			optlambda_a = *lambda_armatus_in // armatus score weight
		} else { optlambda_a = 50.0 }
		if *lambda_is_in > -1 {
			optlambda_is = *lambda_is_in // IS score weight
		} else { optlambda_is = 6.23 } 
		if *lambda_di_in > -1 {
			optlambda_di = *lambda_di_in // DI score weight
		} else { optlambda_di = 6.0 } 

	}

	//compute TADs with opt params and write them to file
	for chrnum := 1; chrnum < 23; chrnum += 1 {
		chrkey := "chr"+strconv.Itoa(chrnum)
		filename := *outfilename+chrkey+".txt"

		fmt.Println(chrkey) //, filename)


		opttads := ftutil.RunDynamicProgram(newhicmaps[chrkey], newchrlengths[chrkey], allgdscores[chrkey], meangdscores[chrkey], optgamma_a, allisscores[chrkey], mindists[chrkey], optalpha_is, optlambda_a, optlambda_is, optlambda_di, discores[chrkey], didists[chrkey])

		//fmt.Println(opttads)
		//ftutil.WriteScoresToFile(alltadscores,*outfilename+chrkey+"_tadscores.txt")
		//os.Exit(1)

		cleantads := hicutil.PostprocessTADs(opttads, allgdscores[chrkey], centromerelocs[chrkey])

		fmt.Println(cleantads)

		ftutil.WriteTADsToFile(cleantads,*res,optalpha_is,optgamma_a,optlambda_a,optlambda_is,optlambda_di,filename)
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
