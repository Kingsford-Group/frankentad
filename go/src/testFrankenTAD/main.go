package main

import(
	"ftutil"
	"hicutil"
	"flag"
	"fmt"
	"os"
)

// this script should take in a file with TADs and return the score in each of the component subfunctions
// inputs: tadset, corresponding hic data, res

// output: print score from each function


func main() {

	tadfile := flag.String("tads","","File with TAD set")
        hicfile := flag.String("hic","","Hi-C file(s)")
        res := flag.Int("res",1,"resolution of Hi-C data")

	//internal params for Armatus, IS, TADTree
	alpha_is := flag.Float64("alpha_is",1.0,"Boundary insulation internal parameter alpha")
	gamma_a := flag.Float64("gamma_a",1.0,"Graph density internal parameter gamma")
	//gamma_tt := flag.Float64("gamma_tt",1.0,"Linear enrichment internal parameter gamma")

	flag.Parse()

	windowsize := 5
        p := 3
        q := 12
	diwindowsize := 3

	// read in files
	orighicmap,origchrlength,_ := hicutil.ReadHiCFile(*hicfile,*res, []float64{})
	//orighicmap,origchrlength,centromereloc := hicutil.ReadHiCFile(*hicfile,*res, []float64{})
	//newhicmap,newchrlength := hicutil.RemoveCentromere(orighicmap,centromereloc,origchrlength)
	tadset := hicutil.ReadTADFile(*tadfile,*res)

	// precompute needed values
	allgdscores, meangdscores, allisscores, mindists, _, _, _, tadmeans, nontadmean, discores, didists := ftutil.PrecomputeAll_nott(orighicmap, origchrlength, windowsize, p, q, *res)
	//allgdscores, meangdscores, allisscores, mindists, diridxtt, avgfreq, exptadfreq, tadmeans, nontadmean, discores := ftutil.PrecomputeAll(newhicmap, newchrlength, windowsize, p, q, *res)

	// score TAD set on each function and print results
	gdscore,gdtime := ftutil.ComputeGraphDensScore(allgdscores, &meangdscores, &tadset, *gamma_a)
	fmt.Println("graph density score =",gdscore)
	fmt.Println("graph density computation took",gdtime,"\n")

        biscore,bitime := ftutil.ComputeBdyInsScore(&allisscores, &mindists, &tadset, *alpha_is)
	fmt.Println("boundary insulation score =",biscore)
	fmt.Println("boundary insulation computation took",bitime,"\n")

	//lescore,_,letime := ftutil.ComputeLinEnrichScore(&tadset, newchrlength, *gamma_tt, diridxtt, avgfreq, exptadfreq, newhicmap, map[hicutil.Pair]float64{})
	//fmt.Println("linear enrichment score =",lescore)
	//fmt.Println("linear enrichment computation took",letime,"\n")

        //bcscore,_,bctime := ftutil.ComputeBlockConstScore(&tadset, origchrlength, nontadmean, tadmeans, orighicmap, map[hicutil.Pair]float64{})
	//fmt.Println("block constant score =",bcscore)
	//fmt.Println("block constant computation took",bctime,"\n")

	bdscore,bdtime := ftutil.ComputeBiasDirScore(discores, didists, &tadset, diwindowsize)
	fmt.Println("bias direction score =",bdscore)
	fmt.Println("bias direction computation took",bdtime,"\n")

	os.Exit(1)
	fmt.Println(tadmeans,nontadmean)
	/*gdscorebytad := 0.0 
	biscorebytad := 0.0 
	bcscorebytad := 0.0
	bdscorebytad := 0.0
	for _,tad := range tadset {
		tadscore,_ := ftutil.ComputeGraphDensScorebyTAD(allgdscores, &meangdscores, tad, *gamma_a)
		gdscorebytad += tadscore
		tadscore,_ = ftutil.ComputeBdyInsScorebyTAD(&allisscores, &mindists, tad, *alpha_is)
		biscorebytad += tadscore
		tadscore,_,_ = ftutil.ComputeBlockConstScorebyTAD(tad, origchrlength, nontadmean, tadmeans, orighicmap, map[hicutil.Pair]float64{})
		bcscorebytad += tadscore
		tadscore,_ = ftutil.ComputeBiasDirScorebyTAD(discores, didists, tad, diwindowsize)
		bdscorebytad += tadscore
	}

	fmt.Println("graph density full: ",gdscore,", by TAD: ",gdscorebytad)
	fmt.Println("boundary insulation full: ",biscore,", by TAD: ",biscorebytad)
	fmt.Println("block constant full: ",bcscore,", by TAD: ",bcscorebytad)
	fmt.Println("bias direction full: ",bdscore,", by TAD: ",bdscorebytad)*/
}
