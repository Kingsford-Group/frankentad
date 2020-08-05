package ftutil

import(
	"time"
//	"fmt"
//	"os"
)

func ComputeBiasDirScore(discores []float64, didists []int, tadset *[][]int, windowsize int) (float64,time.Duration) {

// normalize drop by nearby values? mean-center by mean within a window? nope, still doesn't work. need to somehow penalize for bad boundaries
// only add to score if goes from negative to positive DI at boundary? need a way to incorporate sign change.. avg of next n bins (should be positive) - avg of prev n bins (should be negative)? includes baby tads around true boundaries
// only add to score if this difference is the maximum diff within the window? others need to be < 0 or can still be included with no penalty or include actual difference at boundary?

// ideally: scores at start of TAD very positive, scores near end of TAD very negative
// score tads that way? +scores of first n bins, -scores of last n bins? --> working pretty well, but have extra small tads around true boundaries that shouldn't be there - still need to penalize or differentiate tads which are not ideal
// subtract distances from nearest max (but instead location where sign changes), like with IS?



	//windowsize := 3

        starttime := time.Now()
        score := 0.0

	/*
	// first find sign changes in DI, if this works, move this to precompute function - only care when DI goes from negative to positive
	disignchanges := make([]int,len(discores))
	firstloc := 0
	num := 0
	for i,dival := range discores {
		if i == 0 || i == len(discores)-1 { continue }
		if dival < 0 && discores[i+1] > 0 {
			if firstloc == 0 {firstloc = i}
			disignchanges[num] = i
			num++
		} else if dival > 0 && discores[i-1] < 0 {
			if firstloc == 0 {firstloc = i}
			disignchanges[num] = i
			num++
		}
	}
	disignchanges = disignchanges[:num]


	didists := make([]int, len(discores))
	lastmin := 0
	nextmin := firstloc
	minidx := 0
	for i := 0; i < len(discores); i++ {
		didists,lastmin,nextmin,minidx = updateMinDists(disignchanges,didists,i,lastmin,nextmin,minidx,len(discores))
	}*/


	//var bdyweight float64
        for _,tad := range *tadset {
		//for i:=tad[0]; i < tad[1]+1; i++{
		//fmt.Println(tad,discores[tad[0]],discores[tad[1]])
		//}
		tadend := tad[1]
		if tadend >= len(discores) { tadend = len(discores)-1 }
		tadscore := 0.0
		for x := 0; x < windowsize; x++ {
			if tad[0]+x < len(discores) { 
				tadscore += discores[tad[0]+x] 
			}
			if tadend-x > 0 { 
				tadscore -= discores[tadend-x] 
			}
		}
		
		distscore := float64(didists[tad[0]] + didists[tadend])
		score += tadscore - distscore
		//fmt.Println(tad,tadscore,didists[tad[0]],didists[tadend])
		/*for idx,bdy := range tad {
			avgafter := 0.0
			numafter := 0
			avgbefore := 0.0
			numbefore := 0
			maxloc := 0
			maxdiff := -1e10
			/*if bdy+idx < len(discores) && bdy-idx >=0 {
				bdyweight = discores[bdy+idx]-discores[bdy-idx] // if idx == 0, want to subtract bdy - (bdy-1), if idx == 1, want to subtract (bdy+1)-bdy
				if bdy > 0 && bdy+1 < len(discores) {fmt.Println(idx,bdy,bdyweight,discores[bdy-1:bdy+2])}
			} else {
				bdyweight = 0.0
			}*/
			
		/*	for x := 0; x < windowsize; x++ {
				if bdy-x > 0 {
					avgbefore += discores[bdy-x-1+idx] //if idx == 1, include bdy with befores
					numbefore += 1
					
				}
				if bdy+x+idx < len(discores) {
					avgafter += discores[bdy+x+idx] //if idx == 0, include bdy with afters
					numafter += 1
				}
			}
			if numafter > 0 && numbefore > 0 {
				score += avgafter/float64(numafter) - avgbefore/float64(numbefore)
			}
		}
		//fmt.Println(tad,score)
                //tadmidpt := (tad[0]+tad[1])/2
                //for _,i := range tad {
                //        if i > tadmidpt {
                //                score -= discores[i]
		//		} else {
                //                score += discores[i]
		//			}
		//	}
                //fmt.Println(tad,score)
		//for _,bdy := range tad {
		//      score += discores[bdy]
		//}*/
	}
	//fmt.Println("final score",score)
	//os.Exit(1)
        return score,time.Since(starttime)
}
