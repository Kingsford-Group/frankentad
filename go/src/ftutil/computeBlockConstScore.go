package ftutil

import(
	"hicutil"
	"time"
//	"os"
//	"fmt"
	"math"
)

func ComputeBlockConstScore(tadset *[][]int, chrlength int, nontadmean float64, tadmeans map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, blockconstscores map[hicutil.Pair]float64) (float64,map[hicutil.Pair]float64,time.Duration) {
//s map[hicutil.Pair]
        starttime := time.Now()
        score := 0.0

        for _, tad := range *tadset {
                //fmt.Println(tad)
                if tadscore, ok := blockconstscores[hicutil.Pair{tad[0],tad[1]}]; ok {
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
                        blockconstscores[hicutil.Pair{tad[0],tad[1]}] = tadscore
                        score += tadscore
                        //fmt.Println(tad,tadscore)
                }
                //os.Exit(1)
        }
        //fmt.Println(score)
        //os.Exit(1)
        return score,blockconstscores,time.Since(starttime)
}
