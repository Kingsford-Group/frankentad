package ftutil

import(
	"hicutil"
	"math"
	"time"
	"fmt"
	"os"
)

func ComputeLinEnrichScorebyTAD(tad []int, chrlength int, gamma float64, diridxtt []float64, avgfreq []float64, atilde map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, linenrichscores map[hicutil.Pair]float64) (float64, map[hicutil.Pair]float64, time.Duration) {

        starttime := time.Now()
        score := 0.0
	
	if tadscore, ok := linenrichscores[hicutil.Pair{tad[0],tad[1]}]; ok {
		score += tadscore
	} else {
		tadscore := 0.0
		tadscore += gamma * diridxtt[tad[0]] + gamma * diridxtt[tad[1]]
		for a := tad[0]; a <= tad[1]; a++ {
			for b := a; b < chrlength; b++{
				//fmt.Println(a,b)
				if b <= tad[1] {
					tadscore -= (1.0-gamma)*math.Pow((atilde[hicutil.Pair{a,b}] - hicmap[hicutil.Pair{a,b}]),2)
				} else {
                                        tadscore -= (1.0-gamma)*math.Pow(avgfreq[b-a] - hicmap[hicutil.Pair{a,b}],2)
				}
				if math.IsNaN(tadscore) {
					fmt.Println(tad,a,b)
					os.Exit(1)
				}
			}
		}
		score += tadscore
		linenrichscores[hicutil.Pair{tad[0],tad[1]}] = tadscore
	}
        return score,linenrichscores,time.Since(starttime)
}
