package ftutil

import(
	"hicutil"
	"math"
	"time"
	"fmt"
	"os"
)

func ComputeLinEnrichScore(tadset *[][]int, chrlength int, gamma float64, diridxtt []float64, avgfreq []float64, atilde map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, linenrichscores map[hicutil.Pair]float64) (float64, map[hicutil.Pair]float64, time.Duration) {

        starttime := time.Now()
        score := 0.0

        for _, tad := range *tadset {
                if tadscore, ok := linenrichscores[hicutil.Pair{tad[0],tad[1]}]; ok {
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
                        linenrichscores[hicutil.Pair{tad[0],tad[1]}] = tadscore
                }
                //fmt.Println(tad,score)
                //os.Exit(1)
        }
        //fmt.Println(score)
        //os.Exit(1)
        return score,linenrichscores,time.Since(starttime)
}
