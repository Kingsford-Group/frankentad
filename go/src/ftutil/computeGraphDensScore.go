package ftutil

import(
	"time"
	"hicutil"
	"fmt"
	"os"
	"math"
)

func ComputeGraphDensScore(allgdscores map[hicutil.Pair]float64, meanscores *[][]float64, tadset *[][]int, gamma float64) (float64,time.Duration) {
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
                //      normfactor = 1.0
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
}
