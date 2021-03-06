package ftutil

import(
	"time"
	"hicutil"
	"fmt"
	"os"
	"math"
)

func ComputeGraphDensScorebyTAD(allgdscores map[hicutil.Pair]float64, meanscores *[][]float64, tad []int, gamma float64) (float64,time.Duration) {
        // compute graph density (Armatus) score for particular tad set

        //runtime.GOMAXPROCS(1)
        timestart := time.Now()
        //var tadscore,normfactor float64

        totalscore := 0.0
	
	tadscore := allgdscores[hicutil.Pair{tad[0],tad[1]}]
	normfactor := math.Pow(float64(tad[1]-tad[0]+1), gamma) // armatus code has b-a+1
	if tad[1]-tad[0] < 0 {
		fmt.Println("mistake in TAD set detected")
		fmt.Println(tad)
		os.Exit(1)
	}
	if tad[1]-tad[0] < len(*meanscores) {
		totalscore += tadscore/normfactor - ( (*meanscores)[tad[1]-tad[0]][0] / ( (*meanscores)[tad[1]-tad[0]][1] * normfactor ) )
	}
	return totalscore,time.Since(timestart)
}
