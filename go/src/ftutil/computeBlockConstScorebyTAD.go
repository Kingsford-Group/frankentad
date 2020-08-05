package ftutil

import(
	"hicutil"
	"time"
//	"os"
//	"fmt"
	"math"
)

func ComputeBlockConstScorebyTAD(tad []int, chrlength int, nontadmean float64, tadmeans map[hicutil.Pair]float64, hicmap map[hicutil.Pair]float64, blockconstscores map[hicutil.Pair]float64) (float64,map[hicutil.Pair]float64,time.Duration) {
//s map[hicutil.Pair]
        starttime := time.Now()
        score := 0.0

	if tadscore, ok := blockconstscores[hicutil.Pair{tad[0],tad[1]}]; ok {
		score += tadscore
	} else {
		tadscore := 0.0
		//fmt.Println(tad)
		for a := 0; a <= tad[1]; a++ {
			bstart := math.Max(float64(a),float64(tad[0]))
			//fmt.Println("row =",a,"mincol =",bstart)
			for b := int(bstart); b <= tad[1]; b++ {
				if a >= tad[0] {
					tadscore -= math.Pow(hicmap[hicutil.Pair{a,b}] - tadmeans[hicutil.Pair{tad[0],tad[1]}],2)
				} else {
					tadscore -= math.Pow(hicmap[hicutil.Pair{a,b}] - nontadmean,2)
				}
			}
		}
		blockconstscores[hicutil.Pair{tad[0],tad[1]}] = tadscore
		score += tadscore
	}

        return score,blockconstscores,time.Since(starttime)
}
