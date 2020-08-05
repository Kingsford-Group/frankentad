package ftutil

import(
	"time"
//	"fmt"
//	"os"
)

func ComputeBiasDirScorebyTAD(discores []float64, didists []int, tad []int, windowsize int) (float64,time.Duration) {

// ideally: scores at start of TAD very positive, scores near end of TAD very negative
// score tads that way? +scores of first n bins, -scores of last n bins? --> working pretty well, but have extra small tads around true boundaries that shouldn't be there - still need to penalize or differentiate tads which are not ideal
// subtract distances from nearest max (but instead location where sign changes), like with IS?

        starttime := time.Now()

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
	score := tadscore - distscore
	
        return score,time.Since(starttime)
}
