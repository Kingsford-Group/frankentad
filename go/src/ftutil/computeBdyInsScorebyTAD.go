package ftutil

import(
	"time"
	"math"
//	"fmt"
//	"os"
)

func ComputeBdyInsScorebyTAD(allisscores *[]float64, mindists *[]int, tad []int, alpha float64) (float64,time.Duration) {

        starttime := time.Now()

        // compute insulation score for a particular tad
        tadscore := 0.0
        steepscore := 0.0
        chrlength := len(*allisscores)
	distscore := 0

	windowsize := 3

	if tad[1] >= len(*mindists) { 
		distscore += (*mindists)[tad[0]] + (*mindists)[len(*mindists)-1]
		//fmt.Println("distance to nearest min for tad bdy:",tad[0],(*mindists)[tad[0]])
		//fmt.Println("distance to nearest min for tad bdy:",tad[1],(*mindists)[len(*mindists)-1])
	} else {
		distscore += (*mindists)[tad[0]] + (*mindists)[tad[1]]
		//fmt.Println("distance to nearest min for tad bdy:",tad[0],(*mindists)[tad[0]])
		//fmt.Println("distance to nearest min for tad bdy:",tad[1],(*mindists)[tad[1]])
	}
	slopeavg := 0.0
	startbin := tad[0] - windowsize
	if startbin < 0 {startbin = 0}
	endbin := tad[0] + 1 + windowsize
	if endbin > chrlength - 1 {endbin = chrlength-1}
	numbins := 0
	for i := startbin; i < endbin; i++ {
		numbins += 1
		if i < tad[0] {
			slopeavg += (*allisscores)[i] - (*allisscores)[i+1]
		} else {
			slopeavg += (*allisscores)[i+1] - (*allisscores)[i]
		}
	}
	steepscore += math.Abs(slopeavg)/float64(numbins)
	//fmt.Println("avg slope around tad bdy:",tad[0],math.Abs(slopeavg)/float64(numbins))
	numbins = 0
	slopeavg = 0.0
	//endis := (*allisscores)[tad[1]]
	startbin = tad[1] - windowsize
	if startbin < 0 {startbin = 0}
	endbin = tad[1] + 1 + windowsize
	if endbin > chrlength - 1 {endbin = chrlength-1}
	for i := startbin; i < endbin; i++ {
		//fmt.Println(i, (*allisscores)[i])
		//slopeavg += (*allisscores)[i] - endis
		numbins += 1
		if i < tad[1] {
			slopeavg += (*allisscores)[i] - (*allisscores)[i+1]
		} else {
			slopeavg += (*allisscores)[i+1] - (*allisscores)[i]
		}
	}
	steepscore += slopeavg/float64(numbins)
	//fmt.Println("avg slope around tad bdy:",tad[1],math.Abs(slopeavg)/float64(numbins))
	tadscore = -float64(distscore) + alpha*steepscore
        return tadscore,time.Since(starttime)
}
