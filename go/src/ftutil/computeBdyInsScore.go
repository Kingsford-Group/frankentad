package ftutil

import(
	"time"
	"math"
//	"fmt"
//	"os"
)

func ComputeBdyInsScore(allisscores *[]float64, mindists *[]int, tadset *[][]int, alpha float64) (float64,time.Duration) {

	// need to modify this fcn to match actual TopDom results better - scores their TADs surprisingly badly even though they're actually really good


        //runtime.GOMAXPROCS(1)
        starttime := time.Now()

        // compute insulation score for a particular tad set
        totalscore := 0.0
        steepscore := 0.0
        chrlength := len(*allisscores)
	distscore := 0

	windowsize := 3

        //closestbdy := 0
        //tadidx := 0
        //weight := 0.0
        /*for i := 0; i < chrlength; i++ {
                score := math.Abs(math.Abs(float64(i-closestbdy)) - float64((*mindists)[i]))
		fmt.Println("closest bdy",i-closestbdy,(*mindists)[i])
                // add weight for slope around i - worse to be off on strong boundaries than weak ones
                // this doesn't make sense for all points - only want steepness to matter around actual boundary points
                if (i == closestbdy) && i > windowsize-1 && i < chrlength-windowsize-1 {
                        centerval := (*allisscores)[i]
                        for j := i-windowsize; j < i+windowsize+1; j++{
                                if j < i {
                                        steepscore += ((*allisscores)[j] - centerval)/float64(windowsize)
                                        //fmt.Println(i,j,((*allisscores)[j] - centerval)/10)
                                } else if j > i {
                                        steepscore += (centerval - (*allisscores)[j])/float64(windowsize)
                                        //fmt.Println(i,j,(centerval-(*allisscores)[j])/10)
                                }
                        }
                        fmt.Println("steepscore",i,steepscore)
                }
                if math.IsNaN(score) {
                        fmt.Println("insulation score of NaN from dists")
                        fmt.Println(i, closestbdy, (*mindists)[i])
                        os.Exit(1)
                }
                if math.IsNaN(steepscore){
                        fmt.Println("insulation score of NaN from slope")
                        fmt.Println(i, (*allisscores)[i-windowsize:i+windowsize+1])
                        os.Exit(1)
                }
                totalscore -= score
                //fmt.Println(i,closestbdy,math.Abs(float64(i-closestbdy)),(*mindists)[i])
                if closestbdy == (*tadset)[tadidx][0] && i+1 >= ((*tadset)[tadidx][0]+(*tadset)[tadidx][1])/2 {
                        // switching from closest boundary being prev bdy to next bdy
                        closestbdy = (*tadset)[tadidx][1]
                } else if i == (*tadset)[tadidx][1] && tadidx+1 < len(*tadset) {
                        // switching to next tad
                        closestbdy = (*tadset)[tadidx+1][0]
                        tadidx +=1
                }
        }*/
        for _, tad := range *tadset {
		//taddistscore := 0.0
		//tadsteepscore := 0.0
		//fmt.Println(tad,len(*mindists))
		if tad[1] >= len(*mindists) { 
			distscore += (*mindists)[tad[0]] + (*mindists)[len(*mindists)-1]
			//taddistscore = (*mindists)[tad[0]] + (*mindists)[len(*mindists)-1]
		} else {
			distscore += (*mindists)[tad[0]] + (*mindists)[tad[1]]
			//taddistscore = (*mindists)[tad[0]] + (*mindists)[tad[1]]
		}
                //startis := (*allisscores)[tad[0]]
                slopeavg := 0.0
                startbin := tad[0] - windowsize
                if startbin < 0 {startbin = 0}
                endbin := tad[0] + 1 + windowsize
                if endbin > chrlength - 1 {endbin = chrlength-1}
                numbins := 0
                for i := startbin; i < endbin; i++ {
                        //fmt.Println(i, (*allisscores)[i])
                        //slopeavg += (*allisscores)[i] - startis
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
                //fmt.Println("avg slope around tad bdy:",tad[1],slopeavg/float64(numbins))
                //fmt.Println(tad)
                //if tad[1] < len(*mindists)-1 {fmt.Println("distances to nearest local mins =",(*mindists)[tad[0]], (*mindists)[tad[1]])}//, slopeavg/float64(numbins))
                //if tad[0] > 100 { os.Exit(1)}
        }
        //totalscore = 1.0/float64(distscore) + alpha*steepscore
	totalscore = -float64(distscore) + alpha*steepscore
        //totalscore += alpha*steepscore
        //fmt.Println(totalscore)
        //os.Exit(1)
        return totalscore,time.Since(starttime)
}
