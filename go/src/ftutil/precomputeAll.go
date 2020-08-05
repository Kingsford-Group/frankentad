package ftutil

import(
	"hicutil"
	"math"
	"fmt"
	"os"
)

type LinRegParams interface {
	// y = ax+b
	getA() float64
        getB() float64
}

type LinRegIntermed struct {

        n float64
        dx,dy float64
	varX,covXY float64
        meanX,meanY float64
}

func (x LinRegIntermed) getA() float64 {
        if x.varX > 0 {
                return x.covXY/x.varX
        } else {
                return 0
        }
}

func (x LinRegIntermed) getB() float64 {
        return x.meanY - x.getA()*x.meanX
}


func PrecomputeAll(hicmap map[hicutil.Pair]float64, chrlength int, iswindowsize int, p int, q int,res int) (map[hicutil.Pair]float64, [][]float64, []float64, []int, []float64, []float64, map[hicutil.Pair]float64, map[hicutil.Pair]float64, float64, []float64, []int) {

        //res int, diwindow int)

        //2nd to last return: map[hicutil.Pair]

        // loop through Hi-C matrix once to precompute everything needed
        //for Armatus scoring
        gdscores := make(map[hicutil.Pair]float64) // contain scores for each (a,b) interval on HiC matrix NOT NORMALIZED BY LENGTH^GAMMA
        meangdscores := make([][]float64, chrlength) // sum of scores for all tads of given length (length = array index) and number of intervals of that length (to compute mean)

        //for TopDom/Ins Score scoring
        allisscores := make([]float64, chrlength)
        sumisscores := 0.0
        minlocs := make([]int, chrlength-2*iswindowsize)

        //for TADTree scoring
        diridxtt := make([]float64, chrlength)
        avgfreq := make([]float64, chrlength)
        exptadfreq := make(map[hicutil.Pair]float64)
        linregvals :=make(map[hicutil.Pair]LinRegIntermed)

        //for HiCSeg scoring
        nontadmean := 0.0
        tadmeans := make(map[hicutil.Pair]float64)
        nontadmeans := make(map[hicutil.Pair]float64)
        //nontadmeans2 := make(map[hicutil.Pair]float64)
        colsums := make([][]float64,chrlength)

        //for DI scoring
        diridx := make([]float64,chrlength)
        //diridxscores := make([]float64,chrlength)
	disignchanges := make([]int,chrlength)
	dilen := 20//00000/res
        //diwindow := 5

        cornercount := 0
        for a := 0; a < chrlength/4; a++ {
                for b:= 3*chrlength/4; b < chrlength; b++ {
                        nontadmean += hicmap[hicutil.Pair{a,b}]
                        cornercount++
                }
        }
        nontadmean = nontadmean/float64(cornercount)
        //fmt.Println(nontadmean)

	meangdscores[0] = make([]float64,2)
        colsums[0] = make([]float64,chrlength)
        for a := 0; a < chrlength; a++ {
                // initialize armatus scoring
                gdscores[hicutil.Pair{a,a}] = hicmap[hicutil.Pair{a,a}]
                meangdscores[0][0] += gdscores[hicutil.Pair{a,a}]
                meangdscores[0][1] += 1.0

                //tad means for hicseg scoring
                tadmeans[hicutil.Pair{a,a}] = hicmap[hicutil.Pair{a,a}]

                colsums[0][a] = hicmap[hicutil.Pair{0,a}]

                //compute directionality index for TADtree objective
                diridxtt[a] = computeDirIdxTT(hicmap,a,p,q,chrlength)
                avgfreq[0] += hicmap[hicutil.Pair{a,a}] * (1.0/float64(chrlength))
                //if hicmap[hicutil.Pair{a,a}] > 0 {fmt.Println(hicmap[hicutil.Pair{a,a}])}
        }

        colsum := 0.0
        minidx := 0
	numdi := 0
        for b := 1; b < chrlength; b++ {
                colsums[b] = make([]float64,chrlength)
                //directionality index
                if b > dilen && b < chrlength-dilen {
                        diridx[b] = computeDI(b,dilen,hicmap)
                } else {
                        diridx[b] = -1
                }
                //fmt.Println(b,diridxscores[b])
                //fmt.Println(b)
                colsum = hicmap[hicutil.Pair{b,b}]
                if b >= iswindowsize/2 && b < chrlength-iswindowsize/2+1 {
                        for c := b; c < b+iswindowsize/2+1; c++ {
                                allisscores[b] += hicmap[hicutil.Pair{b,c}]
                        }
                }
                colsums[b][chrlength-1] = colsums[b-1][chrlength-1]+hicmap[hicutil.Pair{b,chrlength-1}]
                for a := b-1; a > -1; a-- {
                        //nontadmeans[hicutil.Pair{a,b}] += rowsum/float64((a-1)*b)
                        x := int(math.Abs(float64(a-b+1)))

                        //hicseg nontadmeans computed per rectangle instead of 1 overall
                        if x > 0 {colsums[x][b] = colsums[x-1][b] + hicmap[hicutil.Pair{x,b}]}
                        if b > 2 && a > 1 {nontadmeans[hicutil.Pair{a,b-1}] = (nontadmeans[hicutil.Pair{a,b-2}]*float64((a-1)*(b-2)) + colsums[a-1][b-1])/float64((a-1)*(b-1))}
                        if b == chrlength-1 && a > 0 { nontadmeans[hicutil.Pair{a,b}] = (nontadmeans[hicutil.Pair{a,b-1}]*float64((a-1)*(b-1)) + colsums[a-1][b])/float64((a-1)*b) }
                        /*for row := 0; row < a; row++ {
                                for col := a; col <= b; col++ {
                                        nontadmeans[hicutil.Pair{a,b}] += hicmap[hicutil.Pair{row,col}]//float64((a-1)*b)
                                }
                        }
                        if b > 1 && b < 100 {
                                fmt.Println(a, b-1, nontadmeans[hicutil.Pair{a,b-1}],nontadmeans2[hicutil.Pair{a,b-1}])
                        } else if b > 100 {
                                os.Exit(1)
                        }*/


                        //if b-a > chrlength-100 { continue }
                        avgfreq[a] += hicmap[hicutil.Pair{b,b+a}] * (1.0/float64(chrlength-a))

                        // Armatus scoring
                        colsum += hicmap[hicutil.Pair{a,b}]
                        gdscores[hicutil.Pair{a,b}] = gdscores[hicutil.Pair{a,b-1}] + colsum

                        //fmt.Println("a=",a,"b=",b,"score=",allscores[hicutil.Pair{a,b}])
                        if len(meangdscores[b-a]) == 0 {
                                meangdscores[b-a] = make([]float64,2)
                        }
                        meangdscores[b-a][0] += gdscores[hicutil.Pair{a,b}]
                        meangdscores[b-a][1] += 1.0
                        // hicseg scoring
                        tadmeans[hicutil.Pair{a,b}] = gdscores[hicutil.Pair{a,b}]/triangleNum(b-a+1)
                        //fmt.Println(a,b,tadmeans[hicutil.Pair{a,b}],triangleNum(b-a+1))
                        //fmt.Println(a,b,"hicval =",hicmap[hicutil.Pair{a,b}])

                        // insulation score/ topdom scoring
                        if b >= iswindowsize/2 && b < chrlength-iswindowsize/2+1 {
                                if a >= b-iswindowsize/2 {
                                        for c := b; c < b+iswindowsize/2+1; c++{
                                                allisscores[b] += hicmap[hicutil.Pair{a,c}]
                                                //if b == 20 {fmt.Println(a,c,hicmap[hicutil.Pair{a,c}])}
                                        } // maybe missing the row where a = b because a starts at b-1
                                        allisscores[b] += hicmap[hicutil.Pair{b,a+iswindowsize/2+1}]
                                        //if b == 20 {fmt.Println(b,a+iswindowsize/2)}
                                }
                                if allisscores[b] < 0 {
                                        fmt.Println("got negative ins score")
                                        fmt.Println(allisscores[b],b)
                                        os.Exit(1)
                                }
                                sumisscores += allisscores[b]
                        }
                }
                if b > 1 &&  allisscores[b] > allisscores[b-1] && allisscores[b-2] > allisscores[b-1] {
                        minlocs[minidx] = b-1
                        minidx += 1
                        //fmt.Println("local min found at",b)
                }
		// find where DI goes from <0 to >0
		if b > 1 && (diridx[b] > 0 && diridx[b-1] < 0) {
			disignchanges[numdi] = b-1
			numdi++
			disignchanges[numdi] = b
			numdi++
		}
        }
        //for i:= 0; i < len(avgfreq); i++ {
        //      fmt.Println(i,avgfreq[i])
        //}
        minlocs[minidx] = chrlength-1
        minlocs = minlocs[:minidx+1] // note this just finds all things that are less than the two before and 1 after -- is this the best way to identify minima? maybe consider changing or at least testing to see if it works decently
        disignchanges[numdi] = chrlength-1
	disignchanges = disignchanges[:numdi+1]
	//fmt.Println("done with primary calculations")
        // normalize all scores as log_2 (binscore/meanscores)
        normfactor := sumisscores/float64(chrlength - 2*iswindowsize)
        mindists  := make([]int, chrlength)
        lastmin := 0
        nextmin := minlocs[0]
        minidx = 0

	didists := make([]int, chrlength)
	lastmindi := 0
	nextmindi := disignchanges[0]
	diidx := 0
        //var xvals []float64
        //var yvals []float64
        //disum := 0.0
        for i := 0; i < chrlength; i++ {
                /*// compute differences in DI for DI scoring
                if i > 0 {
                        diridxscores[i] = diridx[i] - diridx[i-1]
                        //for j := 0; j < diwindow; j++{
                                //diridxscores[i] +=  diridx[i+1+j] - diridx[i-j]
                        //}
                        //disum += diridxscores[i]
                } else { diridxscores[i] = 0 }*/
                //fmt.Println(i,diridx[i],diridxscores[i])
                //fmt.Println(i)
                if i >= iswindowsize/2 && i < chrlength-iswindowsize/2+1 {
                        if allisscores[i] != 0 {
                                allisscores[i] = math.Log2(allisscores[i]/normfactor)
                                if math.IsNaN(allisscores[i]) {
                                        fmt.Println("uh oh, an insulation score is NaN")
                                        os.Exit(1)
                                }
                        } else if i+1 < chrlength && allisscores[i+1] != 0 {
                                allisscores[i] = math.Log2(allisscores[i+1]/normfactor)-1
                        } else {
                                allisscores[i] = -10000
                        }
                        if math.IsInf(allisscores[i],0) {
                                fmt.Println("IS score of infinity")
                                fmt.Println(i, normfactor)
                                os.Exit(1)
                        }
                }

                /*xvals = []float64{0.0}
                if avgfreq[0] > 0 {
                        yvals = []float64{hicmap[hicutil.Pair{i,i}] / avgfreq[0]}
                } else {
                        yvals = []float64{0.0}
                }*/
                //fmt.Println("done with ins score")
                mindists,lastmin,nextmin,minidx = updateMinDists(minlocs,mindists,i,lastmin,nextmin,minidx,chrlength)
                didists,lastmindi,nextmindi,diidx = updateMinDists(disignchanges,didists,i,lastmindi,nextmindi,diidx,chrlength)
		//fmt.Println("done updating min")
                // compute Atilde from TADtree - linear regression for every tad is crazy slow, how to make this better??
                //exptadfreq[i],prevxvals,prevyvals = calcExpectedContactFreq(hicmap,prevxvals,prevyvals,i,chrlength,avgfreq)
                for j := i+1; j < chrlength; j++ {
                        //fmt.Println(j)

                        //if j-i > chrlength-100 { continue }
                        //exptadfreq[hicutil.Pair{i,j}] = calcExpectedContactFreq(hicmap,i,j,avgfreq)
                        // update linregvals for all new points
                        linregvals[hicutil.Pair{i,j}] = linregvals[hicutil.Pair{i,j-1}]
                        newvals := linregvals[hicutil.Pair{i,j}]
                        //var avgfreqval float64
                        for k := 1; k <= j - i; k++ {
                                var newy float64
                                //xvals = append(xvals,float64(k))
                                if avgfreq[k] == 0 {
                                        newy = 0
                                } else {
                                        newy = hicmap[hicutil.Pair{j-k,j}] / avgfreq[k]
                                }
                                if math.IsNaN(newy) {
                                        fmt.Println(i,j,k,hicmap[hicutil.Pair{j-k,j}],avgfreq[k])
                                        os.Exit(1)
                                }
                                //avgfreqval = avgfreq[k]
                                //yvals = append(yvals,newy) // change to be 0 if avgfreq == 0?
                                //fmt.Println(k,avgfreq[k],hicmap[hicutil.Pair{j-k,j}])
                                //fmt.Println(hicmap[hicutil.Pair{j-k,j}] / avgfreq[k])

                                //newy := hicmap[hicutil.Pair{j-k,j}]/avgfreq[k]
                                newvals.n += 1
                                newvals.dx = float64(k) - newvals.meanX
                                newvals.dy = newy - newvals.meanY
                                newvals.varX += (((newvals.n-1)/newvals.n)*newvals.dx*newvals.dx - newvals.varX)/newvals.n
                                newvals.covXY += (((newvals.n-1)/newvals.n)*newvals.dx*newvals.dy - newvals.covXY)/newvals.n
                                newvals.meanX += newvals.dx/newvals.n
                                newvals.meanY += newvals.dy/newvals.n
                        }

                        //fmt.Println(xvals,yvals)
                        //if j > 5 {os.Exit(1)}
                        // compute lin reg params
                        delta := newvals.getA()
                        //if math.IsNaN(delta) {
                        //      delta = 0
                        //      newvals
                        //}
                        beta := newvals.getB()
                        if math.IsNaN(delta) || math.IsNaN(beta) {
                                fmt.Println(delta,beta)
                                fmt.Println(newvals.n)
                                fmt.Println(newvals.dx)
                                fmt.Println(newvals.dy)
                                fmt.Println(newvals.varX)
                                fmt.Println(newvals.covXY)
                                fmt.Println(newvals.meanX)
                                fmt.Println(newvals.meanY)
                        }
                        // need to compare to delta, beta from normal linear regression
                        /*beta2, delta2 := stat.LinearRegression(xvals, yvals, nil, false)
                        if math.IsNaN(delta2) || math.IsNaN(beta2) {
                                fmt.Println("linear regression returned nan")
                                fmt.Println(xvals,yvals,avgfreqval,hicmap[hicutil.Pair{i,j}])
                                os.Exit(1)
                        }
                        if delta - delta2 > 0.1 || beta - beta2 > 0.1 {
                                fmt.Println("linear regression results dont match up")
                                fmt.Println(i,j)
                                fmt.Println(delta,delta2)
                                fmt.Println(beta,beta2)
                        }*/
                        linregvals[hicutil.Pair{i,j}] = newvals
                        // calc actual exptadfreq value
                        exptadfreq[hicutil.Pair{i,j}] = (float64(j-i)*delta + beta)*avgfreq[j-i] // should this be i-j or i-j+1???
                        if math.IsNaN(exptadfreq[hicutil.Pair{i,j}]) {
                                fmt.Println("exp tad freq NaN")
                                fmt.Println("HiC coords",i,j)
                                fmt.Println(delta,beta)
                                fmt.Println(avgfreq[j-i])
                                os.Exit(1)
                        }
                }
                // do we need ANOTHER loop through everything to compute full TAD scores from linear regression values? the TADTree scoring function is the absolute worst. 3 quadratic loops -_-

                //if i > 1 {os.Exit(1)}
        }
        return gdscores, meangdscores, allisscores, mindists, diridxtt, avgfreq, exptadfreq, tadmeans, nontadmean, diridx, didists//scores
}


func triangleNum(n int) float64 {
	if n < 2 {
                fmt.Println("cant compute triangle number for n =",n)
                os.Exit(1)
		}
        return float64(n*(n+1))/2.0
}


func updateMinDists(minlocs []int, mindists []int, i int, lastmin int, nextmin int, minidx int, chrlength int) ([]int,int,int,int) {

        if lastmin == 0 {
                mindists[i] = nextmin - i
        }
        if i == nextmin {
                mindists[i] = 0
                lastmin = i
                if len(minlocs) > minidx+1 {
                        nextmin = minlocs[minidx+1]
                } else {
                        nextmin = chrlength*2
                }
                minidx += 1
        } else if (i - lastmin >= nextmin - i) && lastmin != 0 {
                mindists[i] = nextmin - i
        } else if (nextmin - i > i - lastmin) && lastmin != 0 {
                mindists[i] = i - lastmin
        }
        return mindists,lastmin,nextmin,minidx
}

