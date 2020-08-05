package ftutil


import(
	"hicutil"
	"math"
)


func computeDirIdxTT(hicmap  map[hicutil.Pair]float64, a int, p int, q int, chrlength int) float64{

        var kmax int
        var lmin,lmax int

        dival := 0.0
        if a == 0 { return 0.0 }

        if a < p || a+p > chrlength { kmax = a+1 } else { kmax = p+1 }

        if a-q < 0 {lmin = 0} else {lmin = a-q}
        if a+q >= chrlength {lmax = chrlength} else {lmax = a+q}
        if a-p < 0 {kmax = a+1} else if a+p >= chrlength {kmax = chrlength-a-1} else {kmax = p+1}

        for l := lmin; l < lmax; l++ {
                innersum := 0.0
                for k := 1; k < kmax; k++ {
                        innersum += hicmap[hicutil.Pair{l,a+k}] - hicmap[hicutil.Pair{l,a-k}]
                }
                dival += math.Abs(innersum)
        }
        /*tosubtract := 0.0
        toadd := 0.0
        for k := 1; k < kmax; k++ {
                if a > q {
                        tosubtract += (hicmap[hicutil.Pair{a-q,a+k}] - hicmap[hicutil.Pair{a-q, a-k}])
                }
                if a + q < chrlength {
                        toadd += (hicmap[hicutil.Pair{a+q,a+k}] - hicmap[hicutil.Pair{a+q,a-k}])
                }
        }

        dival = prevvalue - math.Abs(tosubtract) + math.Abs(toadd)
        fmt.Println(a,naivedi,dival)*/

        return dival
}
