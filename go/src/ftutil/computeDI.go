package ftutil

import(
	"math"
	"hicutil"
)

func computeDI(loc int, dilen int, hicmap map[hicutil.Pair]float64) float64 {

        upstream := 0.0
        downstream := 0.0
        for i := 1; i <= dilen; i++{
                upstream += hicmap[hicutil.Pair{loc-i,loc}]
                downstream += hicmap[hicutil.Pair{loc,loc+i}]
        }
        if upstream == downstream {
                return 0.0
        }
        expected := (upstream+downstream)/2.0
        diridx := ((downstream-upstream)/math.Abs(downstream-upstream))*(math.Pow(upstream-expected,2)/expected + math.Pow(downstream-expected,2)/expected)
        return diridx
}
