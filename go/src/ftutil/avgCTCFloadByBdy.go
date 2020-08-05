package ftutil

import "runtime"

func AvgCTCFloadByBdy(tadset [][]int, ctcfweights []int) ([]float64) {

	runtime.GOMAXPROCS(1)

	ctcfload := make([]float64,len(tadset)-2)
	var prevbdy int
	//bdycount := 0
	
	for tadidx,tad := range tadset {
		if tadidx > 0 && tadidx < len(ctcfload) {
			ctcfload[tadidx-1] = float64(prevbdy + ctcfweights[tad[0]])
		}
		prevbdy = ctcfweights[tad[1]]
		//bdycount += 2
	}
	return ctcfload
}
