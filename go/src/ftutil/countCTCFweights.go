package ftutil

import ("runtime"
	//"fmt"
	)

func CountCTCFweights(ctcfmap map[string][]int) (map[string][]int) {

	runtime.GOMAXPROCS(1)

	weightedmap := make(map[string][]int)

	for key,clist := range ctcfmap {
		wlist := make([]int, clist[len(clist)-1]+1)
		for _,c := range clist {
			wlist[c] += 1
		}
		weightedmap[key] = wlist
	}
	return weightedmap
}
