package hicutil

import( 
//	"fmt" 
//	"os"
)

func PostprocessTADs(tadset [][]int, tadvals map[Pair]float64, centromereloc []int) ([][]int) {

	var cleantads [][]int

	//fmt.Println("starting with", tadset, centromereloc)
	cmlen := centromereloc[1] - centromereloc[0]
	for _,tad := range tadset {
		tadsum := tadvals[Pair{tad[0],tad[1]}]
		//fmt.Println(tad,tadsum)
		if tadsum < 1e-10 { continue }
		//fmt.Println("nonzero tad values", tad)
		if tad[0] < centromereloc[0] && tad[1] > centromereloc[0] {
			//fmt.Println("first if")
			cleantads = append(cleantads, []int{tad[0],centromereloc[0]-1}) // tad ends just before centromere
			if tad[1] > centromereloc[1] {
				//fmt.Println("sub if")
				cleantads = append(cleantads,[]int{centromereloc[1]+1,tad[1]+cmlen}) 
			}
		} else if tad[0] >= centromereloc[0] {
			//fmt.Println("second if")
			cleantads = append(cleantads, []int{tad[0]+cmlen,tad[1]+cmlen})
		} else {
			//fmt.Println("last")
			cleantads = append(cleantads, tad)
		}
	}
	//fmt.Println("ending with", cleantads)
	//os.Exit(1)
	return cleantads
}
