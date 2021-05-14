package ftutil

import(
//	"runtime"
	"fmt"
	"math"
	"os"
)


func AvgCTCFload(tadset [][]int, ctcfweights []int) (float64) {

	//runtime.GOMAXPROCS(1)

	ctcfload := 0
	bdycount := 0
	
	//, centromereloc []int
	//cmlen := centromereloc[1] - centromereloc[0]

	for _,tad := range tadset {
		if tad[1] < len(ctcfweights) {
			ctcfload += ctcfweights[tad[0]] + ctcfweights[tad[1]]
		} else if ( tad[0] < len(ctcfweights) && tad[1] > len(ctcfweights) ) {
			ctcfload += ctcfweights[tad[0]] + ctcfweights[len(ctcfweights)-1]
		}
		bdycount += 2
		/*if tad[1] < centromereloc[1] || tad[0] == centromereloc[0] { // TAD ends before centromere does, or tad starts at centromere start (and is therefore equivalent to centromere)
			ctcfload += ctcfweights[tad[0]] + ctcfweights[tad[1]]
			bdycount += 2
		} else if tad[0] > centromereloc[0] { // we've removed centromere so have to add back in its location if TAD started past centromere end
			ctcfload += ctcfweights[tad[0]+cmlen] + ctcfweights[tad[1]+cmlen]
			bdycount += 2
		} else if centromereloc[0] > 0 { // if we have a TAD that gets broken up by the centromere
			ctcfload += ctcfweights[tad[0]] + ctcfweights[centromereloc[0]-1]
			ctcfload += ctcfweights[centromereloc[1]+1] + ctcfweights[tad[1]+cmlen]
			bdycount += 4
		}*/ 
		if math.IsNaN(float64(ctcfload)) {
			fmt.Println("got NaN in CTCF load")
			fmt.Println(tadset)
			fmt.Println(tad)//,centromereloc)
			os.Exit(1)
		}
	}
	avgctcf := float64(ctcfload)/float64(bdycount)
	if math.IsNaN(avgctcf) {
		fmt.Println("avg ctcfload is NaN")
		fmt.Println(ctcfload, bdycount)
		fmt.Println(tadset)
		os.Exit(1)
	}
	return avgctcf
}
