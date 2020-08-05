package ftutil

import(
	"math/rand"
	"fmt"
	"os"
)

func MergeRandomTAD(tadset [][]int) [][]int {

	var newtadset [][]int

	tadidx := rand.Intn(len(tadset))
	// need to account for case where tad selected is from start to centromere or centromere to end - can't merge w/ anything
	tad := tadset[tadidx]

	if tadidx != 0 { //&& tad[0] != centromereloc[1]+1 {
		newtad := []int{tadset[tadidx-1][0],tad[1]}
		if newtad[0] > newtad[1] {
			fmt.Println(newtad, tadidx)
			os.Exit(1)
		}
		//fmt.Println(newtad)
		newtadset = make([][]int, len(tadset)-1)
		copy(newtadset[:tadidx-1],tadset[:tadidx-1])
		newtadset[tadidx-1] = newtad
		copy(newtadset[tadidx:],tadset[tadidx+1:])
	} else {
		//fmt.Println(tadidx,tadidx+1,len(currtadset))
		newtad := []int{tad[0],tadset[tadidx+1][1]}
		if newtad[0] > newtad[1] {
			fmt.Println(newtad, tadidx)
			os.Exit(1)
		}
		newtadset = make([][]int, len(tadset)-1)
		copy(newtadset[:tadidx],tadset[:tadidx])
		newtadset[tadidx] = newtad
		copy(newtadset[tadidx+1:],tadset[tadidx+2:])
	}
	return newtadset
}
