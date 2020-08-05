package ftutil

import(
	"math/rand"
	"os"
	"fmt"
)

func SplitRandomTAD(tadset [][]int, chrlength int) [][]int {

	var newloc int
        var newtadset [][]int
        var idx int
        var tad []int
	
	validlocs := make([]int,tadset[len(tadset)-1][1]+1)
	for idx,_ := range validlocs {
		//if idx < centromereloc[0] || idx > centromereloc[1] {
		validlocs[idx] = idx
		//}
	}
	for _,tad := range tadset {
		validlocs[tad[0]] = 0
		if len(validlocs) <= tad[0]+1 {
			fmt.Println("something wrong w/ tad lengths")
			fmt.Println(tadset)
			os.Exit(1)
		}
		validlocs[tad[0]+1] = 0
		if tad[0] > 0 {validlocs[tad[0]-1] = 0}
		validlocs[tad[1]] = 0
		if tad[1] < chrlength-1 {validlocs[tad[1]+1] = 0}
		validlocs[tad[1]-1] = 0
	}
	//fmt.Println(validlocs)
	//need to remove nil values from validlocs
	validlocs_clean := removeZeros(validlocs)

	if len(validlocs_clean) == 0 {
		fmt.Println(tadset)
		fmt.Println(validlocs)
		fmt.Println("no TADs big enough to split")
		os.Exit(1)
	}
	newloc = validlocs_clean[rand.Intn(len(validlocs_clean))]
	// still need to find which tad this newloc is in
	var tadidx int
	var splittad []int
	for idx,tad := range tadset {
		if newloc > tad[0] && newloc < tad[1] {
			tadidx = idx
			splittad = tad
			break
		}
	}
	if len(splittad) == 0 {
		fmt.Println(tadset)
		fmt.Println(newloc,tad)
		os.Exit(1)
	}

	newtads := [][]int{[]int{splittad[0],newloc}, []int{newloc+1,splittad[1]}}
	if newtads[0][0] > newtads[0][1] || newtads[1][0] > newtads[1][1] {
		fmt.Println(newloc, tad, idx)
		os.Exit(1)
	}
	newtadset = make([][]int, len(tadset)+1)
	copy(newtadset[:tadidx],tadset[:tadidx])
	newtadset[tadidx] = newtads[0]
	newtadset[tadidx+1] = newtads[1]
	copy(newtadset[tadidx+2:],tadset[tadidx+1:])

	return newtadset
}

func removeZeros(x []int) []int {

        var xclean []int
        for _,val := range x {
                if val != 0 {
                        xclean = append(xclean,val)
                }
        }
        return xclean
}
