package hicutil

import (
	//"fmt"
)

func RemoveCentromere(orighic map[Pair]float64, centromereloc []int, chrlength int) (map[Pair]float64,int) {

	cmlen := centromereloc[1] - centromereloc[0]
	newhic := make(map[Pair]float64)

	for origidx,val := range orighic {
		origrow := origidx.A.(int)
		origcol := origidx.B.(int)
		newidx := []int{origrow,origcol}
		if origrow > centromereloc[1] {
			newidx[0] -= cmlen
		}
		if origcol > centromereloc[1] {
			newidx[1] -= cmlen
		}
		//fmt.Println(centromereloc, origidx, newidx)
		newhic[Pair{newidx[0],newidx[1]}] = val
	}
	return newhic, chrlength-cmlen
}

//func (p*Pair)getPair() { return p }
