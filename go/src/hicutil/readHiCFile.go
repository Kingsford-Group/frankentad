package hicutil

import (
	"bufio"
//	"log"
	"strings"
	"strconv"
	"os"
	"math"
	"sort"
	"fmt"
)

// define a pair (tuple) type for keys to Hi-C maps
type Pair struct {
     A,B interface{}
}

func ReadHiCFile(filename string, res int, normvals []float64) (map[Pair]float64,int,[]int) {

     hicmap := make(map[Pair]float64)
     nonzerolocmap := make(map[int]struct{})
     var row,col int
     var val float64
     var norm bool
     var chrlength int
     centromereloc := []int{0,0}

     if len(normvals) > 0 {
	     norm = true
     } else {
	     norm = false
     }

     if file, err := os.Open(filename); err == nil {
     	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
	    vals := strings.Fields(scanner.Text())
	    row, err = strconv.Atoi(vals[0])
	    row = row/res
	    col, err = strconv.Atoi(vals[1])
	    col = col/res
	    
	    if row > chrlength {
	       chrlength = row
	    }
	    if col > chrlength {
	       chrlength = col
	    }

	    val, err = strconv.ParseFloat(vals[2],64)
	    if norm {
	       val = float64(val) / (float64(normvals[row])*float64(normvals[col]))
	    }
	hicmap[Pair{row,col}] = math.Log(val)
	//	hicmap[Pair{row,col}] = val

	nonzerolocmap[row] = struct{}{}
	nonzerolocmap[col] = struct{}{}
		//	hicmap[Pair{row,col}] = val
	}

	// find the centromere
	nonzerolocs := make([]int, 0, len(nonzerolocmap))
	for l := range nonzerolocmap {
		nonzerolocs = append(nonzerolocs, l)
	}
	sort.Ints(nonzerolocs)
	//fmt.Println(nonzerolocs)

	locdiffs := make([]int,len(nonzerolocs)-1)
	for i,x := range nonzerolocs[:len(nonzerolocs)-1] {
		locdiffs[i] = nonzerolocs[i+1]-x
	}

	maxdiff := MaxIntIdx(locdiffs,nonzerolocs[0])
	if maxdiff != 0 && nonzerolocs[maxdiff+1] < chrlength {
		centromereloc[0] = nonzerolocs[maxdiff]
		centromereloc[1] = nonzerolocs[maxdiff+1]
	} else {
		centromereloc[1] = nonzerolocs[0]
	}
	/*if err = scanner.Err(); err != nil {
	   log.Fatal(err)
	} else {
	  log.Fatal(err)
	}*/
     } else {
	     fmt.Println("could not locate Hi-C file, please check filepath")
	     os.Exit(1)
     }
	return hicmap, chrlength, centromereloc
}

func MaxIntIdx(x []int, maxval int) int {
	
	var maxidx int

	for idx, v := range x {
		if v > maxval {
			maxidx = idx
			maxval = v
		}
	}

	return maxidx
}
