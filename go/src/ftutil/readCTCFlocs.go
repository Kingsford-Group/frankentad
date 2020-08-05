package ftutil

import(
	"bufio"
	"os"
	"strings"
	"strconv"
	"sort"
//	"runtime"
//	"fmt"
)


func ReadCTCFlocs(filename string, res int) (map[string][]int) {

// currently specific to column format of CTCFBSDB_all_exp_sites_Sept12_2012.txt

	ctcfmap := make(map[string][]int)
	var start,end,midpt int

	if file,err := os.Open(filename); err == nil {
		defer file.Close()
		scanner := bufio.NewScanner(file)
		for scanner.Scan(){
			line := strings.Fields(scanner.Text())
			chrlabel := line[4]
			start,err = strconv.Atoi(line[5])
			end,err = strconv.Atoi(line[6])
			midpt = (start+end)/2/res
			//fmt.Println(midpt)
			if ctcflist,ok := ctcfmap[chrlabel]; ok {
				// extend list w/ new midpt
				ctcfmap[chrlabel] = append(ctcflist, midpt)
			} else {
				ctcfmap[chrlabel] = []int{midpt}
				//fmt.Println(ctcfmap[chrlabel])
			}

		}
	}

	for key,sites := range ctcfmap {
		sort.Ints(sites)
		ctcfmap[key] = sites
	}

	return ctcfmap
}
