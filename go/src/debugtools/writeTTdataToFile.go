package debugtools

import("os"
	"fmt"
	"bufio"
	"strings"
	"strconv"
	"hicutil"
)

func WriteTTdataToFile(tad []int, hicmap map[hicutil.Pair]float64, res int, backdist []float64, slope float64, lgint float64, filename string) {
	//also need linear regression parameters
	
	f,err := os.Create(filename)
	if err != nil { panic(err) }

	w := bufio.NewWriter(f)
	//write linear reg params first
	lgdata := make([]string,2)
	lgdata[0] = fmt.Sprintf("%f",slope)
	lgdata[1] = fmt.Sprintf("%f",lgint)
	lgstring := strings.Join(lgdata,"\t")
	fmt.Fprintf(w,lgstring+"\n")

	//then write tad values
	for a := tad[0]; a <= tad[1]; a++ {
		for b := a; b <= tad[1]; b++ {
			data := make([]string, 2)
			data[0] = strconv.Itoa((b-a)*res)
			data[1] = fmt.Sprintf("%f", hicmap[hicutil.Pair{a,b}]/backdist[b-a])
			dataline := strings.Join(data,"\t")
			fmt.Fprintf(w,dataline+"\n")
		}
	}
	w.Flush()
	//w.Close()
	fmt.Println("wrote TADtree data for",tad,"to",filename)
}
