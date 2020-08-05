package ftutil

import( "os"
	"fmt"
	"bufio"
	"strings"
	"strconv"
	)

func WriteTADsToFile(tadset [][]int, res int, outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
        labelline := []string{"TADstart", "TADend"}
        //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,line1+"\n")

	for _,vals := range tadset {
                strvals := make([]string, 2)
                strvals[0] = strconv.Itoa(vals[0]*res)
                strvals[1] = strconv.Itoa(vals[1]*res + res-1)
                newline := strings.Join(strvals, "\t")
                fmt.Fprintf(w,newline+"\n")
        }
        w.Flush()
        f.Close()
        fmt.Println("Wrote output values to", outfile)
}
