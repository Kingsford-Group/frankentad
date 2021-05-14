package ftutil

import( "os"
	"fmt"
	"bufio"
	"strings"
	"strconv"
	)

func WriteTADsToFile(tadset [][]int, res int, alpha float64, gamma float64, lambda_a float64, lambda_is float64, lambda_di float64, outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
	paramstring := fmt.Sprintf("alpha=%f gamma=%f lambda_a=%f lambda_is=%f lambda_di=%f\n",alpha,gamma,lambda_a,lambda_is,lambda_di)
        labelline := []string{"TADstart", "TADend"}
        //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,paramstring)
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
