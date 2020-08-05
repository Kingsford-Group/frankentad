package ftutil

import( "os"
	"fmt"
	"bufio"
//        "github.com/c-bata/goptuna"
	"strings"
//	"strconv"
	)

func WriteParamsToFile(paramvals []float64, ctcfscores []float64, label string, outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
        labelline := []string{label,"ctcfscore"}
        //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,line1+"\n")

	for i,p := range paramvals {
                strvals := make([]string, 2)
                strvals[0] = fmt.Sprintf("%f",p)
                strvals[1] = fmt.Sprintf("%f",ctcfscores[i])
		newline := strings.Join(strvals, "\t")
		//newline := strconv.Itoa(t.ID)
                fmt.Fprintf(w,newline+"\n")
        }
        w.Flush()
        f.Close()
        fmt.Println("Wrote parameter values with score to", outfile)
}
