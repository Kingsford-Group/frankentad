package ftutil

import( "os"
	"fmt"
	"bufio"
        "github.com/c-bata/goptuna"
	"strings"
	"strconv"
	)

func WriteBOInfoToFile(trials []goptuna.FrozenTrial, outfile string) {

        //write values to file
        f,err := os.Create(outfile)
        if err != nil {
		panic(err)
		}
        //defer f.Close()

        w := bufio.NewWriter(f)
        labelline := []string{"trialnumber", "value","alpha_is","gamma_a","lambda_is","lambda_di"}
        //fmt.Println(strings.Join(labelline, "\t"))
	line1 := strings.Join(labelline, "\t")
        //fmt.Println(line1)
	fmt.Fprintf(w,line1+"\n")

	for _,t := range trials {
                strvals := make([]string, 6)
                strvals[0] = strconv.Itoa(t.ID)
                strvals[1] = fmt.Sprintf("%f",t.Value)
		strvals[2] = fmt.Sprintf("%f", t.Params["alpha_is"])
                strvals[3] = fmt.Sprintf("%f", t.Params["gamma_a"])
		strvals[4] = fmt.Sprintf("%f", t.Params["lambda_is"])
		strvals[5] = fmt.Sprintf("%f", t.Params["lambda_di"])
		newline := strings.Join(strvals, "\t")
		//newline := strconv.Itoa(t.ID)
                fmt.Fprintf(w,newline+"\n")
        }
        w.Flush()
        f.Close()
        fmt.Println("Wrote optimization info to", outfile)
}
