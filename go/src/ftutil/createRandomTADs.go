package ftutil

import(
	"math/rand"
	"os"
	"fmt"
)

func CreateRandomTADs(chrlength int) ([][]int) {

        //runtime.GOMAXPROCS(1)

        //fmt.Println("chr length:",chrlength)

        tadend := 0
	var tadset [][]int

        for tadend < chrlength-1 {
		newtadlen := rand.Intn(14)+2
                if tadend + 1 + newtadlen > chrlength-1 {newtadlen = chrlength - tadend - 2}
                newtad := make([]int,2)
		if len(tadset) == 0 {
                        newtad[0] = 0
                        newtad[1] = newtadlen
                } else if chrlength-1 - tadend < 20 {
			        //fmt.Println("got here")
			newtad[0] = tadend+1
                        newtad[1] = chrlength - 1
                        //fmt.Println(newtad)
			} else {
                        newtad[0] = tadend+1
                        newtad[1] = tadend+1+newtadlen
                }
                tadset = append(tadset,newtad)
                if tadset[0][1] == -1 {
                        fmt.Println(tadset, newtad)
                        os.Exit(1)
                }
                tadend = newtad[1]
                //fmt.Println(newtad)
                for idx,tad := range tadset {
                        if tad[0] > tad[1] {
                                fmt.Println("mistake in randomly generating TADs")
                                fmt.Println(tadset)
                                fmt.Println(newtad,chrlength)
                                fmt.Println(idx,tad)
                                os.Exit(1)
                        } else if tadend >= chrlength {
                                fmt.Println("TAD too far")
                                fmt.Println(tadset)
                                fmt.Println(newtad, chrlength)
                                os.Exit(1)
                        }
                }
        }
        //fmt.Println("randomly generated TADs:",tadset)
        return tadset
}
