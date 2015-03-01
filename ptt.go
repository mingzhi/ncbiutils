package ncbiutils

import (
	"bufio"
	"io"
	"os"
	"strconv"
	"strings"
)

type Ptt struct {
	Loc         Location
	Length      int
	PID         string
	Gene        string
	SynonymCode string
	COG         string
	Product     string
}

type PttFile struct {
	r io.ReadCloser
}

func NewPttFile(fileName string) *PttFile {
	p := PttFile{}
	var err error
	p.r, err = os.Open(fileName)
	if err != nil {
		panic(err)
	}

	return &p
}

func (p *PttFile) Close() {
	p.r.Close()
}

func (p *PttFile) ReadAll() (ptts []Ptt) {
	rd := bufio.NewReader(p.r)
	skipLines := 3
	for i := 0; i < skipLines; i++ {
		rd.ReadString('\n')
	}

	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			} else {
				break
			}
		}
		line = strings.TrimSpace(line)
		fields := strings.Split(line, "\t")
		locTerms := strings.Split(fields[0], "..")
		start, _ := strconv.Atoi(locTerms[0])
		end, _ := strconv.Atoi(locTerms[1])
		strand := fields[2]
		ptt := Ptt{}
		ptt.Loc = Location{From: start, To: end, Strand: strand}
		ptt.PID, ptt.Gene, ptt.SynonymCode, ptt.COG, ptt.Product =
			fields[3], fields[4], fields[5], fields[6], fields[7]
		ptts = append(ptts, ptt)
	}

	return
}
