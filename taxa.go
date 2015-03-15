package ncbiutils

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
)

var GCTABLES map[string]*GeneticCode

type Taxa struct {
	Id               string       // node id in GenBank taxonomy database
	Name             string       // the unique variant of names
	Parent           string       // parent node id in GenBank taxonomy database
	Rank             string       // rank of this node (superkingdom, kingdom ...)
	EmblCode         string       // locus-name prefix
	Division         string       // division
	InheritedDivFlag string       // 1 if node inherits division from parent
	GeneticCode      *GeneticCode // genetic code
	InheriteGCFlag   string       // 1 if node inherits genetic code from parent
	MitochondrialGC  *GeneticCode // mitochondrial genetic code
	InheriteMGCFlag  string       // 1 if node inherits mitochondrial genetic code from parent
	Comments         string       // free-text comments and citations
}

// read taxonomy from NCBI taxonomy database dmp
// returns a map[id]Taxa
func ReadTaxas(dir string) map[string]Taxa {
	taxaMap := make(map[string]Taxa)

	namesFilePath := filepath.Join(dir, "names.dmp")
	f, err := os.Open(namesFilePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	nameMap := getNames(f)

	gcFilePath := filepath.Join(dir, "gencode.dmp")
	f, err = os.Open(gcFilePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	gencodes := ReadGeneticCodes(f)

	nodesFilePath := filepath.Join(dir, "nodes.dmp")
	f, err = os.Open(nodesFilePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	r := bufio.NewReader(f)
	for {
		l, err := r.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		fields := strings.Split(l, "\t|\t")
		id := fields[0]
		parent := fields[1]
		rank := fields[2]
		embl := fields[3]
		division := fields[4]
		idivflag := fields[5]
		gcid := fields[6]
		igcflag := fields[7]
		mgcid := fields[8]
		imgcflag := fields[9]
		comments := fields[12]

		taxa := Taxa{
			Id:               id,
			Parent:           parent,
			Rank:             rank,
			EmblCode:         embl,
			Division:         division,
			InheritedDivFlag: idivflag,
			GeneticCode:      gencodes[gcid],
			InheriteGCFlag:   igcflag,
			MitochondrialGC:  gencodes[mgcid],
			InheriteMGCFlag:  imgcflag,
			Comments:         comments,
			Name:             nameMap[id],
		}

		taxaMap[id] = taxa
	}

	return taxaMap
}

type GeneticCode struct {
	Id           string          // GenBank genetic code id
	Abbreviation string          // genetic code name abbreviation
	Name         string          // genetic code name
	Table        map[string]byte // translate table for this genetic code
	Starts       []string        // start codon for this genetic code
	FFCodons     map[string]bool // four-fold codons
}

func GeneticCodes() (gcMap map[string]*GeneticCode) {
	buf := bytes.NewBufferString(GCSTRING)
	gcMap = ReadGeneticCodes(buf)
	return
}

func ReadGeneticCodes(f io.Reader) (gcMap map[string]*GeneticCode) {
	gcMap = make(map[string]*GeneticCode)
	rd := bufio.NewReader(f)
	for {
		l, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		fields := strings.Split(l, "\t|\t")
		id := fields[0]
		abb := fields[1]
		name := fields[2]
		cde := fields[3]
		starts := fields[4]
		table, ffs := getTables(cde[:64])
		startCodons := getStartCodons(starts[:64])
		gc := GeneticCode{
			Id:           id,
			Abbreviation: abb,
			Name:         name,
			Table:        table,
			Starts:       startCodons,
			FFCodons:     ffs,
		}
		gcMap[id] = &gc
	}

	return
}

func getTables(s string) (table map[string]byte, ffCodons map[string]bool) {

	table = make(map[string]byte)
	nn := "TCAG"
	l := 0
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				c := string([]byte{nn[i], nn[j], nn[k]})
				table[c] = s[l]
				l++
			}
		}
	}

	ffCodons = make(map[string]bool)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			ff := true
			aa := table[string([]byte{nn[i], nn[j], nn[0]})]
			for k := 0; k < 4; k++ {
				c := string([]byte{nn[i], nn[j], nn[k]})
				if table[c] != aa {
					ff = false
					break
				}
			}
			for k := 0; k < 4; k++ {
				c := string([]byte{nn[i], nn[j], nn[k]})
				ffCodons[c] = ff
			}
		}
	}

	return
}

func getStartCodons(s string) (starts []string) {

	nn := "TCAG"
	l := 0
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				c := string([]byte{nn[i], nn[j], nn[k]})
				if s[l] == 'M' {
					starts = append(starts, c)
				}
				l++
			}
		}
	}
	return
}

// returns a map[id][scientific name]
func getNames(f io.Reader) map[string]string {
	nameMap := make(map[string]string)
	r := bufio.NewReader(f)
	for {
		l, err := r.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		fields := strings.Split(l, "\t|\t")
		id := fields[0]
		name := fields[1]
		class := strings.Split(strings.TrimSpace(fields[3]), "\t|")[0]
		if class == "scientific name" {
			nameMap[id] = name
		}
	}
	return nameMap
}
