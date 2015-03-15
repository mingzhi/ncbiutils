package ncbiutils

import (
	"github.com/mingzhi/biogo/seq"
	"log"
	"os"
	"path/filepath"
	"strings"
)

// A SeqRecord contains the protein sequence
// and its corresponding nucleotide sequence.
type SeqRecord struct {
	Id     string   // protein id
	Name   string   // gene name
	Prot   []byte   // protein sequence
	Nucl   []byte   // nucelotide sequence
	Genome string   // genome accession
	Code   string   // codon table id
	Loc    Location // genome location
}

type Location struct {
	From, To int
	Strand   string
}

type SeqRecords []SeqRecord

// Implement functions for sort.Interface
func (r SeqRecords) Len() int      { return len(r) }
func (r SeqRecords) Swap(i, j int) { r[i], r[j] = r[j], r[i] }

type ByGenome struct{ SeqRecords }

func (b ByGenome) Less(i, j int) bool {
	return b.SeqRecords[i].Genome < b.SeqRecords[j].Genome
}

func ReadSeqRecords(genome, dir string, gc string) SeqRecords {
	// Clean genome accession.
	acc := strings.Split(genome, ".")[0]
	// Read protein sequences in a .faa file.
	faaPath := filepath.Join(dir, acc+".faa")
	faaSeqs := readFasta(faaPath)
	// Read genome nucleotide sequences in a .fna file.
	fnaPath := filepath.Join(dir, acc+".fna")
	fnaSeqs := readFasta(fnaPath)
	// Read ptt file.
	pttPath := filepath.Join(dir, acc+".ptt")
	pttFile := NewPttFile(pttPath)
	defer pttFile.Close()
	ptts := pttFile.ReadAll()

	g := fnaSeqs[0].Seq // genome sequence.

	// pttMap: {pid: ptt}
	pttMap := make(map[string]Ptt)
	for _, ptt := range ptts {
		pttMap[ptt.PID] = ptt
	}

	records := SeqRecords{}
	for _, prot := range faaSeqs {
		id := prot.Id
		ptt, found := pttMap[id]
		if found {
			if ptt.Loc.From > 0 && ptt.Loc.To > 0 && ptt.Loc.To <= len(g) && ptt.Loc.From <= len(g) {
				var nucl []byte
				if ptt.Loc.To >= ptt.Loc.From {
					nucl = g[ptt.Loc.From-1 : ptt.Loc.To]
				} else {
					nucl = g[ptt.Loc.From-1:]
					nucl = append(nucl, g[:ptt.Loc.To]...)
				}

				if ptt.Loc.Strand == "-" {
					nucl = revcomp(nucl)
				}

				sr := SeqRecord{}
				sr.Id = id
				sr.Name = ptt.Gene
				sr.Genome = acc
				sr.Prot = prot.Seq
				sr.Nucl = nucl
				sr.Loc = ptt.Loc
				sr.Code = gc
				if CheckCodons(sr) {
					records = append(records, sr)
				}
			}
		}
	}

	return records
}

func readFasta(fileName string) []*seq.Sequence {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	fr := seq.NewFastaReader(f)
	fr.DeflineParser = func(x string) string { return strings.Split(x, "|")[1] }
	seqs, err := fr.ReadAll()
	if err != nil {
		panic(err)
	}

	return seqs
}

// reverse and complement a nucleotide sequence
func revcomp(seq []byte) []byte {
	// reverse
	for i, j := 0, len(seq)-1; i < j; i, j = i+1, j-1 {
		seq[i], seq[j] = seq[j], seq[i]
	}
	// complement
	for i := 0; i < len(seq); i++ {
		seq[i] = comp(seq[i])
	}
	return seq
}

// complement a nucleotide sequence (standard version)
func comp(a byte) byte {
	switch a {
	case 'A', 'a':
		return 'T'
	case 'T', 't':
		return 'A'
	case 'G', 'g':
		return 'C'
	case 'C', 'c':
		return 'G'
	}
	return 0
}

// check codons
func CheckCodons(sr SeqRecord) (good bool) {
	if GCTABLES == nil {
		GCTABLES = GeneticCodes()
	}
	na := sr.Nucl
	aa := sr.Prot
	id := sr.Id
	ncodons := len(na) / 3
	gc := GCTABLES[sr.Code]

	if ncodons < len(aa) {
		log.Printf("%s: nucl length %d is smaller than prot sequence %d\n", id, ncodons, len(aa))
		return
	}

	// check start codon
	sc := string(na[0:3])
	ok := false
	for _, c := range gc.Starts {
		if sc == c {
			ok = true
		}
	}
	if !ok {
		log.Printf("%s: the first codon %s is not a start codon in table %s\n", id, sc, gc.Id)
		return
	}

	// check following codons
	for i := 1; i < len(aa); i++ {
		c := string(na[i*3 : (i+1)*3])
		a := aa[i]
		if a != gc.Table[c] {
			log.Printf("%s: the codon at %d (%s) is not corresponed to %s\n", id, i, c, string(a))
			return
		}
	}

	good = true

	return
}
