# genome-rearrangements

These scripts find rearrangements, such as inversions and
translocations, in genome alignments.  They do not find pure deletions
or insertions/duplications.

## genome-rearrangements.py

This script identifies rearrangements in genome alignments.  You need
to give it a *one-to-one* alignment of *two* genomes, in [MAF][]
format.  [Here](https://zenodo.org/record/17436) are some suitable
alignments.  Optionally (but recommended), you can supply the
locations of unsequenced gaps in each genome, as [gap.txt][] files.

[MAF]: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
[gap.txt]: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

Typical usage:

    genome-rearrangements.py -1 top-genome/gap.txt -2 bot-genome/gap.txt in.maf > out.txt

The output looks like this:

    hg19.chr2:3510663] hg19.chr2:3510663[ panTro4.chr2A:5062364] ...
    hg19.chr4:7611136] hg19.chr4:7611136[ panTro4.chr4:7773318] gap889,105 ...

Each line is one rearrangement, listing the alignment endpoints
participating in that rearrangement.  Each endpoint has a genome name
(e.g. `hg19`), chromosome name (e.g. `chr2`), coordinate, and either
`[` (mnemonic for left edge of an alignment) or `]` (mnemonic for
right edge of an alignment).  Coordinates are zero-based: if an
alignment covers the first 100 bases of a sequence, its left
coordinate is 0 and its right coordinate is 100.

The output also indicates unsequenced gaps: e.g. `gap889,105`
indicates two unsequenced gaps, of length 889 and 105, between the
flanking endpoints.

## supported-rearrangements.py

This script reads a set of "query" rearrangements (e.g. human-chimp),
and writes only those whose endpoints match "reference" rearrangements
(e.g. human-orangutan).  Typical usage:

    supported-rearrangements.py ref-rearrangements query-rearrangements > supported-query-rearrangements

## Filtering spliced retrosequences

The preceding scripts are not supposed to find pure
insertions/duplications, such as retrosequence insertions.  But
occasionally they wrongly show a translocation that is actually an
insertion of a spliced retrosequence.  We can deal with this in two
steps:

1. Find spliced retrosequences.  The following example uses
   [LAST](http://last.cbrc.jp/) to find high-similarity alignments
   between the human genome and a set of human mRNA sequences.  Then,
   last-spliced-retroseqs.py gets the subset of alignments that cross
   splice junctions but lack introns.  The files refMrna.fa (mRNA
   sequences) and refSeqAli.txt (splice junctions) can be obtained
   from [UCSC](http://genome.ucsc.edu/).

        lastdb -cR01 -uNEAR my-rna-db refMrna.fa

        lastal -p human-chimp.v2.mat -e3000 -C1 -m50 -f0 my-rna-db human-genome.fa |
        last-spliced-retroseqs.py refSeqAli.txt - > retros.tab

2. Get rearrangements that do not match these retrosequences.

        rearrangement-retrofilter.py retros.tab rearrangements.txt > good-rearrangements.txt
