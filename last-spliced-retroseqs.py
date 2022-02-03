#! /usr/bin/env python
# Copyright 2016 Martin C. Frith

# Read gene annotations in PSL format, and genome-to-RNA alignments in
# LAST tabular format.  Write the alignments that seem to indicate
# spliced retrosequences.

from __future__ import print_function

import bisect, operator, optparse, signal, sys

def myOpen(fileName):  # faster than fileinput
    if fileName == '-': return sys.stdin
    else:               return open(fileName)

def isExtraPslField(fields):
    return fields[9] in "+-"

def pslSplit(field):
    n = field.rstrip(",").split(",")
    return [int(i) for i in n]

def exonJunctionsFromPslBlocks(blockSizes, qStarts, tStarts):
    minIntronLength = 50
    e = len(blockSizes)
    for i in range(e):
        j = i + 1
        if j < e:
            qGap = qStarts[j] - qStarts[i] - blockSizes[i]
            assert qGap >= 0
            tGap = tStarts[j] - tStarts[i] - blockSizes[i]
            assert tGap >= 0
            if qGap == 0 and tGap >= minIntronLength: yield qStarts[j]

def exonsFromPsl(line):
    fields = line.split()
    if isExtraPslField(fields): fields.pop(0)
    strand = fields[8]
    qName = fields[9]
    qSize = fields[10]
    blockSizes = pslSplit(fields[18])
    qStarts = pslSplit(fields[19])
    tStarts = pslSplit(fields[20])
    ej = list(exonJunctionsFromPslBlocks(blockSizes, qStarts, tStarts))
    if strand == "-":
        s = int(qSize)
        ej = [s - i for i in reversed(ej)]
    return qName, qSize, ej

def readGenes(lines):
    genes = {}
    for line in lines:
        qName, qSize, exonJunctions = exonsFromPsl(line)
        if not exonJunctions: continue
        default = qSize, []
        qSize0, ejList = genes.setdefault(qName, default)
        if qSize != qSize0: raise Exception("inconsistent length for " + qName)
        ejList.append(exonJunctions)
    return genes

def isMultiExon(alignmentBeg, alignmentEnd, exonJunctions):
    # does the alignment extend at least X bp either side of any junction?
    minBasesPastJunction = 50
    i = bisect.bisect_left(exonJunctions, alignmentBeg + minBasesPastJunction)
    j = bisect.bisect_right(exonJunctions, alignmentEnd - minBasesPastJunction)
    return i < j

def isBigInsertion(lastAlignmentBlocks):
    maxInsertLength = 20
    for i in lastAlignmentBlocks.split(","):
        if ":" in i:
            delete, insert = i.split(":")
            if int(insert) > maxInsertLength: return True
    return False

def lastSplicedRetroseqs(opts, args):
    genesFile = myOpen(args[0])
    genes = readGenes(genesFile)
    lastFile = myOpen(args[1])
    for line in lastFile:
        if not line[0] == "#":
            fields = line.split()
            name = fields[1]
            if name not in genes: continue
            qSize, ejList = genes[name]
            beg = int(fields[2])
            end = beg + int(fields[3])
            if not any(isMultiExon(beg, end, i) for i in ejList): continue
            seqlen = fields[5]
            if seqlen != qSize: continue  # can happen, sadly
            blocks = fields[11]
            if isBigInsertion(blocks): continue
        print(line, end="")

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message
    usage = "%prog refSeqAli.txt alignments.tab > retros.tab"
    description = "Get spliced retrosequences, from genome-RNA alignments in LAST tabular format."
    op = optparse.OptionParser(usage=usage, description=description)
    opts, args = op.parse_args()
    if len(args) != 2: op.error("I need 2 file names")
    lastSplicedRetroseqs(opts, args)
