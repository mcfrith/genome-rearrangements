#! /usr/bin/env python
# Copyright 2016 Martin C. Frith

# Warning: the output seems to depend on the order of the input.  In
# fact, the identified rearrangements do not vary with input order,
# but the way they are written does.

from __future__ import print_function

import bisect, operator, optparse, signal, sys, warnings

def myOpen(fileName):  # faster than fileinput
    if fileName == '-': return sys.stdin
    else:               return open(fileName)

def chromosomeNameOnly(chromosomeName):  # remove genome name, if any
    return chromosomeName.split(".")[-1]

def chromosomeNameBase(chromosomeName):
    return chromosomeName.split("_")[0]

def isOrderedChromosome(name):  # OK for hg19, panTro4, ponAbe2, mm10, canFam3
    return "_" not in name and "U" not in name

def isKnownChromosome(sequenceName):  # xxx a bit fragile
    if "U" in sequenceName: return False
    return "chr" in sequenceName

def isCompatibleSequenceNames(name1, name2):
    if not isKnownChromosome(name1): return True
    if not isKnownChromosome(name2): return True
    return chromosomeNameBase(name1) == chromosomeNameBase(name2)

def isExtraFirstGapField(fields):
    return fields[4].isdigit()

def readGaps(lines):
    for line in lines:
        fields = line.split()
        if isExtraFirstGapField(fields): fields = fields[1:]
        if fields[4] not in "NU": continue
        seqName = fields[0]
        end = int(fields[2])
        beg = end - int(fields[5])  # zero-based coordinate
        if isOrderedChromosome(seqName) or fields[7] == "yes":
            isOrderedGap = True
        else:
            isOrderedGap = False
        yield seqName, beg, end, isOrderedGap

def getGaps(fileName):
    if fileName: gaps = list(readGaps(myOpen(fileName)))
    else:        gaps = []
    gaps.sort()
    unorderedGaps = [i for i in gaps if not i[3]]
    return gaps, unorderedGaps

def gapLength(gap):
    return gap[2] - gap[1]

def isReliableMaf(aLine, maxMismap):
    for i in aLine.split():
        x = i.split("=")
        if len(x) > 1 and x[0] == "mismap":
            return float(x[1]) <= maxMismap
    return True

# This reads pair-wise local alignments, and returns 4 edges per
# alignment (2 ends x 2 sequences).  Each edge has these fields:
# 0 serial number
# 1 genome number
# 2 sequence name
# 3 coordinate
# 4 whether it's a start or an end
# 5 serial number of the edge that is aligned to this one
def alignmentEdgesFromMaf(lines, maxMismap):
    i = 0  # serial number for alignment ends
    for line in lines:
        if line[0] == "a":
            genomeNumber = 0
            isWanted = isReliableMaf(line, maxMismap)
        if line[0] == "s" and isWanted:
            genomeNumber += 1
            s, seqName, beg, span, strand, seqLen, aln = line.split(None, 6)
            if strand == "+":
                beg = int(beg)
                end = beg + int(span)
                begType = "start"
                endType = "end"
            else:
                beg = int(seqLen) - int(beg)
                end = beg - int(span)
                begType = "end"
                endType = "start"
            if genomeNumber == 1:
                yield [i+0, 1, seqName, beg, begType, i+1]
                yield [i+2, 1, seqName, end, endType, i+3]
            if genomeNumber == 2:
                yield [i+1, 2, seqName, beg, begType, i+0]
                yield [i+3, 2, seqName, end, endType, i+2]
                i += 4

def gapsBetween(edgeA, edgeB, bothGenomeGaps):
    genome1gaps, genome2gaps = bothGenomeGaps
    assert edgeA[1:3] == edgeB[1:3]  # same genome and sequence
    if edgeA[1] == 1: gaps = genome1gaps
    else:             gaps = genome2gaps
    seqName = chromosomeNameOnly(edgeA[2])
    coordinateA = edgeA[3]
    coordinateB = edgeB[3]
    beg = min(coordinateA, coordinateB)
    end = max(coordinateA, coordinateB)
    fakeBegGap = seqName, beg, beg, False
    fakeEndGap = seqName, end, end, False
    i = bisect.bisect(gaps, fakeBegGap)
    gapList = []
    while i < len(gaps) and gaps[i] < fakeEndGap:
        if gaps[i][2] > end:
            warnings.warn("a gap overlaps an alignment")
        gapList.append(gaps[i])
        i += 1
    return gapList

def isFacing(edgeA, edgeB, maxDistance, unorderedGaps):
    if edgeA[1:3] != edgeB[1:3]: return False
    assert edgeA[3] <= edgeB[3]
    if edgeB[3] - edgeA[3] > maxDistance: return False
    if gapsBetween(edgeA, edgeB, unorderedGaps): return False
    return True

def appendFacingSerialNumbers(edges, maxDistance, unorderedGaps):
    edges.sort(key=operator.itemgetter(1, 2, 3, 4))
    for i, x in enumerate(edges):
        facing = -1
        if i % 2:
            assert x[4] == "end"
            if i+1 < len(edges):
                y = edges[i+1]
                if isFacing(x, y, maxDistance, unorderedGaps):
                    facing = y[0]
        else:
            assert x[4] == "start"
            if i > 0:
                y = edges[i-1]
                if isFacing(y, x, maxDistance, unorderedGaps):
                    facing = y[0]
        x.append(facing)

def isClosedLoop(linkedEdges):
    return linkedEdges[0][6] >= 0

def linkedEdgesAndGaps(linkedEdges, bothGenomeGaps):
    for i, x in enumerate(linkedEdges):
        yield x
        j = i + 1
        if j == len(linkedEdges): continue
        y = linkedEdges[j]
        if x[1] != y[1]: continue
        gapList = gapsBetween(x, y, bothGenomeGaps)
        if gapList: yield gapList

def edgeOrGapsToString(e):
    try:
        prefix = e[2]
        if e[4] == "start": suffix = "["
        else:               suffix = "]"
        return prefix + ":" + str(e[3]) + suffix
    except:
        return "gap" + ",".join(str(gapLength(i)) for i in e)

def printMe(edgesAndGaps):
    out = map(edgeOrGapsToString, edgesAndGaps)
    print(*out)

def isCompatibleEdges(x, y, unorderedGaps):
    xSeqName = x[2]
    ySeqName = y[2]
    if xSeqName == ySeqName:
        return gapsBetween(x, y, unorderedGaps)
    else:
        return isCompatibleSequenceNames(xSeqName, ySeqName)

def isEndJoin(linkedEdges, unorderedGaps):
    if isClosedLoop(linkedEdges): return False
    if len(linkedEdges) != 4: return False
    return isCompatibleEdges(linkedEdges[0], linkedEdges[3], unorderedGaps)

def isGapFill(linkedEdges, gaps, unorderedGaps):
    if isClosedLoop(linkedEdges): return False
    if len(linkedEdges) != 8: return False
    e0, e1, e2, e3, e4, e5, e6, e7 = linkedEdges
    #if e0[2] != e7[2]: return False
    if not isCompatibleSequenceNames(e0[2], e7[2]): return False
    # xxx what if e0[2] or e7[2] is not a fragment/"random" chromosome?
    if e1[2] != e6[2]: return False
    if e1[4] == e6[4]: return False
    if e1[4] > e6[4] and e1[3] >= e6[3]: return False
    if e1[4] < e6[4] and e1[3] <= e6[3]: return False
    if not gapsBetween(e3, e4, gaps): return False
    if not isCompatibleEdges(e0, e3, unorderedGaps): return False
    if not isCompatibleEdges(e7, e3, unorderedGaps): return False
    return True

def isRearranged(linkedEdges, gaps, unorderedGaps):
    n = len(linkedEdges)
    assert not n % 2
    if n == 4 and isClosedLoop(linkedEdges): return False
    if n < 4: return False
    if isEndJoin(linkedEdges, unorderedGaps): return False
    if isGapFill(linkedEdges, gaps, unorderedGaps): return False
    return True

def getLinkedEdges(edges, gaps, unorderedGaps):
    edges.sort()
    for x in edges:
        linkedEdges = []
        y = x
        while 1:
            if y[0] < 0: break
            linkedEdges.append(y)
            y[0] = -1
            facing = y[6]
            if facing < 0: break
            y = edges[facing]
            linkedEdges.append(y)
            y[0] = -1
            aligned = y[5]
            y = edges[aligned]
        y = x
        while 1:
            aligned = y[5]
            y = edges[aligned]
            if y[0] < 0: break
            linkedEdges.insert(0, y)
            y[0] = -1
            facing = y[6]
            if facing < 0: break
            y = edges[facing]
            linkedEdges.insert(0, y)
            y[0] = -1
        if isRearranged(linkedEdges, gaps, unorderedGaps):
            yield linkedEdges

def sortKey(linkedEdges):
    return linkedEdges[0][1:5]

def genomeRearrangements(opts, args):
    gaps1, unorderedGaps1 = getGaps(opts.gap1)
    gaps2, unorderedGaps2 = getGaps(opts.gap2)
    gaps = gaps1, gaps2
    unorderedGaps = unorderedGaps1, unorderedGaps2
    edges = list(alignmentEdgesFromMaf(myOpen(args[0]), opts.mismap))
    appendFacingSerialNumbers(edges, opts.distance, unorderedGaps)
    e = getLinkedEdges(edges, gaps, unorderedGaps)
    s = sorted(e, key=sortKey)
    for i in s:
        j = linkedEdgesAndGaps(i, gaps)
        printMe(j)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message
    usage = "%prog [options] pairwise-alignment-maf-file"
    description = "Find rearrangements in a one-to-one alignment of 2 genomes."
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-1", "--gap1", metavar="FILE",
                  help="read genome1 assembly gaps from agp or gap file")
    op.add_option("-2", "--gap2", metavar="FILE",
                  help="read genome2 assembly gaps from agp or gap file")
    op.add_option("-m", "--mismap", metavar="PROB", type="float", default=1e-5,
                  help="omit alignments with mismap probability > PROB (default: %default)")
    opts, args = op.parse_args()
    if len(args) != 1: op.error("I need 1 file name")
    #opts.distance = 1000  # xxx ???
    opts.distance = 1e9
    genomeRearrangements(opts, args)
