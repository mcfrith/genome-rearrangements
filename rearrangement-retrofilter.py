#! /usr/bin/env python
# Copyright 2016 Martin C. Frith

import bisect, optparse, signal, sys

def myOpen(fileName):  # faster than fileinput
    if fileName == '-': return sys.stdin
    else:               return open(fileName)

def retroseqFromLine(line):
    fields = line.split()
    gene = fields[1]
    chrom = fields[6]
    strand = fields[9]
    if strand == "+":
        beg = int(fields[7])
        end = beg + int(fields[8])
    else:
        end = int(fields[10]) - int(fields[7])
        beg = end - int(fields[8])
    return chrom, beg, end, strand, gene

def retroseqText(r):
    return r[0] + ":" + str(r[1]) + "-" + str(r[2]) + "|" + r[4] + r[3]

def retroseqLength(r):
    return r[2] - r[1]

def edgeFromText(t):
    head, tail = t.split(":")
    genome, chrom = head.split(".")
    tail = tail.rstrip("$")
    pos = int(tail[:-1])
    bracket = tail[-1]
    return chrom, pos, bracket

def overlapJaccardIndex(beg1, end1, beg2, end2):
    assert beg1 <= end1
    assert beg2 <= end2
    num = min(end1, end2) - max(beg1, beg2)
    if num < 0: return 0.0
    den = max(end1, end2) - min(beg1, beg2)
    return 1.0 * num / den

def bestRetroOverlap(queryChrom, queryBeg, queryEnd, retros, maxRetroLength):
    bestOverlap = 0.0
    bestRetro = None
    fakeRetro = queryChrom, queryBeg, queryEnd, "", ""
    k = bisect.bisect(retros, fakeRetro)
    ka = k
    while ka < len(retros):
        r = retros[ka]
        chrom, beg, end, strand, gene = r
        if chrom > queryChrom: break
        if beg >= queryEnd: break
        overlap = overlapJaccardIndex(beg, end, queryBeg, queryEnd)
        if overlap > bestOverlap: bestOverlap, bestRetro = overlap, r
        ka += 1
    kb = k
    while kb > 0:
        kb -= 1
        r = retros[kb]
        chrom, beg, end, strand, gene = r
        if chrom < queryChrom: break
        if beg + maxRetroLength <= queryBeg: break
        overlap = overlapJaccardIndex(beg, end, queryBeg, queryEnd)
        if overlap > bestOverlap: bestOverlap, bestRetro = overlap, r
    return bestOverlap, bestRetro

def bestOverlaps(edges, retros, maxRetroLength, opts):
    overlaps = {}
    for j in range(len(edges)):
        for i in range(j):
            iChrom, iPos, iBracket = edges[i]
            jChrom, jPos, jBracket = edges[j]
            if iChrom != jChrom: continue
            if iPos == jPos: continue
            if iBracket == "[": continue
            if jBracket == "]": continue
            if i + 1 == j: continue  # xxx
            overlap, retro = bestRetroOverlap(iChrom, iPos, jPos, retros,
                                              maxRetroLength)
            if overlap < opts.min_overlap: continue
            if retro not in overlaps or overlaps[retro] < overlap:
                overlaps[retro] = overlap
    return overlaps

def rearrangementRetrofilter(opts, args):
    retroFile = myOpen(args[0])
    retroLines = (i for i in retroFile if i[0] != "#")
    retros = sorted(map(retroseqFromLine, retroLines))
    maxRetroLength = max(map(retroseqLength, retros))

    for line in myOpen(args[1]):
        fields = line.split()
        edgeFields = (i for i in fields if opts.genome in i)
        edges = sorted(map(edgeFromText, edgeFields))
        overlaps = bestOverlaps(edges, retros, maxRetroLength, opts)
        if opts.show:
            for k, v in overlaps.iteritems():
                print "%#.3g" % v, retroseqText(k), line,
        else:
            if not overlaps:
                print line,

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message
    usage = "%prog retros.tab rearrangements"
    description = "Filter out rearrangements that match retrosequences."
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-m", "--min-overlap", metavar="FRAC",
                  type="float", default=0.9,
                  help="check for overlaps with Jaccard index >= FRAC (default: %default)")
    op.add_option("-s", "--show", action="store_true",
                  help="show the overlaps (default: write rearrangements without overlaps)")
    op.add_option("-g", "--genome", metavar="NAME", default="hg19",
                  help="genome name (default: %default)")
    opts, args = op.parse_args()
    if len(args) != 2: op.error("I need 2 file names")
    rearrangementRetrofilter(opts, args)
