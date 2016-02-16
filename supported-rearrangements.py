#! /usr/bin/env python
# Copyright 2016 Martin C. Frith

# Read 2 files of genome rearrangements: write the 2nd file,
# indicating which breaks are supported by the 1st file.

import bisect, collections, optparse, signal

def myOpen(fileName):  # faster than fileinput
    if fileName == '-': return sys.stdin
    else:               return open(fileName)

def isEdgeString(s):
    return s.count(":") == 1  # xxx ???

def edgeFromString(s):
    chromosomeName, tail = s.split(":")
    coordinate = int(tail[:-1])
    isStart = tail[-1] == "["
    return chromosomeName, coordinate, isStart

def edgesFromLines(lines):
    edges = collections.defaultdict(list)
    for line in lines:
        w = line.split()
        for i in w:
            if isEdgeString(i):
                chromosomeName, coordinate, isStart = edgeFromString(i)
                v = [coordinate, isStart, -1]
                edges[chromosomeName].append(v)
    for i in edges.itervalues():
        i.sort()
    return edges

def findEquivalentEdges(edges1, edges2, maxDistance):
    minDistance = maxDistance + 1
    i1best = -1
    i2best = -1
    beg2 = 0
    for i1, e1 in enumerate(edges1):
        coordinate1, isStart1, junk1 = e1
        i2 = beg2
        while i2 < len(edges2):
            e2 = edges2[i2]
            coordinate2, isStart2, junk2 = e2
            if isStart1 == isStart2:
                distance = abs(coordinate1 - coordinate2)
                if distance < minDistance:
                    minDistance = distance
                    i1best = i1
                    i2best = i2
                if coordinate2 <= coordinate1: beg2 = i2 + 1
                if coordinate2 >= coordinate1: break
            i2 += 1
    if minDistance <= maxDistance:
        edges1[i1best][2] = i2best
        edges2[i2best][2] = i1best
        findEquivalentEdges(edges1[:i1best], edges2[:i2best], maxDistance)
        i1next = i1best + 1
        i2next = i2best + 1
        findEquivalentEdges(edges1[i1next:], edges2[i2next:], maxDistance)

def findAllEquivalentEdges(edgeDict1, edgeDict2, maxDistance):
    for k in edgeDict1:
        if k in edgeDict2:
            findEquivalentEdges(edgeDict1[k], edgeDict2[k], maxDistance)

def isSupported(edgeString, edgeDict):
    chromosomeName, coordinate, isStart = edgeFromString(edgeString)
    e = edgeDict[chromosomeName]
    v = [coordinate, isStart, -2]
    k = bisect.bisect(e, v)
    assert k < len(e)
    assert e[k][0] == coordinate
    assert e[k][1] == isStart
    return e[k][2] >= 0

def annotatedFields(fields, edgeDict):
    for i in fields:
        if isEdgeString(i) and isSupported(i, edgeDict):
            i += "$"
        yield i

def supportedRearrangements(opts, args):
    refInput = myOpen(args[0])
    refEdges = edgesFromLines(refInput)
    queryLines = list(myOpen(args[1]))
    queryEdges = edgesFromLines(queryLines)
    findAllEquivalentEdges(refEdges, queryEdges, opts.distance)
    for line in queryLines:
        fields = line.split()
        if opts.all:
            w = annotatedFields(fields, queryEdges)
            print " ".join(w)
        else:
            edgeStrings = filter(isEdgeString, fields)
            n = len(edgeStrings)
            assert n % 2 == 0  # xxx
            s = sum(1 for i in edgeStrings if isSupported(i, queryEdges))
            assert s <= n // 2
            if s == n // 2: print line,

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message
    usage = "%prog ref-rearrangements query-rearrangements"
    description = 'Write "query" rearrangements that are supported by "reference" rearrangements.'
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-a", "--all", action="store_true",
                  help="write all query rearrangements, with '$' at supported edges")
    op.add_option("-d", "--distance", metavar="BASES",
                  type="int", default=1000, help=
                  "maximum distance to supporting edge (default: %default)")
    opts, args = op.parse_args()
    if len(args) != 2: op.error("I need 2 file names")
    supportedRearrangements(opts, args)
