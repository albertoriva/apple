#!/usr/bin/env python

import sys

class CXwriter():
    adjfile = None
    nodesdb = {}
    nedges = {}                 # Number of edges for each node (degree)
    nodecnt = 0                 # Counter for @id fields in nodes
    edgecnt = 0                 # Counter for @id fields in edges
    visibleNodes = {}           # Used for filtering
    mindegree = None            # If set, show only nodes with this degree or higher
    maxnodes = None             # If set, show this number of genes in order of decreasing degree
    nodesWritten = 0            # Number of nodes actually written
    edgesWritten = 0            # Number of edges actually written

    def __init__(self, adjfile):
        self.adjfile = adjfile
        self.nodesdb = {}
        self.nedges = {}
        self.nodecnt = 0
        self.edgecnt = 0
        self.visibleNodes = {}
        self.nodesWritten = 0
        self.edgesWritten = 0

    def writePreamble(self, stream):
        stream.write("""
[ {
  "numberVerification" : [ {
    "longNumber" : 281474976710655
  } ]
}, """)

    def writeNetworkAttributes(self, stream, fields):
        comma = False
        stream.write("""{
  "networkAttributes" : [ """)
        for name, value in fields.iteritems():
            if comma:
                stream.write(", ")
            stream.write('{\n    "n" : "' + name + '",\n    "v" : "' + value + '"\n  }')
            comma = True
        stream.write('  ]\n}')

    def writeMetadata(self, stream):
        stream.write("""{
  "metaData" : [ {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.nodesWritten) + """,
    "idCounter" : 1,
    "name" : "nodes",
    "properties" : [ ],
    "version" : "1.0"
  }, {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.edgesWritten) + """,
    "idCounter" : 1,
    "name" : "edges",
    "properties" : [ ],
    "version" : "1.0"
  } ]
}""")

    def addNode(self, name):
        """Adds a node to the network. Updates the node count."""
        if name not in self.nodesdb:
            self.nodecnt += 1
            self.nodesdb[name] = self.nodecnt
            self.visibleNodes[name] = True         # All nodes are initially visible
        return self.nodesdb[name]
    
    def countEdge(self, name):
        """Increment the edge counter for node `name'."""
        if name in self.nedges:
            self.nedges[name] += 1
        else:
            self.nedges[name] = 1

    def getMaxEdges(self):
        """Return the highest number of edges for any node in this network."""
        m = 0
        for e in self.nedges.values():
            if e > m:
                m = e
        return m

    def collectNodes(self):
        self.edgecnt = 0
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    self.addNode(hub)
                    for i in range(1, len(parsed), 2):
                        self.edgecnt += 1
                        self.addNode(parsed[i])
                        self.countEdge(hub)
                        self.countEdge(parsed[i])
        sys.stderr.write("{} nodes and {} edges read from ADJ file.\n".format(len(self.nodesdb), self.edgecnt))
        sys.stderr.write("Highest degree: {}.\n".format(self.getMaxEdges()))

    def filterNodesByDegree(self):
        """Set to invisible all nodes with a degree smaller than self.mindegree."""
        nvisible = 0
        for (node, deg) in self.nedges.iteritems():
            if deg >= self.mindegree:
                self.visibleNodes[node] = True
                nvisible += 1
            else:
                self.visibleNodes[node] = False
        sys.stderr.write("{} nodes with degree >= {}.\n".format(nvisible, self.mindegree))

    def filterHighestDegree(self):
        """Set as visible the first self.maxnodes nodes in order of decreasing degree, and hide the rest."""
        lowdegree = 0
        degs = [(deg, node) for (node, deg) in self.nedges.iteritems() ]
        degs.sort(key=lambda d: d[0], reverse=True)
        good = degs[:self.maxnodes]
        bad = degs[self.maxnodes:]
        for g in good:
            self.visibleNodes[g[1]] = True
            lowdegree = g[0]
        for b in bad:
            self.visibleNodes[b[1]] = False
        sys.stderr.write("Top {} nodes visible, lowest degree={}.\n".format(self.maxnodes, lowdegree))

    def writeNodes(self, stream):
        self.nodesWritten = 0
        comma = False
        stream.write('{  "nodes" : [ ')
        for node, idx in self.nodesdb.iteritems():
            if self.visibleNodes[node]:
                if comma:
                    stream.write(", ")
                stream.write('{\n    "@id" : ' + str(idx) + ',\n    "n" : "' + node + '"\n  }')
                comma = True
                self.nodesWritten += 1
        stream.write(' ]\n}')
        sys.stderr.write("{} nodes written to CX stream.\n".format(self.nodesWritten))

    def writeEdges(self, stream):
        self.edgesWritten = 0
        comma = False
        mis = []
        stream.write('{  "edges" : [ ')
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    hubi = self.nodesdb[hub]
                    for i in range(1, len(parsed), 2):
                        gene = parsed[i]
                        if self.visibleNodes[hub] or self.visibleNodes[gene]:
                            genei = self.nodesdb[gene]
                            mi = parsed[i+1]
                            if comma:
                                stream.write(", ")
                            self.edgecnt += 1
                            mis.append((self.edgecnt, mi))
                            stream.write('{\n    "@id" : ' + str(self.edgecnt) + ',\n    "s" : ' + str(hubi) + ',\n    "t" : ' + str(genei) + '\n  }')
                            comma = True
                            self.edgesWritten += 1
        stream.write('  ]\n}')

        # Now write edge attributes
        stream.write(', {  "edgeAttributes" : [ ')
        comma = False
        for m in mis:
            if comma:
                stream.write(", ")
            stream.write('{\n    "po" : ' + str(m[0]) + ',\n    "n" : "MI",\n    "v" : "' + m[1] + '"\n  }')
            comma = True
        stream.write('  ]\n}')
        sys.stderr.write("{} edges written to CX stream.\n".format(self.edgesWritten))

    def writeClosing(self, stream):
        stream.write(' ]\n')

    def writeNetwork(self, stream=sys.stdout, attributes={}):
        self.writePreamble(stream)
        self.writeNetworkAttributes(stream, attributes)
        stream.write(", ")
        self.collectNodes()
        if self.mindegree:
            self.filterNodesByDegree()
        elif self.maxnodes:
            self.filterHighestDegree()
        self.writeNodes(stream)
        stream.write(", ")
        self.writeEdges(stream)
        stream.write(", ")
        self.writeMetadata(stream)
        self.writeClosing(stream)
        
if __name__ == "__main__":
    C = CXwriter(sys.argv[1])
    C.maxnodes = 4
    C.writeNetwork(stream=sys.stdout, attributes={'name': 'TestNet', 'description': 'A test network converted from ADJ'})
