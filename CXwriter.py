#!/usr/bin/env python

import sys

class CXwriter():
    nodesdb = {}
    nedges = {}
    nodecnt = 0                 # Counter for @id fields in nodes
    edgecnt = 0                 # Counter for @id fields in edges
    adjfile = None
    
    def __init__(self, adjfile):
        self.nodesdb = {}
        self.nodecnt = 0
        self.edgecnt = 0
        self.adjfile = adjfile

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
    "elementCount" : """ + str(self.nodecnt) + """,
    "idCounter" : 1,
    "name" : "nodes",
    "properties" : [ ],
    "version" : "1.0"
  }, {
    "consistencyGroup" : 1,
    "elementCount" : """ + str(self.edgecnt) + """,
    "idCounter" : 1,
    "name" : "edges",
    "properties" : [ ],
    "version" : "1.0"
  } ]
}""")

    def addNode(self, name, nedges):
        if name not in self.nodesdb:
            self.nodecnt += 1
            self.nodesdb[name] = self.nodecnt
        return self.nodesdb[name]
        
    def collectNodes(self):
        with open(self.adjfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    self.addNode(hub, nedges)
                    for i in range(1, len(parsed), 2):
                        self.addNode(parsed[i])
        sys.stderr.write("{} nodes read from ADJ file.\n".format(len(self.nodesdb)))
        
    def writeNodes(self, stream):
        comma = False
        stream.write('{  "nodes" : [ ')
        for node, idx in self.nodesdb.iteritems():
            if comma:
                stream.write(", ")
            stream.write('{\n    "@id" : ' + str(idx) + ',\n    "n" : "' + node + '"\n  }')
            comma = True
        stream.write(' ]\n}')

    def writeEdges(self, stream):
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
                        genei = self.nodesdb[gene]
                        mi = parsed[i+1]
                        if comma:
                            stream.write(", ")
                        self.edgecnt += 1
                        mis.append((self.edgecnt, mi))
                        stream.write('{\n    "@id" : ' + str(self.edgecnt) + ',\n    "s" : ' + str(hubi) + ',\n    "t" : ' + str(genei) + '\n  }')
                        comma = True
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
        sys.stderr.write("{} edges read from ADJ file.\n".format(self.edgecnt))

    def writeClosing(self, stream):
        stream.write(' ]\n')

    def writeNetwork(self, stream=sys.stdout, attributes={}):
        C.writePreamble(stream)
        C.writeNetworkAttributes(stream, attributes)
        stream.write(", ")
        C.collectNodes()
        C.writeNodes(stream)
        stream(", ")
        C.writeEdges(stream)
        stream(", ")
        C.writeMetadata(stream)
        C.writeClosing(stream)
        
if __name__ == "__main__":
    C = CXwriter(sys.argv[1])
    C.writeNetwork(stream=sys.stdout, attributes={'name': 'TestNet', 'description': 'A test network converted from ADJ'})
