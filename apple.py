#!/usr/bin/env python

#####################################################
#                                                   #
# apple.py - Aracne Processing PipeLine Extensions. #
#                                                   #
# (c) 2015, Alberto Riva, Son Le                    #
#           University of Florida                   #
#####################################################

import sys
import glob
import math
import os.path
import random
from scipy.stats import norm

## ToDo:
## Check for division by zero in normalization (sigma=0) - Done
## Don't output hub genes with no edges - Done
## Finish -f and -c arguments
## Add command for table stats
## Validate arguments for top-level commands

# Utils

def readLines(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return [ s.rstrip("\n") for s in lines ]

def changeExtension(filename, newext):
    (name, ext) = os.path.splitext(filename)
    return name + newext

def countSamples(filename):
    """Returns the number of data columns in file `filename', minus two (to
account for columns containing gene identifier and name)."""
    with open(filename, "r") as f:
        line = f.readline().split("\t")
        return len(line) - 2

def samplingWithReplacement(m):
    """Returns a list of m elements in the range 0..m-1 using
random sampling with replacement."""
    return [ random.randrange(m) for i in range(m) ]

def bootstrapData(filename, nrounds):
    """Given a data file `filename', generate `nrounds' new files
obtained by bootstrapping its columns. Returns the list of new filenames."""
    nsamples = countSamples(filename)
    outfiles = []

    for i in range(nrounds):
        columns = samplingWithReplacement(nsamples)
        (name, ext) = os.path.splitext(filename)
        outfilename = name + "-{}".format(i) + ext
        outfiles.append(outfilename)

        print "Writing {} ".format(outfilename)
        with open(outfilename, "w") as o:
            with open(filename, "r") as f:
                for line in f:
                    parsed = line.rstrip("\n").split("\t")
                    new = parsed[0:2]
                    for c in columns:
                        new.append(parsed[c+2])
                    o.write("\t".join(new) + "\n")
        
    return outfiles

def writeSubmitScript(outfile, datafiles, outdir="results/", pval="1e-7", dpi=0.1):
    with open(outfile, "w") as o:
        o.write("#!/bin/bash\n\n")
        o.write("module load dibig_tools\n\n")
        o.write("mkdir -p {}\n\n".format(outdir))
        for df in datafiles:
            of = changeExtension(df, ".out")
            oof = changeExtension(df,".adj")
            o.write("PREV1=`submit aracne-step1.qsub {} {}.aa IDSaa`\n".format(df, of))
            o.write("PREV2=`submit aracne-step1.qsub {} {}.ab IDSab`\n".format(df, of))
            o.write("PREV3=`submit aracne-step1.qsub {} {}.ac IDSac`\n".format(df, of))
            o.write("PREV4=`submit aracne-step1.qsub {} {}.ad IDSad`\n".format(df, of))
            o.write("submit -after $PREV1 -after $PREV2 -after $PREV3 -after $PREV4 aracne-step2.qsub {} {} {}/{} {} {}\n".format(df, of, outdir, oof, pval, dpi))
            o.write("\n")

def writeSubmitScriptSplit(outfile, datafiles, split, outdir="results/", pval="1e-7", dpi=0.1):
    spl = readLines(split)
    after = ""
    cnt = 1
    with open(outfile, "w") as o:
        o.write("#!/bin/bash\n\n")
        o.write("module load dibig_tools\n\n")
        o.write("mkdir -p {}\n\n".format(outdir))
        for df in datafiles:
            of = changeExtension(df, ".out")
            oof = changeExtension(df,".adj")
            for s in spl:
                o.write("PREV{}=`submit aracne-step1.qsub {} {}.{} SIDS{}`\n".format(cnt, df, of, s, s))
                after += " -after $PREV{}".format(cnt)
                cnt += 1

            o.write("submit {} aracne-step2.qsub {} {} {}/{} {} {}\n".format(after, df, of, outdir, oof, pval, dpi))
            o.write("\n")

def doAracneBootstrap(filename, outdir="results/", nrounds=100, pval="1e-7", dpi=0.1):
    files = bootstrapData(filename, nrounds)
    writeSubmitScript("run.sh", files, outdir=outdir, pval=pval, dpi=dpi)
    return files
    
"""
STEPS:

1. Call doAracneBootstrap() to write run.sh file;
2. Execute run.sh to submit jobs to queue;
3. When jobs are done, delete or compress all .out files (optional);
4. Call getconsensus.qsub to generate consensus network from all .adj files.
"""

# Consensus reconstruction

class bstable():
    totsupport = None
    totedge = None
    totbs = 0                   # Number of bootstrap files loaded
    mu = 0.0
    sigma = 0.0
    pvalues = None
    names = None
    invnames = None
    nameidx = 0
    supportDist = {}

    def __init__(self):
        self.totsupport = {}
        self.totedge = {}
        self.pvalues = {}
        self.names = {}
        self.invnames = {}
        self.supportDist = {}

    def internName(self, name):
        """Return the index associated with `name', creating
a new one if necessary."""
        if name in self.names:
            return self.names[name]
        else:
            new = self.nameidx
            self.names[name] = new
            self.invnames[new] = name
            self.nameidx += 1
            return new

    def decodeName(self, nameid):
        """Return the name associated with index `nameid' (this
is the inverse of internName)."""
        return self.invnames[nameid]

    def addEdge(self, hub, gene, mi):
        """Add an edge from `hub' to `gene', with
mutual information `mi'."""
        hubid = self.internName(hub)
        geneid = self.internName(gene)
        if hubid in self.totsupport:
            hubtable = self.totsupport[hubid]
        else:
            hubtable = {}
            self.totsupport[hubid] = hubtable

        if geneid in hubtable:
            edge = hubtable[geneid]
            edge[0] += 1
            edge[1] += mi
        else:
            edge = [1, mi]
            hubtable[geneid] = edge
        return edge

    def countEdges(self):
        """Returns the total number of edges currently in this bstable."""
        n = 0
        for (hub, table) in self.totsupport.iteritems():
            n += len(table)
        return n

    def getTableHeader(self, filename):
        """Returns the header of file `filename' as a string."""
        hdr = ""
        with open(filename, "r") as f:
            for line in f:
                if line[0] == ">":
                    hdr += line
                else:
                    return hdr

    def loadBootstrap(self, filename, idx):
        """Load a bootstrap file `filename' into the current
bstable, giving it index `idx'."""
        numedges = 0
        with open(filename, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    for i in range(1, len(parsed), 2):
                        gene = parsed[i]
                        mi = float(parsed[i+1])
                        self.addEdge(hub, gene, mi)
                        numedges += 1
        self.totedge[idx] = numedges
        self.totbs += 1
        return numedges

    def loadAllBootstrap(self, filenames):
        """Load all the files in the list `filenames' into the current bstable."""
        idx = 1
        prevedges = 0
        for filename in filenames:
            if os.path.isfile(filename):
                print "Parsing {}...".format(filename),
                self.loadBootstrap(filename, idx)
                newedges = self.countEdges()
                if prevedges == 0:
                    print "{} edges.".format(newedges)
                else:
                    print "{} edges ({}% increase).".format(newedges, int(round(100 * (newedges-prevedges)/prevedges)))
                prevedges = newedges
                idx += 1
            else:
                print "File {} does not exist, skipping.".format(filename)

    def computeMuSigma(self):
        """Compute mu and sigma for the current bstable."""
        totedge = self.countEdges()
        for i in range(self.totbs):
            prob = 1.0 * self.totedge[i+1] / totedge
            self.mu += prob
            self.sigma += prob * (1 - prob)
        self.sigma = math.sqrt(self.sigma)
        print "Mu = {}, Sigma = {}".format(self.mu, self.sigma)

    def getpval(self, support):
        """Returns the P-value corresponding to `support'. Saves computed
values in the pvalues dictionary for speed."""
        if support in self.pvalues:
            return self.pvalues[support]
        if self.sigma == 0:
            return 1
        else:
            z = (support - self.mu) / self.sigma
            p = 1 - norm.cdf(z)
            self.pvalues[support] = p
            return p

    def supportDistribution(self):
        """Returns a dictionary mapping each support count (ie,
number of bootstrap files an edge appears in) to the number of
times it occurs in totsupport."""
        result = {}
        for (hub, table) in self.totsupport.iteritems():
            for (gene, edge) in table.iteritems():
                support = edge[0]
                if support in result:
                    result[support] += 1
                else:
                    result[support] = 1
        result

    def saveNetworkData(datafile):
        print "Saving network data to {}...".format(datafile)
        with open(datafile, "w") as out:
            for (hub, table) in self.totsupport.iteritems():
                for (gene, edge) in table.iteritems():
                    out.write("{}\t{}\t{}\t{}\t{}\n".format(self.decodeName(hub), self.decodeName(gene), edge[0], edge[1], self.getpval(edge[0])))

    def saveConsensusByPval(self, filename, pval, header, bonferroni=True):
        if bonferroni:
            pval = pval / self.countEdges()
        print "Effective P-value threshold: {}".format(pval)
        print "Saving consensus network to {}...".format(filename)
        nwritten = 0
        with open(filename, "w") as out:
            out.write(header)
            for (hub, table) in self.totsupport.iteritems():
                if len(table) > 0:
                    out.write(self.decodeName(hub))
                    for (gene, edge) in table.iteritems():
                        support = edge[0]
                        mi = edge[1]
                        p = self.getpval(support)
                        if p < pval:
                            nwritten += 1
                            out.write("\t{}\t{}".format(self.decodeName(gene), mi / support))
                    out.write("\n")
        print "{} edges written to consensus network.".format(nwritten)

        return nwritten

    def saveConsensusBySupport(self, filename, minsupport, header):
        print "Saving consensus network to {}...".format(filename)
        nwritten = 0
        with open(filename, "w") as out:
            out.write(header)
            for (hub, table) in self.totsupport.iteritems():
                if len(table) > 0:
                    out.write(self.decodeName(hub))
                    for (gene, edge) in table.iteritems():
                        support = edge[0]
                        mi = edge[1]
                        if support >= minsupport:
                            nwritten += 1
                            out.write("\t{}\t{}".format(self.decodeName(gene), mi / support))
                    out.write("\n")
        print "{} edges written to consensus network.".format(nwritten)

        return nwritten

# Extract a subset of genes from an .adj file

def addToSeen(seen, g1, g2):
    if g1 in seen:
        seen[g1].append(g2)
    else:
        seen[g1] = [g2]

def wasNotSeen(seen, g1, g2):
    if g1 in seen:
        return not (g2 in seen[g1])
    else:
        return True

def extractGeneSubset(infile, outfile, genes, both=False):
    """Extract only the edges where one of the two genes belongs to the `genes' list
from adj file `infile' and write them to `outfile' in tab-delimited format (g1, g2, mi).
Also filters out symmetrical edges (ie, each pair is written only once regardless of order).
If `both' is True, only output edges where both genes belong to the list."""
    print "Extracting {} genes from {} (both={})".format(len(genes), infile, both)
    seen = {}
    with open(outfile, "w") as out:
        out.write("#Gene1\tGene2\tMI\n")
        with open(infile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    if hub in genes:
                        for i in range(1, len(parsed), 2):
                            gene = parsed[i]
                            if (both == False) or (gene in genes):
                                if wasNotSeen(seen, gene, hub):
                                    out.write("{}\t{}\t{}\n".format(hub, gene, parsed[i+1]))
                                    addToSeen(seen, hub, gene)
                    else:
                        if both == False:
                            for i in range(1, len(parsed), 2):
                                gene = parsed[i]
                                if gene in genes:
                                    if wasNotSeen(seen, gene, hub):
                                        out.write("{}\t{}\t{}\n".format(hub, gene, parsed[i+1]))
                                        addToSeen(seen, hub, gene)

# Some random utils

def writeReconstruct(outfile, start, end):
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        for i in range(start, end):
            out.write("echo -n {} \n".format(i))
            out.write("submit ../reconstruct.qsub ../SPLIT all.gene_RPKM.aracne.nozero-{}.adj {}/all.gene_RPKM.aracne.nozero-{}\n".format(i, i, i))

def writeStep2(outfile, start, end):
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        for i in range(start, end):
            out.write("submit aracne.qsub ../all.gene_RPKM.aracne.nozero.txt all.gene_RPKM.aracne.nozero-{}.dpi.adj dpi=0.1 adj=all.gene_RPKM.aracne.nozero-{}.adj\n".format(i, i))

def writeStep2b(outfile, start, end):
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        for i in range(start, end):
            dpi = "results/all.gene_RPKM.aracne.nozero-{}.dpi.adj".format(i)
            if not os.path.isfile(dpi):
                out.write("submit aracne.qsub ../all.gene_RPKM.aracne.nozero.txt all.gene_RPKM.aracne.nozero-{}.dpi.adj dpi=0.1 adj=all.gene_RPKM.aracne.nozero-{}.adj\n".format(i, i))

def aracneTableStats(filename):
    """Returns a tuple containing: the number of rows in the file, the sum of the
number of elements of each row, and the ratio between the sum and the number of
rows (ie, the average number of elements in each row)."""
    nrows = 0
    nfields = 0
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                nrows += 1
                nfields += (line.count("\t") / 2)
    return (nrows, nfields, 1.0 * nfields / nrows)

def aracneAllStats(outfile):
    with open(outfile, "w") as out:
        for i in range(1, 100):
            print i
            dpi = "results/all.gene_RPKM.aracne.nozero-{}.dpi.adj".format(i)
            if os.path.isfile(dpi):
                tp1 = aracneTableStats("results/all.gene_RPKM.aracne.nozero-{}.adj".format(i))
                tp2 = aracneTableStats(dpi)
                out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i, tp1[0], tp1[1], tp1[2], tp2[1], tp2[2]))

# Main function and top-level commands

def usage():
    print "Usage: apple.py command-and-arguments...\n"
    print "Where command-and-arguments can be:\n"
    print "  bootstrap filename rounds - bootstrap a file into `rounds' new files.\n"
    print "  extract [-a] adj outfile genesfile - extract edges for the genes in file genesfile from the input adj file and write them to outfile in tab-delimited format.\n"
    print "  consensus [options] outfile infiles... - generate a consensus network from multiple .adj files."
    print "    Options:"
    print "    [-c countsfile] - write a tab-delimited file with two columns: support, number of occurrences of support."
    print "    [-d datafile]   - write a tab-delimited file with five columns: hub, gene, support, sum of MI, P-value."
    print "    [-p pval] [-nb] - Use the specified P-value to filter edges in output; -nb disables Bonferroni correction."
    print "    [-s support]    - Only output edges found in at least `support' bootstrap files."
    print "    [-f fraction]   - Like -s, but determines the support in order to have the specified fraction of edges in output."
    
class consensusArgs():
    outfile = None
    csvfile = None
    datafile = None             # File for full network data
    infiles = []
    pval = None
    bonf = True
    support = None
    fraction = None
    
    def __init__(self, args):
        self.outfile = None
        self.infiles = []
        next = ""
        for arg in args: 
            if next == "-p":
                self.pval = float(arg)
                next = ""
            elif next == "-d":
                self.datafile = arg
                next = ""
            elif next == "-c":
                self.csvfile = arg
                next = ""
            elif next == "-s":
                self.support = int(arg)
                next = ""
            elif next == "-f":
                self.fraction = float(arg)
                next = ""
            elif arg in ["-p", "-c", "-d", "-s", "-f"]:
                next = arg
            elif arg == "-nb":
                bonf=False
            elif self.outfile == None:
                self.outfile = arg
            else:
                self.infiles.append(arg)

def runConsensus(outfile, infiles, params): # pval=None, bonf=True, csvfile=None):
    print "Consensus reconstruction from {} input files.".format(len(infiles))
    print "  Output to: {}".format(outfile)
    if params.pval:
        print "  P-value threshold: {} (Bonferroni={})".format(params.pval, params.bonf)
    if params.support:
        print "  Support threshold: {}".format(params.support)
    if params.datafile:
        print "  Network data to: {}".format(params.datafile)
    if params.csvfile:
        print "  Support counts to: {}".format(params.csvfile)
    bs = bstable()
    bs.loadAllBootstrap(infiles)
    bs.computeMuSigma()
    bs.supportDist = bs.supportDistribution()
    if params.pval:
        bs.saveConsensusByPval(outfile, params.pval, bs.getTableHeader(infiles[0]), bonferroni=params.bonf)
    elif params.support:
        bs.saveConsensusBySupport(outfile, params.support, bs.getTableHeader(infiles[0]))
    if params.datafile:
        bs.saveNetworkData(params.datafile)

class extractArgs():
    adjfile = None
    outfile = None
    genesfile = None
    geneslist = None
    both = False

    def __init__(self, args):
        for arg in args:
            if arg == "-a":
                self.both = True
            elif self.adjfile == None:
                self.adjfile = arg
            elif self.outfile == None:
                self.outfile = arg
            else:
                self.genesfile = arg

        if (self.adjfile == None) or (self.outfile == None) or (self.genesfile == None):
            usage()
            exit()

        self.geneslist = readLines(self.genesfile)
            

def main():
    if len(sys.argv) == 1:
        return usage()

    command = sys.argv[1]
    if command == "consensus":
        parsed = consensusArgs(sys.argv[2:])
        runConsensus(parsed.outfile, parsed.infiles, parsed)
    elif command == "bootstrap":
        bootstrapData(sys.argv[2], int(sys.argv[3]))
    elif command == "extract":
        parsed = extractArgs(sys.argv[2:])
        extractGeneSubset(parsed.adjfile, parsed.outfile, parsed.geneslist, both=parsed.both)
    else:
        usage()

if __name__ == "__main__":
    print "apple.py - Aracne Processing PipeLine Extensions.\n"
    main()
    

