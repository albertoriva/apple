#!/usr/bin/env python

#####################################################
#                                                   #
# apple.py - Aracne Processing PipeLine Extensions. #
#                                                   #
# (c) 2015, Alberto Riva (ariva@ufl.edu)            #
#           Son Le (nuibang@gmail.com               #
#           University of Florida                   #
#####################################################

import sys
import glob
import math
import gzip
import os.path
import random
import sqlite3
from scipy.stats import norm
import numpy.random

## ToDo:
## Check for division by zero in normalization (sigma=0) - Done
## Don't output hub genes with no edges - Done
## Finish -f and -c arguments
## Add command for table stats - Done
## Validate arguments for top-level commands
## Complete documentation
## Multiple files for translate command
## Command for COMPARE-MI-HISTOGRAMS
## Histogram command should return range to stdout
## Clean up output (messages to stderr)

# COMMANDS = []

# Utils

def message(string, *args):
    sys.stderr.write(string.format(*args) + "\n")

def genOpen(filename, mode):
    """Generalized open() function - works on both regular files and .gz files."""
    (name, ext) = os.path.splitext(filename)
    if ext == ".gz":
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def readLines(filename):
    """Returns the contents of `filename' as a list of strings, removing
the trailing return character."""
    with open(filename, "r") as f:
        lines = f.readlines()
    return [ s.rstrip("\n\r") for s in lines ]

def changeExtension(filename, newext):
    """Returns a new filename obtained by replacing the extension in `filename' with `newext'."""
    (name, ext) = os.path.splitext(filename)
    return name + newext

def countSamples(filename):
    """Returns the number of data columns in file `filename', minus two (to
account for columns containing gene identifier and name)."""
    with open(filename, "r") as f:
        line = f.readline().split("\t")
        return len(line) - 2

def ensureInt(s):
    try:
        return int(s)
    except ValueError:
        return None

def ensureFloat(s):
    try:
        return float(s)
    except ValueError:
        return False

def ensureFile(s):
    if os.path.isfile(s):
        return s
    else:
        return False

# Top level commands and usage message

class toplevelCall():
    name = ""
    argdesc = ""                # short description of arguments
    shortdesc = ""              # short description of command
    longdesc = ""               # long description of command (multi-line)

    def parse(self, args):
        """Parse the command-line arguments `args' for this command."""
        pass

    def run(self):
        """Run this command."""
        pass



    # print "Where command-and-arguments can be:\n"
    # print "  bootstrap filename rounds - bootstrap a file into `rounds' new files.\n"
    # print "  extract [-a] adj outfile genesfile - extract edges for the genes in file genesfile from the input adj file and write them to outfile in tab-delimited format.\n"
    # print "  stats filenames... - print statistics on all supplied filenames (in adj format).\n"
    # print "  random outfile genesfile nsamples - generate random expression data for the genes in `genesfile' on `nsamples' samples using a negative binomial distribution.\n"
    # print "  consensus [options] outfile infiles... - generate a consensus network from multiple .adj files."
    # print "    Options:"
    # print "    [-c countsfile] - write a tab-delimited file with two columns: support, number of occurrences of support."
    # print "    [-d datafile]   - write a tab-delimited file with five columns: hub, gene, support, sum of MI, P-value."
    # print "    [-p pval] [-nb] - Use the specified P-value to filter edges in output; -nb disables Bonferroni correction."
    # print "    [-s support]    - Only output edges found in at least `support' bootstrap files."
    # print "    [-f fraction]   - Like -s, but determines the support in order to have the specified fraction of edges in output."
    
#
# Bootstrap command
#

class bootstrapCommand(toplevelCall):
    filename = ""
    nrounds = 0

    name = "bootstrap"
    argdesc = "filename rounds"
    shortdesc= "bootstrap a file into `rounds' new files."
    longdesc = """
This command takes as input a file containing gene expression values, and generates `rounds' new files through a bootstrap procedure.

The input file is assumed to have genes in the rows and samples in the columns. The first two columns are reserved for gene identifiers. All remaining columns contain data for different samples.

Each output file will have the same number of columns as the input file, chosen at random from the input file, with replacement. Therefore a column from the input file may appear more than once (or not at all) in the output file.
"""

    def parse(self, args):
        if len(args) != 2:
            return False
        self.filename = args[0]
        if not ensureFile(self.filename):
            print "File {} does not exist.".format(self.filename)
            return False
        self.nrounds = ensureInt(args[1])
        if not self.nrounds:
            print "Argument `{}' should be a number.".format(args[1])
            return False
        return True

    def run(self):
        bootstrapData(self.filename, self.nrounds)

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

#
# Consensus command
#
# Consensus reconstruction
# This part of the program replaces the getconsensusnet.pl script in the
# Aracne distribution, which has a number of problems.

class consensusCommand(toplevelCall):
    outfile = None
    csvfile = None
    datafile = None             # File for full network data
    infiles = []
    pval = None
    bonf = True
    support = None
    fraction = None

    name = "consensus"
    argdesc = "[options] outfile infiles..."
    shortdesc = "generate a consensus network from multiple .adj files."
    longdesc = """
This command generates a consensus network from a list of .adj files (usually generated through a bootstrap procedure). The following options are available:
  [-c countsfile] - write a tab-delimited file with two columns: support, number of occurrences of support.
  [-d datafile]   - write a tab-delimited file with five columns: hub, gene, support, sum of MI, P-value.
  [-p pval]       - Use the specified P-value to filter edges in output (with Bonferroni correction).
  [-nb]           - If specified, disables Bonferroni correction (used with -p).
  [-s support]    - Only output edges found in at least `support' bootstrap files.
  [-f fraction]   - Like -s, but determines the support in order to have the specified fraction of edges in output.
"""

    def parse(self, args):
        self.outfile = None
        self.infiles = []
        next = ""
        for arg in args: 
            if next == "-p":
                self.pval = ensureFloat(arg)
                if not self.pval:
                    print "The value of -p should be a P-value."
                    return False
                next = ""
            elif next == "-d":
                self.datafile = arg
                next = ""
            elif next == "-c":
                self.csvfile = arg
                next = ""
            elif next == "-s":
                self.support = ensureInt(arg)
                if self.support == None:
                    print "The value of -s, `{}', should be an integer number.".format(arg)
                    return False
                next = ""
            elif next == "-f":
                self.fraction = ensureFloat(arg)
                if not self.fraction:
                    print "The value of -f should be a number between 0 and 1."
                    return False
                next = ""
            elif arg in ["-p", "-c", "-d", "-s", "-f"]:
                next = arg
            elif arg == "-nb":
                bonf=False
            elif self.outfile == None:
                self.outfile = arg
            else:
                self.infiles.append(arg)

        if self.outfile == None or len(self.infiles) == 0:
            print "The consensus command requires an output file and one or more input files."
            return False

        return True

    def run(self):
        outfile = self.outfile
        infiles = self.infiles
        print "Consensus reconstruction from {} input files.".format(len(infiles))
        print "  Output to: {}".format(outfile)
        if self.pval:
            print "  P-value threshold: {} (Bonferroni={})".format(self.pval, self.bonf)
        if self.support:
            print "  Support threshold: {}".format(self.support)
        if self.datafile:
            print "  Network data to: {}".format(self.datafile)
        if self.csvfile:
            print "  Support counts to: {}".format(self.csvfile)
        bs = bstable()
        bs.loadAllBootstrap(infiles)
        bs.computeMuSigma()
        bs.supportDist = bs.supportDistribution()
        print bs.supportDist
        if self.pval:
            bs.saveConsensusByPval(outfile, self.pval, bs.getTableHeader(infiles[0]), bonferroni=self.bonf)
        elif self.support:
            bs.saveConsensusBySupport(outfile, self.support, bs.getTableHeader(infiles[0]))
        if self.datafile:
            bs.saveNetworkData(self.datafile)
        if self.csvfile:
            bs.saveSupportDist(self.csvfile)

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
        """Return the index associated with `name', creating a new one if necessary."""
        if name in self.names:
            return self.names[name]
        else:
            new = self.nameidx
            self.names[name] = new
            self.invnames[new] = name
            self.nameidx += 1
            return new

    def decodeName(self, nameid):
        """Return the name associated with index `nameid' (this is the inverse of internName)."""
        return self.invnames[nameid]

    def addEdge(self, hub, gene, mi):
        """Add an edge from `hub' to `gene', with mutual information `mi'."""
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
        with genOpen(filename, "r") as f:
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
        return result

    def saveNetworkData(self, datafile):
        print "Saving network data to {}...".format(datafile)
        with open(datafile, "w") as out:
            for (hub, table) in self.totsupport.iteritems():
                for (gene, edge) in table.iteritems():
                    out.write("{}\t{}\t{}\t{}\t{}\n".format(self.decodeName(hub), self.decodeName(gene), edge[0], edge[1], self.getpval(edge[0])))

    def saveSupportDist(self, csvfile):
        print "Saving support distribution to {}...".format(csvfile)
        with open(csvfile, "w") as out:
            out.write("#Support\tCount\n")
            for (supp, cnt) in self.supportDist.iteritems():
                out.write("{}\t{}\n".format(supp, cnt))

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

class bstableDB(bstable):
    DBconn = None
    DBcurs = None

    def openDBconn(self, name):
        """Open a connection to database `name'."""
        self.DBconn = sqlite3.connect(name)
        self.DBcurs = curs = self.DBconn.cursor()
        try:
            curs.execute("DROP TABLE edges")
        except sqlite3.OperationalError:
            pass                # ignore error if table does not exist
        curs.execute("CREATE TABLE edges(gene1 int, gene2 int, support int, totmi float);")
        curs.execute("CREATE INDEX edges1idx ON edges(gene1);")
        curs.execute("CREATE INDEX edges2idx ON edges(gene1);")

    def closeDBconn(self):
        self.DBconn.commit()
        self.DBconn.close()

    def addEdge(self, hub, gene, mi):
        """Add an edge from `hub' to `gene', with mutual information `mi',
or update it if it already exists."""
        hubid = self.internName(hub)
        geneid = self.internName(gene)
        s1 = self.DBcurs.execute("SELECT support, totmi FROM edges WHERE gene1=? AND gene2=?", (hubid, geneid))
        data = s1.fetchone()
        if data == None:
            self.DBcurs.execute("INSERT INTO edges(gene1, gene2, support, totmi) VALUES (?, ?, ?, ?)", (hubid, geneid, 1, mi))
        else:
            newsupport = data[0] + 1
            newmi = data[1] + mi
            self.DBcurs.execute("UPDATE edges SET support=?, totmi=? WHERE gene1=? AND gene2=?", (newsupport, newmi, hubid, geneid))
            
    def countEdges(self):
        """Returns the total number of edges currently in this bstable."""
        s1 = self.DBcurs.execute("SELECT count(*) FROM edges")
        data = s1.fetchone()
        if data:
            return data[0]
        else:
            return 0

# Extract a subset of genes from an .adj file

class extractCommand(toplevelCall):
    adjfile = None
    outfile = None
    genesfile = None
    geneslist = None
    genesfile2 = None
    geneslist2 = None
    both = False

    name = "extract"
    argdesc = "[-a] [-o outfile] adj outfile genesfile [genesfile2]"
    shortdesc = "extract edges for the genes in file genesfile from the input adj file and write them in tab-delimited format."
    longdesc = """
This command extracts the genes specified in `genesfile' from the .adj file in input and writes their edges to `outfile'.

File `genesfile' should have a single column containing gene identifiers (one per line).

The output (sent to standard output, or to a file specified with the -o option) is tab-delimited with three columns: 
hub gene, target gene, MI. The hub gene is always one of the genes specified in `genesfile', while MI is the mutual 
information of the edge connecting it to the target gene.

If the -a option is specified, both the hub gene and the target gene are required to be in `genesfile'.
"""

    def parse(self, args):
        next = ""
        for arg in args:
            if next == "-o":
                self.outfile = arg
                next = ""
            elif arg == "-o":
                next = "-o"
            elif arg == "-a":
                self.both = True
            elif self.adjfile == None:
                self.adjfile = arg
            elif self.genesfile == None:
                self.genesfile = arg
            else:
                self.genesfile2 = arg

        if self.adjfile == None:
            print "The adjfile argument is required."
            return False
        if self.genesfile == None:
            print "The genesfile argument is required."
            return False
        self.geneslist = readLines(self.genesfile)
        if self.genesfile2 != None:
            self.geneslist2 = readLines(self.genesfile2)
        return True

    def run(self):
        extractGeneSubset(self.adjfile, self.outfile, self.geneslist, both=self.both, genes2=self.geneslist2)

def addToSeen(seen, g1, g2):
    """Record the fact that an edge connecting `g1' and `g2' was seen in the current network."""
    if g1 in seen:
        seen[g1].append(g2)
    else:
        seen[g1] = [g2]

def wasNotSeen(seen, g1, g2):
    """Returns True if an edge connecting `g1' to `g2' has NOT been seen yet in the current network."""
    if g1 in seen:
        return not (g2 in seen[g1])
    else:
        return True

def extractGeneSubset(infile, outfile, genes, both=False, genes2=None):
    """Extract only the edges where one of the two genes belongs to the `genes' list
from adj file `infile' and write them to `outfile' in tab-delimited format (g1, g2, mi).
Also filters out symmetrical edges (ie, each pair is written only once regardless of order).
If `both' is True, only output edges where both genes belong to the list."""

    seen = {}
    extracted = 0
    paired = (genes2 != None)
    if paired:
        message("Extracting {} and {} genes from {} (paired mode)", len(genes), len(genes2), infile)
    else:
        message("Extracting {} genes from {} (both={})", len(genes), infile, both)

    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, "w")
    try:
        out.write("#Gene1\tGene2\tMI\n")
        with genOpen(infile, "r") as f:
            for line in f:
                if line[0] != ">":
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]

                    if paired:
                        found = False
                        if hub in genes:
                            found = 1
                        elif hub in genes2:
                            found = 2
                        if found:
                            for i in range(1, len(parsed), 2):
                                gene = parsed[i]
                                if ((found == 1) and (gene in genes2)) or ((found == 2) and (gene in genes)):
                                    if wasNotSeen(seen, gene, hub):
                                        out.write("{}\t{}\t{}\n".format(hub, gene, parsed[i+1]))
                                        addToSeen(seen, hub, gene)
                                        extracted += 1
                            
                    elif hub in genes:
                        for i in range(1, len(parsed), 2):
                            gene = parsed[i]
                            if (both == False) or (gene in genes):
                                if wasNotSeen(seen, gene, hub):
                                    out.write("{}\t{}\t{}\n".format(hub, gene, parsed[i+1]))
                                    addToSeen(seen, hub, gene)
                                    extracted += 1

                    else:
                        if both == False:
                            for i in range(1, len(parsed), 2):
                                gene = parsed[i]
                                if gene in genes:
                                    if wasNotSeen(seen, gene, hub):
                                        out.write("{}\t{}\t{}\n".format(hub, gene, parsed[i+1]))
                                        addToSeen(seen, hub, gene)
                                        extracted += 1
    finally:
        if outfile != None:
            out.close()
    message("{} edges extracted.".extracted)


#
# Stats command
#

class statsCommand(toplevelCall):
    filenames = None
    outfile = None

    name = "stats"
    argdesc = "[-o outfile] filenames..."
    shortdesc = "print statistics on all supplied filenames (in adj format)."
    longdesc = """
This command prints statistics on all the .adj files supplied as arguments. For each file, the command
prints: the number of hub genes, the total number of edges, and the average number of edges per hub.

Output is in tab-delimited format, and is written to standard output or to the file specified with the
-o option.
"""

    def parse(self, args):
        self.filenames = []
        next = ""

        for a in args:
            if next == "-o":
                self.outfile = a
                next = ""
            elif a == "-o":
                next = a
            else:
                self.filenames.append(a)
        return True

    def run(self):
        doAracneStats(self.filenames, self.outfile)

def aracneTableStats(filename):
    """Returns a tuple containing: the number of rows in the file, the sum of the
number of elements of each row, and the ratio between the sum and the number of
rows (ie, the average number of elements in each row)."""
    nrows = 0
    nfields = 0
    with genOpen(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                nrows += 1
                nfields += (line.count("\t") / 2)
    return (nrows, nfields, 1.0 * nfields / nrows)

def doAracneStats(filenames, outfile=None):
    """Call aracneTableStats on all files in `filenames' printing the 
results to standard output (or to `outfile' if provided) in tab-delimited format."""
    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, "w")

    out.write("Filename\tRows\tTotEdges\tAvgEdges\n")
    for f in filenames:
        if os.path.isfile(f):
            stats = aracneTableStats(f)
            out.write("{}\t{}\t{}\t{}\n".format(f, stats[0], stats[1], stats[2]))
    if outfile:
        out.close()

def aracneAllStats(outfile):
    with open(outfile, "w") as out:
        for i in range(1, 100):
            print i
            dpi = "results/all.gene_RPKM.aracne.nozero-{}.dpi.adj".format(i)
            if os.path.isfile(dpi):
                tp1 = aracneTableStats("results/all.gene_RPKM.aracne.nozero-{}.adj".format(i))
                tp2 = aracneTableStats(dpi)
                out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i, tp1[0], tp1[1], tp1[2], tp2[1], tp2[2]))

# Histogram of MI values

class histogramCommand(toplevelCall):
    infile = None
    outfile = None
    mifile = None
    summi = False
    nbins = 100
    low = None
    high = None
    overflow = False

    name = "histogram"
    argdesc = "[options] infile"
    shortdesc = "generate histogram of MI values from adj files."
    longdesc = """This command computes the histogram of MI values for the edges in the specified .adj file.
The following options are available:
  [-o outfile] - write output to `outfile' instead of standard output.
  [-n nbins]   - Specifiy number of bins to use (100 by default).
  [-r min max] - Only consider values between `min' and `max' (by default, the whole range of MIs is used).
  [-v]         - If specified, values higher than `max' are added to the last bin.
  [-s]         - If specified, the histogram is computed on the sum of the MIs of each row.
  [-m mifile]  - Write all distinct MI values to `mifile'.
"""

    def parse(self, args):
        next = None

        for a in args:
            if a in ["-m", "-n", "-o", "-r"]:
                next = a
            elif next == "-m":
                self.mifile = a
                next = None
            elif next == "-n":
                n = ensureInt(a)
                if n == None:
                    print "The value of -n, {}, should be a number.".format(a)
                    return False
                self.nbins = n
                next = None
            elif next == "-o":
                self.outfile = a
                next = None
            elif next == "-r":
                if self.low == None:
                    self.low = ensureFloat(a)
                else:
                    self.high = ensureFloat(a)
                    next = None
            elif a == "-s":
                self.summi = True
            elif a == "-v":
                self.overflow = True
            else:
                self.infile = a

        if self.infile == None or not os.path.isfile(self.infile):
            print "This command requires an input file."
            return False

        return True

    def run(self):
        doMIhistogram(self.infile, self.outfile, mifile=self.mifile, nbins=self.nbins, summi=self.summi,
                      low=self.low, high=self.high, overflow=self.overflow)

def doMIhistogram(filename, outfile, mifile=None, nbins=100, summi=False, low=None, high=None, overflow=False):
    mis = []
    bins = [None for i in range(nbins+1)]
    outstream = sys.stdout

    if summi:
        message("Computing histogram of sum of MI values for each row.")
    else:
        message("Computing histogram of all MI values.")
    if outfile:
        message("Writing histogram to file {}.", outfile)

    with genOpen(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                parsed = line.rstrip("\n").split("\t")
                if summi:
                    tot = 0
                    for i in range(2, len(parsed), 2):
                        mi = float(parsed[i])
                        tot += mi
                    mis.append(tot)
                else:
                    for i in range(2, len(parsed), 2):
                        mi = float(parsed[i])
                        mis.append(mi)
    mis.sort()
    message("Actual MI range: {} - {}", mis[0], mis[-1])
    if mifile != None:
        with open(mifile, "w") as out:
            for m in mis:
                out.write(str(m) + "\n")

    if low == None:
        minmi = mis[0]
    else:
        minmi = low
    if high == None:
        maxmi = mis[-1]
    else:
        maxmi = high
    message("Histogram range: {} - {}", minmi, maxmi)
    step = (maxmi - minmi) / nbins
    limit = minmi + step
    for i in range(nbins+1):
        bins[i] = [ limit, 0]
        limit += step

    this = 0
    bin = bins[this]
    terminate = False

    for m in mis:
        while m >= bin[0]:       # if we are outside the current bin
            if this == nbins:
                if overflow:
                    pass
                else:
                    terminate = True
                break
            else:
                this += 1           # move to next bin
            bin = bins[this]    # until we find the one that contains this mi
        if terminate:
            break
        bin[1] += 1             # and increment its counter

    if outfile != None:
        outstream = open(outfile, "w")
    try:
        for bin in bins:
            outstream.write("{}\t{}\n".format(bin[0], bin[1]))
    finally:
        if outfile != None:
            outstream.close()
    return bins

# Main function and top-level commands

class translateCommand(toplevelCall):
    table = None
    infiles = []
    outfiles = []

    name = "translate"
    argdesc = "table infile outfile ..."
    shortdesc = "translate identifiers in `infile' writing them to `outfile'."
    longdesc = """
This command converts the gene identifiers in `infile' according to a supplied translation table.
Gene identifiers are looked for in the first two columns of each input file.
File `table' should contain two columns: original gene name, converted name. 

Multiple pairs of input and output files may be specified on the command line. E.g.:

  apple.py translate table.txt in1.csv out1.csv in2.csv out2.csv in3.csv out3.csv
"""

    def parse(self, args):
        nargs = len(args)
        if nargs > 1 and (nargs % 2 == 1):
            self.table = args[0]
            for i in range(1, len(args), 2):
                self.infiles.append(args[i])
                self.outfiles.append(args[i+1])
        else:
            message("The translate command requires at least three filenames as arguments: table file, input file, output file.")
            return False

        if not os.path.isfile(self.table):
            message("Translation table `{}' does not exist.", self.table)
            return False

        return True
    
    def run(self):
        table = {}
        with open(tablefile, "r") as f:
            for line in f:
                parsed = line.rstrip("\n").split("\t")
                if parsed[0] != "\\N":
                    table[parsed[0]] = parsed[1]
        message("Translation table loaded, {} genes.", len(table))
        for i in range(len(self.infiles)):
            translateFile(table, self.infiles[i], self.outfiles[i])

def translateFile(table, infile, outfile):
    nseen = 0
    nfound = 0
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            for line in f:
                nseen += 1
                parsed = line.rstrip("\n").split("\t")
                if parsed[0] in table:
                    parsed[0] = table[parsed[0]]
                    nfound += 1
                if parsed[1] in table:
                    parsed[1] = table[parsed[1]]
                    nfound +=1
                out.write("\t".join(parsed) + "\n")
    message("{}: {} rows, {} identifiers translated", infile, nseen, nfound)

#
# Generate random dataset
#

class randomCommand(toplevelCall):
    outfile = None
    genesfile = None
    nsamples = None

    name = "random"
    argdesc = "[-o outfile] [-nb nsamples] genesfile"
    shortdesc = "generate random expression data for the genes in `genesfile' on `nsamples' samples using a negative binomial distribution."
    longdesc = """
This command generates a simulated gene expression dataset, using one of two different methods:

  If -nb is specified, the program will generate `nsamples' values for each gene listed in the first 
  column of file `genesfile', using a negative binomial distribution. The output file will have a
  number of columns equal to nsamples+2, with the first two columns containing the gene name.

  Otherwise, the expression values in `genesfile' (all values in the row except for the first two)
  will be shuffled.

Output will be written to standard output or to the file specified with the -o option.
"""
    
# apple.py random infile -o outfile -nb nsamples   - random data using negative binomial
# apple.py random infile -o outfile                - shuffling of original data

    def parse(self, args):
        nargs = 0
        next = ""

        for a in args:
            if a in ["-o", "-nb"]:
                next = a
            elif next == "-o":
                self.outfile = a
                next = ""
            elif next == "-nb":
                n = ensureInt(a)
                if n == None:
                    print "The argument {} should be a number.".format(a)
                    return False
                else:
                    self.nsamples = n
            else:
                self.genesfile = a

        if self.genesfile == None:
            print "This command requires an input file."
            return False

        if not os.path.isfile(self.genesfile):
            print "File {} does not exist.".format(self.genesfile)
            return False

        return True
    
    def run(self):
        if self.nsamples:
            doRandomDataset(self.genesfile, self.nsamples, outfile=self.outfile)
        else:
            doShuffledDataset(self.genesfile, outfile=self.outfile)

def doRandomDataset(genesfile, nsamples, outfile=None):
    print "Generating random dataset with {} samples.".format(nsamples)
    with open(outfile, "w") as out:
        out.write("GeneId\tGeneName\n")
        for i in range(0, nsamples):
            out.write("\tSample{}".format(i+1))
        out.write("\n")

        with open(genesfile, "r") as gf:
            for gene in gf:
                gene = gene.rstrip("\n\r")
                r = random.randint(1, nsamples)
                p = random.uniform(0.0, 1.0)
                out.write("{}\t{}".format(gene, gene))
                for i in range(0, nsamples):
                    out.write("\t{}".format(numpy.random.negative_binomial(r, p)))
                out.write("\n")

def doShuffledDataset(infile, outfile=None):
    print "Shuffling dataset {}.".format(infile)
    out = sys.stdout
    if outfile != None:
        print "Writing dataset to file {}.".format(outfile)
        out = open(outfile, "w")

    try:
        with open(infile, "r") as f:
            for line in f:
                parsed = line.rstrip("\n").split("\t")
                values = parsed[2:]
                random.shuffle(values)
                out.write("{}\t{}".format(parsed[0], parsed[1]))
                for x in values:
                    out.write("\t{}".format(x))
                out.write("\n")
    finally:
        if outfile != None:
            out.close()

# Filter command

class filterCommand(toplevelCall):
    threshold = None
    total = False
    infile = None
    outfile = None

    name = "filter"
    argdesc = "[options] infile threshold"
    shortdesc = "filter an adj file keeping only edges with MI over the threshold."
    longdesc = """
This command writes a new .adj file containing only the edges with an MI value over
the specified `threshold'. The following options are available:

  [-o outfile] - write output to `outfile' instead of standard output.
  [-t]         - apply the threshold to the sum of all MI values for a hub, and write
                 the whole line if successful.
"""

    def parse(self, args):
        next = ""
        for a in args:
            if a in ["-o"]:
                next = a
            elif next == "-o":
                self.outfile =a
                next = ""
            elif a == "-t":
                self.total = True
            elif self.infile == None:
                self.infile = a
            elif self.threshold == None:
                self.threshold = ensureFloat(a)
        if self.infile == None:
            print "This command requires an input file."
            return False
        if self.threshold == None:
            print "This command requires an MI threshold."
            return False
        return True

    def run(self):
        doFilter(self.infile, self.threshold, outfile=self.outfile, total=self.total)

def doFilter(infile, threshold, outfile=None, total=False):
    out = sys.stdout
    nin = 0
    nout = 0
    nrows = 0

    if outfile != None:
        out = open(outfile, "w")
    try:
        with genOpen(infile, "r") as f:
            for line in f:
                if line[0] == ">":
                    out.write(line)
                else:
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]

                    if total:
                        nin += 1
                        tot = 0
                        for i in range(2, len(parsed), 2):
                            tot += float(parsed[i])
                        if tot >= threshold:
                            out.write(line)
                            nout += 1
                    else:
                        towrite = True
                        for i in range(1, len(parsed), 2):
                            gene = parsed[i]
                            mi = float(parsed[i+1])
                            nin += 1
                            if mi >= threshold:
                                if towrite:
                                    out.write(hub)
                                    nrows += 1
                                    towrite = False
                                out.write("\t{}\t{}".format(gene, mi))
                                nout += 1
                        if not towrite:
                            out.write("\n")
    finally:
        if outfile != None:
            out.close()
    if total:
        print "{} hub genes seen, {} written.".format(nin, nout)
    else:
        print "{} edges seen, {} edges written for {} hub genes.".format(nin, nout, nrows)

# Network -> adj conversion

class convertCommand(toplevelCall):
    infile = None
    outfile = None
    operator = None
    support = 0

    name = "convert"
    argdesc = "[options] op infile outfile"
    shortdesc = "convert `infile' to a different format according to operator `op' and write the results to `outfile'."
    longdesc = """
This command converts between different file formats, according to the specified `op'. `op' can be one of:

  na - convert from networkData format to adj 
  nc - convert from networkData format to cytoscape
  ca - convert from cytoscape format to adj
"""

    def parse(self, args):
        next = ""
        nargs = 0
        for a in args:
            if next == "-s":
                self.support = ensureInt(a)
                next = ""
            elif a in ["-s"]:
                next = a
            elif nargs == 0:
                self.operator = a
                nargs += 1
            elif nargs == 1:
                self.infile = a
                nargs += 1
            elif nargs == 2:
                self.outfile = a
                nargs += 1
        if nargs < 3:
            print "This command requires three arguments."
            return False
        return True

    def run(self):
        op = self.operator
        if op == "na":
            networkToAdj(self.infile, self.outfile, support=self.support)
        elif op == "nc":
            networkToCytoscape(self.infile, self.outfile, support=self.support)
        elif op == "ca":
            cytoscapeToAdj(self.infile, self.outfile)
        else:
            print "Operator should be one of: na, nc, ca."

def networkToAdj(infile, outfile, support=0):
    current = ""
    genes = []
    with open(outfile, "w") as out:
        with genOpen(infile, "r") as f:
            for line in f:
                parsed = line.rstrip("\n").split("\t")
                hub = parsed[0]
                if hub != current:
                    if len(genes) > 0:
                        out.write("{}".format(hub))
                        for g in genes:
                            out.write("\t{}\t{}".format(g[0], g[1]))
                        out.write("\n")
                        current = hub
                        genes = []
                supp = int(parsed[2])
                if supp >= support:
                    genes.append((parsed[1], float(parsed[3]) / supp))
        if len(genes) > 0:
            out.write("{}".format(hub))
            for g in genes:
                out.write("\t{}\t{}".format(g[0], g[1]))
            out.write("\n")

def cytoscapeToAdj(infile, outfile):
    current = ""
    hdr = True
    with open(outfile, "w") as out:
        with genOpen(infile, "r") as f:
            for line in f:
                if hdr:
                    hdr = False
                else:
                    parsed = line.rstrip("\n").split("\t")
                    hub = parsed[0]
                    if hub != current:
                        if current != "":
                            out.write("\n")
                        out.write("{}".format(hub))
                        current = hub
                    out.write("\t{}\t{}".format(parsed[1], parsed[2]))
        out.write("\n")

def networkToCytoscape(infile, outfile, support=0):
    with open(outfile, "w") as out:
        with genOpen(infile, "r") as f:
            for line in f:
                parsed = line.rstrip("\n").split("\t")
                supp = int(parsed[2])
                if supp >= support:
                    out.write("{}\t{}\t{}\n".format(parsed[0], parsed[1], float(parsed[3]) / supp))
                
# Main

def main():
    argc = len(sys.argv)
    if argc == 1:
        return usage()
    elif argc == 2:
        return usage(sys.argv[1])

    name = sys.argv[1]
    command = findCommand(name)

    if command == None:
        return usage()

    if command.parse(sys.argv[2:]):
        command.run()
    else:
        return usage(name)

COMMANDS = {"bootstrap": bootstrapCommand(),
            "consensus": consensusCommand(),
            "extract": extractCommand(), 
            "stats": statsCommand(),
            "histogram": histogramCommand(),
            "random": randomCommand(),
            "translate": translateCommand(),
            "filter": filterCommand(),
            "convert": convertCommand()}

def findCommand(name):
    global COMMANDS
    if name in COMMANDS:
        return COMMANDS[name]
    else:
        return None

def usage(what=None):
    global COMMANDS

    if what == None:
        message("Usage: apple.py command arguments...\n")
        message("The following commands are available:\n")
        for (name, c) in COMMANDS.iteritems():
            message("  {} - {}", name, c.shortdesc)
        message("\nUse 'apple.py command' to get a description of each command and its arguments.\n")

    else:
        c = findCommand(what)
        if c == None:
            usage()
        else:
            message("Usage: apple.py {} {} - {}", c.name, c.argdesc, c.shortdesc)
            if c.longdesc != "":
                message(c.longdesc)
            
if __name__ == "__main__":
    sys.stderr.write("apple.py - Aracne Processing PipeLine Extensions.\n\n")
    main()
    
# Everything from this point on is not part of the apple.py program. These functions
# are leftovers from tests and one-off scripts.

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

"""
STEPS:

1. Call doAracneBootstrap() to write run.sh file;
2. Execute run.sh to submit jobs to queue;
3. When jobs are done, delete or compress all .out files (optional);
4. Call getconsensus.qsub to generate consensus network from all .adj files.
"""

def doAracneBootstrap(filename, outdir="results/", nrounds=100, pval="1e-7", dpi=0.1):
    files = bootstrapData(filename, nrounds)
    writeSubmitScript("run.sh", files, outdir=outdir, pval=pval, dpi=dpi)
    return files
    
# Some random utils

def writeReconstruct(outfile, start, end):
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        for i in range(start, end):
            out.write("echo -n {} \n".format(i))
            # out.write("submit ../reconstruct.qsub ../SPLIT all.gene_RPKM.aracne.nozero-{}.adj {}/all.gene_RPKM.aracne.nozero-{}\n".format(i, i, i))
            out.write("submit reconstruct.qsub ../SPLIT randomdata-{}.adj {}/randomdata-{}\n".format(i, i, i))

def writeStep2(outfile, start, end, tfs=False, source="all.gene_RPKM.aracne.nozero"):
    tfsopt = ""
    if tfs:
        tfsopt = " tfs={}".format(tfs)
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        src = source + ".txt"
        dst = source + "-{}.dpi.adj"
        adj = source + "-{}.adj"
        for i in range(start, end):
            idst = dst.format(i)
            iadj = adj.format(i)
            out.write("submit aracne.qsub {} {} dpi=0.1 adj={} {}\n".format(src, idst, iadj, tfsopt))

def writeStep2b(outfile, start, end):
    with open(outfile, "w") as out:
        out.write("#!/bin/bash\n\n")
        for i in range(start, end):
            dpi = "results/all.gene_RPKM.aracne.nozero-{}.dpi.adj".format(i)
            if not os.path.isfile(dpi):
                out.write("submit aracne.qsub ../all.gene_RPKM.aracne.nozero.txt all.gene_RPKM.aracne.nozero-{}.dpi.adj dpi=0.1 adj=all.gene_RPKM.aracne.nozero-{}.adj\n".format(i, i))

