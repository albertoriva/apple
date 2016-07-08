# apple.py
A command-line tool that facilitates working with the ARACNE program. Usage:

```
apple.py command arguments...
```

The following commands are available:

|Command|Description|
|-------|-----------|
|bootstrap|Bootstrap a file into *rounds* new files, each containing *samplesize* columns.|
|consensus|Generate a consensus network from multiple .adj files.|
|convert|Convert *infile* to a different format according to operator *op* and write the results to *outfile*.|
|extract|Extract edges for the genes in file genesfile from the input adj file and write them in tab-delimited format.|
|filter|Filter an adj file keeping only edges with MI over the threshold.|
|histogram|Generate histogram of MI values from adj files.|
|random|Generate random expression data for the genes in *genesfile* on *nsamples* samples using a negative binomial distribution.|
|stats|Print statistics on all supplied filenames (in adj format).|
|translate|Translate identifiers in *infile* writing them to *outfile*.|

Use 'apple.py command' to get a description of each command and its arguments.

## Command descriptions
Detail usage of each apple.py command is provided below.

### Bootstrap

Usage: 

```
apple.py bootstrap [-z samplesize] filename rounds
```

This command takes as input a file containing gene expression values, and generates *rounds* new files through a bootstrap procedure.

The input file is assumed to have genes in the rows and samples in the columns. The first two columns are reserved for gene identifiers. All remaining columns contain data for different samples.

Each output file will have the same number of columns as the input file (unless a different number is specified with the -z argument), chosen at random from the input file, with replacement. Therefore a column from the input file may appear more than once (or not at all) in the output file.

### Consensus

Usage: 

```
apple.py consensus [options] outfile infiles...
```

This command generates a consensus network from a list of .adj files (usually generated through a bootstrap procedure). The following options are available:

```
  [-c countsfile] - write a tab-delimited file with two columns: support, number of occurrences of support.
  [-d datafile]   - write a tab-delimited file with five columns: hub, gene, support, sum of MI, P-value.
  [-p pval]       - Use the specified P-value to filter edges in output (with Bonferroni correction).
  [-nb]           - If specified, disables Bonferroni correction (used with -p).
  [-s support]    - Only output edges found in at least `support' bootstrap files.
  [-f fraction]   - Like -s, but determines the support in order to have the specified fraction of edges in output.
```

### Convert

Usage: apple.py convert [options] op infile outfile - convert `infile' to a different format according to operator `op' and write the results to `outfile'.

This command converts between different file formats, according to the specified `op'. `op' can be one of:

  na - convert from networkData format to adj
  nc - convert from networkData format to cytoscape
  ca - convert from cytoscape format to adj
  co - convert from cytoscape format to connections
  ac - convert from adj to cytoscape

Formats:
cytocscape: gene1, gene2, mi
connection: hub, number of connected genes, list of connected genes
