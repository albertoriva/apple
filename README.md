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

