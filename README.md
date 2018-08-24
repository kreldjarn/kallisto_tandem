# kallisto_tandem.py

This script serves as a pre- and postprocessing pipeline for using kallisto to align against separate transcriptomes of varying disparity in tandem.

## Preprocessing
The preprocessing step combines one or more `fasta` files and gives the entries in each file an identifier corresponding to the data it contains.
```
./kallisto_tandem.py concatenate -i /path/to/input/file -o /path/to/output/file.fasta
```
Here, `/path/to/input/file` is a config file, where each line corresponds to a single `fasta` file to be concatenated. This file is on the following format:
```
path/to/file1.fasta class variant
path/to/file2.fasta class variant
...
```
`class` and `variant` can be whichever identifiers you feel are descriptive of the fasta files, but usually `class` refers to an organism, and `variant` refers to different versions of its transcriptome, e.g. one of `processed, unprocessed`. Once combined, you can create a `kallisto` index from the resulting concatenated `fasta`file.

# Postprocessing
Postprocessing consists mainly of splitting the `abundance.tsv` file output from a tandem kallisto run into separate matrices for each `class` and for each `class-variant` combination, and then renormalizing the respective TPMs w.r.t. the size of the corresponding subset of the sample. Note that the `-o` flag should point to a folder in which the results are to be stored.

```
./kallisto_tandem.py renormalize -i /path/to/abundance.tsv -o /path/to/output_folder
```
