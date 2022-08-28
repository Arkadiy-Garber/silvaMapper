Script for querying 16S reads from a metagenomic library against the SILVA database.

This script requires the SILVA database file (in FASTA format) as input. This database can be found here: https://www.arb-silva.de/no_cache/download/archive/current/Exports/.

The second required input is a FASTA file of raw metagenomic reads. To speed things up, I recommend providing only the 16S reads, pulled out from the metagenomic library using a program like SortMeRNA: https://bioinfo.lifl.fr/RNA/sortmerna/

Alternatively, instead of the reads, you can provide BLAST output (outfmt 6). This will skip the most time-intensive step of the pipeline.

This script is pre-compiled and does not need to be installed. Simply put it in your $PATH and run like 

    ./silvaMapper.py -h
