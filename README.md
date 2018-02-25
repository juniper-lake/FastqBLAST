# FastqBLAST
***
> **Author**: Danielle Novick  
> **Last Update**: February 15, 2017  
> **Creation Date**: October 24, 2017

> **Modified Date**: Febuary 24, 2018 by M. Joseph Tomlinson IV  
Summary Report file, Ability to Change Databases and Perform Organism Searches

## Program Description

This program takes a sample of sequences from a fastq file, trims the low quality ends, BLASTs them, fetches additional info from NCBI, and produces a tabular report.

## System Requirements
**Python 3** + the following packages

* biopython (download here: [http://biopython.org/wiki/Download](http://biopython.org/wiki/Download))
* urllib
* sys
* random
* argparse
* time
* collections
* numpy

## Running the Program
To run the program with all default settings, use the following command.

`python <scriptname> <filename> --email <email>`

For example,

`python FastqBLAST.py demo_data.fastq --email myemail@udel.edu` 


## Help File & Parameters
Running `python FastqBLAST.py --help` should return the following help file with information about all of the parameters. 


<pre>
usage: FastqBLAST.py [-h] --email EMAIL [--ascii64 ASCII64] [--nPercent NPERCENT]
                 [--nAbsolute NABSOLUTE] [--leadingQ LEADINGQ]
                 [--trailingQ TRAILINGQ] [--hitlistSize HITLISTSIZE] [--database DATABASE]
                 
                 filename

This script takes a sample of sequences from a fastq file, trims the low
quality ends, BLASTs them, fetches additional info from NCBI, and produces a
report.

positional arguments:
  filename              The filename of the fastq file you wish to sample and
                        BLAST

optional arguments:
  -h, --help            show this help message and exit
  --ascii64 ASCII64, -a ASCII64
                        Select true if Phred quality scores are encoded as
                        ASCII 64 (most are ASCII 33), default is False
  --nPercent NPERCENT, -np NPERCENT
                        A float between 0 and 100, this argument takes
                        precedence over nAbsolute, default is 0
  --nAbsolute NABSOLUTE, -na NABSOLUTE
                        An integer between 0 and the number of sequences in
                        your fastq file, this argument is superseded by
                        nPercent, default is 100
  --leadingQ LEADINGQ, -lq LEADINGQ
                        The minimum quality required to keep a base at the
                        leading end of a read, default is 20
  --trailingQ TRAILINGQ, -tq TRAILINGQ
                        The minimum quality required to keep a base at the
                        trailing end of a read, default is 20
  --hitlistSize HITLISTSIZE, -hs HITLISTSIZE
                        The number of blast hits to keep for the final report,
                        default is 1
  --taxID TAXID, -ID
                        The specific entrez species ID, to limit the blast search
                        to that organism, default is None
                        Example Species Entrez IDs:
                        chickens - 9031
                        humans - 9606
                        mouse - 10090
                        giant panda - 9646
  --database DATABASE, -db DATABASE 
                        NCBI database used in blast search,
                        default is nt 
                        Example Database Names from NCBI:
                        nt - nucleotide database (sequences from tons of databases)
                        est_human - human subset of EST database
                        est_mouse - mouse subset of EST database
                        Note: EST stands for "Expressed Sequence TAG"

                        There are many more databases that can be blasted and a
                        full breakdown of NCBI databases can be found at the following link:
                        https://www.ncbi.nlm.nih.gov/books/NBK62345/#_blast_ftp_site_The_blastdb_subdirectory/

required named arguments:
  --email EMAIL, -e EMAIL
                        A valid email address is required to use NCBI tools
                        and will be used if NCBI observes requests that
                        violate their policies.
</pre>




## Future Work
1. Wrap the BLAST function in a try/except block to make 2-3 attempts before giving up and throwing the error about the query being rejected.  
2. Use the final tabular report to calculate some summary statistics and maybe make some charts based on that output
3. Include info in tabular report for sequences that result in "NO HIT", such as trimmed sequence length. Right now it's just a catch-all.



