# FastqBLAST 
***
> **Author**: Danielle Novick  
> **Last Update**: February 15, 2017  
> **Creation Date**: October 24, 2017


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
                 [--trailingQ TRAILINGQ] [--hitlistSize HITLISTSIZE]
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

## Versions
**v1.0.1** Fixed program failure from key error when sequence was trimmed to length of 0 bases  
 
**v1.0.0** First release

