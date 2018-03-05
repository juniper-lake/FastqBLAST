#!/usr/bin/python3

#Version FastqBLAST_iv4.0.1.py

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"""
AUTHOR:         Danielle Novick
DATE CREATED:   October 24, 2017
LAST UPDATE:    February 15, 2018
LAST MODIFIED WITH PERMISSION: March 2, 2018 by M. Joseph Tomlinson IV
MODIFICATIONS: Created Summary Report files, Built in ability to Change Databases and Perform Organism Searches
               Split the Blast results into two files (Hit vs. No Hits), Fixed some minor issues, Fix a big bug with 
               gene duplicates in data not being reported properly

OBJECTIVE:      This script takes a sample of sequences from a fastq file, trims the low quality ends, BLASTs them,
                fetches additional info from NCBI, and produces a report.
NCBI's BLAST Usage Guidelines
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
NCBI's Databases 
https://www.ncbi.nlm.nih.gov/books/NBK62345/#_blast_ftp_site_The_blastdb_subdirectory_
BioPython's Manual
http://biopython.org/DIST/docs/tutorial/Tutorial.html
"""

import sys
import random
import argparse
import time
import numpy as np
import time
import os
from collections import defaultdict
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
from urllib.error import HTTPError


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - G L O B A L  D E C L A R A T I O N S  - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

entryLine = "\n\n* * * * * alpha snail on the trail * * * * *\n\n"
middleLine = "\n\n* * * * * gamma salamander in the alley * * * * *\n\n"
exitLine = "\n\n* * * * * beta toad on the road * * * * *\n\n"


parser = argparse.ArgumentParser(description='This script takes a sample of sequences from a fastq file, trims the \
                                low quality ends, BLASTs them, fetches additional info from NCBI, and produces a report.')

# positional arguments
parser.add_argument('filename', action="store", type=str,
                     help='The filename of the fastq file you wish to sample and BLAST')

# required named arguments
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('--email','-e', action="store", type=str, required=True,
                           help='A valid email address is required to use NCBI tools and will be used if NCBI observes \
                           requests that violate their policies.')

# optional arguments
parser.add_argument('--ascii64','-a', action="store", default=False, type=bool,
                     help='Select true if Phred quality scores are encoded as ASCII 64 (most are ASCII 33), '
                          'default is False')
parser.add_argument('--nPercent','-np', action="store", default=0, type=float,
                     help='A float between 0 and 100, this argument takes precedence over nAbsolute, default is 0')
parser.add_argument('--nAbsolute','-na', action="store", default=100, type=int,
                     help='An integer between 0 and the number of sequences in your fastq file, this argument is '
                          'superseded by nPercent, default is 100')
parser.add_argument('--leadingQ','-lq', action="store", default=20, type=int,
                     help='The minimum quality required to keep a base at the leading end of a read, default is 20')
parser.add_argument('--trailingQ','-tq', action="store", default=20, type=int,
                     help='The minimum quality required to keep a base at the trailing end of a read, default is 20')
parser.add_argument('--hitlistSize','-hs', action="store", default=1, type=int,
                     help='The number of blast hits to keep for the final report, default is 1')
parser.add_argument('--database','-db', action="store", default='nt', type=str,
                     help='Please visit NCBIs website for various database names, default is nt (nucleotide sequence database)')
parser.add_argument('--taxID','-ID', action="store", default='', type=str,
                     help='Please visit NCBIs website for various species taxonomy IDs, default is none')

args = parser.parse_args()

Entrez.email = args.email



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def random_sample(fastq_filename, n_absolute, n_percent):
    """
    Counts the number of lines/sequences in a FASTQ file and then takes a sample that can be entered either as a percent
    or an absolute number (the percent supersedes the absolute if both are entered).
    :param fastq_filename: The name of a fastq file
    :param n_absolute: An absolute number of samples that will be taken, superseded by n_percent if that's larger than 0
    :param n_percent: A percent of samples that will be taken, if this is specified it will take precedence over n_absolute
    :return: a list of integers that correspond to the first line of each sequence that is in the sample
    """
    summary_report = open('summary_blast_report.txt', 'w')


    with open(fastq_filename) as file:
        for counter, line in enumerate(file, 1):  # start index at 1 instead of 0
            pass
    num_lines = counter
    num_sequences = int(num_lines/4)
    print("You have %s sequences in your FASTQ file.\n" % "{:,}".format(num_sequences))
    sample_size = n_absolute if n_percent == 0 else int(n_percent*num_sequences/100)
    if sample_size > num_sequences:
        print("Your sample size is greater than your population size. Please select a number less than %s or consider "
              "using the percentage parameter." % "{:,}".format(num_sequences))
        sys.exit()
    elif sample_size == 0:
        print("You're trying to take a sample of size 0. Please adjust either nPercent or nAbsolute to correct the issue.")
        sys.exit()
    else: print("Sampling %s sequences from your file..." % "{:,}".format(sample_size))
    sample_set = [x * 4 for x in random.sample(range(0,int(num_lines/4)), sample_size)]

    
    summary_report.write("You have %s sequences in your FASTQ file.\n" % "{:,}".format(num_sequences))
    summary_report.write("\n")
    summary_report.write("Sampling %s sequences from your file..." % "{:,}".format(sample_size))
    summary_report.write("\n")
    summary_report.close()

    return sample_set


def fastq_to_dict(fastq_filename, sample_list):
    """
    Uses a list of sequence starting lines to pull out sequences from a fastq file and stores them in a dictionary,
    then decodes the phred scores and stores those as well
    :param fastq_filename: The name of a fastq file
    :param sample_list: A list of integers that correspond to the first line of each sequence that is in the sample
    :return: a two-level defaultdict {header:{'sequence': }, {'ascii':}, {'phred'}}
    """
    sample_dict = defaultdict(lambda: defaultdict())
    with open(fastq_filename) as file:
        for counter, line in enumerate(file):
            if counter in sample_list:
                header = line.rstrip().split("\t")[0]
            elif (counter - 1) in sample_list:
                sample_dict[header]['sequence'] = line.rstrip()
            elif (counter - 3) in sample_list:
                sample_dict[header]['ascii'] = line.rstrip()
    base = 33 if args.ascii64 == False else 64
    for key in sample_dict.keys():
        sample_dict[key]['phred'] = [ord(x) - base for x in list(sample_dict[key]['ascii'])]
    return sample_dict


def trim_ends(sample_dictionary, leadingQthreshold, trailingQthreshold):
    """
    Trims low quality bases from the leading and trailing ends of sequences
    :param sample_dictionary:   a two-level defaultdict {header:{'sequence': }, {'ascii':}, {'phred'}}
    :param leadingQthreshold:   the minimum quality required to keep a base at the leading end of a read
    :param trailingQthreshold:  the minimum quality required to keep a base at the trailing end of a read
    :return: the parameter sample_dictionary, but with two additional keys (trimmed_phred and trimmed_sequence)
    """
    print("Trimming the low-quality ends...")
    sample_dict = sample_dictionary.copy()
    below_blasting_threshold = 0
    #introduced verdict when sequence is completely trimmed
    for key in sample_dict.keys():
        for base, Q in enumerate(sample_dict[key]['phred']):
            if Q < leadingQthreshold:
                continue
            else:
                sample_dict[key]['trimmed_phred'] = sample_dict[key]['phred'][base:]
                sample_dict[key]['trimmed_sequence'] = sample_dict[key]['sequence'][base:]
                break
        for base, Q in reversed(list(enumerate(sample_dict[key]['phred']))):
            if Q < trailingQthreshold:
                continue
            else:
                sample_dict[key]['trimmed_phred'] = sample_dict[key]['trimmed_phred'][:base+1]
                sample_dict[key]['trimmed_sequence'] = sample_dict[key]['trimmed_sequence'][:base+1]
                break

    return (sample_dict)


def write_fasta(sample_dictionary):
    """
    Writes a FASTA file with the sequence IDs and trimmed sequences from the sample dictionary
    :param sample_dictionary:  a two-level defaultdict with information about the sequences to be BLASTed
    :return: blast_queries.fasta
    """
    OUT = open('Log_Directory/blast_queries.fasta', 'w')
    for key in sample_dictionary:
        # failsafe for sequences that are trimmed to be 0 bases long
        if 'trimmed_sequence' in  sample_dictionary[key].keys():
            OUT.write('>' + key[1:] + '\n' + sample_dictionary[key]['trimmed_sequence'] + '\n')
    OUT.close()


def blast_reads(number_hits, ncbi_database, organism):
    #blast_reads(number_hits, ncbi_database, entrez_query)
    """
    Uses Biopython's qblast() to BLAST sequences from a FASTA file, then write the blast results to a file
    :param number_hits: The maximum number of hits to return for each BLAST query sequence
    :return: blast_results.xml
    """
    print("Searching for BLAST hits...")
    fasta_string = open("Log_Directory/blast_queries.fasta").read()
    print ("The ncbi database being searched is:", ncbi_database)
    if len(organism) > 0:
        print ("The organism being searched is: ", organism)
        query ='"txid'+str(organism)+'"'
        result_handle = NCBIWWW.qblast("blastn", ncbi_database, fasta_string, entrez_query=query, hitlist_size=number_hits,
            expect=10.0, nucl_penalty=-2, nucl_reward=1, megablast=True, word_size=28, expect_low=True, gapcosts='0 2')
    else:
        print ("No organism is designated")
        result_handle = NCBIWWW.qblast("blastn", ncbi_database, fasta_string, hitlist_size=number_hits)
    blast_result = open("Log_Directory/blast_results.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()


def blast_to_dict():
    """
    Parses BLAST results and stores useful information in a dictionary. Throws an error if the blast_results.xml file
    only contains the queries with no results, which is an indicator that the BLAST was rejected by NCBI
    :return: a two-level defaultdict with information from the BLAST results and a flat list of the genes identified by BLAST
    """
    print("Parsing the BLAST results...")
    GeneIDs = []
    blast_dict = defaultdict(lambda: defaultdict())
    for record in NCBIXML.parse(open("Log_Directory/blast_results.xml")):
        for align in record.alignments:
            for hsp in align.hsps:
                percent_identity = round(100 * float(hsp.identities) / float(hsp.align_length),2)  # https://www.dnastar.com/megalign_help/index.html#!Documents/calculationofpercentidentity.htm
                hit_id = align.title.split('|')
                # this uses NCBI's gi number (GenInfo Identifier) which is reliable now but getting phased out, so might
                # need to change to hit_id[3] at some point
                GeneIDs.append(hit_id[1])
                blast_dict[record.query]['Hit_ID'] = align.title
                blast_dict[record.query]['Gene_ID'] = hit_id[1]
                blast_dict[record.query]['Sequence'] = hsp.query
                blast_dict[record.query]['SeqLength'] = len(hsp.query)
                blast_dict[record.query]['Description'] = hit_id[4]
                blast_dict[record.query]['Accession'] = hit_id[3]
                blast_dict[record.query]['Db'] = hit_id[2]
                blast_dict[record.query]['Score'] = hsp.score
                blast_dict[record.query]['E_value'] = hsp.expect
                blast_dict[record.query]['Percent_Identity'] = percent_identity
    GeneIDs = list(set(GeneIDs))
    if not GeneIDs:
        print('\nYour BLAST query was rejected. Please enter a smaller sample size or try running this script \
              at a better time.\nNCBI asks that you run scripts on weekends or between 9pm and 5am Eastern \
              time on weekdays if more than 50 searches will be submitted.')
        sys.exit()

    return blast_dict, GeneIDs,


def fetch_gene_info(gene_list, batch_size=100):
    """
    Uses an NCBI tool called efetch to look up more information about the genes identified by BLAST, then writes
    the results to a file. Epost is used here as good practice for large submissions to efetch
    :param gene_list: a list of NCBI gi's that will be submitted to efetch
    :param batch_size: the size of batches of gi's that get submitted to efetch to prevent overloading it
    :return: fetch_results.txt
    """
    print("Looking up additional information about the genes identified by BLAST...")
    post_handle = Entrez.epost(db="nucleotide", id=",".join(gene_list))
    result = Entrez.read(post_handle)
    post_handle.close()
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    count = len(gene_list)
    OUT = open("Log_Directory/fetch_results.txt", "w")
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Fetching records %i through %i" % (start + 1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start, retmax=batch_size,
                                            webenv=webenv, query_key=query_key)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise
        OUT.write(fetch_handle.read())
        fetch_handle.close()
    OUT.close()


def fetch_to_dict(blast_dictionary):
    """
    Stores information from an efetch results file into a dictionary and then merges that dictionary with the BLAST
    results dictionary
    :param blast_dictionary: a two level dictionary containing results from a BLAST search
    :return: the blast_dictionary parameter, but with additional keys (Organism, Source, Domain, Taxonomy)
    """
    blast_dict = blast_dictionary.copy()
    fetch_dict = defaultdict(lambda: defaultdict())
    for record in SeqIO.parse("Log_Directory/fetch_results.txt", "genbank"):
        fetch_dict[record.id]['Organism'] = record.annotations['organism']
        fetch_dict[record.id]['Source'] = record.annotations['source']
        fetch_dict[record.id]['Domain'] = record.annotations['taxonomy'][0]
        fetch_dict[record.id]['Taxonomy'] = record.annotations['taxonomy']
    for record in blast_dict.keys():
        for accession in fetch_dict.keys():
            if accession in blast_dict[record]["Accession"]:
                for accession_item in next(iter(fetch_dict.values())).keys():
                    blast_dict[record][accession_item] = fetch_dict[accession][accession_item]

    return blast_dict


def tabular_report(sample_dictionary, blast_dictionary):
    """
    Writes a report about the sample set using information from BLAST and eFetch
    :param sample_dictionary: a two-level defaultdict with information about the sequences from a fastq file
    :param blast_dictionary: a two level dictionary containing results from a BLAST search and eFetch
    :return: blast_hits_report.txt and blast_no_hits_report.txt
    """
    print("Writing the report...")
    sample_dict = sample_dictionary.copy()
    blast_dict = blast_dictionary.copy()

    #creating quick dictionary to pull out trimmed sequence
    trimmed_data_dict={}
    for key in sample_dict.keys():
        try:
            trimmed_sequence=(sample_dict[key]['trimmed_sequence'])
            key = key.strip('@')
            trimmed_data_dict.update({key:trimmed_sequence})
        except KeyError:
            continue

    samples = []
    for sequenceID in sample_dict:
        samples.append(sequenceID[1:])
    
    records = blast_dict.keys()
    
    columns = ["SeqID", "Trimmed Sequence", "Trimmed Sequence Length","BLAST Sequence", 
    "BLAST SeqLength", "Description", "Accession", "Db", 
    "Score", "E_value", "Percent_Identity", "Organism", 
    "Source", "Domain", "Taxonomy"]
    
    #Writing
    OUT = open("blast_hits_report.txt", "w")
    OUT.write('\t'.join(columns) + '\n')
    
    NO_HITS_OUT = open("blast_no_hits_report.txt", "w")
    NO_HITS_OUT.write("SeqID\tOriginal Seq\tOriginal Seq Length"
        "\tTrimmed Seq\tTrimmed Seq Length\tResult\n")

    for record in blast_dict.keys():

        trimmed_sequence = trimmed_data_dict[record]
        length_trimmed_sequence=len(trimmed_data_dict[record])
        
        #Used Brute force coding to be able to manipulate and add new columns to output
        try: 
            OUT.write(str(record)
                +'\t'+str(trimmed_sequence)
                +'\t'+str(length_trimmed_sequence)
                +'\t'+str(blast_dict[record]['Sequence'])
                +'\t'+str(blast_dict[record]['SeqLength'])
                +'\t'+str(blast_dict[record]['Description'])
                +'\t'+str(blast_dict[record]['Accession'])
                +'\t'+str(blast_dict[record]['Db'])
                +'\t'+str(blast_dict[record]['Score'])
                +'\t'+str(blast_dict[record]['E_value'])
                +'\t'+str(blast_dict[record]['Percent_Identity'])
                +'\t'+str(blast_dict[record]['Organism'])
                +'\t'+str(blast_dict[record]['Source'])
                +'\t'+str(blast_dict[record]['Domain'])
                +'\t'+str(blast_dict[record]['Taxonomy'])+'\n')
        except KeyError:
            continue

    for sample in samples:

        if sample not in records:
            sample_stripped = sample.split("\t")[0]

            #Get original trimmed sequence for reference
            try:
                trimmed_sequence = trimmed_data_dict[sample_stripped]
                length_trimmed_sequence=len(trimmed_sequence)
            except KeyError:
                trimmed_sequence = '-'
                length_trimmed_sequence=0

            NO_HITS_OUT.write(sample_stripped
                #Commend out original data being reported
                + '\t' + sample_dict['@'+sample]['sequence'] 
                + '\t' + str(len(sample_dict['@'+sample]['sequence'])) 
                + '\t' + str(trimmed_sequence)
                + '\t' + str(length_trimmed_sequence)

                + '\t' + 'NO HIT OR SEQUENCE QUALITY BELOW THRESHOLD\n')
    OUT.close()
    NO_HITS_OUT.close()


def parsing_no_hits_data(global_avg_trimmed_length):
    """ Analyzes the blast_no_hits_report.txt and creates summary reports of the final data
     :param start_time (see how long program takes to run)
     :return: counters (no_hit_counter, totally_trimmed_counter) and lists (average_trimmed_no_hit_lenght, 
        global_avg_trimmed_lenght)
     """

    #No Hit Counter
    no_hit_counter = 0

    #Totally trimmed counter
    totally_trimmed_counter = 0

    #No hits results
    average_trimmed_no_hit_length=[]

    #Opening and Parsing blast_no_hits_report.txt
    no_hit_results = open('blast_no_hits_report.txt', 'r')
    for line in no_hit_results:
        data = line.split("\t")
        
        if line.startswith('SeqID'):
            continue
        else:
            average_trimmed_no_hit_length.append(float(data[4]))
            global_avg_trimmed_length.append(float(data[4]))
            
            no_hit_counter +=1
            
            if float(data[4]) == 0:
                totally_trimmed_counter +=1
            continue
    no_hit_results.close


    return {'no_hit_counter':no_hit_counter, 'totally_trimmed_counter':totally_trimmed_counter, 
    'average_trimmed_no_hit_length':average_trimmed_no_hit_length, 'global_avg_trimmed_length':global_avg_trimmed_length}

def parsing_hits_data():
    """ Analyzes the blast_hits_report.txt and creates summary stats of the final data
     :param none (reads file directly)
     :return: numerous counters and lists for final analysis
     """
    #counters
    blast_hit_counter = 0
    
    #ALL DATA Results
    global_avg_trimmed_length=[]
    
    #Only hits results
    hits_avg_trimmed_length=[]
    hits_avg_blast_length=[]
    hits_avg_score=[]
    hits_avg_percent_identity=[]

    #Key word counters
    predicted_counter=0

    #creating parsing dictionary for program (hits only)
    blast_hit_dict = {}
    
    #Opening and Parsing blast_report.txt
    blast_hit_results = open('blast_hits_report.txt', 'r')

    for line in blast_hit_results:
        data = line.split("\t")
   
        if line.startswith('SeqID'):
            continue

        else:
            blast_hit_counter +=1

            #See How Many Genees are Predicted
            gene_description=(data[5]).lstrip(' ')
            if gene_description.startswith('PREDICTED'):
                predicted_counter += 1
            
            #Trimmed Sequence Stats
            global_avg_trimmed_length.append(float(data[2]))

            #Hits Stats
            hits_avg_trimmed_length.append(float(data[2]))
            hits_avg_blast_length.append(float(data[4]))
            hits_avg_score.append(float(data[8]))
            hits_avg_percent_identity.append(float(data[10]))
   
            #Test to see if organism in dictionary
            verdict = blast_hit_dict.get(data[11])
            
            #If not in 
            if str(verdict) == "None":
                #creating new entry
                key = data[11]
                #Value[Counts, Trimmed_Length, Blast Length, Blast_Score, Blast_Percent_Identity]
                value=[1, [float(data[2])], [float(data[4])], [float(data[8])], [float(data[10])] ]
                blast_hit_dict.update({key:value})
            else:
                #Fills dictionary based on organism name
                (blast_hit_dict[data[11]][0])+=1
                (blast_hit_dict[data[11]][1]).append(float(data[2]))
                (blast_hit_dict[data[11]][2]).append(float(data[4]))
                (blast_hit_dict[data[11]][3]).append(float(data[8]))
                (blast_hit_dict[data[11]][4]).append(float(data[10]))

    blast_hit_results.close()

    return {'blast_hit_counter': blast_hit_counter,  'global_avg_trimmed_length': global_avg_trimmed_length,
    'hits_avg_trimmed_length': hits_avg_trimmed_length, 'hits_avg_blast_length': hits_avg_blast_length,
    'hits_avg_score':hits_avg_score, 'hits_avg_percent_identity': hits_avg_percent_identity,
    'predicted_counter': predicted_counter, 'blast_hit_dict':blast_hit_dict}


def summary_blast_report(start_time):
    """ Prints the final analysis reports from  the blast_report.txt and blast_no_hits_report.txt 
     :param start_time (see how long program takes to run)
     :return: summary_blast_report.txt
     """

    #Getting Hit Results from File
    hit_data=parsing_hits_data()
    blast_hit_counter=hit_data['blast_hit_counter']
    global_avg_trimmed_length=hit_data['global_avg_trimmed_length']
    hits_avg_trimmed_length=hit_data['hits_avg_trimmed_length']
    hits_avg_blast_length=hit_data['hits_avg_blast_length']
    hits_avg_score=hit_data['hits_avg_score']
    hits_avg_percent_identity=hit_data['hits_avg_percent_identity']
    predicted_counter=hit_data['predicted_counter']
    blast_hit_dict=hit_data['blast_hit_dict']

    #Getting No Hit Results from File
    no_hit_data= parsing_no_hits_data(global_avg_trimmed_length)
    no_hit_counter=no_hit_data['no_hit_counter']
    totally_trimmed_counter=no_hit_data['totally_trimmed_counter']
    average_trimmed_no_hit_length=no_hit_data['average_trimmed_no_hit_length']
    global_avg_trimmed_length=no_hit_data['global_avg_trimmed_length']

    total_counts= blast_hit_counter + no_hit_counter
   
    #Printing all final results to output files
    summary_report = open('summary_blast_report.txt', 'a')
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("The number of sequences analyzed was: " + str(total_counts)+"\n")
    summary_report.write("The average trimmed length of ALL sequences: " 
        + str(round((np.average(global_avg_trimmed_length)),2))+"\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("Sequences with NO BLAST Results\n")
    summary_report.write("The number of sequences with no results: " + str(no_hit_counter)+ "\n")
    #Else/if for trimmed sequences when none exist in the file
    if no_hit_counter == 0:
        summary_report.write("The average sequence length was: N.A." + "\n")
        summary_report.write("Number of sequences totally trimmed: N.A." + "\n")

    else:
        summary_report.write("The average sequence length was: " 
            + str(round((np.average(average_trimmed_no_hit_length)),2)) +"\n")
        summary_report.write("Number of sequences totally trimmed: " + str(totally_trimmed_counter) + "\n")
    
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("Sequences with BLAST Hit Results\n")
    summary_report.write("The number of sequences with BLAST hit results: " + str(blast_hit_counter)+"\n")
    summary_report.write("The number of genes described as PREDICTED were: " + str(predicted_counter) + "\n")
    summary_report.write("The average trimmed sequence length (pre-BLASTING): " 
        + str(round((np.average(hits_avg_trimmed_length)),2))+"\n")
    summary_report.write("The average sequence BLAST hit length was: " 
        + str(round((np.average(hits_avg_blast_length)),2)) + "\n")
    summary_report.write("The average blast hits score was: "
                         + str(round((np.average(hits_avg_score)),2))+ "\n")
    summary_report.write("The average percent identity was: "
                         + str(round((np.average(hits_avg_percent_identity)),2))+ "\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("Organism\tCounts\tTrimmed_Seq_Length\tHits_Avg_Seq_Length\t"
        "Hits_Avg_Score\tHits_Avg_Percent_Identity\n")
                    
    for data in blast_hit_dict:
        organism =(data)
        counts_of_records = (blast_hit_dict[data][0])
        overall_trimmed_length = round((np.average(blast_hit_dict[data][1])),2) 
        overall_avg_length = round((np.average(blast_hit_dict[data][2])),2)
        overall_avg_score = round((np.average(blast_hit_dict[data][3])),2)
        overall_avg_identity = round((np.average(blast_hit_dict[data][4])),2)

        summary_report.write(str(organism)+"\t"+str(counts_of_records)+"\t"
            +str(overall_trimmed_length)+"\t" 
            +str(overall_avg_length)+"\t"
            +str(overall_avg_score)+"\t" 
            +str(overall_avg_identity)+"\n")
    

    summary_report.write("\n")
    summary_report.write("\n")
    summary_report.write("\n")

    #Started a time in the program to see how long runs take
    total_time = time.clock() - start_time
    summary_report.write("Program ran for a total of " + str(round(total_time, 2))+" seconds")

    summary_report.close()

def tallying_genes():
    """ Parses the blast_report.txt to tally and count genes that showed up more than once.
        Specifically searches on accession id
        param: none
        return: gene dict (dictionary based on accession ids with corresponding information
     """
    #Creating a tallying Mechanism of genes with multiple sequences in file and
    # an output file for future alignment of sequences 
    blast_hit_results = open('blast_hits_report.txt', 'r')
    gene_dict={}

    for line in blast_hit_results:
        data = line.split("\t")
       
        if line.startswith('SeqID'):
            continue
        else:
            #Test to see if organism in dictionary
            verdict = gene_dict.get(data[6])
                
            if str(verdict) == "None":
                #creating new entry
                key = data[6]
                seq_info=str(data[0])+"|"+str(data[1])
                counter = 1
                #Value[Counts, Trimmed_Length, Blast Length, Blast_Score, Blast_Percent_Identity]
                value=[data[5], counter, [seq_info]]
                gene_dict.update({key:value})
            else:
                #Fills dictionary based on organism name
                seq_info=str(data[0])+"|"+str(data[1])
                gene_dict[data[6]][1]+=1
                gene_dict[data[6]][2].append(seq_info)
    blast_hit_results.close()
    return(gene_dict)


def printing_summary_gene_report(gene_dict):
    """ Parses the gene dict and prints the final summary report
    param: gene_dict
    return: none
     """
    #Creating a summary report for data
    summary_gene_report = open('summary_gene_report.txt', 'w')
    summary_gene_report.write("Accesssion ID\tDescription\tCounts\n")
    for key in gene_dict:
        accession_ID=str(key)    
        gene_description=str(gene_dict[key][0])
        gene_counts=str(gene_dict[key][1])
        summary_gene_report.write(accession_ID+"\t"+gene_description+"\t"+gene_counts+"\n")
    summary_gene_report.close()


def printing_blat_searchable_data(gene_dict):
    """ Parses the gene dict and prints a report of genes and corresponding sequences for that gene
    that can blasted/blat to un-cover further information about the sequences and genes
    param: gene_dict
    return: none
     """
    #Creating a report of all sequences that can searched and then blasted
    blat_gene_report = open('Log_Directory/blat_gene_seq_report.txt', 'w')
    blat_gene_report.write("This report was created to allow a user to search specific groups of sequences\n")
    blat_gene_report.write("for a gene using either BLAST or UCSC Genome Browser to try and possibly identify\n")
    blat_gene_report.write("a feature that caused enrichment for that gene in the data (length, CNV, highly expressed etc)\n")
    blat_gene_report.write("\n")
    blat_gene_report.write("\n")
    for key in gene_dict:
        blat_gene_report.write("Gene\tDescription\tCounts\n")
        gene_description=str(gene_dict[key][0])
        gene_counts=str(gene_dict[key][1])
        accession_ID=str(key)
        gene_description=str(gene_dict[key][0])
        blat_gene_report.write(accession_ID+"\t"+gene_description+"\t"+gene_counts+"\n")

        blasted_sequences=gene_dict[key][2]
        for sequences in blasted_sequences:
            sequences=sequences.split("|")
            seq_ID=str(sequences[0])
            sequence=str(sequences[1])
            blat_gene_report.write(">"+seq_ID+"\n")
            blat_gene_report.write(sequence+"\n")

        blat_gene_report.write("\n")
        blat_gene_report.write("\n")
    blat_gene_report.close()

def main():
    print(entryLine)
    print ("")

    #Setting up program
    start_time=time.clock()
    home_directory = os.getcwd()

    #Making a directory to put all important, but non-essential results in
        # See if directory exists otherwise make it
    verdict = os.path.exists('Log_Directory')
    if str(verdict) == 'False':
        os.makedirs('Log_Directory')
    else:
        print("Log Directory already exists")
        print("")

    
    sample_set = random_sample(fastq_filename=args.filename, n_absolute=args.nAbsolute, n_percent=args.nPercent)
    sample_dict = fastq_to_dict(fastq_filename=args.filename, sample_list=sample_set)
    sample_dict = trim_ends(sample_dictionary=sample_dict,leadingQthreshold=args.leadingQ, trailingQthreshold=args.trailingQ)
    write_fasta(sample_dictionary=sample_dict)
    blast_reads(number_hits=args.hitlistSize, ncbi_database=args.database, organism=args.taxID)
    blast_dict, GeneIDs = blast_to_dict()
    fetch_gene_info(gene_list=GeneIDs)
    blast_dict = fetch_to_dict(blast_dictionary=blast_dict)
    tabular_report(sample_dictionary=sample_dict, blast_dictionary=blast_dict)

    print(middleLine)
    
    summary_blast_report(start_time)
    #Getting the results from blast_hits_report.txt
    gene_dict=tallying_genes()
    #Printing the summary reports of analysis
    printing_summary_gene_report(gene_dict)
    printing_blat_searchable_data(gene_dict)
    

    total_time = time.clock() - start_time
    
    print ("The program ran for ", total_time, " seconds")
    print ("")
    print(exitLine)


if __name__ == "__main__":
    main()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - E n d   o f   F i l e - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
