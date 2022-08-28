#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber
from collections import defaultdict
import statistics
import re
import os
import textwrap
import argparse


parser = argparse.ArgumentParser(
    prog="ssuSilva.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    '''))

parser.add_argument('-blast_result', help='BLAST results of reads against the SILVA database. This is the output '
                                          '(format 6) of a BLAST of your raw reads against the SILVA nucleotide database. '
                                          'This input is optional. You can provide ssuSilva the raw reads and this program '
                                          'will run BLAST. To save time, I recommend providing only the subset of reads '
                                          'that you know are derived from the 16S rRNA gene. You can pull these out using a '
                                          'program like SortMeRNA: https://github.com/biocore/sortmerna', default="NA")
parser.add_argument('-reads', help="location of your raw reads (QC'ed, of course). Provide this file only if you did"
                                   " not provide the blast results file. If you would like ssuSilva to run BLAST, I "
                                   "recommend running this program using \'screen\' or \'nohup\', as it may take a while. "
                                   "To speed this up, you can allow up to 16 threads to be used for the BLAST command.", default="NA")
parser.add_argument('-silva_DB', help="location of SILVA database")
parser.add_argument('-t', help="number of threads to use for BLAST (default = 1)", default=1)
parser.add_argument('-qcov_hsp_perc', help="BLAST option: percent alignment coverage of query sequence (default = 100.00)", default=100.00)
parser.add_argument('-perc_identity', help="BLAST option: percent identity between query and target sequence (default = 99)", default=99.00)
parser.add_argument('-min_aln', help="minimum length of base pairs over which there is alignment (default = 50). Change this if your read length is different.", default=50)
parser.add_argument('-out', help="base name of output file (default = out)", default="out")

args = parser.parse_args()


def extender(listOfcoords):
    counter = 0
    max = 0
    outlist = []
    count = 0
    for i in listOfcoords:
        if counter == 0:
            counter += 1
            start = int(i[0])
            end = int(i[1])
        else:
            if int(i[0]) < end:
                if int(i[1]) > end:
                    end = int(i[1])
                    count = 0
            else:
                count += 1
                outlist.append([start, end])
                start = int(i[0])
                end = int(i[1])
    outlist.append([start, end])
    return outlist


def topRange(listOfRanges):
    max = 0
    for i in listOfRanges:
        length = i[1] - i[0]
        if length > max:
            max = length
            bestPair = [i[0], i[1]]
    return [max, bestPair]


def readCompile(listOfcoords, range):
    outList = []
    start = range[0]
    end = range[1]
    count = 0
    for i in listOfcoords:
        if i[0] >= start and i[1] <= end:
            outList.append(count)
            count += 1
        else:
            count += 1
    return outList


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


vList = [[0, 10000], [33, 99], [137, 242], [433, 497], [576, 682],
         [822, 879], [986, 1043], [1117, 1173], [1243, 1294], [1435, 1465]]

print("\nssuSilva up and running...")

if args.reads != "NA":
    print("Beginning BLAST of provided reads: " + args.reads + " against database: " + args.silva_DB +
          " with " + str(args.t) + " threads")
    os.system("blastn -query " + args.reads + " -db " + args.silva_DB + " -outfmt 6 -out tmp.blast -num_threads " +
              args.t + " -qcov_hsp_perc " + args.qcov_hsp_perc  + " -perc_identity " + args.perc_identity + " -max_target_seqs 1")
    print("BLAST finished")
else:
    print("Please provide the reads!")
print("\n")

count = 0
print("Reading database file: " + args.silva_DB)
silva = open(args.silva_DB, "r")
silvaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in silva:
    if re.match(r'>', i):
        count += 1
        print(count)
        ls = i.rstrip().split(" ")
        head = ls[0][1:]
        try:
            silvaDict[head] = ls[1]
        except IndexError:
            silvaDict[head] = head

SILVA = fasta(silva)


if args.blast_result == "NA":
    ssumap = open("tmp.blast", "r")
    print("Reading BLAST results")
else:
    ssumap = open(args.blast_result, "r")
    print("Reading BLAST results file: " + args.blast_result)

print("Creating a BLAST-generated alignment map\n")
MapDict = defaultdict(lambda: defaultdict(list))
matchedReads = 0
queryList = []
for i in ssumap:
    ls = i.rstrip().split("\t")
    query = ls[0]
    try:
        coords = [int(ls[8]), int(ls[9])]
        coords = sorted(coords)
        id = float(ls[2])
        match = ls[1]
        if query not in queryList:
            queryList.append(query)
            MapDict[match]["coords"].append(coords)
            MapDict[match]["id"].append(id)
            matchedReads += 1
    except IndexError:
        pass


outfile = open(args.out + ".csv", "w")
outfile.write("DB match" + "," + "overall_percent_match" + "," + "proportion_of_reads_mapped" + "," +
              "total_reads_mapped" + "," + "range_contiguous_16s_coverage" + "," + "hypervariable_regions_cov" +
              "," + "total_aln_length" + "," + "mean_read_depth" + "," + "stdev_read_depth" + "\n")
count = 0
for i in MapDict.keys():
    bpList = []
    totalRange = 0
    Allrange = (extender(sorted(MapDict[i]["coords"])))
    best = topRange(Allrange)
    aln = best[0]
    bestRange = best[1]
    print("Processing " + silvaDict[i])
    VarList = []
    if int(bestRange[1]-bestRange[0]) > int(args.min_aln):
        outfile.write(silvaDict[i] + "," + str(statistics.mean(MapDict[i]["id"])) + "," + str(len(MapDict[i]["id"])/matchedReads) + "," + str(len(MapDict[i]["id"])) + ",")
        if len(Allrange) > 1:
            for j in Allrange:
                totalRange += j[1] - j[0] + 1
                outfile.write(str(j[0]) + "-" + str(j[1]) + "; ")
                vregions = readCompile(vList, j)
                for v in vregions:
                    VarList.append(v)
            outfile.write(",")
            for var in VarList:
                outfile.write("V" + str(var) + "; ")
            outfile.write(",")
        else:
            totalRange += (Allrange[0][1] - Allrange[0][0] + 1)
            outfile.write(str(Allrange[0][0]) + "-" + str(Allrange[0][1]) + ",")
            vregions = readCompile(vList, Allrange[0])
            for v in vregions:
                outfile.write("V" + str(v) + "; ")
            outfile.write(",")
        outfile.write(str(totalRange) + ",")
        ran = (MapDict[i]["coords"])
        depthList = []
        for Ran in Allrange:
            for bp in range(Ran[0], Ran[1]):
                counter = 0
                for read in ran:
                    if int(bp) >= int(read[0]) and int(bp) <= int(read[1]):
                        counter += 1
                depthList.append(counter)
        outfile.write(str(statistics.mean(depthList[int(len(depthList)*0.1):int(len(depthList)*0.9)])) + ",")
        outfile.write(str(statistics.stdev(depthList[int(len(depthList)*0.1):int(len(depthList)*0.9)])) + ",")
        outfile.write(SILVA[i] + "\n")
        count += 1


inFile = open(args.out + ".csv", "r")



print("Done!")