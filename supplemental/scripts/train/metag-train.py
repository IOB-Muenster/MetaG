#!/usr/bin/env python3

# AUTHORS

# Felix Manske, felix.manske@uni-muenster.de
# Marc-Nicolas Bendisch
# Norbert Grundmann
# Author for correspondence: Wojciech Makalowski, wojmak@uni-muenster.de


# COPYRIGHT

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO,  PROCUREMENT  OF  SUBSTITUTE GOODS  OR  SERVICES;
# LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


#==============================================================#
# Analyze the training results of MetaG
#--------------------------------------------------------------#
# 
# DESCRIPTION
#
# Analyze the raw output of MetaG training and calculate
# performance metrics: Bray-Curtis similarity index (BC)
# or Rk (multiclass generalization of Matthews's correlation
# coefficient).
#
#
# USAGE
# 
# Calculate performance metric:
# ./metag-train.py --query QUERY FILE --exp EXPECTED TAXONOMY \
#    --train TRAIN RESULTS --rank RANK --metric [rk, bc] \
#    --name ANALYSIS NAME
#
# QUERY FILE is the FASTA file used to generate the raw
# training results. EXPECTED TAXONOMY is a text file that
# provides the expected classifications. Its format depends on
# the metric (see examples).
#
#
# OUTPUT
#
# Settings and their Rk or BC values will be written to STDOUT.
#
#
# EXAMPLES
#
# The expectation file consists of three sections. A line starting
# with "@@ranks" followed by a tab and the ranks of the taxonomy
# divided by semicolon. The next section starts with the line
# "@@negative" and contains information about reads that should
# be unclassified. The next line just contains a number (BC metric)
# specifying the amount of unclassified reads or patterns that appear
# in the read headers and uniquely identify these reads (Rk metric).
# A similar structure applies to the "@@positive" section, that
# contains the expected classifications of classified reads.
# For the BC metric, each line starts with the expected count,
# followed by a tab and the expected taxonomy (divided by
# semicolon). For the Rk metric, each line starts with a pattern
# that appears in the read headers and uniquely identifies
# classifiable reads with a common classification. After a
# tab, the expected taxonomy is provided (see above).
# No empty lines are allowed in the files. Mock examples for both
# expectation files are provided below:
# BC metric:
# 
# @@ranks domain;phylum;class;order;family;genus;species
# @@negative
# 10
# @@positive
# 100    d1;p1;c1;o1;f1;g1;s1
# 100    d2;p2;c2;o2;f2;g2;s2
#
# Rk metric:
#
# @@ranks domain;phylum;class;order;family;genus;species
# @@negative
# shuffled
# @@positive
# id1    d1;p1;c1;o1;f1;g1;s1
# id2    d2;p2;c2;o2;f2;g2;s2
#
#
# DEPENDENCIES
#
# python3
#    matplotlib
#    numpy
#    pandas
#    scipy
#    sklearn
#==============================================================#


import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import pathlib as pl
import re
from scipy.spatial.distance import braycurtis
from sklearn.metrics import matthews_corrcoef
import sys


#-----------------------------------------------------------#
# Read CLI arguments
#-----------------------------------------------------------#
parser = argparse.ArgumentParser(description='Analyze MetaG training')
parser.add_argument('--query', type=str, required=True, help="File with input reads (required)")
parser.add_argument('--exp', type=str, required=True, help="File with expected classifications (required)")
parser.add_argument('--train', type=str, required=True, help="File with temporary training results (required)")
parser.add_argument('--rank', type=str, required=True, help="Rank to analyze (required)")
parser.add_argument('--metric', type=str, required=True, choices=["rk", "bc"], help="Analysis metric (required)")
parser.add_argument('--name', type=str, required=True, help="Name for the run (required)")
args = parser.parse_args()

# Check if files from arguments exist
for arg in (args.query, args.exp, args.train):
    path = pl.Path(arg)
    if not path.exists() or path.stat().st_size == 0:
        print("ERROR: The file " + arg + " is empty or does not exist.\n")
        sys.exit(1)
        break
            

#-----------------------------------------------------------#
# Read expected classifications of sequence IDs. Get read
# names for expected classifications.
#-----------------------------------------------------------#
isNegative = False
negPatterns = {}
posPatterns = {}
reads = {}
abunds = {}
totalAbund = 0


if (args.metric == "rk"):
    # Read patterns and expected taxonomy from expected file.
    # Patterns are used to get the read IDs from the query file.
    with open(args.exp, "r") as expec:
        for line in expec:
            # Ranks of the taxonomy
            if re.match('^@@ranks', line):
                line = line.strip("\n")
                rank = line.split("\t")[1]
                ranks = rank.split(";")
                
                # Index of rank to analyze
                idx = ranks.index(args.rank)
                continue
            # Unalignable reads
            elif re.match('^@@negative', line):
                isNegative = True
                continue
            # Identifiable reads
            elif re.match('^@@positive', line):
                isNegative = False
                continue
                
            if isNegative == True:
                line = line.strip("\n")
                negPatterns[line] = ""
            else:
                line = line.strip("\n")
                (pattern, lineage) = line.split("\t")
                lineage = lineage.split(";")[idx]
                if lineage == "0":
                    lineage = np.NaN
                    
                posPatterns[pattern] = lineage
    
    # Use patterns to get read IDs from query file.
    with open(args.query, "r") as query:
        isMatch = False
        for line in query:
            if (re.match('^>', line)):
                isMatch = False
                line = line.strip(">")
                line = line.strip("\n")
                
                for pattern in negPatterns.keys():
                    if (re.search(pattern, line)):
                        reads[line] = {'exp': np.NaN}
                        isMatch = True
                        break
                if isMatch == False:
                    for pattern in posPatterns.keys():
                        if (re.search(pattern, line)):
                            reads[line] = {'exp': posPatterns[pattern]}
                            isMatch = True
                            break    
                if isMatch == False:
                    print("ERROR: Unexpected read ID: " + line + "\n")
                    sys.exit(1)
elif args.metric == "bc" :
     with open(args.exp, "r") as expec:
        for line in expec:
            # Ranks of the taxonomy
            if re.match('^@@ranks', line):
                line = line.strip("\n")
                rank = line.split("\t")[1]
                ranks = rank.split(";")
                
                # Index of rank to analyze
                idx = ranks.index(args.rank)
                continue
            # Unalignable reads
            elif re.match('^@@negative', line):
                isNegative = True
                continue
            # Identifiable reads
            elif re.match('^@@positive', line):
                isNegative = False
                continue
                
            if isNegative == True:
                line = line.strip("\n")
                totalAbund += int(line)
                if "UNMATCHED" in abunds.keys():
                    abunds["UNMATCHED"]["exp"] += int(line)
                else:
                    abunds["UNMATCHED"] = {"exp": int(line)}
            else:
                line = line.strip("\n")
                (abundance, lineage) = line.split("\t")
                abundance = int(abundance)
                lineage = lineage.split(";")[idx]
                totalAbund += abundance
                if lineage in abunds.keys():
                    abunds[lineage]["exp"] += abundance
                else:
                    abunds[lineage] = {"exp": abundance}   
                    

#-----------------------------------------------------------#
# Check, how reads were classified by MetaG with the
# different parameter combinations.
#-----------------------------------------------------------#
settings = []
if (args.metric == "rk"):
    with open(args.train, "r") as trainRes:
        readId = ""
        setting = ""
        foundC = 0
        for line in trainRes:
            line = line.strip("\n")
            if re.match('^@@', line):
                setting = re.sub('^@@_', '', line, 1)
                settings.append(setting)
                continue
                
            if setting != "":
                if re.match('^>', line):
                    line = line.strip(">")
                    readId = line
                # Empty values in df (unmatched reads) are automatically entered
                # as NaN by pandas
                elif re.match('^No match', line):
                    readId = ""
                else:
                    if readId != "" and re.match('^' + args.rank + ':', line):
                        if readId in reads:
                            if setting in reads[readId]:
                                print("ERROR: The read " + readId + "has been classified twice in the same run\n")
                                sys.exit(1)
                            else:
                                lineage = line.split(": ")[1]
                                reads[readId][setting] = lineage
                                foundC += 1
                        else:
                            print("ERROR: Unexpected read ID in MetaG output:" + readId + "\n")
                            sys.exit(1)
                
        if foundC == 0 and setting == "":
            print ("ERROR: No setting found. Does your file have @@_* lines?\n")
            sys.exit(1)      
elif args.metric == "bc" :
    with open(args.train, "r") as trainRes:
        taxon = ""
        readId = ""
        setting = ""
        foundC = 0
        for line in trainRes:
            line = line.strip("\n")
            if re.match('^@@', line):
                setting = re.sub('^@@_', '', line, 1)
                settings.append(setting)
                # Init unmatched taxa with the total amount of reads.
                # Subtract indentified reads later.
                if "UNMATCHED" in abunds.keys():
                    abunds["UNMATCHED"][setting] = totalAbund
                else:
                    abunds["UNMATCHED"]= {setting: totalAbund}
                continue
                
            if setting != "":
                if re.match('^' + args.rank + ':', line):
                    lineage = line.split(": ")[1]
                    
                    # This taxon is expected
                    if lineage in abunds.keys():
                        if setting in abunds[lineage].keys():
                            abunds[lineage][setting] += 1
                        else:
                            abunds[lineage][setting] = 1
                    # An unexpected taxon. Pandas will
                    # fill abunds[lineage]["exp"] with NaN
                    else:
                        abunds[lineage] = {setting: 1}
                        
                    abunds["UNMATCHED"][setting] -= 1
                    foundC += 1
                    
        if foundC == 0 and setting == "":
            print ("ERROR: No setting found. Does your file have @@_* lines?\n")
            sys.exit(1)
                      
#-----------------------------------------------------------#
# Calculate performance metrics
#-----------------------------------------------------------#
if args.metric == "rk":
    # df: Row = Expected classification or parameter setting
    # df: Column = Read ID
    df = pd.DataFrame.from_dict(reads)
    
    # sklearn does not like NaN
    df = df.fillna("unclassified")

    expec = df.loc["exp"]
    
    rks = []
    for setting in settings:
        rk = matthews_corrcoef(expec, df.loc[setting])
        rk = round(rk, 4)
        rks.append((setting, rk))

    # Write output
    # Sort based on Rk (reverse: "-"(!) rk[1]) and setting
    rks = sorted(rks, key = lambda rk: (-rk[1], rk[0]), reverse = False)
    for tuple in rks:
        print (str(tuple[0]) + "\t" + str(tuple[1]))
        
elif args.metric == "bc" :
    # df: Row = Expected classification or parameter setting
    # df: Column = Taxon
    df = pd.DataFrame.from_dict(abunds)
    
    # Bray-Curtis cannot handle NaN
    df = df.fillna(0)
    
    expec = df.loc["exp"]
    bcs = []
    for setting in settings:
        bc = braycurtis(expec, df.loc[setting])
        # dissimilarity --> similarity index
        bc = 1 - round(bc, 4)
        bcs.append((setting, bc))
    
    # Write output
    # Sort based on Bray-Curtis (reverse: "-"(!) bc[1]) and setting
    bcs = sorted(bcs, key = lambda bc: (-bc[1], bc[0]), reverse = False)
    for tuple in bcs:
        print (str(tuple[0]) + "\t" + str(tuple[1]))
