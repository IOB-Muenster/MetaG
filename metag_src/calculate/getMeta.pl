#!/usr/bin/env perl

# AUTHORS

# Felix Manske, felix.manske@uni-muenster.de
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


# Give perl the locations of modules. Will look in production directory or cwd.
if (-d '/bioinf/projects/metag/calculate') {
	BEGIN {push @INC , '/bioinf/projects/metag/calculate'}
}
else {
	use FindBin;
	FindBin::again();
	use lib $FindBin::Bin;
}

use strict;
use warnings;

# Custom modules
use modules::ReadQuery;
use modules::ReadDB;
use modules::ProcessAlign qw(readAlign processAlign);
use modules::GetSpec qw(getSpec visOut getStat getUnmatchedSeq);
use modules::GetPatho;

# Use unbuffered stdout
select(STDOUT);
$| = 1;

# Get runtime
sub runTime {
	my $prevTime = $_[0];
	
	my $now = time;
	my $runTime = $now - $prevTime;
	
	return $runTime;
}

#===============================================================================================================#
# Grep query file, alignment file and parameters from STDIN
#---------------------------------------------------------------------------------------------------------------#

my $argLen = @ARGV;
print "\nStarted $0 with parameters: @ARGV\n";
my $queryFile = undef;
my $alignFile = undef;
my $eCut = undef;
my $alignThresh = undef;
my $confCut = undef;
my $type = undef;
my $dbPath = undef;
my $pdbPath = undef;
my $krona = undef;
my $isVirus = 0;

for (my $i=0; $i <= $#ARGV; $i++) {
	if ($argLen == 9 or $argLen == 10) {
		$queryFile = $ARGV[0];
		$alignFile = $ARGV[1];
		$eCut = $ARGV[2];
		$alignThresh = $ARGV[3];
		$confCut = $ARGV[4];
		$type = $ARGV[5];
		$dbPath = $ARGV[6];
		$pdbPath = $ARGV[7];
		$krona = $ARGV[8];
		$isVirus = $ARGV[9];
		last;
	}
	elsif ($argLen == 18 or $argLen == 19) {
		if ($ARGV[$i] eq "-q") {
			$queryFile = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-a") {
			$alignFile = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-e") {
			$eCut = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-ac") {
			$alignThresh = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-cc") {
			$confCut = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-m") {
			$type = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-db") {
			$dbPath = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-pdb") {
			$pdbPath = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-kro") {
			$krona = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-vir") {
			$isVirus = 1;
		}
	}
	else {
		if ($ARGV[$i] eq "-h") {
			print "A script to determine the species of a metagenomic sample\n";
			print "Arguments needed:\n";
			print "(-q) Query fasta file\n";
			print "(-a) Alignment maf file\n";
			print "(-e) E-value threshold\n";
			print "(-ac) Alignment score cutoff\n";
			print "(-cc) Confidence cutoff\n";
			print "(-m) Method for average confidence calculation\n\twilliams: Williams' mean\n\tgeometric: Geometric mean\n\tarithmetic: Arithmetic mean\n\tharmonic: Harmonic mean (sugg.)\n";
			print "(-db) Path to restructured database file\n";
			print "(-pdb) Path to pathogen database file\n";
			print "(-kro) Path to the Krona main folder\n";
			print "(-vir) Tag to indicate analysis using viral database\n\n";
			exit(0);
		}
		else {
			die "Expected 9 or 19 arguments, $argLen given. Use -h for help.";
		}
	}
}
if (defined $queryFile and defined $alignFile and defined $eCut and defined $alignThresh and defined $confCut and defined $type and defined $dbPath and defined $pdbPath and defined $krona) {
	print "\nQuery file: $queryFile";
	print "\nAlignment file: $alignFile\n";
	print "\nE-value cutoff: $eCut";
	print "\nAlignm. score cutoff: $alignThresh";
	print "\nConfidence cutoff: $confCut";
	print "\nMethod for average confidence calculation: $type\n";
	print "\nDatabase taxonomy file: $dbPath\n";
	print "\nPathogen database file: $pdbPath\n";
	
	if ($isVirus == 1) {
		print "\nPath to the Krona main folder: $krona\n";
		print "\nViral database: True\n\n";
	}
	else {
		print "\nPath to the Krona main folder: $krona\n\n";
	}
}
else {
	die "Undefined essential parameters (-q, -a, -e, -ac, -cc, -m, -db, -pdb, -kro). Use -h for help.\n"
}

#Measure runtime
my $start = time;


#my $outPrefix = (split('\.', $queryFile))[0]."_E".$eCut."_AC".$alignThresh."_CC".$confCut;
my $outPrefix = $queryFile =~ s/[^\/]*$/calc./r;
$eCut = 1*exp($eCut);


#Read queries
print "Reading queries: ".runTime($start)."\n";
my ($queriesR, $allQ, $queryHref) = readQuery($queryFile);

#Read 16S/18S DB
print "Reading database: ".runTime($start)."\n";
my ($dbHashR, $ranksR) = readDB($dbPath);

#Read MAF alignment
print "Reading alignment: ".runTime($start)."\n";
my $seqHashR = readAlign($alignFile);

#Apply e-value and alignment score cutoffs
print "Processing alignment: ".runTime($start)."\n";
($seqHashR, my $lostEval) = processAlign($seqHashR, $eCut, $alignThresh);

#Assign species to query reads and write visual output file and assignment file for each read
print "Assigning species: ".runTime($start)."\n";
my ($resultsR, $pathosR, $unmatchedR) = getSpec($ranksR, $outPrefix, $queriesR, $seqHashR, $dbHashR, $confCut, $isVirus);

#Print unmatched read IDs per rank
print "Getting unmatched sequences: ".runTime($start)."\n";
getUnmatchedSeq($unmatchedR, $queryHref, $outPrefix);

#Write total statistics file
print "Writing statistics: ".runTime($start)."\n";
getStat($ranksR, $outPrefix, $allQ, $resultsR, $lostEval, $type, $confCut);

#Lookup viruses, if necessary
if ($isVirus == 1) {
	print "Looking for viral hosts: ".runTime($start)."\n";
	getHostVir($pathosR, $outPrefix."PATHO.txt", $pdbPath);
}
#Lookup strains in BV-BRC
else{
	print "Looking for pathogens: ".runTime($start)."\n";
	getPatho($pathosR, $outPrefix."PATHO.txt", $pdbPath);
}

#Create html visualizations with krona
print "Creating diagrams: ".runTime($start)."\n";
visOut($krona, $outPrefix."VIS.txt");
visOut($krona, $outPrefix."PATHO.txt");

print "DONE: ".runTime($start)."\n\n";