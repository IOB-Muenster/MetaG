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


#=====================================================================================================#
# This script builds a PATRIC database in MetaG format.
# First you need to download genome_lineage and genome_metadata from
# ftp://ftp.patricbrc.org/RELEASE_NOTES/ to the directory of this script.
# 
# Then run
# awk -F '\t' '{if (tolower($46) ~ /^human|^homo/) {print $1"\tHomo sapiens\t"$53}}' genome_metadata > patricHuman.txt
# To get only information about human hosts.
#
# This script will produce the patho.PATRIC.txt file which can be supplied as -pdbPath to MetaG.
# The file has the following format:
#
# PATRICgenomeID;lineage;#host;#resistance
#------------------------------------------------------------------------------------------------------#


use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $dir = dirname(__FILE__);


# PATRICS own taxonomy
my %lins = ();
open(LIN, "<", "$dir/genome_lineage") or die "Can't open PATRIC lineage";
while(<LIN>) {
	
	# Loose header
	next if ($_ =~ m/^genome_id/);
	chomp($_);
	
	my @splits = split("\t", $_);
	my $id = $splits[0];
	my $species = $splits[1];
	my ($kingdom, $phylum, $class, $order, $family, $genus,) = @splits[3..6];
	
	foreach my $rank ($kingdom, $phylum, $class, $order, $family, $genus, $species) {
		$rank = "0" if (not $rank);
	}
	
	#Separate species and strain
	my $strain = 0;
	if ($species ne "0") {
		
		# Improves matching across dbs
		$species =~ s/"//g;
		
		my @specSplits = split(" ", $species);
		my $lenSpec = @specSplits;
		
		if ($lenSpec >= 3) {
			
			if ($specSplits[1] =~ m/\.$/ or $specSplits[1] eq "bacterium" or $specSplits[1] eq "undefined") {
				$species = "unclassified";
				$strain = join(" ", @specSplits[2..$#specSplits]);
			}
			elsif($specSplits[0] eq "uncultured") {
				$species = "uncultured";
			}
			elsif ($specSplits[0] =~ m/^[a-z0-9]/) {
				$species = "unclassified";
				$strain = join(" ", @specSplits);
			}
			else {
				$species = join(" ", @specSplits[0..1]);
				$strain = join(" ", @specSplits[2..$#specSplits]);
			}
		}
		else {
			
			if ($lenSpec > 1 and $specSplits[1] =~ /\.$/ or $specSplits[1] eq "bacterium" or $specSplits[1] eq "undefined") {
				$species = "unclassified";
			}
			elsif($specSplits[0] eq "uncultured") {
				$species = "uncultured";
			}
			elsif ($specSplits[0] =~ m/^[a-z0-9]/) {
				$species = "unclassified";
				$strain = join(" ", @specSplits);
			}
			else {
				$species = join(" ", @specSplits);
			}
		}
	}
	
	$strain =~ s/;/,/g;
	$lins{$id} = "$kingdom;$phylum;$class;0;$order;0;$family;$genus;$species;$strain";
    
}


open(OUT, ">", "$dir/patho.PATRIC.txt") or die "Can't open OUTFILE";
print OUT "#PATRICgenomeID;lineage;#host;#resistance";


# Get PATRIC metadata and connect with taxonomy. Metadata was processed previously with
# awk -F '\t' '{if ($46 ~ /^Human|^Homo/) {print $1"\t"$46"\t"$53}}' genome_metadata > patricHuman.txt
# to get host name and resistances for human hosts, only.
open(META, "<", "$dir/patricHuman.txt") or die "Can't open PATRIC metadata";
while(<META>) {
	
	next if ($_ =~ m/^#/);
	chomp($_);
	
	my @splits = split("\t", $_);
	my $id = $splits[0];
	my $host = $splits[1];
	
	# Some people can't spell Homo sapiens :)
	$host = "Human" if ($host =~ m/^Human/ or $host =~ m/^Homo/);
	$host =~ s/;/,/g;
	
	my $resist = $splits[2];
	$resist = "0" if (not $resist);
	$resist =~ s/;/+/g;
	
	my $lin = undef;
	
	if (exists $lins{$id}) {
		$lin = $lins{$id};
	}
	else {
		die "#$id# not in lineage file."
	}
	
	print OUT "\n>$id;$lin;#$host;#$resist";

}
print "DONE\n"

