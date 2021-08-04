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


use strict;
use warnings;
use Getopt::Long;


#=========================================================================#
# Filter reads
#=========================================================================#
#
# DESCRIPTION
#
# Extracts reads from an alignment file.
#
#
# USAGE
#
# ./filterReads.pl -alignment ALIGN.MAF -sequences SEQ.FA [--unaligned]
#		[--include PATTERN | --exclude PATTERN]
#
# Given an alignment and a query fasta file, extract all sequences that
# were aligned (default) or unaligned (--unaligned). It is also possible
# to filter the output by choosing inclusion or exclusion patterns for
# the read IDs.
# Results are printed to STDOUT.
#=========================================================================#


# Extract Reads from MAF
my $inclPat = "";
my $exclPat = "";

# Search in fasta
my $isUnaligned = 0;

my $mafF = "";
my $seqF = "";

GetOptions ("alignment=s" => \$mafF,
            "sequences=s"   => \$seqF,
            "unaligned"  => \$isUnaligned,
	    	"include=s"	=> \$inclPat,
	    	"exclude=s"	=> \$exclPat,
) or die("ERROR: Unrecognized argument\n");

die "ERROR: Missing essential argument" if (not $mafF or not $seqF);
die "ERROR: Can only include or exclude, not both" if ($inclPat and $exclPat);

die "ERROR: Invalid or empty alignment file" if (! -s $mafF);
die "ERROR: Invalid or empty sequence file" if (! -s $seqF);

my %headers = ();


#------------------------------------------------------------------#
# Read headers from MAF which match or don't match the patterns
#------------------------------------------------------------------#
my $sCounter = 0;

open(MAF, "<", $mafF) or die "Could not open MAF file, $mafF";
while(<MAF>) {
	if ($_ =~ m/^s\s/) {
		my $id = (split(" ", $_))[1];

		$sCounter++;		

		# DB header is 1st s-line. Query header in 2nd.		
		if ($sCounter == 2) {
			if ($inclPat) {
				next if ($id !~ m/$inclPat/);
				$headers{$id} = ""; 
			}
			elsif ($exclPat) {
				next if ($id =~ m/$exclPat/);
				$headers{$id} = "";
			}
			else {
				$headers{$id} = "";
			}
			$sCounter = 0;
		}
	}
}
close(MAF);


#------------------------------------------------------------------#
# Get the sequences for all aligned or unaligned reads 
#------------------------------------------------------------------#
my $isMatch = 0;

open(SEQ, "<", $seqF) or die "Could not open sequence file, $seqF";
while (<SEQ>) {
	if ($_ =~ m/^>/) {
		chomp ($_);

		my $id = (split(" ", $_))[0];
		$id =~ s/^>//;

		# Get all aligned reads --> reads that appear in MAF
		if ($isUnaligned == 0) {
			if (exists $headers{$id}) {
				# Full header, not just ID
				print $_, "\n";
				$isMatch = 1;
			}
			else {
				$isMatch = 0;
			}
		}
		# Get all unaligned reads --> reads that don't appear in MAF
		else {
			if (exists $headers{$id}) {
            	$isMatch = 0;
            }
            else {
				print $_, "\n";
                $isMatch = 1;
            }
		}
	}
	else {
		print $_ if ($isMatch == 1);
	}
}
close(SEQ);
