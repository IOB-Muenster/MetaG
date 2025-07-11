#!/usr/bin/env perl


# AUTHORS

# Felix Manske, felix.manske@uni-muenster.de
# Norbert Grundmann
# Author for correspondence: Wojciech Makalowski, wojmak@uni-muenster.de


# COPYRIGHT

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
#	disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
#	disclaimer in the documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# POSSIBILITY OF SUCH DAMAGE.


#======================================================================================================================#
# Build a MetaG compatible database from BV-BRC (formely PATRIC)
#======================================================================================================================#
# DESCRIPTION
#
#	Build a MetaG compatible database from BV-BRC's genome_lineage and genome_metadata files.
#
#
# USAGE
# 
# 	Download the genome_lineage and genome_metadata files from ftp://ftp.bvbrc.org/RELEASE_NOTES/ and run makeBVBRC.pl.
#		
#	makeBVBRC.pl --metadata genome_metadata --lineage genome_lineage
#
#	OR
#
#	makeBVBRC.pl -m genome_metadata -l genome_lineage
#
#
# OUTPUT
#
#	A host file called patho.BVBRC.txt in the directory of the genome_lineage file. This can be used as a pathogen file
#	for MetaG. The file only contains pathogens from human hosts and has the following format:
#
#	BV-BRCgenomeID;lineage;#host;#resistance
#======================================================================================================================#


use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#----------------------------------------------------------------------------------------------------------------------#
# Parse and check arguments
#----------------------------------------------------------------------------------------------------------------------#
my $argC		=	scalar(@ARGV);
my $metaF		=	"";
my $lineageF	=	"";
my $help		=	0;
my $usage		=	<<'EOF';
#======================================================================================================================#
# Build a MetaG compatible database from BV-BRC (formely PATRIC)
#======================================================================================================================#
DESCRIPTION

	Build a MetaG compatible database from BV-BRC's genome_lineage and genome_metadata files.


USAGE
 
 	Download the genome_lineage and genome_metadata files from ftp://ftp.bvbrc.org/RELEASE_NOTES/ and run makeBVBRC.pl.
		
	makeBVBRC.pl --metadata genome_metadata --lineage genome_lineage

	OR

	makeBVBRC.pl -m genome_metadata -l genome_lineage


OUTPUT

	A host file called patho.BVBRC.txt in the directory of the genome_lineage file. This can be used as a pathogen file
	for MetaG. The file only contains pathogens from human hosts and has the following format:

	BV-BRCgenomeID;lineage;#host;#resistance
#======================================================================================================================#
EOF
;

GetOptions (
	"metadata:s"	=>	\$metaF,
	"lineage:s"		=>	\$lineageF,
	'help|?'		=>	\$help
) or die $usage;
           	
if ($help > 0 or not $metaF or not $lineageF) {
	print $usage;
	exit 0;
}
my $outP		=	dirname($lineageF);


#----------------------------------------------------------------------------------------------------------------------#
# Parse BV-BRC taxonomy
#----------------------------------------------------------------------------------------------------------------------#
my %lins		=	();
open(LIN, "<", "$lineageF") or die "ERROR: Can't open BV-BRC lineage ->$lineageF<-";
while(<LIN>) {
	# Loose header
	next if ($_ =~ m/^genome_id/);
	chomp($_);
	
	my @splits		=	split("\t", $_, -1);
	my $id			=	$splits[0] // "";
	my $species		=	$splits[1] // "";
	die "ERROR: No ID\n$_" if (not $id);
	
	my (
		$kingdom,
		$phylum,
		$class,
		$order,
		$family,
		$genus
	)			=	@splits[3..8];
	foreach my $rank ($kingdom, $phylum, $class, $order, $family, $genus, $species) {
		$rank	=	"0" if (not $rank);
	}
	
	# Separate species and strain
	my $strain		=	0;
		
	# Improves matching across dbs
	$species		=~	s/"//g;
	
	my @specSplits	=	split(" ", $species);
	my $lenSpec 	=	scalar(@specSplits);
	
	if ($lenSpec < 2) {
		$species	=	"unclassified";
	}
	else {
		if ($specSplits[1] =~ m/\.$/ or $specSplits[1] eq "bacterium" or $specSplits[1] eq "undefined") {
			$species	=	"unclassified";
			if ($lenSpec >= 3) {
				$strain	=	join(" ", @specSplits[2..$#specSplits]);
			}
		}
		elsif($specSplits[0] eq "uncultured") {
			$species	=	"uncultured";
		}
		elsif ($specSplits[0] =~ m/^[a-z0-9]/) {
			$species	=	"unclassified";
			$strain		=	join(" ", @specSplits);
		}
		else {
			$species	=	join(" ", @specSplits[0..1]);
			if ($lenSpec >= 3) {
				$strain	=	join(" ", @specSplits[2..$#specSplits]);
			}
		}
	}
	$strain 		=~	s/;/,/g;
	# Insert placeholder "0" for subclass and suborder
	$lins{$id}		=	"$kingdom;$phylum;$class;0;$order;0;$family;$genus;$species;$strain";
    
}


#----------------------------------------------------------------------------------------------------------------------#
# Parse BV-BRC metadata and connect it to taxonomy
#----------------------------------------------------------------------------------------------------------------------#
open(OUT, ">", "$outP/patho.BVBRC.txt") or die "ERROR: Can't open output file";
print OUT "#BV-BRCgenomeID;lineage;#host;#resistance";
open(META, "<", "$metaF") or die "ERROR: Can't open BV-BRC metadata ->$metaF<-";
while(<META>) {
	# Loose header
	next if ($_ =~ m/^genome_id/);
	chomp($_);
	
	my @splits		=	split("\t", $_, -1);
	my $id			=	$splits[0] // "";
	my $host		=	$splits[45] // "";
	die "ERROR: No ID\n$_" if (not $id);
	
	# Retain only species from human host
	if (lc($host) =~ m/^human/ or lc($host) =~ m/^homo/) {
		$host	=	"Human";
	}
	else {
		next
	}
	my $resist		=	$splits[52] // "";
	$resist			=	"0" if (not $resist);
	$resist			=~	s/;/+/g;
	
	my $lin			=	undef;
	if (exists $lins{$id}) {
		$lin		=	$lins{$id};
	}
	else {
		die "ERROR: ID ->$id<- not in lineage file."
	}
	print OUT "\n>$id;$lin;#$host;#$resist";
}
print "DONE\n"