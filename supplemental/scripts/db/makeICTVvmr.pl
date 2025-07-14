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
# Build a MetaG compatible database from ICTV VMR
#======================================================================================================================#
# DESCRIPTION
# 
# 	Given the ICTV VMR spreadsheet, write MetaG sequences, taxonomy, and host database files.
#
#
# USAGE
# 
# 	Download the ICTV VMR spreadsheet from https://ictv.global/vmr . Save the sheet containing the taxonomy as a CSV
#	file, replace empty cells with "NA," and remove	line breaks in fields. The field separator must be tab. Then run
#	makeICTVvmr.pl on the CSV file. If --allow-dups is provided and there are duplicate GenBank IDs, the script will
#	keep one of the copies. By default, both are removed. If a region is provided for a sequence in the format
#	"ID (min.max)," (1-based) this script automatically attempts to slice the sequence to the relevant regions.
#		
#	makeICTVvmr.pl --input vmr.csv [ --allow-dups ]
#
#	OR
#
#	makeICTVvmr.pl -i vmr.csv [ -a ]
#
#
# OUTPUT
# 
# 	The output is printed to three separate files in the directory of the input file: a taxonomy file, a file with host
#	information and one FASTA file with the sequences.
#	
#	tax.ICTV_vmr.txt:
#		#RANK_NAMES
# 		genbankID1;tax
# 		genbankID2;tax
#		...
#
#	patho.ICTV_vmr.txt:
#		genbankID1;host
#		genbankID2;host
#		...
#
#	
# LIMITATIONS
#	
#	GenBank IDs must only be used once. Otherwise, only the first (--allow-dups) or no (default) entry is kept. You may
#	want to report duplicate IDs to the ICTV team. IDs with no reply from the NCBI API are present in the
#	tax.ICTV_vmr.txt and patho.ICTV_vmr.txt, but not in the ICTV_vmr.fa. Additional taxa in these files do not influence
#	downstream analyses. You may want to check, if the GenBank IDs are valid and contact NCBI or ICTV, depending on the
#	situation.
#
#
# DEPENDENCIES
#
#	LWP:UserAgent
#======================================================================================================================#


use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use LWP::UserAgent();


#----------------------------------------------------------------------------------------------------------------------#
# Parse and check arguments
#----------------------------------------------------------------------------------------------------------------------#
my $argC		=	@ARGV;
my $inF			=	"";
my $help		=	0;
my $allowDups	=	0;
my $usage		=	<<'EOF';
#======================================================================================================================#
# Build a MetaG compatible database from ICTV VMR
#======================================================================================================================#
DESCRIPTION

	Given the ICTV VMR spreadsheet, write MetaG sequences, taxonomy, and host database files.


USAGE

	Download the ICTV VMR spreadsheet from https://ictv.global/vmr . Save the sheet containing the taxonomy as a CSV
	file, replace empty cells with "NA," and remove	line breaks in fields. The field separator must be tab. Then run
	makeICTVvmr.pl on the CSV file. If --allow-dups is provided and there are duplicate GenBank IDs, the script will
	keep one of the copies. By default, both are removed. If a region is provided for a sequence in the format
	"ID (min.max)," (1-based) this script automatically attempts to slice the sequence to the relevant regions.
	
	makeICTVvmr.pl --input vmr.csv [ --allow-dups ]
	
	OR

	makeICTVvmr.pl -i vmr.csv [ -a ]


OUTPUT

	The output is printed to three separate files in the directory of the input file: a taxonomy file, a file with host
	information and one FASTA file with the sequences.

	tax.ICTV_vmr.txt:
		#RANK_NAMES
		genbankID1;tax
		genbankID2;tax
		...

	patho.ICTV_vmr.txt:
		genbankID1;host
		genbankID2;host
		...


LIMITATIONS

	GenBank IDs must only be used once. Otherwise, only the first (--allow-dups) or no (default) entry is kept. You may
	want to report duplicate IDs to the ICTV team. IDs with no reply from the NCBI API are present in the
	tax.ICTV_vmr.txt and patho.ICTV_vmr.txt, but not in the ICTV_vmr.fa. Additional taxa in these files do not influence
	downstream analyses. You may want to check, if the GenBank IDs are valid and contact NCBI or ICTV, depending on the
	situation.


DEPENDENCIES

	LWP:UserAgent
#======================================================================================================================#
EOF
;

GetOptions (
	'input:s'		=>	\$inF,
	'allow-dups|a'	=>	\$allowDups,
    'help|?'		=>	\$help
) or die $usage;
           	
if ($help > 0 or not $inF) {
	print $usage;
	exit 0;
}

my $outP		=	dirname($inF);
my $outTaxP		=	"$outP/tax.ICTV_vmr.txt";
my $outHostP	=	"$outP/patho.ICTV_vmr.txt";
my $outSeqP		=	"$outP/ICTV_vmr.fa";


#----------------------------------------------------------------------------------------------------------------------#
# Parse input taxonomy and write output host and taxonomy files.
#----------------------------------------------------------------------------------------------------------------------#
print "INFO: Reading input taxonomy from $inF\n";
my $rankNames	=	"realm;subrealm;kingdom;subkingdom;phylum;subphylum;class;subclass;order;suborder;family;" .
						"subfamily;genus;subgenus;species;isolate_seqName";
my %ids			=	();
my %dups		=	();
my %regions		=	();

open(OUTTAX, ">", $outTaxP) or die "ERROR: Cannot open output taxonomy file ->$outTaxP<-";
print OUTTAX "#".$rankNames."\n";
open(OUTHOST, ">", $outHostP) or die "ERROR: Cannot open output host file ->$outHostP<-"; 

open(INTAX, "<", $inF) or die "ERROR: Cannot open input taxonomy table ->$inF<-";
while(<INTAX>) {
	#Skip header
	next if ($_	=~	m/GENBANK accession/);
	chomp($_);
	$_			=~	s/"//g;
	my @splits	=	split("\t", $_);
		
	# Filter strange undef entries
	die "ERROR: ->" . join("\t", @splits) . "<-\nFound only ->" . scalar(@splits) . "<- entries. Expected 27\n"
		if (scalar(@splits) < 27);
	
	my $genbID	=	$splits[23];
	my $host	=	$splits[$#splits - 1];
	my @hosts	=	split(',', $host);
	
	# Use only entries having a GenBank ID
	next if ($genbID eq "NA" or $genbID eq "");

	# Get and reformat taxonomy
	my $tax		=	join(";", @splits[3..17]);
	$tax		.=	";".$splits[22];
	$tax		=~	s/;NA/;0/g;
	$tax		=~	s/^NA;/0;/g;
	
	# Sometimes one entry has multiple GenBank IDs. Create separate entries per ID.
	$genbID		=~	s/,/;/g;
	$genbID		=~	s/; /;/g;
	my @genbIDs	=	split(';', $genbID);
	foreach my $id (@genbIDs) {
		$id						=~	s/;$//;
		my $seqName 			=	"";
		
		# Again: Sometimes different GenBank IDs are noted with seqName of the nucleotides
		if ($id =~ m/:/) {
			my @seqNames	=	split(':', $id);
			$seqName		=	$seqNames[0];			
			$id				=	$seqNames[1];
			$id				=~	s/^ //;
		}
		# Remove trailing and leading whitespaces in ids
		$id						=~	s/^\s+|\s+$//; 
		
		# Only a certain region of the record should be kept
		my ($tmp, $min, $max)	=	($1, $2, $3) if ($genbID =~ m/^(.*)\((\d+)\.(\d+)\)$/);
		if (defined $min and not defined $max) {
			die "ERROR: Invalid region ->" . ($min // "") . "<- ->" . ($max // "") . "<-";
		}
		elsif (defined $max and not defined $min) {
			die "ERROR: Invalid region ->" . ($min // "") . "<- ->" . ($max // "") . "<-";
		}
		elsif (defined $max and defined $min) {
			die "ERROR: Invalid region ->" . $min . "<- ->" . $max . "<-" if ($min > $max);
			die "ERROR: Invalid ID" if (not defined $tmp);
			$tmp			=~	s/\s+$//;
			$id				=	$tmp;
			$regions{$id}	=	[$min, $max]
		}


		if (exists $ids{$id}) {
			# Allow one of the duplicate IDs in the FASTA file
			if ($allowDups == 1) {
				print "ERROR: GenBank id ->$id<- is not unique. Removed duplicate entry.\n";
				print "\tCURRENT: ->$tax<-: ->$seqName<-\n\tFIRST SEEN: ->" . $ids{$id} . "<-\n";
				next;
			}
			# Remove all entries belonging to a duplicate ID
			else {
				print "ERROR: GenBank id ->$id<- is not unique. Removing all occurences.\n";
				print "\tCURRENT: ->$tax<-: ->$seqName<-\n\tFIRST SEEN: ->" . $ids{$id} . "<-\n";
				$dups{$id}	=	undef;
				next;
			}
		}
		else {
			$ids{$id}		=	$tax . ": $seqName";
		}

		
		if ($seqName ne "") {
			print OUTTAX $id.";".$tax."_".$seqName."\n";
		}
		else{
			print OUTTAX $id.";".$tax."\n";
		}
		
		foreach my $host (@hosts) {
			$host =~ s/^ //;
			# Host may be a species name or isolation source
			if ($host		=~	m/ \(S\)$/) {
				$host		=~	s/ \(S\)$//;
				$host		=	"source: " . $host;
			}
			print OUTHOST $id.";".$host."\n";
		}
	}
}
close(INTAX);
close(OUTTAX);
close(OUTHOST);

# Delete duplicate IDs
if ($allowDups == 0 and %dups) {
	foreach my $dup (keys(%dups)) {
		delete $ids{$dup};
	}
}


#----------------------------------------------------------------------------------------------------------------------#
# Request nucleotide sequences by GenBank accession from NCBI and write output sequence file.
#----------------------------------------------------------------------------------------------------------------------#
my @ids			=	keys(%ids);
die "ERROR: No IDs remaining" if (not @ids);
my $i			=	0;
my $chunkSize	=	1000;
print "INFO: Web request to NCBI: " . scalar(@ids) . " sequences. Sending chunks of ->$chunkSize<-. " . 
	"This may take a while.\n";

# Obtain sequences and IDs and write as FASTA file
open(OUTSEQ ,">", $outSeqP)  or die "ERROR: Cannot open output sequence file ->$outSeqP<-";
my @entries		=	();
my $chunkC		=	1;
while (@ids) {
	my @queries	=	splice (@ids, 0, $chunkSize);
	
	# POST request to NCBI eFetch
	my $ua		=	LWP::UserAgent->new;
	$ua->timeout(45);
	$ua->env_proxy;
	my $query	=	join(",", @queries);
	my %query	=	("db" => "nuccore", "id" => $query, "rettype" => "fasta", "retmode" => "text");

	# Response from NCBI. No "www.", causes error.
	my $seq = $ua->post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", \%query);

	# Crash, if there was an error
	die "ERROR: Error querying NCBI with chunk ->$chunkC<-. Error code: ->" . $seq->status_line .
		"<-\n" unless ($seq->is_success);

	# Response body
	$seq	=	$seq->content;
	$seq	=~	s/^\s+|\s+$//;
	$seq	=~	s/^>//;
	my @res	=	split('>', $seq);
	chomp(@res);
	push(@entries, @res);

	# Maximum three request per second
	select(undef, undef, undef, 0.34);
	$chunkC++;
}

print "INFO: Retrieved ->" . scalar(@entries) . "<- sequences from NCBI.\n";
foreach my $entry (@entries) {
	next if ($entry eq "");
	my @lines	=	split("\n", $entry);
	my $id 		=	(split(" ", $lines[0]))[0];

	# Some GenBank IDs were sent without a version tag, but I received a response with a tag. In this case, I need to
	# save the sequence data without the tag. If I queried GenBank with a tag, I want to keep it. This keeps the naming
	# consistent between the output files.
	my $idNoTag	=	(split('\.', $id))[0];
	if (exists $ids{$id}) {
		delete $ids{$id}
	}
	elsif (exists $ids{$idNoTag}) {
		delete $ids{$idNoTag};
		$id = $idNoTag
	}
	else {
		next;
	}
	
	my $seq		=	"";
	foreach my $line (@lines[1..$#lines]) {
		$seq	.=	$line;          
	}  
	# Optionally slice to relevant region (1-based coordinates)
	if (exists $regions{$id}) {
		my ($min, $max)	=	@{$regions{$id}};
		# 1-based notation to 0-based for slicing
		$min--;
		$seq			=	substr($seq, $min, $max - $min)
	}
	$seq .= "\n";
	print OUTSEQ ">".$id."\n".$seq;

	$i++;
}
close(OUTSEQ);

print "WARNING: Could not find GenBank IDs: ->" . join(", ", keys(%ids)) . "<-\n" if (%ids);
print "\nDONE\n";
