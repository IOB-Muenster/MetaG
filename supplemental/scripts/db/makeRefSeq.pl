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
# Build a MetaG compatible database from RefSeq's Targeted Loci
#======================================================================================================================#
# DESCRIPTION
#
#	Build a MetaG compatible database from a RefSeq FASTA file.
#
#
# USAGE
#
#	Choose and download the refseq database files (*fna.gz) from ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/ . You
#	may concatenate multiple files of interest into one input file and supply it via --input to create a merged MetaG
#	database. Note that the input files	need to be extracted. Download and extract the taxdump files for the NCBI
#	Taxonomy Toolkit (see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
#	directory of your input file, you need to provide the path using --kitdata. Use --debug to enable debug messages.
#
#	makeRefSeq.pl --input seqfile.fna [ --kitdata directory ] [--debug]
#
#	OR
#
#	makeRefSeq.pl -i seqfile.fna [ -k directory ] [-d]
#
#
# OUTPUT
#
#	A FASTA file with "metag." prefix and a taxonomy file with "tax.metag." prefix. Both will be created in the
#	directory of your input file.
#
#
# DEPENDENCIES
#
# 	*	NCBI Taxonomy Toolkit v0.20.0 (https://doi.org/10.1016/j.jgg.2021.03.006)
#		https://github.com/shenwei356/taxonkit
#		Must be located in the PATH!
#	*	JSON::Tiny
#		https://metacpan.org/pod/distribution/JSON-Tiny/lib/JSON/Tiny.pod
# 	*	LWP:UserAgent
#		https://metacpan.org/pod/LWP::UserAgent
#
#
# MISCELLANEOUS
#
#	*	The first rank is either "domain," "cellular root," or "superkingdom."
#	*	The last rank after species is either "strain" or "subspecies."
#	*	If taxonkit indicates that a certain taxonomy ID could not be found, make sure to use the current version of the
#		NCBI taxdump files (see USAGE).
#
#	  
# KNOWN BUGS
# 
#	*	The API may occasionally produce "wrong UID" errors. Checking the UID on the NCBI website and, if valid,
#		rerunning the script may fix the error. True errors would always appear for the same ID and in the same data
#		chunk (check debug messages), since the order of IDs will not change between runs of this program, unless the
#		input data was changed. Expert users may want to check additionally (replace YOUR_UID with the erroneous UID).
#		https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=YOUR_UID
#======================================================================================================================#


use strict;
use warnings;
use LWP::UserAgent;
use JSON::Tiny qw(decode_json);
use File::Basename;
use Getopt::Long;


#----------------------------------------------------------------------------------------------------------------------#
# Parse and check arguments
#----------------------------------------------------------------------------------------------------------------------#
my $help		=	0;
my $debug		=	0;
my $seqF		=	"";
my $taxkitData	=	"";
my $url 		=	"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json";
my $usage		=	<<'EOF';
#======================================================================================================================#
# Build a MetaG compatible database from RefSeq's Targeted Loci
#======================================================================================================================#
DESCRIPTION

	Build a MetaG compatible database from a RefSeq FASTA file.


USAGE

	Choose and download the refseq database files (*fna.gz) from ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/ . You
	may concatenate multiple files of interest into one input file and supply it via --input to create a merged MetaG
	database. Note that the input files	need to be extracted. Download and extract the taxdump files for the NCBI
	Taxonomy Toolkit (see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
	directory of your input file, you need to provide the path using --kitdata. Use --debug to enable debug messages.

	makeRefSeq.pl --input seqfile.fna [ --kitdata directory ] [--debug]

	OR

	makeRefSeq.pl -i seqfile.fna [ -k directory ] [-d]


OUTPUT

	A FASTA file with "metag." prefix and a taxonomy file with "tax.metag." prefix. Both will be created in the
	directory of your input file.


DEPENDENCIES

 	*	NCBI Taxonomy Toolkit v0.20.0 (https://doi.org/10.1016/j.jgg.2021.03.006)
		https://github.com/shenwei356/taxonkit
		Must be located in the PATH!
	*	JSON::Tiny
		https://metacpan.org/pod/distribution/JSON-Tiny/lib/JSON/Tiny.pod
 	*	LWP:UserAgent
		https://metacpan.org/pod/LWP::UserAgent


MISCELLANEOUS

	*	The first rank is either "domain," "cellular root," or "superkingdom."
	*	The last rank after species is either "strain" or "subspecies."
	*	If taxonkit indicates that a certain taxonomy ID could not be found, make sure to use the current version of the
		NCBI taxdump files (see USAGE).

	  
KNOWN BUGS
 
	*	The API may occasionally produce "wrong UID" errors. Checking the UID on the NCBI website and, if valid,
		rerunning the script may fix the error. True errors would always appear for the same ID and in the same data
		chunk (check debug messages), since the order of IDs will not change between runs of this program, unless the
		input data was changed. Expert users may want to check additionally (replace YOUR_UID with the erroneous UID).
		https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=YOUR_UID
#======================================================================================================================#
EOF
;

my $argC		=	scalar(@ARGV);
GetOptions (
	"input:s"	=>	\$seqF,
	"kitdata:s"	=>	\$taxkitData,
	'help|?'	=>	\$help,
	'debug'		=>	\$debug
) or die $usage;       	
if ($help > 0 or not $seqF) {
	print $usage;
	exit 0;
}

if (not $taxkitData) {
	$taxkitData		=	dirname($seqF);
	print "INFO: --kitdata omitted. Assuming: ->$taxkitData<-\n";
}
# Quick check to see, if taxkitData is OK.
system("echo \"9606\" | taxonkit lineage --data-dir $taxkitData >/dev/null");
my $rc			=	$? >> 8;
die "ERROR: TaxonKit failed. Check --kitdata or the installation." if ($rc != 0);

my $outFasta	=	dirname($seqF) . "/metag." . basename($seqF);
my $outTax		=	dirname($seqF) . "/tax.metag." . basename($seqF);
$outTax			=~	s/\.[^.]*$/\.txt/;
my $header		=	"";
my $id			=	"";
my @ids			=	();
my %ids			=	();
my %seqs		=	();
my %taxids		=	();


#----------------------------------------------------------------------------------------------------------------------#
# Reformat the FASTA file
#----------------------------------------------------------------------------------------------------------------------#
print "INFO: Writing FASTA file\n";
open(SEQ, "<", $seqF) or die "Could not open sequence FASTA file ->$seqF<-";
while(<SEQ>) {
	chomp($_);
	if ($_ =~ m/^>/) {
		$header			=	(split(" ", $_))[0];
		$id				=	$header	=~	s/^>//r;
				
		# ID without version needed to check, if all IDs were returned by the API. API looses version tag at the end,
		# but I also need the original sequence IDs lateron.
		my $idNoVers	=	$id		=~	s/\..*$//r;
		if (exists $ids{$idNoVers}) {
			die "ERROR: RefSeq ID must be unique, but is not: ->$idNoVers<-"
		}
		else {
			$ids{$idNoVers}	=	$id;
			push (@ids, $id);
		}
	}
	else {
		die "ERROR: Unexpected FASTA: Sequence not preceeded by header in line ->" . $. . "<-\n" if (not $id);
		if (exists $seqs{$id}) {
			$seqs{$id}	.=	$_;
		}
		else {
			$seqs{$id}	=	$_;
		}
	}
}
close(SEQ);


#----------------------------------------------------------------------------------------------------------------------#
# Get tax IDs from eUTILs using sequence IDs
#----------------------------------------------------------------------------------------------------------------------#
# Max is 500 for JSON
my $chunkSize	=	450;
my $chunkC		=	1;
# Maximum number of retries for API requests
my $maxRetry	=	3;
@ids			=	sort{$a cmp $b} @ids;
print "INFO: Web request to NCBI: ->" . scalar(@ids) . "<- sequences. Sending chunks of ->$chunkSize<-. " .
		"This may take a while.\n";
	
while (@ids) {
	print "DEBUG: Querying chunk ->$chunkC<-\n" if ($debug == 1);
    my @queries		=	splice (@ids, 0, $chunkSize);	
	my %queries		=	('id' => join(",", @queries));
	my $ua			=	LWP::UserAgent->new(timeout => 40);
	
	# Try to query a chunk again, if there was any error. Sometimes NCBI API returns "random" error for valid ID.
	my $retryC		=	0;
	my $isSuccess	=	0;
	while ($isSuccess == 0 and $retryC < $maxRetry) {
		my $response	=	$ua->post($url, \%queries);
		if (not $response->is_success) {
			print "ERROR: API error in chunk ->$chunkC<-: ->" . $response->status_line . "<-\n"
		}
		else {
			if (not defined $response->content or not $response->content) {
				print "WARNING: Unexpected API response\n";
			}
			else {
				$response	=	decode_json ($response->content);
				if (exists $response->{'error'}) {
					print "API ERROR: ->" . $response->{'error'} . "<-\n";
				}
				else {
					if (exists $response->{'result'}) {
						if (exists $response->{'result'}->{'uids'}) {
							my @uids		=	@{$response->{'result'}->{'uids'}};
							my %seqidApis	=	();
							if (scalar(@queries) == scalar(@uids)) {
								for (my $i = 0; $i <= $#uids; $i++) {
									my $uid		=	$uids[$i];
									if (exists $response->{'result'}->{$uid}) {
										my $taxid		=	$response->{'result'}->{$uid}->{'taxid'} // "";
										my $seqidApi	=	$response->{'result'}->{$uid}->{'caption'} // "";
										if ($taxid ne "" and $seqidApi ne "" and exists $ids{$seqidApi}) {
											$seqidApis{$seqidApi}	=	undef;
											my $seqid				=	$ids{$seqidApi};
											$taxids{$seqid} 		=	$taxid;
										}
									}
									else {
										print "WARNING: Unexpected API response.\n"
									}
								}
								# Only remove uids from candidate list after all in batch were successfully retrieved
								if (scalar(@uids) == scalar(keys(%seqidApis))) {
									# Flag to break out of the while loop
									$isSuccess	=	1;
									foreach my $seqidApi (keys(%seqidApis)) {
										delete $ids{$seqidApi};
									}
								}
								else {
									print "WARNING: Unexpected API response. Some results missing.\n"
								}
							}
							else {
								print "WARNING: Too few IDs returned by API\n";
							}
						}
						else {
							print "WARNING: Unexpected API response\n";
						}
					}
					else {
						print "WARNING: Unexpected API response\n";
					}
				}
			}
		}
		if ($isSuccess == 0) {
			$retryC++;
			if ($retryC < $maxRetry) {
				print "INFO: Error in API response. Trying again (->$retryC<- of ->$maxRetry<-).\n";
				sleep(8);
			}
			else {
				die "ERROR: No valid API response within retry limit ->$maxRetry<-\n";
			}
		} 
	}
	# Only 2 requests per second
	select(undef, undef, undef, 0.5);
	$chunkC++;
}
if (keys(%ids)) {
	die "ERROR: Could not get tax IDs for sequence ids: ->" . join("<-; ->", keys(%ids)) . "<-\n";
}

open(OUTTAX, ">", $outTax) or die "Could not open output taxonomy file ->$outTax<-";
foreach my $seqid (keys(%taxids)) {
	print OUTTAX $seqid, "\t", $taxids{$seqid}, "\n";
}


#----------------------------------------------------------------------------------------------------------------------#
# Use the NCBI Taxonomy Toolkit to get the lineage for each sequence ID using its tax ID
#----------------------------------------------------------------------------------------------------------------------#
print "INFO: Getting lineage info from tax IDs\n";
my $tmpTax		=	dirname($outTax) . "/tmp." . basename($outTax);
# Use the tax ID in the second column to write a temporary taxonomy with 8 ranks, missing ranks are "0"
my $command 	=	(
	"taxonkit reformat2 -I 2 --data-dir $taxkitData --miss-rank-repl \"0\" -f " .
		"\"{domain|cellular root|superkingdom}\t{phylum}\t{class}\t{order}\t{family}\t{genus}\t{species}\t" .
		"{strain|subspecies}\" $outTax > $tmpTax"
);
system($command);
$rc				=	$? >> 8;
die "ERROR: TaxonKit failed." if ($rc != 0);


#----------------------------------------------------------------------------------------------------------------------#
# Reformat output
#----------------------------------------------------------------------------------------------------------------------#
print "INFO: Writing taxonomy file ->$outTax<-\n";
open(OUTSEQ, ">",$outFasta) or die "ERROR: Could not open output FASTA file ->$outFasta<-";
open(OUTTAX, ">", $outTax) or die "ERROR: Could not open final output taxonomy file ->$outTax<-";
open(TMPTAX, "<", $tmpTax) or die "ERROR: Could not open temporary taxonomy file ->$tmpTax<-";
while(<TMPTAX>) {
	chomp($_);
	my @splits		=	split("\t", $_);
	my $seqid		=	$splits[0];
	
	# Get the whole tab-delimited taxonomy
	my @taxa		=	@splits[2..$#splits];
	if (not @taxa) {
		print "DEBUG: Lost ->$seqid<-. No lineage.\n" if ($debug == 1);
		next;
	}
	# Loose entries that have unknown ("0") taxon at all ranks. This can happen, if the associated taxonomy ID was
	# deleted by NCBI.
	my %tmps		=	map {$_ => undef} @taxa;
	if (scalar(keys(%tmps)) == 1 and exists $tmps{"0"}) {
		print "DEBUG: Lost ->$seqid<-. Only unknown taxa in lineage.\n" if ($debug == 1);
		next;
	}
	# Reformatted taxonomy has 8 ranks. MetaG needs 10. Insert subclass and suborder as "0" = unknown
	splice(@taxa, 3, 0, "0");
	splice(@taxa, 5, 0, "0");
	print OUTTAX $seqid, ";", join(";", @taxa), "\n";
	print OUTSEQ ">" . $seqid . "\n" . $seqs{$seqid} ."\n";
}
close(TMPTAX);
close(OUTTAX);
close(OUTSEQ);

# Cleanup
system("rm $tmpTax");
$rc				=	$? >> 8;
die "ERROR: Could not remove temporary taxonomy file ->$tmpTax<-" if ($rc != 0);

print "\nDONE\n";
