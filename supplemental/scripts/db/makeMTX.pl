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


#=======================================================================================================#
# Build a MetaG compatible database from MTX
#-------------------------------------------------------------------------------------------------------#
#
# USAGE
#	
#	Download and extract the taxdump files for the NCBI Taxonomy Toolkit
#	(see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
#	directory of your input file, you need to provide the path using -kitdata. Use --debug to enable
#	debug messages.
# 	Locate the database directory of Metaxa2 (https://doi.org/10.1111/1755-0998.12399) and run
# 	blastdbcmd (https://doi.org/10.1186/1471-2105-10-421) on each of your BLAST databases of interest
# 	(see example below for LSU).
#
#	blastdbcmd ‑db LSU/blast ‑entry all ‑out lsu.fa
#
#	Concatenate the resulting FASTA files to create MetaG compatible database files using makeMTX.pl. 
#		
#	makeMTX.pl -input mtx.fa [ -kitdata directory ] [ --debug ]
#
#	OR
#
#	makeMTX.pl -i mtx.fa [ -k directory ] [ -d ]
#
# OUTPUT
#
#	A fasta file with "metag." prefix and a taxonomy file with "tax.metag." prefix. Both will be created
#	in the directory of your input FASTA file.
#
# DEPENDENCIES
#
# 	* NCBI Taxonomy Toolkit v0.10.1 (https://doi.org/10.1101/513523)
#		https://github.com/shenwei356/taxonkit
#		Must be located in the PATH!
#	* JSON::Tiny
#		https://metacpan.org/pod/distribution/JSON-Tiny/lib/JSON/Tiny.pod
# 	* LWP:UserAgent
#		https://metacpan.org/pod/LWP::UserAgent
#
# MISCELLANEOUS
#
#	* The last rank after species is either strain, subspecies or pseudo-strain.
#	* If taxonkit indicates that a certain taxonomy ID could not be found, make sure to use the current
#	  version of the NCBI taxdump files (see USAGE).
#	  
# KNOWN BUGS
# 
#	* The API may occasionally produce "wrong UID" errors. Checking the UID on the NCBI website and, if
#	  valid, rerunning the script may fix the error. True errors would always appear for the same ID
#	  and in the same data chunk (check debug messages), since the order of IDs will not change between
#	  runs of this program, unless the input data was changed.
#	  Expert users may want to check additionally
#	  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=YOUR_UID
#	  (replace YOUR_UID with the erroneous UID).
#	  
#=========================================================================================================#


use strict;
use warnings;
use LWP::UserAgent;
use JSON::Tiny qw(decode_json);
use File::Basename;
use Getopt::Long;


my $help = 0;
my $debug = 0;
my $seqF = "";
my $taxkitData = "";
my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json";

my $usage = <<'EOF';
#=======================================================================================================#
# Build a MetaG compatible database from MTX
#-------------------------------------------------------------------------------------------------------#

 USAGE
	
	Download and extract the taxdump files for the NCBI Taxonomy Toolkit
	(see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
	directory of your input file, you need to provide the path using -kitdata. Use --debug to enable
	debug messages.
 	Locate the database directory of Metaxa2 (https://doi.org/10.1111/1755-0998.12399) and run
 	blastdbcmd (https://doi.org/10.1186/1471-2105-10-421) on each of your BLAST databases of interest
 	(see example below for LSU).

	blastdbcmd ‑db LSU/blast ‑entry all ‑out lsu.fa

	Concatenate the resulting FASTA files to create MetaG compatible database files using makeMTX.pl. 
		
	makeMTX.pl -input mtx.fa [ -kitdata directory ] [ --debug ]

	OR

	makeMTX.pl -i mtx.fa [ -k directory ] [ -d ]

 OUTPUT

	A fasta file with "metag." prefix and a taxonomy file with "tax.metag." prefix. Both will be created
	in the directory of your input FASTA file.

 DEPENDENCIES

 	* NCBI Taxonomy Toolkit v0.10.1 (https://doi.org/10.1101/513523)
		https://github.com/shenwei356/taxonkit
		Must be located in the PATH!
	* JSON::Tiny
		https://metacpan.org/pod/distribution/JSON-Tiny/lib/JSON/Tiny.pod
 	* LWP:UserAgent
		https://metacpan.org/pod/LWP::UserAgent

 MISCELLANEOUS

	* The last rank after species is either strain, subspecies or pseudo-strain.
	* If taxonkit indicates that a certain taxonomy ID could not be found, make sure to use the current
	  version of the NCBI taxdump files (see USAGE).
	  
 KNOWN BUGS
 
	* The API may occasionally produce "wrong UID" errors. Checking the UID on the NCBI website and, if
	  valid, rerunning the script may fix the error. True errors would always appear for the same ID
	  and in the same data chunk (check debug messages), since the order of IDs will not change between
	  runs of this program, unless the input data was changed.
	  Expert users may want to check additionally
	  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=YOUR_UID
	  (replace YOUR_UID with the erroneous UID).
	  
#=========================================================================================================#
EOF
;


#-----------------------------------------------------------------------------------------------#
# Parse and check arguments
#-----------------------------------------------------------------------------------------------#
my $argC = @ARGV;

GetOptions ('input:s'	=> \$seqF,
            'kitdata:s' => \$taxkitData,
           	'help|?' => \$help,
           	'debug' => \$debug) or die $usage;
           	
if ($help > 0 or not $seqF) {
	print $usage;
	exit 0;
}

if (not $taxkitData) {
	$taxkitData = dirname($seqF);
	print "INFO: -kitdata omitted. Assuming: $taxkitData\n";
}

# Quick check to see, if taxkitData is OK. Exit 0 is error in system().
system("echo \"9606\" | taxonkit lineage --data-dir $taxkitData >/dev/null") and
	die "ERROR: TaxonKit failed. Check -kitdata or the installation.";

my $outFa = dirname($seqF) . "/metag." . basename($seqF);
my $outTax = dirname($seqF) . "/tax.metag." . basename($seqF);
$outTax =~ s/\.[^.]*$/\.txt/;


#-----------------------------------------------------------------------------------------------#
# Load fasta sequences and get NCBI accessions from headers
#-----------------------------------------------------------------------------------------------#
print "INFO: Parsing FASTA file\n";

my $id = "";
my %seqs = ();
my %ids = ();

open(FASTA, "<", $seqF) or die "ERROR: Could not open input FASTA file $seqF.";
while (<FASTA>) {
	if ($_ =~ m/^>/) {
		$id = (split(' ', $_))[0];
		$id =~ s/^>//;
		
		my $acc = "";
		if ($id =~ m/^([A-Z]+_{0,1}[0-9]+)/) {
			$acc = $1
		}
		else {
			die "ERROR: Unexpected ID: ->" . $id ."<-";
		}
		
		# Save raw sequence ids
		if (exists $ids{$acc}) {
			push(@{$ids{$acc}}, $id);
		}
		else {
			$ids{$acc} = [$id];
		}
	}
	else {
		chomp($_);
		
		if (not $id) {
			die "ERROR: Unexpected FASTA: Sequence not preceeded by header in line " . $. . "\n";
		}
		
		if (exists $seqs{$id}) {
			$seqs{$id} .= $_;
		}
		else {
			$seqs{$id} = $_;
		}
	}
}
close(FASTA);


#-----------------------------------------------------------------------------------------------#
# Get taxids from eUTILs using sequence accessions
#-----------------------------------------------------------------------------------------------#
my %taxids = ();
my @accs = keys(%ids);

# Max is 500 for JSON
my $chunkSize = 450;

print "INFO: Web request to NCBI: " . @accs .
	" sequences. Sending chunks of $chunkSize. This may take a while.\n";

my $chunkC = 1;

while (@accs) {
	print "DEBUG: Querying chunk $chunkC\n" if ($debug == 1);
	
    my @queries = splice (@accs, 0, $chunkSize);	
	my %queries = ('id' => join(",", @queries));
	
	my $ua = LWP::UserAgent->new(timeout => 40);
	
	# Try to query a chunk again, if there was any error.
	# Sometimes NCBI API returns "random" error for valid ID.
	# Retry at most 3 times.
	my $retryC = 0;
	my $maxRetry = 3;
	my $isSuccess = 0;
	while ($isSuccess == 0 and $retryC < $maxRetry) {
		my $response = $ua->post($url, \%queries);
	
		if (not $response->is_success) {
			print "ERROR: API error with chunk $chunkC: " . $response->status_line . "\n"
		}
		else {
			if (not defined $response->content or not $response->content) {
				print "WARNING: Unexpected API response\n";
			}
			else {
				$response = decode_json ($response->content);
				
				if (exists $response->{'error'}) {
					print "API ERROR: " . $response->{'error'} . "\n";
				}
				else {
					if (defined $response->{'result'}->{'uids'}) {
					
						foreach my $uid (@{$response->{'result'}->{'uids'}}) {
							my $taxid = 0;
							
							if (exists $response->{'result'}->{$uid}->{'taxid'}) {
								$taxid = $response->{'result'}->{$uid}->{'taxid'};
								my $seqidApi = $response->{'result'}->{$uid}->{'caption'};
								my $seqidsR = $ids{$seqidApi};
								
								foreach my $seqid (@{$seqidsR}) {
									$taxids{$seqid} = $taxid;
								}
								delete $ids{$seqidApi};
								
								# Flag to break out of the while loop
								$isSuccess = 1;
							}
							else {
								print "ERROR unexpected API response. No taxid returned.\n"
							}
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
				print "INFO: Error in API response. Trying again ($retryC of $maxRetry).\n";
				sleep(8);
			}
			else {
				die "ERROR: No valid API response in retry limit\n";
			}
		} 
	}
	# Only 2 requests per second
	select(undef, undef, undef, 0.5);
	
	$chunkC++;
}

if (keys(%ids)) {
	die "ERROR: Could not get taxids for sequence ids: " . join("; ", keys(%ids)) . "\n";
}

open(OUTTAX, ">", $outTax) or die "Could not open output taxonomy file $outTax";
foreach my $seqid (keys(%taxids)) {
	print OUTTAX $seqid, "\t", $taxids{$seqid}, "\n";
}


#-----------------------------------------------------------------------------------------------#
# Use the NCBI Taxonomy Toolkit to get the lineage for each sequence id using its taxid
#-----------------------------------------------------------------------------------------------#
print "INFO: Getting lineage info from taxids\n";

my $tmpTax = dirname($outTax) . "/tmp." . basename($outTax);

# Use the taxid in the second column to write a temporary taxonomy with 8 ranks, missing ranks are "0"
my $command = ("taxonkit reformat -I 2 --data-dir $taxkitData --miss-rank-repl \"0\" -S -f \"{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}\" $outTax > $tmpTax");

# For system, exit code 0 is an error.
system($command) and die "ERROR: Taxonkit failed.";


#-----------------------------------------------------------------------------------------------#
# Reformat taxonomy output and write sequences for entries in taxonomy file
#-----------------------------------------------------------------------------------------------#
print "INFO: Writing reformatted output to $outTax\n";

open(OUTSEQ, ">",$outFa) or die "ERROR: Could not open output FASTA file $outFa.";
open(OUTTAX, ">", $outTax) or die "ERROR: Could not open final output taxonomy file $outTax";
open(TMPTAX, "<", $tmpTax) or die "ERROR: Could not open temporary taxonomy file $tmpTax";
while(<TMPTAX>) {
	chomp($_);
	my @splits = split("\t", $_);
	my $seqid = $splits[0];
	
	# Get the whole tab-delimited and reformated taxonomy
	my @taxa = @splits[2..$#splits];
	
	if (not @taxa) {
		print "DEBUG: Lost $seqid. No lineage.\n" if ($debug == 1);
		next;
	}
	
	# Reformatted taxonomy has 8 ranks. MetaG needs 10.
	# Insert subclass and suborder as "0" = unknown
	splice(@taxa, 3, 0, "0");
	splice(@taxa, 5, 0, "0");
	
	print OUTTAX $seqid, ";", join(";", @taxa), "\n";
	print OUTSEQ ">" . $seqid . "\n" . $seqs{$seqid} ."\n";
}
close(TMPTAX);
close(OUTTAX);
close(OUTSEQ);

# Cleanup
# For system, exit code 0 is an error.
system("rm $tmpTax") and die "ERROR: Could not remove temporary taxonomy file ->$tmpTax<-";

print "\nDONE\n";