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
# Build a MetaG compatible database from RefSeq's Targeted Loci
#-------------------------------------------------------------------------------------------------------#
#
# USAGE
#
#	makeNCBI.pl -input seqfile.fna [ -kitdata directory ] [--debug]
#
#	OR
#
#	makeNCBI.pl -i seqfile.fna [ -k directory ] [-d]
#
#	Choose and download the refseq database files (*fna.gz) from
#	ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/ . You may concatenate multiple files of interest into
#	one input file and supply this via -input to create a merged MetaG database. Note that the input files
#	need to be extracted.
#	Download and extract the taxdump files for the NCBI Taxonomy Toolkit
#	(see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
#	directory of your input file, you need to provide the path using -kitdata.
#	Use --debug to enable debug messages.
#
# OUTPUT
#
#	A fasta file with "metag." prefix and a taxonomy file with "tax." prefix. Both will be created in the
#	directory of your input file.
#
# DEPENDENCIES
#
# 	* NCBI Taxonomy Toolkit v0.7.2 (https://doi.org/10.1101/513523)
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
# Build a MetaG compatible database from RefSeq's Targeted Loci
#-------------------------------------------------------------------------------------------------------#

 USAGE

	makeNCBI.pl -input seqfile.fna [ -kitdata directory ] [--debug]

	OR

	makeNCBI.pl -i seqfile.fna [ -k directory ] [-d]

	Choose and download the refseq database files (*fna.gz) from
	ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/ . You may concatenate multiple files of interest into
	one input file and supply this via -input to create a merged MetaG database. Note that the input files
	need to be extracted.
	Download and extract the taxdump files for the NCBI Taxonomy Toolkit
	(see https://bioinf.shenwei.me/taxonkit/usage/#before-use). If the extracted files are not in the
	directory of your input file, you need to provide the path using -kitdata.
	Use --debug to enable debug messages.

 OUTPUT

	A fasta file with "metag." prefix and a taxonomy file with "tax." prefix. Both will be created in the
	directory of your input file.

 DEPENDENCIES

 	* NCBI Taxonomy Toolkit v0.7.2 (https://doi.org/10.1101/513523)
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

GetOptions ("input:s"   => \$seqF,
            "kitdata:s" => \$taxkitData,
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

my $outFasta = dirname($seqF) . "/metag." . basename($seqF);
my $outTax = dirname($seqF) . "/tax.metag." . basename($seqF);
my $header = "";
my $seq = "";
my @ids = ();
my %ids = ();
my %taxids = ();


#-----------------------------------------------------------------------------------------------#
# Reformat the fasta file
#-----------------------------------------------------------------------------------------------#
print "INFO: Making fasta file\n";

open(OUTFA, ">", $outFasta) or die "Could not open output fasta file";
open(SEQ, "<", $seqF) or die "Could not open sequence fasta file $seqF";
while(<SEQ>) {
	chomp($_);
	
	if ($_ =~ m/^>/) {
		$header = (split(" ", $_))[0];
		my $id = $header =~ s/^>//r;
				
		# No version needed to check, if all ids were returned by the API
		# API looses version tag at the end, but I need the original
		# sequence ids lateron.
		my $idNoVers = $id =~ s/\..*$//r;
		
		# Save sequence ids for lineage extraction
		if (exists $ids{$idNoVers}) {
			die "ERROR: RefSeq id must be unique, but is not: ->$idNoVers<-"
		}
		else {
			$ids{$idNoVers} = $id;
			push (@ids, $id);
		}
		
		print OUTFA $seq, "\n" if ($seq);
		$seq = "";
		
		print OUTFA $header, "\n";
	}
	else {
		$seq .= $_;
	}
}
print OUTFA $seq, "\n" if ($seq);
$seq = "";

close(SEQ);
close(OUTFA);


#-----------------------------------------------------------------------------------------------#
# Get taxids from eUTILs using sequence ids
#-----------------------------------------------------------------------------------------------#
# Max is 500 for JSON
my $chunkSize = 450;

print "INFO: Web request to NCBI: ".@ids.
	" sequences. Sending chunks of $chunkSize. This may take a while.\n";

my $chunkC = 1;

while (@ids) {
	print "DEBUG: Querying chunk $chunkC\n" if ($debug == 1);
	
    my @queries = splice (@ids, 0, $chunkSize);	
	my %queries = ('id' => join(",", @queries));
	
	my $ua = LWP::UserAgent->new(timeout => 40);
	my $response = $ua->post($url, \%queries);

	if (not $response->is_success) {
		die "ERROR: API error with chunk $chunkC: " . $response->status_line . "\n"
	}
	else {
		if (not defined $response->content or not $response->content) {
			print "WARNING: Unexpected API response";
		}
				
		$response = decode_json ($response->content);
		
		if (exists $response->{'error'}) {
			die "API ERROR: " . $response->{'error'};
		}
		else {
			if (defined $response->{'result'}->{'uids'}) {
			
				foreach my $uid (@{$response->{'result'}->{'uids'}}) {
					my $taxid = 0;
					
					if (exists $response->{'result'}->{$uid}->{'taxid'}) {
						$taxid = $response->{'result'}->{$uid}->{'taxid'};
						my $seqidApi = $response->{'result'}->{$uid}->{'caption'};
						my $seqid = $ids{$seqidApi};
						$taxids{$seqid} = $taxid;
						delete $ids{$seqidApi};
					}
					else {
						die "ERROR unexpected API response. No taxid returned."
					}
				}
			}
			else {
				print "WARNING: Unexpected API response\n";
			}
		}
	}
	# Only 2 requests per second
	select(undef, undef, undef, 0.5);
	
	$chunkC++;
}

if (keys(%ids)) {
	die "ERROR: Could not get taxid for sequence ids: " . join("; ", keys(%ids)) . "\n";
}

open(OUTTAX, ">", $outTax) or die "Could not open output taxonomy file $outTax";
foreach my $seqid (keys(%taxids)) {
	print OUTTAX $seqid, "\t", $taxids{$seqid}, "\n";
}


#-----------------------------------------------------------------------------------------------#
# Use the NCBI Taxonomy Toolkit to get the lineage for each sequence id using its taxid
#-----------------------------------------------------------------------------------------------#
print "INFO: Getting lineage info from taxid\n";

my $tmpTax = dirname($outTax) . "/tmp." . basename($outTax);

# Use the taxid in the second column to write a temporary taxonomy
# Then, the taxonomy is formatted to 8 ranks, missing ranks are "0"
my $command = ("taxonkit lineage -i 2 $outTax --data-dir $taxkitData | taxonkit reformat -i 3 --data-dir $taxkitData --miss-rank-repl \"0\" -S -f \"{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}\" > $tmpTax");

# For system, exit code 0 is an error.
system($command) and die "ERROR: Taxonkit failed.";


#-----------------------------------------------------------------------------------------------#
# Reformat output
#-----------------------------------------------------------------------------------------------#
print "INFO: Writing reformatted output to $outTax\n";

open(OUTTAX, ">", $outTax) or die "Could not open final output taxonomy file $outTax";
open(TMPTAX, "<", $tmpTax) or die "Could not open temporary taxonomy file $tmpTax";
while(<TMPTAX>) {
	chomp($_);
	my @splits = split("\t", $_);
	my $seqid = $splits[0];
	
	# Get the whole tab-delimited and reformated taxonomy
	my @taxa = @splits[3..$#splits];
	
	# Reformatted taxonomy has 8 ranks. MetaG needs 10.
	# Insert subclass and suborder as "0" = unknown
	splice(@taxa, 3, 0, "0");
	splice(@taxa, 5, 0, "0");
	
	print OUTTAX $seqid, ";", join(";", @taxa), "\n";
	
}
close(TMPTAX);
close(OUTTAX);

# Cleanup
# For system, exit code 0 is an error.
system("rm $tmpTax") and die "ERROR: Could not remove temporary taxonomy file";

print "\nDONE\n";