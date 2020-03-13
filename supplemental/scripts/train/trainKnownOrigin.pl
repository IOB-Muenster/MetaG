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


#===========================================================================#
# This script trains MetaG on a sample in which the origin of each
# fasta read is known. It will return results to STDOUT in the format:
# -e\t-ac\t-cc\tMCC\tsensitivity\tspecificity\tprecision\tTP\tFP\tTN\tFN
#
# REQUIREMENTS
# You need the to install the CPAN modules:
# Parallel::ForkManager and Algorithm::Loops
# MetaG's "modules" folder must be located in the same directory as this
# script. Additionally, a precalculated alignment in MAF format is mandatory.
# You can choose the analyzed rank by modifying $reqRank.
#
# Reads from one source should all start with a common prefix, which is
# separated from the rest of the read ID by a "_", e.g. taxonA_123. Reads
# which were simulated to be unalignable must contain "unaligned" within
# the ID. The file with the expected taxonomy must use the common prefix to
# define a taxonomy.
#
# >taxonA
# domain: A
# pyhlum: AA
# ...
# >taxonB
# ...
#
# Additionally, the total number of reads and the total number of unaligned
# reads must be specified in the file with the expected taxonomy:
#
# >#total 100
# >#unaligned 1
#
# KNOWN BUGS
# May cause core dump on CentOS
#---------------------------------------------------------------------------#

use FindBin;
FindBin::again();
use lib $FindBin::Bin;

use strict;
use warnings;
use Algorithm::Loops qw( NestedLoops );
use modules::ProcessAlign qw(readAlign processAlign);
use modules::ReadDB;
use modules::ReadQuery;
use Parallel::ForkManager;
use File::Basename;
use Cwd;

# Use unbuffered stdout
select(STDOUT);
$| = 1;


my $cwd= getcwd;
my $dir = dirname(__FILE__);


#==========================================================================#
# MODIFY HERE
#--------------------------------------------------------------------------#
my $reqRank = "genus";



#======================================================#
# Calculate confidence
#------------------------------------------------------#

sub getConf {
	#$taxaStat{$taxon} = [$eValue, $aScore, $taxCount]
	my %taxaStat = %{$_[0]};
	
	#One taxon means highest confidence, exit sub
	my $taxNum = keys(%taxaStat);
	return 1 if ($taxNum == 1);
	
	my $maxTaxon = $_[1];
	my $totalScoring = 0;
	my $maxScore = 4;
	
	# Check, if the difference between taxon with most and second most count is too high,
	# i.e., second most count is equal to or less than 50% of most count.
	# Exit directly.
	my @sorted = sort {$taxaStat{$b}->[2] <=> $taxaStat{$a}->[2]} keys(%taxaStat);
	return 1 if  ($taxaStat{$sorted[0]}->[2]*0.5 >=  $taxaStat{$sorted[1]}->[2]);
	
	my $sort = sub {
		$taxaStat{$a}->[0] <=> $taxaStat{$b}->[0]  # by e-value
		or
	    $taxaStat{$b}->[1] <=> $taxaStat{$a}->[1]  # if equal: By alignment score
	    or 
	    $taxaStat{$b}->[2] <=> $taxaStat{$a}->[2]  # if equal: By taxCount
	};
	
	@sorted = sort $sort keys(%taxaStat);
	
	# Difference between first and last taxon in sorted
	
	# Do e-values differ?
	my $index = 0; 
	
	my $diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
	
	# E-value identical, but what about alignment score?
	if ($diff == 0) {
		$index = 1;
		$diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
	}
	# Alignment score identical, but what about taxCount?
	if ($diff == 0) {
		$index = 2;
		$diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
	}
	
	# All the same: Should have been filtered out, previously by maxTaxon = "-"
	die "Identical inputs in confidence calculation" if ($diff == 0);
	
	my $bestVal = $taxaStat{$sorted[0]}->[$index];
	my $conf = 0;
	
	# Make bins and distribute weights
	foreach my $taxon (@sorted) {
		my $value = $taxaStat{$taxon}->[$index];
		my $weight = undef;
		my $firstThird = $bestVal - ($diff / 3);
		my $secThird = $bestVal - (2 * $diff / 3);
			
		# For e-value
		if ($index == 0) {
			
			# ++
			if ($value == $bestVal) {
				$weight = $maxScore 
			}
			# +
			elsif ($value > $bestVal and $value <= $firstThird) {
				$weight = 2
			}
			# -
			elsif ($value > $firstThird and $value <= $secThird) {
				$weight = 1
			}
			# --
			else {
				$weight = 0
			}
			
		}
		# For others
		else {
			
			# ++
			if($value == $bestVal) {
				$weight = $maxScore
			}
			# +
			elsif($value < $bestVal and $value >= $firstThird) {
				$weight = 2
			}
			# -
			elsif($value < $firstThird and $value >= $secThird) {
				$weight = 1
			}
			# --
			else{
				$weight = 0
			}
		}
		
		# Count number of taxa contributing to weight
		$totalScoring += $taxaStat{$taxon}->[2] if ($weight > 0);
		
		if ($taxon eq $maxTaxon) {
			$conf += $weight * $taxaStat{$taxon}->[2]
		}
		else {
			$conf -= $weight * $taxaStat{$taxon}->[2]
		}
	}
	# Calculate realtive score
	$conf = $conf / $totalScoring;
	
	# Rescale to a value between 0 and 1
	# Conf can ranke from -maxScore to + maxScore
	$conf = ($conf + $maxScore) / ($maxScore*2);
	
	return $conf;
}


my $queryFile = undef;
my $alignFile = undef;
my $forks = 1;
my $eRange = undef;
my $dbPath = undef;
my $expecTax = undef;
my @acRange = undef;
my @ccRange = undef;
my $totalUnalignable = 0;
my $totalReads = 0;


my$argLen = @ARGV;
for (my $i=0; $i <= $#ARGV; $i++) {
	if ($argLen == 16) {
		if ($ARGV[$i] eq "-q") {
			$queryFile = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-a") {
			$alignFile = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-e") {
			my @eRange = ((split("_", $ARGV[$i+1]))[0]..(split("_", $ARGV[$i+1]))[-1]); #lower_upper
			$eRange = \@eRange;
		}
		elsif ($ARGV[$i] eq "-ac") {
			# Get decimal steps for nested loops
			@acRange = ((split("_", $ARGV[$i+1]))[0] * 10 .. (split("_", $ARGV[$i+1]))[-1] * 10); #lower_upper
		}
		elsif ($ARGV[$i] eq "-cc") {
			# Get decimal steps for nested loops
			@ccRange = ((split("_", $ARGV[$i+1]))[0] * 10 .. (split("_", $ARGV[$i+1]))[-1] * 10); #lower_upper
		}
		elsif ($ARGV[$i] eq "-db") {
			$dbPath = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-expec") {
			$expecTax = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-forks") {
			$forks = $ARGV[$i+1];
		}
	}
	else {
		if ($ARGV[$i] eq "-h") {
			print "\n#======================================================================================#";
			print "\n#=================================== TRAIN MetaG ======================================#";
			print "\n#======================================================================================#";
			print "\nThis script trains MetaG on a sample in which the origin of each fasta read is known\n";
			print "\n#ARGUMENTS\n";
			print "(-q) Query fasta file\n";
			print "(-a) Alignment maf file\n";
			print "(-e) E-value exponent range: lower_upper\n";
			print "(-ac) Alignment score range: lower_upper\n";
			print "(-cc) Confidence range: lower_upper\n";
			print "(-db) Database Taxonomy\n";
			print "(-expec) Expected taxonomy file\n";
			print "(-forks) OPTIONAL: Specify number of parallel processes. Default 1.\n";
			print "
# OUTPUT
The script will return results to STDOUT in the format:
-e\t-ac\t-cc\tMCC\tsensitivity\tspecificity\tprecision\tTP\tFP\tTN\tFN

# REQUIREMENTS
You need the to install the CPAN modules: Parallel::ForkManager and Algorithm::Loops
MetaG's \"modules\" folder must be located in the same directory as this
script. Additionally, a precalculated alignment in MAF format is mandatory.
You can choose the analyzed rank by modifying \$reqRank.

Reads from one source should all start with a common prefix, which is
separated from the rest of the read ID by a \"_\", e.g. taxonA_123. Reads
which were simulated to be unalignable must contain \"unaligned\" within
the ID. The file with the expected taxonomy must use the common prefix to
define a taxonomy.

>taxonA
domain: A
pyhlum: AA
...
>taxonB
...

Additionally, the total number of reads and the total number of unaligned
reads must be specified in the file with the expected taxonomy:

>#total 100
>#unaligned 1

# KNOWN BUGS:
May cause core dump on CentOS with Slurm workload manager.		
";
			exit(0);
		}
		else {
			die "Expected 16 arguments, $argLen given. Use -h for help.";
		}
	}
}

if (defined $queryFile and defined $alignFile and defined $eRange and @acRange and @ccRange and defined $dbPath and defined $expecTax) {
	print "\nQuery file: $queryFile";
	print "\nAlignment file: $alignFile";
	print "\nExpected taxonomy: $expecTax";
	print "\nE-value range: @{$eRange}[0] to @{$eRange}[-1]";
	print "\nAlignment score range: ".$acRange[0]/10.0." to ".$acRange[-1]/10.0;
	print "\nConfidence range: ".$ccRange[0]/10.0." to ".$ccRange[-1]/10.0;
	print "\nForks: $forks";
	print "\nDatabase taxonomy: $dbPath\n";
}
else {
	die "Undefined essential parameters (-q, -a, -e, -ac, -cc, -db, -expec). Use -h for help.\n"
}

#===========================================================================#
# Read expected taxonomy
#---------------------------------------------------------------------------#

my %exp =();
my $id = undef;
open(EXP, "<", $expecTax) or die "Could not open expected taxonomy";
while(<EXP>) {
	
	# Get ID
	if ($_ =~ m/^>/) {
		chomp($_);
		$id = $_;
		$id =~ s/>//g;
		
		if ( $id =~ m/^#unaligned/ ){
			$totalUnalignable = (split(" ", $id))[1];
			next;
		}

		if ( $id =~ m/^#total/ ){
			$totalReads = (split(" ", $id))[1];
			next;
		}
	}
	# Get taxon
	else {
		if ($_ =~ m/^$reqRank/) {
			chomp($_);
			my $taxon = (split(":", $_))[1];
			$taxon =~ s/\W//g;
			$taxon =~ s/_//g;
			
			$exp{$id} = $taxon
			
		}
	}
}
close(EXP);

print "\nTotal reads: $totalReads";
print "\nTotal unalignable reads: $totalUnalignable\n\n";

#==========================================================================#
# Read ids from 16S input sequences
#--------------------------------------------------------------------------#
my ($queriesR, $allQ) = readQuery ($queryFile);

#=============================================================================================#
# Load database in fasta format into memory
#---------------------------------------------------------------------------------------------#
my ($db, $ranks) = readDB($dbPath);

my @ranks = ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species", "strain");


#==========================#
# Read maf alignment once.
#--------------------------#
my $seqHashR = readAlign($alignFile);


#========================================================================================================================#
# Launch forks with variable parameters
#------------------------------------------------------------------------------------------------------------------------#
# Write to logfile

print "#-e\t-ac\t-cc\tMCC\tsensitivity\tspecificity\tprecision\tTP\tFP\tTN\tFN\n";

#my $tmpdir = "$dir/temp/";

# Spawn max n processes
my $pm = Parallel::ForkManager->new($forks);

# Get data from forks
$pm -> run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
	 
	    # Retrieve data structure from child
	    if (defined($data)) {  # children are not forced to send anything
	      my $result = ${$data};  # child passed a string reference
	      print "$result";
	    }
	    else {  # problems occurring during storage or retrieval will throw a warning
	      die qq|No message received from child process $pid!\n|;
	    }
	}
);


# Get combinations of three parameters. E; AC; CC
my @limits = ($eRange, [map{$_/10} @acRange], [map{$_/10} @ccRange]);
my @combs; NestedLoops(\@limits, sub { push @combs, [ @_ ] });
my $combNumb = @combs;
COMBS: foreach my $paramCom (@combs) {
	
	$pm->start and next COMBS; # do the fork
	
	my ($eCut, $alignThresh, $confCut) = @$paramCom;
	my $ePrint = $eCut;
	$eCut = 1*exp($eCut);
	
	my ($cutSeqHash, $lostEval) = processAlign ($seqHashR, $eCut, $alignThresh);
	
	my ($sensi, $specifi, $preci, $TP, $FP, $TN, $FN, $MCC) = getSpec($cutSeqHash, $confCut, $queriesR, $db, $reqRank, \%exp, $totalUnalignable, $totalReads);
	
	my $result =  "$ePrint\t$alignThresh\t$confCut\t$MCC\t$sensi\t$specifi\t$preci\t$TP\t$FP\t$TN\t$FN\n";
	
	# Pass string ref to run_on_finish
	$pm -> finish(0, \$result);
}

$pm->wait_all_children;

#=========================================================================================================================#
# Assign queries to taxa
#-------------------------------------------------------------------------------------------------------------------------#
sub getSpec {
	my $totalRanks = 0;
	my ($seqHashR, $confCut, $queriesR, $db, $reqRank, $expR, $totalUnalignable, $totalReads) = (@_);
	my %seqHash = %{$seqHashR};
	my @queries = @$queriesR;
	my %dbHash = %$db;
	my @dbKeys = keys(%dbHash);
	
	my %res = ();
	my %exp = %{$expR};
	
	
	my @seqHashKeys = keys(%seqHash);
	
	my $allQ = @queries;
	
	
	#======================================================================================================#
	# Find the lca (lowest common ancestor; most broad shared ancestor) of all LAST matches of a sample
	#------------------------------------------------------------------------------------------------------#
	
	my @ranks = ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species"); #,"strain"
	
	
	# Check the ranks in descending order: From broad to specific
	foreach my $query (@queries) { #seqHashKeys
		
		# Remaining ranks for majority analysis
		if (exists $seqHash{$query}) {
			
			my $matchCount = keys(%{$seqHash{$query}});
						
			##
			my $id = (split("_", $query))[0];
			if ($query =~ m/unaligned/) {
				$id .= "_unaligned"
			}
			$id =~ s/>//g;
			$id =~ s/-//g;
			##
			
			#====================================================================================#
			# Majority decision for ranks = lca, if confCut = 1.0.
			#-------------------------------------------------------------------------------------#
			
			my @DBids = keys(%{$seqHash{$query}});
			my $confidence = undef;
			RANK: for (my $i=0; $i <= $#ranks; $i++) {
				
				$totalRanks++;
				
				my $rank = $ranks[$i];
				my %ranks=();	
			
				# Get each LAST match dbID
				ID: foreach my $dbID (@DBids) {
					
					# Check information in db
					if (exists ($dbHash{$dbID})) { #dbHash: $dbID => {$rank => $taxon}
								
							# Is the rank noted for the id?
							if (exists $dbHash{$dbID}{$rank}) {
								
								my $taxon = $dbHash{$dbID}{$rank};
								
								if (exists $ranks{$taxon}) { 
									
									# Count numbers of taxon-supporting matches from the same DBid.
									# However, similar DBids merged to average before. Count should always be 1.
									if (exists $ranks{$taxon}{$dbID}) {
										$ranks{$taxon}{$dbID} += 1;
									}
									else {
										$ranks{$taxon}{$dbID} = 1
									}
								}
								else {
									$ranks{$taxon} = {$dbID => 1}
								}
							}
					}
				}
				my @taxa = keys(%ranks);
				my $maxCount = undef;
				my $totalCount = 0;
				my $maxTaxon = undef;
				my %taxaStat = ();
				foreach my $taxon (@taxa) {
					
					my @IDs = keys(%{$ranks{$taxon}});
					my $taxCount = 0;
					my $eValue = 0;
					my $aScore = 0;
					my $lenIDs = @IDs;
					
					foreach my $seqID (@IDs) {
						
						$eValue += $seqHash{$query}{$seqID}->[0];
						$aScore += $seqHash{$query}{$seqID}->[1];
						
						# Total count for one rank
						$totalCount += $ranks{$taxon}{$seqID};
						
						# Number of matches supporting the taxon over all matching DBids
						$taxCount += $ranks{$taxon}{$seqID};
						
						
					}
					# Mean eValue, aScore per taxon
					$eValue = $eValue / $lenIDs;
					$aScore = $aScore / $lenIDs;
					
					if ($taxon eq "0") {
						$taxon = "unclassified"
					}
					
					$taxaStat{$taxon} = [$eValue, $aScore, $taxCount];
						
					if (defined $maxCount) {
						if ($taxCount > $maxCount) {
							$maxCount = $taxCount;
							$maxTaxon = $taxon;
							
							# Only dbIDs from maxTaxon should be considered for lower taxa search
							@DBids = @IDs;
							
						}
						elsif ($taxCount == $maxCount) {
							
							# Matches have an equal count. How about the alignment quality?
							# E-value
							if ($taxaStat{$taxon}->[0] < $taxaStat{$maxTaxon}->[0]) {
								$maxTaxon = $taxon;
								@DBids = @IDs;
							}
							elsif ($taxaStat{$taxon}->[0] == $taxaStat{$maxTaxon}->[0]){
								
								# Alignment score
								if ($taxaStat{$taxon}->[1] > $taxaStat{$maxTaxon}->[1]) {
									$maxTaxon = $taxon;
									@DBids = @IDs;
								}
								# Identical in count and alignment quality: Set maxTaxon to special character "-".
								# If local maxTaxon "-" becomes global. Break, as I cannot decide between to equally likely taxa.
								elsif ($taxaStat{$taxon}->[1] == $taxaStat{$maxTaxon}->[1]) {
									$maxTaxon = "-";
									
									# Write "-" to %taxaStat to avoid bugs with following $taxCount == $maxCount
									$taxaStat{"-"} = [$eValue, $aScore, $taxCount];
									
									# Will not be considered in the end. Just if algoirithm will be changed,
									push (@DBids, @IDs);
								}
							}
						}
					}
					else {
						$maxCount = $taxCount;
						$maxTaxon = $taxon;
						@DBids = @IDs;
					}
				}
				
				# Evaluate global maxTaxon for a rank. If it is "-", break and 
				# set this and following ranks to "UNMATCHED"
				if ($maxTaxon eq "-") {
					##
					if ($rank eq $reqRank) {
						if (exists $res{$id}) {
							if (exists $res{$id}{"UNMATCHED"}) {
								$res{$id}{"UNMATCHED"} += 1;
							}
							else {
								$res{$id}{"UNMATCHED"} = 1;
							}
						}
						else {
							$res{$id} = {"UNMATCHED" => 1}
						}
					}
					##
					last RANK;
				}
				# If we have a meaningful maxTaxon, try to remove local "-" from taxaStat
				else {
					delete $taxaStat{"-"} if (exists $taxaStat{"-"});
				}
				
				#================================================================#
				# Calculate confidences of correct assignment
				#----------------------------------------------------------------#
				
				# Probability of correct assignment
				if (defined $confidence) {
					my $lastConf = $confidence;
					# Weighted confidence
					$confidence = getConf(\%taxaStat, $maxTaxon);
					$confidence = $confidence * $lastConf;
				}
				else {
					# Weighted confidence
					$confidence = getConf(\%taxaStat, $maxTaxon);
				}
					
				if ($confidence < $confCut) {
					##
					if ($rank eq $reqRank) {
						if (exists $res{$id}) {
							if (exists $res{$id}{"UNMATCHED"}) {
								$res{$id}{"UNMATCHED"} += 1;
							}
							else {
								$res{$id}{"UNMATCHED"} = 1;
							}
						}
						else {
							$res{$id} = {"UNMATCHED" => 1}
						}
					}
					##
					last RANK;
				}
				
				
				#====================================================================================#
				# Calculate score based on absolute abundance of expected taxa and false postive taxa
				#------------------------------------------------------------------------------------#
				
				##
				
				$maxTaxon =~ s/\W//g; ##
				$maxTaxon =~ s/_//g; ##
			
				if ($rank eq $reqRank) {
					if (exists $res{$id}) {
						if (exists $res{$id}{$maxTaxon}) {
							$res{$id}{$maxTaxon} += 1;
						}
						else {
							$res{$id}{$maxTaxon} = 1;
						}
					}
					else {
						$res{$id} = {$maxTaxon => 1}
					}
				}
				##
				#next RANK if ($maxTaxon eq "unclassified");
				#next RANK if ($rank eq "strain");
			}
		}
	}
	
	
	my $TP = 0;
	my $FP = 0;
	my $TN = $totalUnalignable;
	my $FN = 0;
	
	foreach my $ID (keys(%res)) {
		
		foreach my $taxon (keys(%{$res{$ID}})) {
			
			if ($ID =~ m/unaligned$/ and $taxon ne "UNMATCHED") {
				$FP += $res{$ID}{$taxon};
				$TN -= $res{$ID}{$taxon};
				next;
			}
			# Found the expected taxon
			# Less strict matching because of nomenclature clashes
			
			next if ($ID =~ m/unaligned/);
			if (not exists $exp{$ID}) {
				die "#$ID#" 
			}
			# Found taxon expression is part of expected expression or vice versa.
			if ($exp{$ID} =~ m/$taxon/i or $taxon =~ m/$exp{$ID}/i) {
				$TP += $res{$ID}{$taxon};
			}
			# Wrong taxon or unmatched
			else {
				if ($taxon eq "UNMATCHED") {
					$FN += $res{$ID}{$taxon};
				}
				else {
					$FP += $res{$ID}{$taxon};
				}
			}
			
		}
		delete $exp{$ID};
	}
	
	
	my $total = $TP + $FP + $TN + $FN;
	$FN += ($totalReads - $total);
	$total = $TP + $FP + $TN + $FN;
	my $sensi = "NA";
	my $specifi = "NA";
	my $preci = "NA";
	my $MCC = "NA";
	eval{
		$sensi = sprintf("%.4f", $TP / ($TP +$FN));
		$specifi = sprintf("%.4f", $TN / ($TN + $FP));
		$preci = sprintf("%.4f", $TP / ($TP + $FP));
		$MCC = ($TP*$TN -$FP*$FN)/((($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))**0.5);
	};
	return ($sensi, $specifi, $preci, $TP, $FP, $TN, $FN, $MCC);
}	
