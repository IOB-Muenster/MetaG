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


#========================================================================#
# This script simulates all possible combinations of -e, -ac, -cc
# in a custom range. It rates the results for a chosen rank based on
# the expected total abundances of taxa.
#
# OUTPUT
#
# E\tAC\tCC\ttaxon\tfound_expec\t%correct
# 
# found_expec: total abundances for the taxon
# %correct: foundAbund(taxon) / expecAbund(taxon)
# mock taxon "average": arithmetic mean of all %correct
#
# Rating for the UNMATCHED class also acknowledges reads which were
# removed to the LAST settings of e-value cutoff.
#
# REQUIREMENTS
# You need the to install the CPAN modules:
# Parallel::ForkManager and Algorithm::Loops
# MetaG must be installed. This script must be in the same folder as
# metag.sh and you need to provide the environment variables GETMETA and
# KRONA (see documentation of MetaG). Additionally, a precalculated
# alignment in MAF format is mandatory. You can choose the analyzed rank
# for the fasta query by modifying $reqRank.
#
# The file with the expected abundances must be in this format:
# >rank1
# taxonA:countA
# taxonB:countB
# >rank2
# ...
#
#------------------------------------------------------------------------#

use strict;
use warnings;
use Parallel::ForkManager;
use Algorithm::Loops qw( NestedLoops );
use File::Basename;


# Use unbuffered stdout
select(STDOUT);
$| = 1;

#==========================================================================#
# MODIFY HERE
#--------------------------------------------------------------------------#
my $reqRank = "species";



#=========================================================================#
# Subroutines
#-------------------------------------------------------------------------#

sub getStat{
	
	my $linkPath = $_[0];
	my $paramCom = $_[1];
	my $reqRank = $_[2];
	my $exp = $_[3];
	my %exp = %$exp;
	my @expected = keys(%exp);
	
	my %prints = ();
	
	my $isReq = 0;
	
	my ($e, $ac, $cc) = @$paramCom;
	
	my $removedC = 0;
	
	# Found taxa are in same directory as link to query file
	my $foundPath = $linkPath =~ s/query.txt$/calc.LOG.txt/r;
	
	# Look at found taxa at requested rank. Calculate ranking stats
	open (FOUND, "<", $foundPath) or die "No file with found taxa, $!";
	while (<FOUND>) {
		if ($_ =~ m/^[a-zA-Z]/) {
			
			# Get number of reads removed by LAST or e-value cutoff.
			# These need to be added to the unmatched class at the 
			# required taxonomic level to get the total amount of 
			# unmatched reads.
			if ($_ =~ m/^Removed/) {
				chomp($_);
				$removedC = (split(" ", $_))[1];
			}
			
			
			# Switch to look for taxa only at requested rank
			$isReq = 0;
			
			if ($_ =~ m/^$reqRank/) {
				$isReq = 1;
			}
		}
		else {
			if ($isReq == 1) {
				chomp($_);
				
				foreach my $expected (@expected) {
					if ($_ =~ m/$expected/) {
						my $abund = (split("\t", $_))[2];
						
						# Add number of reads removed by LAST and ec
						# to the UNMATCHED taxa.
						if ("$expected" eq "UNMATCHED") {
							$abund += $removedC
						}
						
						# Taxon => [found abund_expected abund, %correct]
						my $perc = sprintf("%.3f",$abund/$exp{$expected});
						$prints {$expected} = [$abund."_".$exp{$expected}, $perc];
					}
				}
			}
		}
	}
	close(FOUND);

	my $average = 0;
	my $expecNumb = @expected;
	my $result = "";
	
	# Calculate average and record ranking stats per taxon
	foreach my $expected (@expected) {
		
		# Expected taxon was found
		if (defined $prints{$expected}) {
			$result.= $e."\t". $ac."\t". $cc."\t".$expected."\t".@{$prints{$expected}}[0]."\t".@{$prints{$expected}}[1]."\n";
			$average += @{$prints{$expected}}[1];
		}
		# Expected taxon has not been found
		else {
			
			# Consider number of reads removed by LAST and ec
			# for UNMATCHED taxa.
			if ("$expected" eq "UNMATCHED") {
				my $perc = sprintf("%.3f",$removedC/$exp{$expected});
				$result.= $e."\t". $ac."\t". $cc."\t".$expected."\t".$removedC."_".$exp{$expected}."\t".$perc."\n";
				$average += $perc;
			}
			else {
				$result.= $e."\t". $ac."\t". $cc."\t".$expected."\t"."0_".$exp{$expected}."\t0.000\n";
				$average += 0.000;
			}
		}
	}
	$average = $average / $expecNumb;
	$average = sprintf("%.3f",$average);
	
	# Write average
	$result.= $e."\t". $ac."\t". $cc."\taverage\tNA\t".$average."\n";
	
	return \$result;
}

#==================================================================================#
# Get command line arguments
#----------------------------------------------------------------------------------#

my $queryFile = undef;
my $alignFile = undef;
my $eRange = undef;
my $dbPath = undef;
my $expecTax = undef;
my $pdbPath = undef;
my $forks = 0;
my $dir = dirname(__FILE__);

my @acRange = undef;
my @ccRange = undef;


my$argLen = @ARGV;
for (my $i=0; $i <= $#ARGV; $i++) {
	if ($argLen == 16 or $argLen == 18) {
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
		elsif ($ARGV[$i] eq "-pdbPath") {
			$pdbPath = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] eq "-forks") {
			$forks = $ARGV[$i+1];
		}
		elsif ($ARGV[$i] =~ m/^-[a-zA-Z]/) {
			die "ERROR: Unknown option $ARGV[$i]. Use -h for help."
		}
	}
	else {
		if ($ARGV[$i] eq "-h") {
			print "\n#======================================================================================#";
			print "\n#=================================== TRAIN MetaG ======================================#";
			print "\n#======================================================================================#";
			print "\nThis script simulates all possible combinations of -e, -ac, -cc
in a custom range. It rates the results for a chosen rank based on
the expected total abundances of taxa.\n\n";
			print "#ARGUMENTS\n";
			print "(-q) Query fasta file\n";
			print "(-a) Alignment maf file\n";
			print "(-e) E-value exponent range: lower_upper\n";
			print "(-ac) Alignment score range: lower_upper\n";
			print "(-cc) Confidence range: lower_upper\n";
			print "(-db) Path to database taxonomy\n";
			print "(-expec) Expected taxonomy file\n";
			print "(-pdbPath) Path to pathogen database\n";
			print "(-forks) OPTIONAL: Specify number of parallel processes. Default 1.\n";
			print "
# OUTPUT
E\tAC\tCC\ttaxon\tfound_expec\t%correct

found_expec: total abundances for the taxon
%correct: foundAbund(taxon) / expecAbund(taxon)
mock taxon \"average\": arithmetic mean of all %correct

Rating for the UNMATCHED class also acknowledges reads which were
removed to the LAST settings of e-value cutoff.

# REQUIREMENTS
You need the to install the CPAN modules: Parallel::ForkManager and Algorithm::Loops
MetaG must be installed. This script must be in the same folder as
metag.sh and you need to provide the environment variables GETMETA and
KRONA (see documentation of MetaG). Additionally, a precalculated
alignment in MAF format is mandatory. You can choose the analyzed rank
for the fasta query by modifying \$reqRank.

The file with the expected abundances must be in this format:
>rank1
taxonA:countA
taxonB:countB
>rank2
...
";
			exit(0);
		}
		else {
			die "Expected 16 or 18 arguments, $argLen given. Use -h for help.\n";
		}
	}
}

if (defined $queryFile and defined $alignFile and defined $eRange and @acRange and @ccRange and defined $dbPath and defined $expecTax and defined $pdbPath) {
	print "\n# Query file: $queryFile";
	print "\n# Alignment file: $alignFile";
	print "\n# E-value range: @{$eRange}[0] to @{$eRange}[-1]";
	print "\n# Alignment score range: ".$acRange[0]/10.0." to ".$acRange[-1]/10.0;
	print "\n# Confidence range: ".$ccRange[0]/10.0." to ".$ccRange[-1]/10.0;
	print "\n# Path to restructured DB: $dbPath";
	print "\n# Path to expected taxonomy file: $expecTax";
	print "\n# Path to pathogen taxonomy: $pdbPath";
	print "\n# Number of parallel processes: $forks\n";
}
else {
	die "Undefined essential parameters (-q, -a, -e, -ac, -cc, -db, -expec, pdbPath). Use -h for help.\n"
}

#===========================================================================#
# Read expected taxonomy
#
# >Rank
# taxon1:expecAbund
# taxon2:expecAbund
# >Rank
# ...
#---------------------------------------------------------------------------#

my %exp =();
my $rank = undef;

open(EXP, "<", $expecTax) or die "Could not open expected taxonomy";
while(<EXP>) {
	
	# Get taxon
	if ($_ =~ m/^>/) {
		chomp($_);
		$_ =~ s/>//;
		$rank = $_;
	}
	# Get expected abundance
	else{
		if (defined $rank and $rank eq $reqRank) {
			chomp($_);
			my ($taxon, $abund) = (split(":", $_));
			if (exists $exp{$taxon}) {
				$exp{$taxon}+= $abund;
			}
			else {
				$exp{$taxon} = $abund;
			}
		}
	}
}
close(EXP);

#========================================================================================================================#
# Launch metag with variable parameters. Use forking.
#------------------------------------------------------------------------------------------------------------------------#

# Spawn max n processes
my $pm = Parallel::ForkManager->new($forks);

# Run this to retrieve data from children
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


# Get all combinations of three parameters. E; AC; CC
my @limits = ($eRange, [map{$_/10} @acRange], [map{$_/10} @ccRange]);
my @combs; NestedLoops(\@limits, sub { push @combs, [ @_ ] });
my $progrC = 0;
my $combNumb = @combs;

print "# Number of parameter combinations: $combNumb\n";
print "# Checking at taxonomic level: $reqRank\n\n";

# Write all output to directories starting at query file location
my $targetDir = dirname($queryFile);

# Prepare for forking
# Temporary directory containing single metag runs. Will contain subdirectories for each run
my $tempPath = $targetDir."/temp.results";

if (-d $tempPath) {
		die "ERROR: Temporary directory $tempPath should be unique, but exists. Please delete it manually.";
}
else {
	`mkdir $tempPath`;
}

print "# E\tAC\tCC\ttaxon\tfound_expec\t%correct\n";

# Loop over all possible parameter combinations using forked processes
COMBS: for (my $i=0; $i<=$#combs; $i++) {
	
	#============================================#
	# Forked process in single temp directories
	#--------------------------------------------#
	$pm->start and next COMBS; # do the fork
	my $paramCom = $combs[$i];
	
	`mkdir $tempPath/$i`;
	
	# Link query file to temporary dir of process
	# Metag writes output to location of link
	my $linkPath = $tempPath."/".$i."/query.txt";
	$queryFile = (fileparse($queryFile))[0];
	$queryFile = "../../$queryFile";
	
	
	symlink($queryFile, $linkPath);
	
	# Run metag
	my ($e, $ac, $cc) = @$paramCom;
	`$dir/metag.sh -q $linkPath --no_align $alignFile -ltax $dbPath -pdbPath $pdbPath -vir -e $e -ac $ac -cc $cc -m williams`;
	
	# Retrieve scores
	my $result = getStat($linkPath, $paramCom, $reqRank, \%exp);
	
	`rm -r $tempPath/$i`;
	
	# Pass string ref to run_on_finish
	$pm -> finish(0, $result);
	
	#=========================================#

}
$pm->wait_all_children;

# Cleanup temp dir
`rmdir $tempPath`