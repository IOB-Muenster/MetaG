package modules::ProcessAlign;

use modules::Misc qw(searchStarts);
use strict;
use warnings;
use Exporter;

our @ISA       = qw(Exporter);
our @EXPORT    = ();
our @EXPORT_OK = qw(readAlign processAlign);


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


#=========================#
# Read maf alignment file.
#-------------------------#

sub readAlign {

	my ($alignPath) = @_;

	open(DATA, "<", $alignPath) or die "Couldn't open maf file, $!";

	my $DBid     = undef;
	my $counter  = 0;
	my $eValue   = "";
	my $aScore   = "";
	my %missedE  = ();
	my %seqHash  = ();
	my $remCeval = 0;

	my $sCounter = 0;

	while (<DATA>) {

		# Skip comments
		if ($_ !~ m/^#/) {

			if ($_ =~ m/^a/) {

				$sCounter = 0;

				my @line_data = split(' ', $_);
				$aScore = (split('=', $line_data[1]))[1];

				# Get eValue as float not string
				$eValue = (split("=", $line_data[3]))[1];

				# Already noted as float?
				if ($eValue =~ m/e/) {
					my ($factor, $expon) = split("e", $eValue);
					$eValue = $factor * exp($expon);
				}
			}

			if ($_ =~ m/^s/) {

				my @line_data = split ' ', $_;

				# First s-line contains DBid, second one contains query id
				if ($sCounter == 0) {
					$DBid = $line_data[1];
					$sCounter++;
					next;
				}
				my $qID = $line_data[1];

				if (exists $seqHash{$qID}) {
					push(@{$seqHash{$qID}}, [$eValue, $aScore, $DBid]);
				}
				else {
					$seqHash{$qID} = [[$eValue, $aScore, $DBid]];

				}
			}
		}
	}
	close(DATA);

	return (\%seqHash);
}

#=====================================================================#
# Check alignment score threshold and e-value. And restructure seqHash.
#---------------------------------------------------------------------#

sub processAlign {

	my ($seqHashR, $eCut, $alignThresh) = @_;
	my %seqHash  = %$seqHashR;
	my $lostEval = 0;
	my @qIDs     = keys(%seqHash);

	foreach my $qID (@qIDs) {

		# Each element looks like this. [[$eValue, $aScore, $DBid], [$eValue, $aScore, $DBid]]
		# Sort by e-value
		my @sorted = sort {$a->[0] <=> $b->[0]} @{$seqHash{$qID}};

		# Return an array of arrays with e-values smaller than or equal to the cutoff.
		@sorted = searchStarts(\$eCut, \@sorted);

		if (@sorted) {

			#Apply alignment score
			@sorted = sort {$b->[1] <=> $a->[1]} @sorted;

			# Get the maximum $aScore
			my $max = $sorted[0]->[1];

			# Get the index of the last element above the alignment score threshold
			my $lastInd = $#sorted;

			for (my $index = 1 ; $index <= $#sorted ; $index++) {

				if ($sorted[$index]->[1] < $max * $alignThresh) {
					$lastInd = $index - 1;
					last;
				}
			}

			# Replace the hash of array of arrays with a hash of hashes of arrays.
			# The filtered DBids are the keys of the inner hash. The array references contain the eValue and aScore.

			$seqHash{$qID} = undef;

			foreach my $arrayRef (@sorted[0 .. $lastInd]) {

				my $DBid    = $arrayRef->[2];
				my $eValue  = $arrayRef->[0];
				my $aScore  = $arrayRef->[1];
				my $counter = 0;

				if (exists $seqHash{$qID}{$DBid}) {
					$seqHash{$qID}{$DBid}->[0] += $eValue;
					$seqHash{$qID}{$DBid}->[1] += $aScore;

					# DBid counter
					$seqHash{$qID}{$DBid}->[2] += 1;
				}
				else {
					$counter = 1;
					$seqHash{$qID}{$DBid} = [$eValue, $aScore, $counter];
				}
			}

			# Calculate mean e-value and alignment score per query ID
			my @DBids = keys(%{$seqHash{$qID}});

			foreach my $DBid (@DBids) {
				my $eValue  = $seqHash{$qID}{$DBid}->[0];
				my $aScore  = $seqHash{$qID}{$DBid}->[1];
				my $counter = $seqHash{$qID}{$DBid}->[2];
				$seqHash{$qID}{$DBid} = [$eValue / $counter, $aScore / $counter];
			}

			# Should not happen, because at least the highest scoring element per query should be processed.
			die "Error: Alignment threshold lost query element" if (not keys(%{$seqHash{$qID}}));
		}
		else {
			delete $seqHash{$qID};
			$lostEval++;
		}
	}

	return (\%seqHash, $lostEval);
}

1;
