package modules::Misc;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(searchStarts);


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
# Takes a sorted array of arrays and returns an array of arrays with inner
# arrays smaller than or equal to query.
#---------------------------------------------------------------------------#

sub searchStarts {
	
	# Input array has to be sorted
	my $listRef = $_[1];
	
	# Query
	my $query = ${$_[0]};
	
	# Variables
	my @result = (); #undef won't work
	my $lenList = @{$listRef};
	$lenList = sprintf("%.0f", $lenList);
	
	# Position of putative match
	my $posList = $lenList/2;
	$posList = sprintf("%.0f", $posList);
	
	# Variables for controlling the amount of searches
	my $isRep = 1;
	my $lastValue = undef;
	my $isNeighbor = 0;
	
	
	# Check once if the query could be found in the list (greater than or equal to first)
	if ($query >= $listRef -> [0][0]) {
	}
	else {
		$isRep =0;
		@result =();
	}
	# If the item has not been found yet...
	while ($isRep == 1) {
		
		# Check if the array index is still valid.
		if ($posList >= 0 and $posList <= $#$listRef) {
			
			# The query is greater or equal to entry at that position.
			if ($listRef -> [$posList][0] <= $query) {
							
				if ($posList+1 < $#$listRef) {
					
					# The next item should be bigger than the query.
					if ($listRef -> [$posList+1][0] > $query) {
						$isRep = 0;
						@result = @{$listRef}[0..$posList];
					}
					# Else: Increase the index by 50 %.
					else {
						$isNeighbor = 1;
						$lastValue = $posList;
						$posList = $posList + 0.5 * $posList;
						$posList = sprintf("%.0f", $posList);
						
					}
				}
				# Avoids infinite loops for query also bigger last element
				else {
					
					if ($listRef -> [$#$listRef][0] <= $query) {
						$isRep = 0;
						@result = @{$listRef};
					}
					else{
						$isRep = 0;
						@result = @{$listRef}[0..$#$listRef-1];
						
					}
				}
			}
			# The query is smaller than the item at the position. Decrease the index by 50 %.
			else {
				
					# A match was found, but the index was increased to find matches at further positions. 
					# Now, the element at the list index is too big.
					if ($isNeighbor == 1) {
						for my $a (($lastValue..$posList)) {
							if ($listRef -> [$a][0] > $query) {
								@result = @{$listRef}[0..$a-1];
								$isRep = 0;
								last;
							}
						}
						if (not @result) {
							$isRep = 0;
							@result =();
						}
					}
					$lastValue = $posList;
					$posList = $posList - 0.5 * $posList;
					$posList = sprintf("%.0f", $posList);			
			}
		}
		# The array index is invalid...
		else {
			
			# ... because it is smaller than 0.
			if ($posList < 0) {
				for my $a ((0..$lastValue)) {
					if ($listRef -> [$a][0] > $query) {
						if (not $a == 0) {
							@result = @{$listRef}[0..$a-1];
							$isRep = 0;
							last;
						}
					}
				}
			}
			# ... because it is bigger than the last index.
			else {
				if ($listRef -> [$#$listRef][0] <= $query) {
					@result = @{$listRef};
					$isRep = 0;
				}
				else {
					for my $a (($lastValue..$#$listRef)) {
						
						# The query is smaller, break and save previous
						if ($listRef -> [$a][0] > $query) {
							@result = @{$listRef}[0..$a-1];
							$isRep = 0;
							last;
						}
					}
				}
			}
			
			if (not @result) {
				$isRep = 0;
				@result =();
			}
		}
	}
	return @result
}
