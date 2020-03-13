package modules::ReadQuery;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(readQuery);
our @EXPORT_OK = qw(readQuery);


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


#==========================================================================#
# Read ids from 16S input sequences
#--------------------------------------------------------------------------#
sub readQuery {
	
	my $queryPath = $_[0];
	
	open(QUERY, "<", $queryPath) or die "Couldn't open query file $queryPath, $!"; 
	
	# Save file to array
	my @queries = <QUERY>;
	close(QUERY);
	
	my %queryH = ();
	my $qID = "";
	chomp(@queries);
	
	foreach my $query (@queries) {
		
		# Take fasta headers
		if ($query =~ m/^>/) {
			
			
			# Remove ">" and take only first part
			$query =~ s/>//;
			$query = (split(' ', $query))[0];
			$qID = $query;
			
			# Save each header only once
			$queryH{$qID} ="";
		}
		else {
			chomp ($query);
			$queryH{$qID}.=$query;
		}
	}
	@queries = keys(%queryH);
	my $allQ = @queries;
	
	return (\@queries, $allQ, \%queryH);
}

1;