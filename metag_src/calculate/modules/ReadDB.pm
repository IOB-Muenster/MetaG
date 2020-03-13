package modules::ReadDB;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(readDB);
our @EXPORT_OK = qw(readDB);


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


#=============================================================================================#
# Load database taxonomy into memory
#---------------------------------------------------------------------------------------------#

sub readDB {
	
	my $dbPath = $_[0];
	
	my %dbHash = ();
	
	open(DB, "<", $dbPath) or die "Couldn't open DB file, $!";
	
	# Fall back rank structure
	my @ranks = ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species", "strain");
	
	while(<DB>) {
		
		chomp($_);
		
		# Allow for other than standard rank system
		if ($. == 1 and $_ =~ m/^#/) {
			$_ =~ s/^#//;
			@ranks = split(";", $_);
			next;
		}
		
		
		my @line_data = split ';', $_;
		my $encode = $line_data[0];
		my $counter = 1;
		
		foreach my $rank (@ranks) {
			
			my $rankInf = $line_data[$counter];
			
			# Remove newline
			if ($counter ==10) {
				chomp($rankInf);
			}
			$dbHash{$encode}{$rank} = $rankInf;	
			$counter++;
		}
		
	}
	
	close (DB);
	return (\%dbHash, \@ranks);
}



1;