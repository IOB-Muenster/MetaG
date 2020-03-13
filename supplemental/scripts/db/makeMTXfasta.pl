#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;


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


#===================================================================#
# Remove entries from the MTX fasta file, that do not have
# a taxonomic description in the taxonomy file, anymore.
# 
# Taxonomy file tax.MTX.txt and fasta file MTX.fa must be in the
# script directory. The output is out.MTX.fa. You may adapt these
# names in the "MODIFY HERE" section.
#-------------------------------------------------------------------#

my $cwd = getcwd;

#========================= MODIFY HERE =============================#
my $tax = $cwd."/"."tax.MTX.txt";
my $inFasta = $cwd."/"."MTX.fa";
my $outFasta = $cwd."/"."out.MTX.fa";
#-------------------------------------------------------------------#


# Read tax file.
my %taxonomy = ();

open(TAX, "<", $tax) or die "Can't open taxonomy file $tax, $!";
while(<TAX>) {
	
	chomp($_);
	my @splits = split(";", $_);
	my $id = $splits[0];
	
	if (exists $taxonomy{$id}) {
		die "#$id# Already in taxonomy, but should be unique."
	}
	else {
		$taxonomy{$id} = undef;
	}
}
close (TAX);

# Read db.fa, compare with %taxonomy and write matches directly.

# Switch to print sequence, if header matches
my $print=0;

open(OUT, ">",$outFasta) or die "Cant open db outfile $outFasta, $!";
open(FASTA, "<", $inFasta) or die "Can't open input db $inFasta , $!";

while (<FASTA>) {
	
	if ($_ =~ m/^>/) {
		my $id = (split(' ', $_))[0];
		$id =~ s/^>//;
		
		print OUT "\n" if ($print == 1);
		
		if (exists $taxonomy{$id}) {
			print OUT ">".$id."\n";
			$print = 1;
		}
		else {
			$print = 0;
		}
	}
	else {
		if ($print == 1) {
			chomp($_);
			print OUT $_;
		}
	}
}
close(FASTA);
close(OUT);