#!/usr/bin/env perl

use strict;
use warnings;

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


my $helpMSG = <<'END_MSG';

#======================================================#
# Translates fastq to fasta files.
#------------------------------------------------------#


# Expected file structure:
    @id1
    seq1
    +
    quality1
    @id2
    [...]


# Usage:
    fastq2a.pl fastq.fq
    
    fastq2a.pl -h or --help shows this message.

END_MSG


my $fastqF = $ARGV[0];


# Help?
die $helpMSG if (not defined $fastqF or $fastqF eq "-h" or $fastqF eq "--help" or $fastqF eq "-help");


# Check fastq format

my $isFQ = 1;
my $lineC = `wc -l < $fastqF`;

die "\nERROR: Unexpected fastq format.\nLine count is no multiple of 4.\n\n" if ($lineC % 4 != 0);

open(FASTQ, "<", $fastqF) or die "Can't open fastq file";
while (<FASTQ>) {
        $isFQ = 0 if ($_ !~ m/^@/ and $. == 1);
        $isFQ = 0 if ($_ !~ m/^[a,t,g,c,n,A,T,G,C,N]/ and $. == 2);
        $isFQ = 0 if ($_ !~ m/\+/ and $. == 3);
        
        last if ($. > 4 or $isFQ == 0);
}
close (FASTQ);

die "\nERROR: Unexpected fastq format.\n" if ($isFQ == 0);



open(FASTQ, "<", $fastqF) or die "No such fastq file, $!";
while(<FASTQ>) {
        next if ($. % 4 != 1 and $. % 4 != 2);
        if ($_ =~ m/^@/) {
                $_ =~ s/^@/>/;
                $_ = (split(" ", $_))[0];
                $_ .= "\n";
        }
        print $_;
}