#!/usr/bin/env sh


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
# This script is run on the merged taxonomy of Metaxa2 (MTX.txt)
# to create a taxonomy file for MetaG.
#
# It expects MTX.txt to be in the script directory. The output will
# be tax.MTX.txt.
#
# This should be followed by makeMTXfasta.pl to delete records from
# the fasta file which have are no longer present in the tax.MTX.txt
# file.
#-------------------------------------------------------------------#

DIR=$(dirname $0)

# Remove entries with ;; at the end, but species was often given before
grep -v ';;$' "$DIR/MTX.txt" |\

# Remove trailing ;
perl -ne '$_ =~ s/;$//; print $_' |\

# MTX entries can have 7, 8 or 9 ranks. Force taxonomy to 7 ranks. Similar to SILVA-QIIME-version, I take the first 3 and last 4 ranks. Besides, I add "0" for subclass, -order.
perl -ne '@splits=split(";", $_); $len=@splits; if ($len != 7) {@res = (); push(@res, @splits[0..2]); push(@res, @splits[-4..-1]); @splits = @res}; splice(@splits,3,0,0); splice(@splits,5,0,0); $
print = join(";", @splits); print $print' |\

# Set unknown taxa to 0 not ""
perl -ne 'chomp($_); @splits=split(";", $_);$print = ""; foreach $split (@splits) {$split = 0 if ($split eq ""); $print.=$split.";"} $print =~ s/;$//; print $print."\n" '|\

# Especially mitochondria and chloroplast are missing the genus. Use first word of species as genus.
perl -ne '$_ =~ s/0;0;([a-zA-Z]+)/0;$1;$1/; print $_' |\

# Split species and strain
perl -ne 'chomp($_); @splits=split(";", $_); $spec = $splits[-1]; @specs = split(" ", $spec); $spec = join(" ",@specs[0..1]); $strain = join(" ", @specs[2..$#specs]); $spec.=";".$strain; splice(@
splits,-1,1,$spec); $print=join(";", @splits); print $print."\n"' |\

# Add 0 for unknown strains
perl -ne '$_ =~ s/;$/;0/; print $_' |\

# Replace species entries with genus sp. by 0
perl -ne '$_ =~ s/([A-Za-z]+)\s+(sp.);/0;/; print $_' |\

perl -ne '$_=~ s/\t/;/; print $_' > "$DIR/tax.MTX.txt"