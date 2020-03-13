#!/usr/bin/env perl

use strict;
use warnings;

#==========================================================================================#
# Read and process config file giving parameters for analysis and environment variables
#------------------------------------------------------------------------------------------#

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


my $config = $ARGV[0];

my @lineSplits = ();	
my $command = "";
my $lparam = "";
my $matrix = "";
my $env = "";


open(CONF, "<", $config) or die "Could not open config file $config, $!";
while (<CONF>) {
	
	chomp($_);
	@lineSplits = split (" ", $_);
	next if (not @lineSplits);
	
	if ($lineSplits[0] =~ m/^#/) {
		
		if ($lineSplits[0] =~ m/^#env$/) {
			$env .= $lineSplits[1]."=".$lineSplits[2]." "
		}
		elsif ($lineSplits[0] =~ m/^#last$/) {
			$lparam .= join("_", @lineSplits[1..$#lineSplits])."_"
		}
		elsif ($lineSplits[0] =~ m/^#lastsplit$/) {
			$command .= "-lsplit"." ".$lineSplits[2]." ";
		}
		else {
			if ($lineSplits[2] ne "+") {
				$command .= $lineSplits[1]." ".$lineSplits[2]." ";
			}
			else {
				$command .= $lineSplits[1]." "
			}
		}
	}
	else {
		$matrix .= $_."_";
	}
}
close(CONF);


if ($lparam) {
	$lparam =~ s/_$//;
	$command .= "-lparam ".$lparam
}

if ($env) {
	$env =~ s/ $//;
}

chomp($command);

print $command."#".$env."#".$matrix