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


#=========================================================================================================#
# Build a MetaG compatible database from RDP
#---------------------------------------------------------------------------------------------------------#
# 
# DESCRIPTION
# 
#	Using a RDP FASTA file, this script will create a taxonomy and a FASTA file in MetaG format. It aims
#	to improve the taxonomic resulotion of RDP, by setting uninformative taxon names to a common
#	expression (e.g.: uncultured * --> uncultured). Besides, species and strain names are separated.
#
# USAGE
# 
#	Download the RDP FASTA file(s) from https://rdp.cme.msu.edu/misc/resources.jsp. If you downloaded
#	multiple files, concatenate the extracted FASTA files into a single file. Replace all ";" in FASTA
#	headers with ",". Then run:	
#		
#	makeRDP.pl -input rdp.fa
#
#	OR
#
#	makeRDP.pl -i rdp.fa
#
# OUTPUT
# 
#	The following output files will be created in the directory of the input FASTA file:
#
#	tax.rdp.txt		Taxonomy file in MetaG format
#	rdp.fa			Fasta file in MetaG format
#
#
# LIMITATIONS
#
#	Expects the database to contain only the following ranks: "domain", "phylum", "class", "subclass",
#	"order", "suborder", "family", "genus", "species" and "strain".
#	Due to the size of RDP we expect that there are special cases whose correction is beyond the scope
#	of this script. These should be fixed manually.
#
#=========================================================================================================#


use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#-----------------------------------------------------------------------------------------------#
# Parse and check arguments
#-----------------------------------------------------------------------------------------------#
my $argC = @ARGV;

my $inF = "";
my $help = 0;

my $usage = <<'EOF';
#=========================================================================================================#
# Build a MetaG compatible database from RDP
#---------------------------------------------------------------------------------------------------------#
 
 DESCRIPTION
 
	Using a RDP FASTA file, this script will create a taxonomy and a FASTA file in MetaG format. It aims
	to improve the taxonomic resulotion of RDP, by setting uninformative taxon names to a common
	expression (e.g.: uncultured * --> uncultured). Besides, species and strain names are separated.

 USAGE
 
	Download the RDP FASTA file(s) from https://rdp.cme.msu.edu/misc/resources.jsp. If you downloaded
	multiple files, concatenate the extracted FASTA files into a single file. Replace all ";" in FASTA
	headers with ",". Then run:	
		
	makeRDP.pl -input rdp.fa

	OR

	makeRDP.pl -i rdp.fa

 OUTPUT
 
	The following output files will be created in the directory of the input FASTA file:

	tax.rdp.txt		Taxonomy file in MetaG format
	rdp.fa			Fasta file in MetaG format


 LIMITATIONS

	Expects the database to contain only the following ranks: "domain", "phylum", "class", "subclass",
	"order", "suborder", "family", "genus", "species" and "strain".
	Due to the size of RDP we expect that there are special cases whose correction is beyond the scope
	of this script. These should be fixed manually.

#=========================================================================================================#
EOF
;

GetOptions ("input:s"   => \$inF,
           	'help|?' => \$help) or die $usage;
           	
if ($help > 0 or not $inF) {
	print $usage;
	exit 0;
}

my $outP = dirname($inF);
my $outTax = $outP . "/tax.rdp.txt";
my $outFa = $outP . "/rdp.fa";


#-----------------------------------------------------------------------------------------------#
# Parse RDP taxonomy and write output FASTA.
#-----------------------------------------------------------------------------------------------#
my %dbHash =();

print "INFO: Loading and analyzing. This will take a while...\n";

open(FA, ">", $outFa) or die "ERROR: Couldn't open output fasta file $outFa, $!";
open(DB, "<", $inF) or die "ERROR: Couldn't open db file $inF, $!";
while(<DB>) {
	
	# Get fasta headers...
	if ($_ =~m"^\>") {
		# Write to new fasta file
		print FA $_;
		
		my @line_data = split ' ', $_;
		my $encode = $line_data[0];
		
		# Remove greater sign
		$encode = substr ($encode,1);
		
		# Save identifier and systematic nomenclature to hash
		my %lineage = ();
		
		foreach my $a (0..$#line_data) {
			
			# Get information about ancestors
			if (index($line_data[$a], "Lineage")!=-1) {
				
				# Sometimes lineage information is devided by blanc. Reunite the parts.
				# Lineage should be the last part in each row.
				my $lineage = join(' ',@line_data[$a..$#line_data]);     
				
				# Devide it in single ranks
				my @lineages = split ',', $lineage;
				
				# Loose noninformative root
				@lineages = @lineages[2..$#lineages];
				
				# Save rank and rank details
				for (my $index = 0; $index < @lineages; $index += 2) {
					
					# Loose non-informative lineage data and skip rank
					if ($lineages[$index] =~ m"incertae_sedis") {
						$index+=2
					}
					elsif ($lineages[$index] =~ m"^unclassified_") {
						$index+=2
					} 
					else {
						# Some entries finish with a taxon that is not assigned to
						# any rank. --> Skip.
						last if (not $lineages[$index+1]);
						
						if (exists $dbHash{$encode}) {
							$dbHash{$encode}{$lineages[$index+1]} = $lineages[$index]
						}
						else {
							$dbHash{$encode} = {$lineages[$index+1] => $lineages[$index]};
						}
					
					}
				}
				
				# Get species and strain information
				my $species = "";
				my $strain = "";
				my $uncCount = 0;

				# The information is located between the ID and lineage information
				for (my $i = 1; $i < $a; $i++) {
					
					# Replace all uncultured species comments with uncultured
					if ($line_data[$i] eq "uncultured") {
						$species = "uncultured";
					}
					
					# Remove "," from species data and concatenate the data
					else {
						my $info = $line_data[$i];
						
						if ($line_data[$i] =~ m"[,]$") {
							$info = substr($line_data[$i], 0, -1);
						}

						# Special characters, numbers, capital letters inside a word indicate a strain description
						if ($info =~ m"[^A-Za-z,]" or $info =~ m".[A-Z]"){
							
							# Information about the strain will be after sp.
							if ($info eq "sp.") {
								my $s = undef;	
									
								if ($i+1 < $a-1) {
									$s = join(" ", @line_data[$i+1..$a-1]);
								}
								elsif ($i+1 < $a) {
									$s = join(" ", $line_data[$i+1])
								}
								
								$strain = $s;
								last;
							}
							else {
								my $s = undef;
									
								if ($i < $a-1) {
									$s = join(" ", @line_data[$i..$a-1]);
								}
								else {
									$s = join(" ", $line_data[$i])
								}
								
								$strain = $s;
								last;
							}


						}
						else {
							# The species was set to uncultured, now the species information needs to be skipped.
							# Just the strain information is needed.
							if ($species =~ m"^uncultured") {
									my $count = 0;
									my $iLoop = $i; # Only for while loop
									
									# Check elements between topical and lineage element
									while ($count < $a - $iLoop - 1) {

										# Remove genus, but not strains with big letters. Removes information like bacteria or organism. Even if preceeded by whitespace or ",".
										if ($line_data[$i] =~ m"^[\s,]?[A-Z][a-z]" or $line_data[$i] =~ m"^[\s,]?bact" or $line_data[$i] =~ m"^[\s,]?orga") {
											
											# It's a strain!
											last if ($line_data[$i] =~ m"[^A-Za-z,]");
											$i += 1
										}
										# Check one condition per count. Remove sp. or spec.
										elsif ($line_data[$i] =~ m"^sp.*[.]") {
											$i += 1
										}

										$count++
									}
									
									# Extract strain information
									my $s = undef;
											
									if ($i < $a-1) {
										$s = join(" ", @line_data[$i..$a-1]);
									}
									else {
										$s = join(" ", $line_data[$i])
									}
									
									$strain = $s;
									last;
							}
							else {
								# Save species information
								$species = $species." ".$info;
							}
							
						}
						
					}					
				
				}				
				# Post-process species entry
				# Remove front and end whitespaces
				$species =~ s/^\s+|\s+$//g;
				my $specLen= split(' ', $species);
				
				# genus is repeated as spec., we already know that --> unclassified bacterium
				if ($specLen == 1 and $species !~ m/^uncultured/) {
					$species = "unclassified"
				}
				# remove acidophilic... or bacterium...., but leave unclassified, uncultured etc.
				if ($species =~ m"^[a-tv-z]") {
					$strain = "" if (not $strain);
					$strain = $species." ".$strain;
					$species = undef;
				}
				$strain = "0" if (not $strain);
				
				# Save strain
				$dbHash{$encode}{"strain"} = $strain;
				
				# No entry for species
				if (not $species) {
					$species = "unclassified";
				}
				#print "$encode\tSpecies: $species\tStrain: $strain\n\n";
				$dbHash{$encode}{"species"} = $species;
				
				last;
			}
			
		}
		
	}
	else {
		# Print to new fasta file
		# Upper case to avoid issues with repeat masking
		print FA uc($_);
	}
	
}
close (DB);
close (FA);


#-----------------------------------------------------------------------------------------------#
# Write taxonomy
#-----------------------------------------------------------------------------------------------#
print "INFO: Writing taxonomy\n";

my @encList = keys(%dbHash);

# Write new db file
open(TAX, ">", $outTax) or die "ERROR: Couldn't open output taxonomy file $outTax, $!";

foreach my $enc (@encList) {
	my $printstr = "";
	$printstr = $printstr.$enc;
	foreach my $rank ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species", "strain") {
		if (exists $dbHash{$enc}{$rank}) {
			$printstr = $printstr.";".$dbHash{$enc}{$rank};
		}
		else {
			$printstr = $printstr.";"."0"
		}
	}
	
	$printstr = $printstr."\n";
	print TAX $printstr;
}
close (TAX);

print "\nDONE\n";