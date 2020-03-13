package modules::GetPatho;

use strict;
use warnings;
use Exporter;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT = qw(getPatho getHostVir);
our @EXPORT_OK = qw(getPatho getHostVir);


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


# Internal subroutine
sub readPatho {
	my $pdbPath = $_[0];
	my %pathogens = ();
	
	open (PATRIC, "<", $pdbPath) or die "Couldn't open PATRIC file, $!";
	while(<PATRIC>) {
		my $line = $_;
		
		next if ($line =~ m/^#/);
		chomp($line);
		
		my @lineSplits = split (';', $line);
		
		# Prepare less sensitive matching, because of possible nomenclature clashes
		my $spec = $lineSplits[-4];
		$spec = lc($spec);
		$spec =~ s/\W//g;
		$spec =~ s/_//g;
	
		my $strain = $lineSplits[-3];
		
		# Get only meaningful species, strain combination.
		next if ($spec =~ m/^un/ and $strain eq "0");
		
		my $host = $lineSplits[-2];
		$host = substr($host, 1);
		
		my $resist = $lineSplits[-1];
		$resist = substr($resist, 1);
		
		$resist = "No information" if ($resist eq "0");
		
		if (exists $pathogens{$spec}) {
			$pathogens{$spec}{$strain} =  [$host, $resist];
		}
		else {
			$pathogens{$spec} = {$strain => [$host, $resist]};
		}		
	}
	close (PATRIC);
	return \%pathogens;	
}

# External suroutines

sub getPatho {
	my ($candidatesR, $outP, $pdbPath) = @_;
	my %candidates = %$candidatesR; # could also contain (species => {strain => count}), now only (species => count)
	
	my $pathogensR = readPatho($pdbPath);
	my %pathogens = %$pathogensR;
	
	open (PATHO, ">", $outP) or die "Could not open patho out file, $_!";
	# I want to keep the proportions of called species, but I also want to report the hosts/resistances of all strains
	# for a species. This would make my printed pathogen calls: #calls from getMeta * #host/resistances.
	# To correct the bias, I multiply the abundance of each species-host-resistance phenotype by a correction factor.
	# The factor is #calls from get Meta / # total species-host-resistance phenotypes in local PATRIC.
	# NOTE: Although the proportions of species are not biased, the percentage of single species-host-resistance phenotypes
	# 		is influenced by the abundance of the phenotypes in PATRIC and is not a biological pattern. It is a probability.

	
	my $correctFact = 0;
	my $expecHarmlC = 0;
	my %pathoPrints = ();
	
	foreach my $species (keys(%candidates)) {
		
		# Adjusted species for better comparison between amplicon db and PATRIC, but save unadjusted amplicon DB name
		my $adjustedSpec = $species =~ s/\W//gr;
		$adjustedSpec =~ s/_//g;
		$adjustedSpec = lc($adjustedSpec);

		
		# Loose non-informative uncultured/ unclassified. Without further rank, those could be anything
		if ($species =~ m/^un/) {
			$expecHarmlC += $candidates{$species};
			next;
		}
		
		if (exists $pathogens{$adjustedSpec}) {
			
			foreach my $strain (keys(%{$pathogens{$adjustedSpec}})) { 	
					my $host = $pathogens{$adjustedSpec}{$strain}->[0];
					my $resist = $pathogens{$adjustedSpec}{$strain}->[1];
					my @resists = split('\+', $resist);
					for (my $i = 0; $i<= $candidates{$species}; $i++) {
						
						foreach my $res (@resists) {
							
							if (exists $pathoPrints{$species}) {
								if (exists $pathoPrints{$species}{$host}) {
									if (exists $pathoPrints{$species}{$host}{$res}) {
										$pathoPrints{$species}{$host}{$res} += 1;
									}
									else {
										$pathoPrints{$species}{$host}{$res} = 1;
									}
								}
								else {
									$pathoPrints{$species}{$host} = {$res => 1};
								}
							}
							else {
								$pathoPrints{$species} = {$host => {$res => 1}};
							}
							
							#Count number of prints
							if (exists $pathoPrints{$species}{"COUNT"}) {
								$pathoPrints{$species}{"COUNT"} += 1
							}
							else{
								$pathoPrints{$species}{"COUNT"} = 1
							}
						}
					}
			}
			

		}
		else {
			$expecHarmlC += $candidates{$species};
		}
	}
	foreach my $pathoSpec (keys(%pathoPrints)) {
		
		# Escape " to help Krona
		my $speciesPrint = $pathoSpec =~ s/"/&quot;/gr;
		
		# Number of species calls in getMeta / number of total prints
		my $correctFact = $candidates{$pathoSpec} / $pathoPrints{$pathoSpec}{"COUNT"};
		
		foreach my $host (keys %{$pathoPrints{$pathoSpec}}) {
			next if ($host eq "COUNT");
			foreach my $res (keys %{$pathoPrints{$pathoSpec}{$host}}) {
				my $count = $pathoPrints{$pathoSpec}{$host}{$res} * $correctFact;
				print PATHO "$count\t$host\t$res\t$speciesPrint\n" #\t$strain
			}
		}
	}
	
	my $printedHarml = $expecHarmlC;
	print PATHO "$printedHarml\tNonpathogenic\n";
	
	close(PATHO);
	return;
}


sub getHostVir {
	my ($virsR, $outP, $dbF) = @_;
	
	my %virs = %$virsR; # {id => [spec_strain, count]}
	my %hosts = ();
	
	# Read host database
	open(DB, "<", $dbF) or die "Could not open host db file";
	while (<DB>) {
		
		chomp($_);
		my @splits = split(/;/, $_);
		my $id = $splits[0];
		my $host = $splits[1];
		
		if (exists $hosts{$id}) {
			push(@{$hosts{$id}}, $host)
		}
		else {
			$hosts{$id} = [$host];
		}
	}
	close(DB);
	
	# Find hosts for detected species and print results
	open (PATHO, ">", $outP) or die "Could not open patho out file!";
	
	foreach my $id (keys(%virs)) {
		
		my $detectNumb = $virs{$id} -> [1];
		my $name = $virs{$id} -> [0];
		
		if (exists $hosts{$id}) {
			
			my $hostNumb = @{$hosts{$id}};			
			$detectNumb = $detectNumb / $hostNumb;
			
			foreach my $host (@{$hosts{$id}}) {
				
				print PATHO "$detectNumb\t$host\t$name\n"
			}
		}
		else {
			print PATHO "$detectNumb\tNonpathogenic\n"
		}
	}
	close(PATHO);
}

1;
