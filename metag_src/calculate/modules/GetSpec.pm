package modules::GetSpec;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(getSpec visOut getStat getUnmatchedSeq);


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


#===========================================#
# Internal subroutines
#-------------------------------------------#

	#=====================================================================================#
	# Calculate geometric, harmonic, arithmetic and Williams' mean of confidence.
	#-------------------------------------------------------------------------------------#
	sub getMean {
		my @confidences = @{$_[0]};
		my $type = $_[1];
		my $mean = undef;
		if ($type eq "geometric") {
			
			$mean = 1.0;
			my $width = 500;
			my $i = 0;
			
			# Calculate the geometric mean
			# (c1*c2*c3....*cn)**(1/n)
			# in an efficient way even for large arrays:
			# (c1*c2*cw)**(1/n)
			# Multiply the above result for all steps of
			# width w in n.
			# This may cause some rounding issues.
			
			while (1) {
			    my $tmp = 1.0;
			    my $count = 0;
			    
			    while ($i + $count < @confidences && $count < $width) {
			    	my $conf = $confidences[$i + $count];
			    	
			    	if ($conf == 0) {
						$mean = undef;
						last;
					}
			        $tmp *= $conf;
			        $count++;
			    }
			    last if ($count == 0);
			    $mean *= $tmp ** (1.0 / @confidences);
			    $i += $count;
			}
		}
		elsif ($type eq "williams") {
			$mean = 1.0;
			my $width = 500;
			my $i = 0;
			
			# Calculate the William's mean
			# ((c1+1)*(c2+1)*(c3+1)....*(cn+1))**(1/n) -1
			# in an efficient way even for large arrays:
			# ((c1+1)*(c2+1)*(cw+1))**(1/n)
			# Multiply the above result for all steps of
			# width w in n and subtract 1
			# This may cause some rounding issues.
			
			while (1) {
			    my $tmp = 1.0;
			    my $count = 0;
			    
			    while ($i + $count < @confidences && $count < $width) {
			        $tmp *= ($confidences[$i + $count] + 1);
			        $count++;
			    }
			    
			    last if ($count == 0);
			    $mean *= $tmp ** (1.0 / @confidences);
			    $i += $count;
			}
			$mean--;
		}
		elsif ($type eq "arithmetic") {
			foreach my $conf (@confidences) {
				$mean += $conf;
			}
			
			$mean = $mean / @confidences
		}
		elsif ($type eq "harmonic") {
			foreach my $conf (@confidences) {
				
				# Should be 0, if calculated for values of 0.
				if ($conf == 0) {
					$mean = 0;
					last;
				}
				
				if (defined $mean) {
					$mean = $mean + (1 / $conf);
				}
				else {
					$mean = 1 / $conf;
				} 
			}
			
			if ($mean != 0) {
				$mean = @confidences / $mean;
			}
		}
		return $mean;
	}
	
	
	#========================================================================================================================#
	# Sort taxa by e-value, alignment score and abundance. Devide sorted in 4 categories: ++, +, -, --. Apply differential relative
	# weight for confidence to these groups. Taxa chosen by the assignment algorithm increase the score; taxa left behind 
	# decrease it. The score is used to weigh the confidence of the assignment algorithm.
	#------------------------------------------------------------------------------------------------------------------------#
	sub getConf {
		#$taxaStat{$taxon} = [$eValue, $aScore, $taxCount]
		my %taxaStat = %{$_[0]};
		
		#One taxon means highest confidence, exit sub
		my $taxNum = keys(%taxaStat);
		return 1 if ($taxNum == 1);
		
		my $maxTaxon = $_[1];
		my $totalScoring = 0;
		my $maxScore = 4;
		
		# Check, if the difference between taxon with most and second most count is too high,
		# i.e., second most count is equal to or less than 50% of most count.
		# Exit directly.
		my @sorted = sort {$taxaStat{$b}->[2] <=> $taxaStat{$a}->[2]} keys(%taxaStat);
		return 1 if  ($taxaStat{$sorted[0]}->[2]*0.5 >=  $taxaStat{$sorted[1]}->[2]);
		
		my $sort = sub {
			$taxaStat{$a}->[0] <=> $taxaStat{$b}->[0]  # by e-value
			or
		    $taxaStat{$b}->[1] <=> $taxaStat{$a}->[1]  # if equal: By alignment score
		    or 
		    $taxaStat{$b}->[2] <=> $taxaStat{$a}->[2]  # if equal: By taxCount
		};
		
		@sorted = sort $sort keys(%taxaStat);
		
		# Difference between first and last taxon in sorted
		
		# Do e-values differ?
		my $index = 0; 
		
		my $diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
		
		# E-value identical, but what about alignment score?
		if ($diff == 0) {
			$index = 1;
			$diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
		}
		# Alignment score identical, but what about taxCount?
		if ($diff == 0) {
			$index = 2;
			$diff = $taxaStat{$sorted[0]}->[$index] - $taxaStat{$sorted[-1]}->[$index];
		}
		
		# All the same: Should have been filtered out, previously by maxTaxon = "-"
		die "Identical inputs in confidence calculation" if ($diff == 0);
		
		my $bestVal = $taxaStat{$sorted[0]}->[$index];
		my $conf = 0;
		
		# Make bins and distribute weights
		foreach my $taxon (@sorted) {
			my $value = $taxaStat{$taxon}->[$index];
			my $weight = undef;
			my $firstThird = $bestVal - ($diff / 3);
			my $secThird = $bestVal - (2 * $diff / 3);
				
			# For e-value
			if ($index == 0) {
				
				# ++
				if ($value == $bestVal) {
					$weight = $maxScore 
				}
				# +
				elsif ($value > $bestVal and $value <= $firstThird) {
					$weight = 2
				}
				# -
				elsif ($value > $firstThird and $value <= $secThird) {
					$weight = 1
				}
				# --
				else {
					$weight = 0
				}
				
			}
			# For others
			else {
				
				# ++
				if($value == $bestVal) {
					$weight = $maxScore
				}
				# +
				elsif($value < $bestVal and $value >= $firstThird) {
					$weight = 2
				}
				# -
				elsif($value < $firstThird and $value >= $secThird) {
					$weight = 1
				}
				# --
				else{
					$weight = 0
				}
			}
			
			# Count number of taxa contributing to weight
			$totalScoring += $taxaStat{$taxon}->[2] if ($weight > 0);
			
			if ($taxon eq $maxTaxon) {
				$conf += $weight * $taxaStat{$taxon}->[2]
			}
			else {
				$conf -= $weight * $taxaStat{$taxon}->[2]
			}
		}
		# Calculate realtive score
		$conf = $conf / $totalScoring;
		
		# Rescale to a value between 0 and 1
		# Conf can ranke from -maxScore to + maxScore
		$conf = ($conf + $maxScore) / ($maxScore*2);
		
		return $conf;
	}
	

#===========================================#
# External subroutines
#-------------------------------------------#

	#===========================================#
	# Identify the species for every query read.
	#-------------------------------------------#
	
	sub getSpec {
		
		my ($ranksR, $outPrefix, $queriesR, $seqHashR, $dbHashR, $confCut, $isVirus) = @_;
		my @ranks = @$ranksR;
		my @queries = @$queriesR;
		my %seqHash = %$seqHashR;
		my %dbHash = %$dbHashR;
		
		my %pathos = (); 
		my %results =();
		my %unmatched = ();
		my %unmatchedIDs = ();
		my $matRank = 0;

		my $lca = "";
		my $visual = "1\t";
		
		# Each read and its assigned taxa
		my $outLinF = $outPrefix."LIN.txt";
		
		# Special format for krona
		my $visOutF = $outPrefix."VIS.txt";
		
		open(LIN, ">", $outLinF) or die "Could not open $outLinF";
		open (VIS, ">", $visOutF) or die "Could not open $visOutF";
		
		#my @ranks = ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species", "strain");
		
		# Check the ranks in descending order: From broad to specific
		foreach my $query (@queries) { #seqHashKeys
			
			# Remaining ranks for majority analysis
			if (exists $seqHash{$query}) {
				
				my $matchCount = keys(%{$seqHash{$query}});
				
				if ($matchCount == 0) {
					print LIN "No matches for $query.\n";
					next;
				}
				
				
				#====================================================================================#
				# Majority decision for ranks = lca, if confCut = 1.0.
				#-------------------------------------------------------------------------------------#
				
				my @DBids = keys(%{$seqHash{$query}});
				my $confidence = undef;
				RANK: for (my $i=0; $i <= $#ranks; $i++) {
					
					my $rank = $ranks[$i];
					my %ranks=();	
				
					# Get each LAST match dbID
					ID: foreach my $dbID (@DBids) {
						
						# Check information in db
						if (exists ($dbHash{$dbID})) { #dbHash: $dbID => {$rank => $taxon}
									
								# Is the rank noted for the id?
								if (exists $dbHash{$dbID}{$rank}) {
									
									my $taxon = $dbHash{$dbID}{$rank};
									
									if (exists $ranks{$taxon}) { 
										
										# Count numbers of taxon-supporting matches from the same DBid.
										# However, similar DBids merged to average before. Count should always be 1.
										if (exists $ranks{$taxon}{$dbID}) {
											$ranks{$taxon}{$dbID} += 1;
										}
										else {
											$ranks{$taxon}{$dbID} = 1
										}
									}
									else {
										$ranks{$taxon} = {$dbID => 1}
									}
								}
						}
					}
					my @taxa = keys(%ranks);
					my $maxCount = undef;
					my $totalCount = 0;
					my $maxTaxon = undef;
					my %taxaStat = ();
					foreach my $taxon (@taxa) {
						
						my @IDs = keys(%{$ranks{$taxon}});
						my $taxCount = 0;
						my $eValue = 0;
						my $aScore = 0;
						my $lenIDs = @IDs;
						
						foreach my $seqID (@IDs) {
							
							$eValue += $seqHash{$query}{$seqID}->[0];
							$aScore += $seqHash{$query}{$seqID}->[1];
							
							# Total count for one rank
							$totalCount += $ranks{$taxon}{$seqID};
							
							# Number of matches supporting the taxon over all matching DBids
							$taxCount += $ranks{$taxon}{$seqID};
							
							
						}
						# Mean eValue, aScore per taxon
						$eValue = $eValue / $lenIDs;
						$aScore = $aScore / $lenIDs;
						
						if ($taxon eq "0") {
							$taxon = "unclassified"
						}
						
						$taxaStat{$taxon} = [$eValue, $aScore, $taxCount];
							
						if (defined $maxCount) {
							if ($taxCount > $maxCount) {
								$maxCount = $taxCount;
								$maxTaxon = $taxon;
								
								# Only dbIDs from maxTaxon should be considered for lower taxa search
								@DBids = @IDs;
								
							}
							elsif ($taxCount == $maxCount) {
								
								# Matches have an equal count. How about the alignment quality?
								# E-value
								if ($taxaStat{$taxon}->[0] < $taxaStat{$maxTaxon}->[0]) {
									$maxTaxon = $taxon;
									@DBids = @IDs;
								}
								elsif ($taxaStat{$taxon}->[0] == $taxaStat{$maxTaxon}->[0]){
									
									# Alignment score
									if ($taxaStat{$taxon}->[1] > $taxaStat{$maxTaxon}->[1]) {
										$maxTaxon = $taxon;
										@DBids = @IDs;
									}
									# Identical in count and alignment quality: Set maxTaxon to special character "-".
									# If local maxTaxon "-" becomes global. Break, as I cannot decide between to equally likely taxa.
									elsif ($taxaStat{$taxon}->[1] == $taxaStat{$maxTaxon}->[1]) {
										$maxTaxon = "-";
										
										# Write "-" to %taxaStat to avoid bugs with following $taxCount == $maxCount
										$taxaStat{"-"} = [$eValue, $aScore, $taxCount];
										
										# Will not be considered in the end. Just if algorithm will be changed,
										push (@DBids, @IDs);
									}
								}
							}
						}
						else {
							$maxCount = $taxCount;
							$maxTaxon = $taxon;
							@DBids = @IDs;
						}
					}
					
					# Evaluate global maxTaxon for a rank. If it is "-", break and 
					# set this and following ranks to "UNMATCHED"
					if ($maxTaxon eq "-") {
						my $rankInfo = "UNMATCHED";
						foreach my $unmRank (@ranks[$i..$#ranks]) {
							if (exists $results{$unmRank}) {
								if (exists $results{$unmRank}{$rankInfo}) {
									$results{$unmRank}{$rankInfo} -> [0] += 1;
								}
								else {
									$results{$unmRank}{$rankInfo} = [1];
								}
							}
							else {
								$results{$unmRank} = {$rankInfo => [1]};
							}
							
							# Write down unmatched read IDs at each rank
							if (exists $unmatched{$unmRank}) {
								push(@{$unmatched{$unmRank}}, $query);
							}
							else {
								$unmatched{$unmRank} = [$query];
							}							
						}
						last RANK;
					}
					# If we have a meaningful maxTaxon, try to remove local "-" from taxaStat
					else {
						delete $taxaStat{"-"} if (exists $taxaStat{"-"});
					}
					
					
					#================================================================#
					# Calculate confidences of correct assignment
					#----------------------------------------------------------------#
					
					# Probability of correct assignment
					if (defined $confidence) {
						my $lastConf = $confidence;
						
						# Weighted confidence
						$confidence = getConf(\%taxaStat, $maxTaxon);
						$confidence = $confidence * $lastConf;
					}
					else {
						
						# Weighted confidence
						$confidence = getConf(\%taxaStat, $maxTaxon);
					}
					if ($confidence < $confCut) {
						
						my $rankInfo = "UNMATCHED";
						foreach my $unmRank (@ranks[$i..$#ranks]) {
							if (exists $results{$unmRank}) {
								if (exists $results{$unmRank}{$rankInfo}) {
									$results{$unmRank}{$rankInfo} -> [0] += 1;
								}
								else {
									$results{$unmRank}{$rankInfo} = [1];
								}
							}
							else {
								$results{$unmRank} = {$rankInfo => [1]};
							}
							
							# Write down unmatched read IDs at each rank
							if (exists $unmatched{$unmRank}) {
								push(@{$unmatched{$unmRank}}, $query);
							}
							else {
								$unmatched{$unmRank} = [$query];
							}	
						}
						last RANK;
					}
					
					if ($isVirus == 0) {
						# For pathogen detection
						if ($rank eq "species") {
							if (exists $pathos{$maxTaxon}) {
								$pathos{$maxTaxon} += 1;
							}
							else {
								$pathos{$maxTaxon} = 1;
							}
							
							 
						}
					}
					elsif ($isVirus == 1) {
						# For host detection in viruses (id => [name, count])
						if ($rank eq $ranks[-1] and $isVirus == 1) {
							
							# Multiple IDs could technically yield the same species and strain with potentially different
							# host entries. So I want to save all....
							my $count = 1 / @DBids;
							
							# ... however, the species and strain name must be the same. So, I arbitrarily look up the first ID in dbHash.
							# It is more convenient to display strains that are not 0 and thus not "unclassified"
							my $strain = $dbHash{$DBids[0]}{$ranks[-1]};
							$strain = "" if ($strain =~ m/^0$/);
							my $name = $dbHash{$DBids[0]}{$ranks[-2]}." ".$strain;
							
							
							
							foreach my $id (@DBids) {
								
								if (exists $pathos{$id}) {
									$pathos{$id}->[1] += $count;
								}
								else {
									$pathos{$id} = [$name, $count]
								}
							}
						}
					}					
					
					#==============================================================#
					# Save resulting annotation
					#--------------------------------------------------------------#
					
					######### EXPERIMENTAL: Remove () in taxa names, potentially critical in strains
					# not for PATRIC: only spec queried and those don't contain () #########
					
					$maxTaxon =~ s/\(//g;
					$maxTaxon =~ s/\)//g;
					
					$lca = $lca."$rank: $maxTaxon: $maxCount ($totalCount)\n";
					$visual = $visual."$maxTaxon\t";
					
					# Save output statistics for all queries
						
						# Save taxa for ranks. For each taxon save count (array element 0) and all confidences.
						if (exists $results{$rank}) {
							if (exists $results{$rank}{$maxTaxon}) {
								$results{$rank}{$maxTaxon} -> [0] += 1;
								
								# If confCut is 1, average confidence will always be 1.
								if ($confCut != 1) {
									push(@{$results{$rank}{$maxTaxon}}, $confidence);
								}
							}
							else {
								$results{$rank}{$maxTaxon} = [1, $confidence];
							}
						}
						else {
							$results{$rank} = {$maxTaxon => [1, $confidence]};
						}
				}
				
				#---------------------------------------------------------------------------#
				$lca = substr($lca, 0, -1);
				$visual = substr($visual, 0, -1);
				
				# Escape " with &quot; for krona.
				$visual =~ s/"/&quot;/g;
				
				# Print results for each read and output for krona
				print LIN ">$query\n$lca\n";
				print VIS "$visual";
				$lca = "";
				$visual = "\n1\t";
			}
			else {
				print LIN "No match for $query and chosen cutoff.\n";
				
				# Write down unmatched read IDs
                if (exists $unmatched{"removed"}) {
                	push(@{$unmatched{"removed"}}, $query);
                }
                else {
                    $unmatched{"removed"} = [$query];
                }
			}
		}
		close(LIN);
		close(VIS);
		
		return (\%results, \%pathos, \%unmatched, \%unmatchedIDs);
	}
	
	
	#======================================#
	# Use krona to visualize species calls
	#--------------------------------------#
	
	sub visOut {
		
		my $krona = $_[0];
		my $visOutP = $_[1];
		
		my $visFileN = (split("/", $visOutP))[-1];
		my $krFileN = $visFileN =~ s/txt$/html/r;
		my $krFileP = $visOutP =~ s/$visFileN$/$krFileN/r;
		qx($krona -o $krFileP $visOutP) or die "Krona failed with $?";
	}
	
	
	#===================================#
	# Print statistic for species calls
	#-----------------------------------#
	sub getStat {
		
		my ($ranksR, $outPrefix, $allQ, $resultsR, $lostEval, $type, $confCut) = @_;
		my @ranks = @$ranksR;
		my %results = %$resultsR;
		
		my $broadRank = $ranks[0];
		
		my $statFile =$outPrefix."LOG.txt";
		open(STAT, ">", $statFile) or die "Could not open $statFile";
		
		# Print overall taxa analysis
		print STAT "################################################\nOverall statistics\n------------------------------------------------\n";
		print STAT "\nTotal read count: $allQ\n";
		
		# Total matches = matches at most unspecific rank, i.e. domain.
		my $totalM = 0;
		if (exists $results{$broadRank}) {
			my @taxa = keys(%{$results{$broadRank}});
			foreach my $taxon (@taxa) {
				$totalM += $results{$broadRank}{$taxon} -> [0];
			}
		}
		
		my $lastRem = $allQ - $totalM - $lostEval;
		my $totalRem = $lostEval + $lastRem;
		print STAT "Matched: $totalM\n";
		print STAT "Removed: $totalRem\n\tDue to LASTal settings: $lastRem\n\tDue to e-value cutoff: $lostEval\n";
		
		#foreach my $rank ("domain", "phylum", "class", "subclass", "order", "suborder", "family", "genus", "species", "strain") {
		foreach my $rank (@ranks) {
			my $sum = 0;
			print STAT "\n$rank";
			my @prints = ();
			if (exists $results{$rank}) {
				my @taxa = keys(%{$results{$rank}});
				foreach my $taxon (@taxa) {
					my $count = $results{$rank}{$taxon} -> [0];
					$sum += $count;
					my @confidences = @{$results{$rank}{$taxon}}[1 ... $#{$results{$rank}{$taxon}}];
					
					# Unmatched does not have confidences
					if (@confidences) {
						my $mean = undef;
						
						if ($confCut == 1) {
							$mean = 1;
						}
						else {
							# Calculate mean of confidences
							$mean = getMean(\@confidences, $type);
							$mean = sprintf("%.2f", $mean);
						}
										
						# Handle errors: Geometric for values containing 0
						if (not defined $mean) {
							print "Error while calculating average confidence. Did you try to calculate geometric mean with 0 in input?\n";
							$mean = "Error";
						}
						
						# Mean and abbreviation of selected mean type
						$mean = $mean."[".substr($type, 0, 1)."]";
						push(@prints, [$count, $taxon, $mean])
					}
					else {
						push(@prints, [$count, $taxon, "NA"])
					}
				}
				# Print high abundant taxa first
				my @sPrints = sort {$b -> [0] <=> $a -> [0] } @prints;
				foreach my $array (@sPrints) {
					my $count = $array -> [0];
					my $taxon = $array -> [1];
					my $mean = $array -> [2];
					if (defined $mean) {
						print STAT "\n\t$taxon(conf: $mean)\t$count";
					}
					else {
						print STAT "\n\t$taxon\t$count"
					}
				}
			}	
			else {
				print STAT "\n$rank not in dataset!";
			}
			
			print STAT "\n\tTotal: $sum ($allQ)\n";
		}
		close(STAT);
	}
	
	
	#===============================================#
	# Write unmatched reads to separate fasta files
	#-----------------------------------------------#
	sub getUnmatchedSeq {
		
		my %unmatched = %{$_[0]};
		my %queryH = %{$_[1]};
		my $outPrefix = $_[2];
		
		foreach my $rank (keys(%unmatched)) {
			
			my $outfile = $outPrefix."unmatched.".$rank;
			open (OUT, ">", $outfile) or die "Could not open $outfile.";
			
			foreach my $ID (@{$unmatched{$rank}}) {
				if (exists $queryH{$ID}){
					print OUT ">".$ID."\n".$queryH{$ID}."\n";
				}
			}
			
			close(OUT)
		}
		
	}


1;