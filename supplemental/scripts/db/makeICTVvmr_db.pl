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


#=====================================================================================#
# This script receives the uptodate taxonomy from ICTV. ICTV also provides genbank
# identifiers to reference sequences. I request the sequences from the NCBI eFetch
# interface. The taxonomy is then assigned to the sequence. The output is printed
# to three files: A taxonomy file, a file with host information and one with
# the sequences.
#
# tax.VMR.txt:
# #realm;subrealm;kingdom;subkingdom;phylum;subphylum;class;subclass;order;suborder;family;subfamily;genus;subgenus;species;isolate_seqName <- rank names
# genbankID1;tax
# genbankID2;tax
# ...
#
# patho.VMR.txt:
# genbankID1;host
# genbankID2;host
#
# VMR.fa:
# >genbankID1
# seq1
# >genbankID2
# seq2
#
# IMPORTANT: The VMR file from ICTV is an excel sheet and can be found here: 
#			 https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository
#			 You will have to save the sheet containing the taxonomy as VMR.csv.
#			 Please first replace empty cells with "NA". Remove line breaks
#			 in fields. After that set the field separator to tab and save
#			 the sheet.
#			 The csv file must be in the same directory as this script. The
#			 output will also be created in the directory of the script.
#
#			 GenBank IDs must only be used once. Otherwise, the duplicate entry/entries
#			 will be lost completely. You may want to report duplicate IDs to the
#			 ICTV team.
#
# KNOWN ISSUES: IDs which are not reported by the NCBI API are present in
#				the tax.VMR.txt and patho.VMR.txt, but not in the VMR.fa.
#				Additional taxa in these files do not influence downstream
#				analyses. You may want to check, if the GenBank IDs are valid
#				and contact NCBI or ICTV, depending on the situation.
#-------------------------------------------------------------------------------------#


use strict;
use warnings;
use LWP::UserAgent();

#=====================================================================================#
# Read input taxonomy.
#-------------------------------------------------------------------------------------#


# Input taxonomy from ICTV.
my $inTaxP = "VMR.csv";
my $outTaxP = "tax.VMR.txt";
my $outHostP = "patho.VMR.txt";
my $outSeqP = "VMR.fa";
my $rankNames = "realm;subrealm;kingdom;subkingdom;phylum;subphylum;class;subclass;order;suborder;family;subfamily;genus;subgenus;species;isolate_seqName";


my %ids = ();

print "INFO: Reading input taxonomy from $inTaxP\n";

open(OUTTAX, ">", $outTaxP);
print OUTTAX "#".$rankNames."\n";

open(OUTHOST, ">", $outHostP); 

open(INTAX, "<", $inTaxP);
while(<INTAX>) {
	#Skip header
	next if ($_ =~ m/GENBANK accession/);
	
	chomp($_);
	$_ =~s/"//g;
	
	my @splits = split("\t", $_);
		
	# Filter strange undef entries
	print "Error: @splits\nFound only ".@splits." entries. Expected 26-\n" if (@splits < 26);
	next if (@splits < 26);
	
	my $genbID = $splits[21];
	my $host = $splits[25];
	my @hosts = split(/,/, $host);
	
	# Use only entries having a genbank ID
	next if ($genbID eq "NA" or $genbID eq "");
	
	# Get and reformat taxonomy
	my $tax = join(";", @splits[2..16]);
	$tax .= ";".$splits[20];
	$tax =~ s/;NA/;0/g;
	$tax =~ s/^NA;/0;/g;
	
	# Sometimes one entry has multiple genbank IDs.
	# Create separate entries per ID.	
	my @genbIDs = split("; ", $genbID);
	
	foreach my $id (@genbIDs) {
		$id =~ s/;$//;
		my $seqName = "";
		
		# Again: Sometimes different genbank IDs are noted with
		# seqName of the nucleotides
		if ($id =~ m/:/) {
			my @seqNames = split(/:/, $id);
			$seqName = $seqNames[0];			
			$id = $seqNames[1];
			$id =~ s/^ //;
		}
		# Remove trailing and leading whitespaces in ids
        $id =~ s/^\s+|\s+$//; 
		
		if (exists $ids{$id}) {
			print "ERROR: GenBank id $id is not unique. Removed duplicate entry.\n";
            print "\tCURRENT: $tax" . ": $seqName\n\tFIRST SEEN:" . $ids{$id} . "\n";
			next;
		}
		else {
			$ids{$id} = $tax . ": $seqName";
		}

		
		if ($seqName ne "") {
			print OUTTAX $id.";".$tax."_".$seqName."\n";
		}
		else{
			print OUTTAX $id.";".$tax."\n";
		}
		
		foreach my $host (@hosts) {
			
			$host =~ s/^ //;
			
			# Host may be a species name or isolation source
			if ($host =~ m/ \(S\)$/) {
					
				$host =~ s/ \(S\)$//;
				$host = "source: ".$host;
			}
			
			print OUTHOST $id.";".$host."\n";
		}
	}
	
}
close(INTAX);
close(OUTTAX);
close(OUTHOST);

#=============================================================================================#
# Requesting nucleotide sequences by genebank accession from NCBI.
#---------------------------------------------------------------------------------------------#
my @ids = keys(%ids);
my $i = 0;
my $chunkSize = 1000;

print "INFO: Web request to NCBI: ".@ids.
	" sequences. Sending chunks of $chunkSize. This may take a while.\n";

# Obtain sequences and IDs and write as fasta file
open(OUTSEQ ,">", $outSeqP);

my @entries = ();

while (@ids) {
        my @queries = splice (@ids, 0, $chunkSize);
        
        # POST request to NCBI eFetch
        my $ua = LWP::UserAgent->new;
        $ua->timeout(30);
        $ua->env_proxy;

        # All genbank IDs in the input taxonomy
        my $query = join(",", @queries);
        my %query = ("db"=>"nuccore", "id"=>$query, "rettype"=>"fasta", "retmode"=>"text");

        # Response from NCBI
        my $seq = $ua->post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", \%query); # no www. causes error!!!

        # Crash, if there was an error
        die "ERROR: Error querying NCBI\n" unless $seq->is_success;

        # Response body
        $seq = $seq->content;
        $seq =~ s/^\s+|\s+$//;
        $seq =~ s/^>//;
        
        my @res = split(">", $seq);
        chomp(@res);
        push(@entries, @res);

        # Maximum three request per second
        select(undef, undef, undef, 0.34);
}

print "INFO: Retrieved " . @entries . " sequences from NCBI.\n";

foreach my $entry (@entries) {
        next if ($entry eq "");
        my @lines = split("\n", $entry);
        
        my $id = (split(" ", $lines[0]))[0];
        
        # Some GenBank ids were sent without a version tag,
        # but I received a response with a tag.
        # In that case, I need to save the sequence data
        # without the tag.
        # If I queried GenBank with a tag, I want to keep it.
        # This keeps the naming consistent between the output files.
        my $idNoTag = (split(/\./, $id))[0];
        
        if (exists $ids{$id}) {
                delete $ids{$id}
        }
        elsif (exists $ids{$idNoTag}) {
                delete $ids{$idNoTag};
                $id = $idNoTag
        }
        else {
                next;
        }

        my $seq = "";
        
        foreach my $line (@lines[1..$#lines]) {
                $seq .= $line;          
        }
        
        $seq .= "\n";
        
        print OUTSEQ ">".$id."\n".$seq;
        
        $i++;
}
close(OUTSEQ);

print "WARNING: Could not find GenBank ids: " . join(", ", keys(%ids)) . "\n" if (%ids);
print "\nFinished!\n";