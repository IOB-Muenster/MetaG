# Access
Due to file size limitations, the databases can be found here:

* GRCh38p13 (filter workflow): http://bioinformatics.uni-muenster.de/tools/metag/download/GRCh38p13.zip
* ICTV_vmr: http://bioinformatics.uni-muenster.de/tools/metag/download/ICTV_vmr.zip
* MTXssuLSU from Metaxa2.2 : http://bioinformatics.uni-muenster.de/tools/metag/download/MTXssuLSU.zip
* PATRIC: http://bioinformatics.uni-muenster.de/tools/metag/download/PATRIC.zip
* RDP16s28s version 11, update 5: http://bioinformatics.uni-muenster.de/tools/metag/download/RDP16s28s.zip
* RefSeq: http://www.bioinformatics.uni-muenster.de/tools/metag/download/refseq.zip

Use the md5sums.txt file to verify your download: http://bioinformatics.uni-muenster.de/tools/metag/download/md5sums.txt

# How to use the databases with MetaG?
Supply "tax." files as -ltax and "patho." files as -pdbPath.
The other files are in LAST format. Use -ldb, but don't give
the file extension. E.g.: -ldb ICTV_vmr/ICTV_vmr
The "version" files indicate the versions of the databases.
For the databases used in the filter workflow, no "tax." and
"patho." files are needed.

# References

## GRCh38p13 (filter workflow)
Human genome used in the filter workflow of MetaG. It was downloaded from RefSeq.

O'Leary NA, Wright MW, Brister JR, Ciufo S, Haddad D, McVeigh R, Rajput B,
Robbertse B, Smith-White B, Ako-Adjei D, Astashyn A, Badretdin A, Bao Y,
Blinkova O, Brover V, Chetvernin V, Choi J, Cox E, Ermolaeva O, Farrell CM,
Goldfarb T, Gupta T, Haft D, Hatcher E, Hlavina W, Joardar VS, Kodali VK, Li W,
Maglott D, Masterson P, McGarvey KM, Murphy MR, O'Neill K, Pujar S, Rangwala SH,
Rausch D, Riddick LD, Schoch C, Shkeda A, Storz SS, Sun H, Thibaud-Nissen F,
Tolstoy I, Tully RE, Vatsan AR, Wallin C, Webb D, Wu W, Landrum MJ, Kimchi A,
Tatusova T, DiCuccio M, Kitts P, Murphy TD, Pruitt KD. 2016.
	Reference sequence (RefSeq) database at NCBI: current status, taxonomic
	expansion, and functional annotation.
	Nucleic Acids Res. 44(D1):D733-45. doi: 10.1093/nar/gkv1189.

## ICTV_vmr
This data base was derived from: https://ictv.global/vmr/
See the ICTV publication.

Walker PJ, Siddell SG, Lefkowitz EJ, Mushegian AR, Dempsey DM,
Dutilh BE, Harrach B, Harrison RL, Hendrickson RC, Junglen S,
et al. 2019.
	Changes to virus taxonomy and the International Code of Virus
	Classification and Nomenclature ratified by the International
	Committee on Taxonomy of Viruses (2019).
	Arch Virol 164: 2417–2429. doi:10.1007/s00705-019-04306-w.

## MTXssuLSU from Metaxa2.2
This database was derived from the Metaxa2.2 program.

Bengtsson-Palme J, Hartmann M, Eriksson KM, Pal C, Thorell K,
Larsson DGJ, Nilsson RH. 2015.
	METAXA2: improved identification and taxonomic classification
	of small and large subunit rRNA in metagenomic data.
	Mol Ecol Resour 15: 1403–1414. doi:10.1111/1755-0998.12399.

It is distributed under the GNU General Public License. 

	This program is free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or(at your option) any later version. This program
	is distributed in the hope that it will be useful, but WITHOUT
	ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	General Public License for more details. You should have received
	a copy of the GNU General Public License along with this program,
	in a file called 'license.txt'. If not, see:
	http://www.gnu.org/licenses/

## PATRIC
See the PATRIC publication.

Wattam AR, Davis JJ, Assaf R, Boisvert S, Brettin T, Bun C,
Conrad N, Dietrich EM, Disz T, Gabbard JL, et al. 2017.
	Improvements to PATRIC, the all-bacterial Bioinformatics Database
	and Analysis Resource Center.
	Nucleic Acids Res 45: D535–D542. doi:10.1093/nar/gkw1017.

## RDP16s28s version 11, update 5
See the RDP publication.

Cole JR, Wang Q, Fish JA, Chai B, McGarrell DM, Sun Y, Brown CT,
Porras-Alfaro A, Kuske CR, Tiedje JM. 2014.
	Ribosomal Database Project: data and tools for high throughput rRNA analysis.
	Nucleic Acids Res 42: D633–D642. doi:10.1093/nar/gkt1244.
	
## RefSeq Targeted Loci
The database contains bacterial (16S, 23S, 5S rRNA), archaeal (16S, 23S, 5S rRNA) and fungal
sequences (18S, 28S rRNA, ITS) downloaded from ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/ .

Sayers EW, Beck J, Bolton EE, Bourexis D, Brister JR, Canese K, Comeau DC,
Funk K, Kim S, Klimke W, et al. 2021.
	Database resources of the National Center for Biotechnology Information.
	Nucleic Acids Res 49: D10–D17. doi:10.1093/nar/gkaa892


