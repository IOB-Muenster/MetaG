# Access
Due to file size limitations, the databases can be found here:

* T2T-CHM13v2 (filter workflow): https://www.bioinformatics.uni-muenster.de/tools/metag/download/T2T-CHM13v2.zip
* ICTV_vmr: https://bioinformatics.uni-muenster.de/tools/metag/download/ICTV_vmr.zip
* MTXssuLSU from Metaxa2.2 : https://bioinformatics.uni-muenster.de/tools/metag/download/MTXssuLSU.zip
* BV-BRC (formerly: PATRIC): https://bioinformatics.uni-muenster.de/tools/metag/download/BVBRC.zip
* RDP16s28s version 11, update 5: https://bioinformatics.uni-muenster.de/tools/metag/download/RDP16s28s.zip
* RefSeq: https://www.bioinformatics.uni-muenster.de/tools/metag/download/refseq.zip

Use the md5sums.txt file to verify your download: https://bioinformatics.uni-muenster.de/tools/metag/download/md5sums.txt

# How to use the databases with MetaG?
Supply "tax." files as -ltax and "patho." files as -pdbPath.
The other files are in LAST format. Use -ldb, but don't give
the file extension. E.g.: -ldb ICTV_vmr/ICTV_vmr
The "version" files indicate the versions of the databases.
For the databases used in the filter workflow, no "tax." and
"patho." files are needed.

# References

## T2T-CHM13v2 (filter workflow)
Human genome used in the filter workflow of MetaG. It was downloaded from RefSeq (GCF_009914755.1).

Nurk S, Koren S, Rhie A, et al. 2022.
	The complete sequence of a human genome.
	Science 376(6588):44-53. doi:10.1126/science.abj6987.

## ICTV_vmr
This data base was derived from: https://ictv.global/vmr/
See the ICTV publication.

Walker PJ, Siddell SG, Lefkowitz EJ, Mushegian AR, Dempsey DM,
Dutilh BE, Harrach B, Harrison RL, Hendrickson RC, Junglen S,
et al. 2019.
	Changes to virus taxonomy and the International Code of Virus
	Classification and Nomenclature ratified by the International
	Committee on Taxonomy of Viruses.
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

## BV-BRC
See the BV-BRC publication.

Olson RD, Assaf R, Brettin T, Conrad N, Cucinell C, Davis JJ,
Dempsey DM, Dickerman A, Dietrich EM, Kenyon RW, et al. 2023.
	Introducing the Bacterial and Viral Bioinformatics Resource
	Center (BV-BRC): a resource combining PATRIC, IRD and ViPR.
	Nucleic Acids Res 51:D678-D689. doi: 10.1093/nar/gkac1003

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


