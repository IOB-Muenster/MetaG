#!/bin/sh

#==========================================================================================#
# Launcher for the MetaG pipeline including LASTAL, LAST-SPLIT and the getMeta script.
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


#=====================================================================#
# Functions
#---------------------------------------------------------------------#


# Test if value is valid or an argument name.
# Is argument, if user forgot to give value to argument before current.
testValue () {
	echo "$1" | grep -qG "^-[a-z]\|^--" && { printf "\nERROR: Argument found when value expected. Did you forget to give a value to the argument before [${1}]?\n\n"; exit 1; }
}

# Increment the argument counter by $2, if variable has not been set
countOpt () {
	if [ ! "$1" ]; then ARGCOUNT=$(expr $ARGCOUNT + $2); fi
}

# Increment the argument counter by 1, if switch has been set to $2
countSwitch () {
	if [ "$1" != "$2" ]; then ARGCOUNT=$(expr $ARGCOUNT + 1); fi
}

# Very basic test for fasta format
# First non-empty line must start with a >
testFasta () {
	ISHEADER=$(perl -ne 'if ($_ !~ m/^\s*$/) { if (m/^>/) {print "0"; last} else {print "1"; last}}' "$1")
	return $ISHEADER	
}

# Check if a number of essential parameters was set and is not empty
checkEssential() {
	for VAR in "$@"; do
		if [ -z $VAR ]; then
			printf "\nERROR: Missing an essential variable. Type -h for help.\n"
			printf "ERROR: Received arguments:\n"
			printf -- "${INPUT}"
			
			# The user gave matrix in config file. Thus, -lmat parameter does not exist in $INPUT,
			# but should be shown here for debugging
			if [ "$LMAT" ]; then printf -- " -lmat ${LMAT}"; fi
			
			printf "\n"
			exit 1
		fi
	done
}

calcAlign () {
	# MAF output of LASTAL | LAST-SPLIT is saved into the same directory as query
	# Name alignment produced by --filter and alignment produced by taxonomic assignment
	# differently.
	MAF=""
	if [ $ISFILTER -eq 0 ]; then
		MAF=$(echo $Q | sed -e 's/[^/]*$/calc.lastalign.maf/')
	else
		MAF=$(echo $Q | sed -e 's/[^/]*$/calc.filter.maf/')
	fi
	
	if [ -e "$MAF" ]; then
		printf "\nERROR: Alignment file ${MAF} already exists. Won't override it.\nDo you want to use this file for the analysis, instead? Set [--no_align].\n\n"
		exit 1
	fi
	
	# Create LASTAL command. Depends on: LAST-paramters: y/n, LAST-subst.matrix: y/n. However, one of the two is needed.
	if [ "$LMAT" = "" -a "$LPARAM" = "" ]; then
		printf "\nERROR: [-lmat] and/or [-lparam] need to be given. Use -h for help.\n\n"
		exit 1
	elif [ "$LMAT" = "" ]; then
		LASTCOMMAND="$LASTAL -P $CORES $LPARAM $LDB $Q"
	elif [ "$LPARAM" = "" ]; then
		LASTCOMMAND="$LASTAL -P $CORES -p $LMAT $LDB $Q"
	else
		LASTCOMMAND="$LASTAL -P $CORES -p $LMAT $LPARAM $LDB $Q"
	fi
	
	#Run LASTAL and optionally LAST-SPLIT
	if [ "$LSPLIT" = "" ]; then
		printf "\n[1/2] Aligning with LAST. This may take a while.\n"
		$LASTCOMMAND > "$MAF" || { printf "\nERROR: Error during alignment\n\n"; exit 1; }
	else
		printf "\n[1/2] Aligning with LAST and LAST-SPLIT. This may take a while.\n"
		$LASTCOMMAND | "$LASTSPLIT" -m "$LSPLIT" -n -fMAF > "$MAF" || { printf "\nERROR: Error during alignment\n\n"; exit 1; }
	fi
}


filterReads () {
	FILTERSCR=$(echo ${GETMETA} | sed 's/getMeta\.pl/filterReads\.pl/')
	FILTERED=$(echo $Q | sed -e 's/[^/]*$/calc\.query\.filtered/')
	"$FILTERSCR" -alignment "$MAF" -sequences "$Q" --unaligned > "$FILTERED"
}


getTaxa () {
	
	# Verbose output?
	if [ $ISVERBOSE -eq 1 ]; then
		
		# Viral database?
		if [ $ISVIR -eq 1 ]; then
			"$GETMETA" -q "$Q" -a "$MAF" -e "$E" -ac "$AC" -cc "$CC" -m "$M" -db "$LTAX" -pdb "$PDBPATH" -kro "$KRONA" -vir
		else
			"$GETMETA" -q "$Q" -a "$MAF" -e "$E" -ac "$AC" -cc "$CC" -m "$M" -db "$LTAX" -pdb "$PDBPATH" -kro "$KRONA"
		fi
	else
		
		# Viral database?
		if [ $ISVIR -eq 1 ]; then
			"$GETMETA" -q "$Q" -a "$MAF" -e "$E" -ac "$AC" -cc "$CC" -m "$M" -db "$LTAX" -pdb "$PDBPATH" -kro "$KRONA" -vir > /dev/null
		else
			"$GETMETA" -q "$Q" -a "$MAF" -e "$E" -ac "$AC" -cc "$CC" -m "$M" -db "$LTAX" -pdb "$PDBPATH" -kro "$KRONA" > /dev/null
		fi
	fi	
}


packUnmatched () {
	
	# Get database ranks from tax file, else take default
	RANKNAMES=$(head -n 1 $LTAX | grep -E "^#")
	RANKTYPE="custom"
	
	if [ $RANKNAMES ];then
		RANKNAMES=$(echo $RANKNAMES | sed -e 's/#//')
		RANKNAMES=$(echo $RANKNAMES | sed -e 's/;/\\n/g')
	else
		RANKNAMES="domain\nphylum\nclass\nsubclass\norder\nsuborder\nfamily\ngenus\nspecies\nstrain"
		RANKTYPE="default"
	fi
		
	
		
	
	UNMATCHEDMSG="
#=================================================================================#
    Readme for the unmatched output of MetaG
#---------------------------------------------------------------------------------#

In this folder you will find at least one calc.unmatched.* fasta file. The * stands 
for multiple possible endings. The endings give you information about the first
unmatched rank of a query read. After a read was declared unmatched for a
given rank, it is also unmatched for all following ranks.

E.g.: Read A is unmatched for the phylum. It is first reported in the
calc.unmatched.phylum file, but also appears in calc.unmatched.class,
calc.unmatched.subclass, calc.unmatched.order,..

Rank endings correspond to the ranks defined in the databases. If you re-run your
analysis with different databases, you may find different endings.
When custom rank definitions are missing from the database, MetaG uses default ranks.
Your current database contains the ${RANKTYPE} ranks:

${RANKNAMES}

The special ending [.removed] indicates that MetaG was either not able to align the
read at all or that it was removed due to your e-value cutoff. Reads in the
[.removed] file are not reported at the lower ranks. This is opposite to the behavior
described above.

The fasta files follow this general format:
>header1
seq1
>header2
seq2
	"
	
	# Check if files with unmatched reads exist and zip them together with a readme
	for FILE in ${DIR}/calc.unmatched.*; do
	    if [ -e "$FILE" ]; then
	    	printf "$UNMATCHEDMSG" > ${DIR}/calc.readme
	    	# Test archive; remove input files; no status report; don't store full path in zip
	    	# Redirect to /dev/null: Gets rid of unzip error message and output in FreeBSD. Does the job
	    	# regardless of unzip error.
	    	zip -Tmqj ${DIR}/calc.unmatched.zip ${DIR}/calc.unmatched.* ${DIR}/calc.readme > /dev/null 2>&1
	    fi
	    break
	done
}

# Print only the values used, if arguments appear more often
printArgs() {
	if [ ! -z "$MAF" ]; then INFO="${INFO}SKIPPING LAST alignment\n\tPATH to MAF file:\t\t${MAF}\n\t"; fi
	if [ ! -z "$Q" ]; then INFO="${INFO}DNA reads:\t\t\t${Q}\n\t"; fi
	if [ "$ISFQ" != 0 ]; then INFO="${INFO}FASTQ input:\t\t\tTrue\n\t"; fi
	if [ "$ISREC" != 0 ]; then ALIGNINFO="${ALIGNINFO}Recursive:\t\t\tTrue\n\t"; fi
	if [ ! -z "$LDB" ]; then ALIGNINFO="${ALIGNINFO}LAST database:\t\t\t${LDB}\n\t"; fi
	if [ ! -z "$LTAX" ]; then ASSIGNINFO="${ASSIGNINFO}Database taxonomy:\t\t${LTAX}\n\t"; fi
	if [ ! -z "$LMAT" ]; then ALIGNINFO="${ALIGNINFO}LAST substitution matrix:\t${LMAT}\n\t"; fi
	if [ ! -z "$LPARAM" ]; then ALIGNINFO="${ALIGNINFO}LAST parameters:\t\t${LPARAM}\n\t"; fi
	if [ ! -z "$LSPLIT" ]; then ALIGNINFO="${ALIGNINFO}LAST-SPLIT -m:\t\t\t${LSPLIT}\n\t"; fi
	if [ ! -z "$PDBPATH" ]; then ASSIGNINFO="${ASSIGNINFO}Pathogen taxonomy:\t\t${PDBPATH}\n\t"; fi
	if [ ! -z "$E" ]; then ASSIGNINFO="${ASSIGNINFO}E-value cutoff:\t\t\t${E}\n\t"; fi
	if [ ! -z "$AC" ]; then ASSIGNINFO="${ASSIGNINFO}Alignment score cutoff:\t\t${AC}\n\t"; fi
	if [ ! -z "$CC" ]; then ASSIGNINFO="${ASSIGNINFO}Confidence cutoff:\t\t${CC}\n\t"; fi
	if [ ! -z "$M" ]; then ASSIGNINFO="${ASSIGNINFO}Method for average confidence:\t${M}\n\t"; fi
	if [ "$ISVIR" != 0 ]; then ASSIGNINFO="${ASSIGNINFO}Viral database:\t\t\tTrue\n\t"; fi
	if [ "$ISVERBOSE" != 0 ]; then INFO="${INFO}Verbose:\t\t\tTrue\n\t"; fi
	if [ "$ISFILTER" != 0 ]; then INFO="${INFO}Filter reads:\t\t\tTrue\n\t"; fi
}


#======================================================================================#
# Read command line argument
#--------------------------------------------------------------------------------------#
readArgs() {
	
	# Number of command line arguments
	ARGCOUNT=0
	
	# Get arguments
	INPUT="$@"
	
	# Default core number for LASTAL
	CORES=1

	INFO="
#==========================================================#
# MetaG -the metagenomics pipeline
#----------------------------------------------------------#

Metagenomic analysis started with parameters:\n\n\t"
	
	if [ -f "$CONFIG" ]; then INFO="${INFO}Reading from config file:\t${CONFIG}\n\t"; fi
	
	while [ "$1" != "" ]; do
		case $1 in
			--no_align)	shift
						testValue $1
						countOpt "$MAF" "2"
						NOALIGN=1
						MAF=$1
						;; 
			-q)		shift
					testValue $1
					countOpt "$Q" "2"
					Q=$1
					;;
			--fastq)	
					countSwitch "$ISFQ" "1"
					ISFQ=1
					;;
			--recursive) ISREC=1
					countSwitch "$ISREC" "1"
					;;
			-ldb)		shift
					testValue $1
					countOpt "$LDB" "2"
					LDB=$1
	
					# Check if user accidentally entered path to database file and not to prefix
					if [ -e "$LDB" ]; then
						printf "\nERROR: [-ldb]: ${LDB} needs to be path to the database prefix and not an existing file.\nUse -h for help.\n\n"
						exit 1
					elif [ ! -e "$LDB"".bck" ]; then
						printf "\nERROR: No such LAST database: ${LDB}. Use -h for help.\n\n"
						exit 1
					fi
					;;
			-ltax)		shift
					testValue $1
					countOpt "$LTAX" "2"
					LTAX=$1		
					;;
			-lmat)		shift
					testValue $1
					countOpt "$LMAT" "2"
					LMAT=$1			
					;;
			-lparam)	shift
					countOpt "$LPARAM" "2"
					LPARAM=$(echo $1 | sed -e 's/_/ /g') #work around, if from config file
					;;
			-lsplit)	shift
					testValue $1
					countOpt "$LSPLIT" "2"
					LSPLIT=$1
					;;
			-cores)		shift
					testValue $1
					countOpt "$CORES" "2"
					CORES=$1
					;;
			-pdbPath)	shift
					testValue $1
					countOpt "$PDBPATH" "2"
					PDBPATH=$1				
					;;
			-e)		shift
					testValue $1
					countOpt "$E" "2"
					E=$1
					;;
			-ac)	shift
					testValue $1
					countOpt "$AC" "2"
					AC=$1
					;;
			-cc)	shift
					testValue $1
					countOpt "$CC" "2"
					CC=$1
					;;
			-m)		shift
					testValue $1
					countOpt "$M" "2"
					M=$1
					;;
			-vir)	countSwitch "$ISVIR" "1"
					ISVIR=1
					;;
			--verbose)	ISVERBOSE=1
					countSwitch "$ISVERBOSE" "1"
					;;
			--filter)	ISFILTER=1
					countSwitch "$ISFILTER" "1"
					;;
			--config) shift
					;;
			-h)		man metag 2>/dev/null || printf "${HELPMSG}" | less
					exit 0
					;;
			*)		printf "\nERROR: Unknown parameter ${1}. Use -h for help.\n\n"
					exit 1
					;;
		esac
		
		# Get error message from shift and replace with more helpful one.
		ERR=$(shift 2>&1)  || { printf "\nERROR: Too few values for given number of arguments. Use -h for help.\n\n" ; exit 1; }
		shift
	done
	
	
	if [ $ARGCOUNT -lt 5 -a  $ARGCOUNT -gt 28 ]; then
		printf "\nERROR: Unexpected number of arguments. Type -h for help\n\n"
		exit 1
	fi
}





# Get environment variables
LASTAL=$(printenv LASTAL)
LASTSPLIT=$(printenv LASTSPLIT)
GETMETA=$(printenv GETMETA)
KRONA=$(printenv KRONA)


# Defaults
NOALIGN=0
ISVIR=0
ISVERBOSE=0
ISFQ=0
ISREC=0
ISFILTER=0
CONFIG=""
ENV=""
CONFREADER="${GETMETA}"
MATRIX=""
INPUT="None"

# This is for information about parameters.
# Only lists those which have an effect on --no_align
# and alignment workflow.
INFO=""

# Same as INFO, but lists parameters which only affect
# alignment workflow, but not --no_align
ALIGNINFO=""

# Same as INFO, but lists parameters which only affect
# metagenomic assignment, but not --filter
ALIGNINFO=""

HELPMSG='
MetaG(1)                                                                                        Documentation                                                                                        MetaG(1)

NAME
       MetaG - performs metagenomic analyses.

SYNOPSIS
       For help

              metag.sh -h

       Filter reads
       		  
       		  metag.sh --filter  -q  query reads -ldb database -lmat subst. matrix | -lparam alignment parameters [-lsplit filtering with LAST-SPLIT] [-cores number of cores for alignment] [--verbose]
       		  [--fastq] [--config  config  file] [--recursive]
       		  
       Filter reads from pre-calculated alignment
       
       		  metag.sh --filter -q  query reads --no_align MAF file [--verbose] [--fastq] [--config  config file]
       
       Metagenomic analysis from reads

              metag.sh  -q  query reads -ldb database -ltax database taxonomy -lmat subst. matrix | -lparam alignment parameters [-lsplit filtering with LAST-SPLIT] [-cores number of cores for alignment]
              -pdbPath pathogen taxonomy [-vir] -e e-value cutoff -ac alignment score cutoff -cc confidence cutoff -m Method to  calc.  average  confidence  [--verbose]  [--fastq]  [--config  config  file]
              [--recursive]

       Metagenomic analysis with pre-calculated alignment

              metag.sh  -q  query  reads  -ltax  database taxonomy -pdbPath pathogen taxonomy [-vir] -e e-value cutoff -ac alignment score cutoff -cc confidence cutoff -m Method to calc. average confidence
              --no_align MAF file [--verbose] [--fastq] [--config  config file]

DESCRIPTION
       MetaG processes input sequencing reads and assigns a taxonomy. This is supplemented with pathogen predictions and, where applicable, antibiotic resistance information.

       The program can be run on fasta/fastq reads from WGS and targeted sequencing and was tested for bacteria, archaea, fungi and viruses. Calculations may start from scratch by alignment of query reads
	   to a database or from a pre-calculated alignment in MAF format. Optionally, it is possible to filter reads first, by using the --filter workflow. This can be used to improve the analysis for species
	   which occur at low read counts, by excluding taxa with very high read counts. It can also be used to exclude contamination which may lead to false positive classifications (e.g. human contamination in
	   a viral analysis). After filtering, the metagenomic analysis can be performed on the retained reads.

       MetaG filters the alignments for each read by a given e-value and alignment score cutoff. If one read can be assigned to different database taxonomies, MetaG assigns the most common taxon. To avoid
       ambiguous assignments, this is also influenced by the quality of the matches and the user given confidence cutoff.

       MetaG needs multiple environment variables. These can also be specified in the config file. See also ENVIRONMENT section.

OPTIONS
       -h     OPT: Show this help message and exit.

       -q     Query reads in fasta or fastq format [FILE] Can be a directory, if --recursive is specified.

       --fastq
              OPT: Indicates fastq reads.

       -ldb   Name of LAST database [NAME].

       -ltax  Database taxonomy [FILE].

       -lmat  LAST substitution matrix. Conflicting with -lparam [FILE].

       -lparam
              LAST alignment parameters, Conflicting with -lmat [STRING].

       -lsplit
              OPT: Filters alignments with LAST-SPLIT.  Translates to last-split -m [FLOAT].

       -cores OPT: Cores to use for LAST alignment [INT].

       -pdbPath
              Pathogen taxonomy [PATH].

       -e     E-value cutoff for alignment filtering [INT].

       -ac    Alignment score cutoff for alignment filtering [INT].

       -cc    Confidence cutoff for taxonomic assignments [INT].

       -m     Method to calculate the average confidence for a taxon [williams, geometric, arithmetic, harmonic].

       --no_align
              Use this pre-calculated MAF alignment [FILE].  Cannot be used with --recursive.

       -vir   OPT: Indicates analysis of viruses.

       --config
              Read parameters from this config file [FILE]. If the same parameters appear in the file and on the command line, the command line overwrites the values in the file.  For exceptions see  EXAM‐
              PLES.

       --verbose
              Print verbose messages.

       --recursive
              Processes multiple query files in the directory defined by -q.  Illegal, if --no_align.
              
       --filter
       		 Indicate the workflow to remove unwanted reads.

EXAMPLES
       All parameters can also be defined in a so-called config file, which may be supplied to MetaG by using the --config parameter. Additional parameters can be specified on the command line.  In case of
       conflicts, parameters on the command line will overwrite those in the config file. This, however, does not apply to the parameters --no_align, --verbose, -vir, --recursive and the environment  vari‐
       ables.  Once set in the config file, they will be used by MetaG.

       Config files can be created interactively, see also metag_setup, or manually. This is the accepted format:

       For all parameters except -lsplit, -lparam and the environment variables:

       #metag param [value]

       For -lsplit:

       #lastsplit -m value

       For -lparam each LAST parameter must be on a single line following this format:

       #last param-1 value-1

       [..]

       #last param-N value-N

       For the environment variables:

       #env param value

ENVIRONMENT
       MetaG depends on the following environment variables to access helper scripts.  These can also be defined in the config file.

       LASTAL Location of LAST aligner [PATH]. OPT, if --no_align.

       LASTSPLIT
              Location of LAST-SPLIT [PATH]. OPT, if --no_align was set or -lsplit was not set.

       GETMETA
              Location of the getMeta.pl script of MetaG [PATH].

       KRONA  Location of Krona [PATH].

FILES
       Output files will be generated in the directory of the query file. The files always start with the calc. prefix with the exception of temp.matrix.txt.
       
   	   calc.query.filtered
   	   		  Contains the reads that were retained in the --filter workflow. Only present, if --filter was used.

       calc.lastalign.maf
              MAF file of LAST alignment. Not present, if --no-align.

       calc.VIS.html
              Interactive graph of taxa abundances.

       calc.VIS.txt
              Raw data for generation of calc.VIS.html.

       calc.PATHO.html
              Interactive graph for pathogen metadata (host and optionally antibiotic resitances).

       calc.PATHO.txt
              Raw data for generation of calc.PATHO.html.

       calc.LOG.txt
              Table  with  read statistics and found abundances of taxa split by taxonomic rank.  Abundances are supplemented with their respective average confidences, followed by the averaging method: w:
              williams, a: arithmetic, g: geometric and h: harmonic.

       calc.LIN.txt
              File displaying the taxonomic assignment per read:
              >readID

              rank1: taxon1.1 supporting matches at current rank (total matches at current rank)

              rank2: taxon1.2 supporting matches at current rank (total matches at current rank)

              [...]

              rankn: taxon1.n supporting matches at current rank (total matches at current rank)

              OR

              No match for readID and chosen cutoff.

       calc.UNMATCHED.zip
              Contains separate files with unmatched reads at each rank or reads which were filtered prior to the taxon assignment. See also calc.readme inside the archive. Not present, if  all  reads  had
              been matched.

       temp.matrix
              Only present, if LAST matrix was specified in the config file of MetaG and not on the command line. This is a temporary file containing the LAST matrix.

       In  case  --recursive  was  specified, folders for each valid query file are created in the input directory.  These contain the individual output files for a single query file (see also above).  The
       directory names start with an abbreviation of the query file and end with an underscore followed by a counter. The counter avoids conflicts with query files of similar names.

DEPENDENCIES
       MetaG depends on perl v5.26.1, LAST and LAST-SPLIT 963 and KronaTools 2.7 or later.  All dependencies, except Perl, are installed with metag_setup.  For more information on a manual installation  of
       the dependencies, see metag_setup.

       MetaG was tested on Ubuntu 18.04.3 LTS, FreeBSD 12.1 and macOS 10.11.6.

AUTHOR
       Felix Manske, felix.manske@uni-muenster.de.

       Norbert Grundmann

       Author for correspondence: Wojciech Makalowski, wojmak@uni-muenster.de.

COPYRIGHT
       Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

       1. Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

       2. Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.

       THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
       ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,  PROCUREMENT  OF  SUBSTITUTE
       GOODS  OR  SERVICES;  LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
       OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

SEE ALSO
       metag_setup(1)

Felix Manske                                                                                         2020                                                                                            MetaG(1)
'

#==============================================================================================================#
# Process arguments from user
#--------------------------------------------------------------------------------------------------------------#

# Pre-process

# User has a config file with additional parameters
# Join command line parameters with those from config file
CONFIG=$(echo "$@" | grep -Eo "\-\-config [^- ]*" | cut -d " " -f2)


# I need GETMETA's value to find readConfig.pl script
# Still user should be able to access help without a crash
HELP=""
HELP=$(echo "$@" | grep "\-h")

if [ "$HELP" ]; then
	man metag 2>/dev/null || printf "${HELPMSG}" | less
	exit 0
fi

# Get path of readConfig.pl
# Mandatory, if config file is given
if [ ! "$CONFREADER" -a -f "$CONFIG" ]; then
	CONFREADER=$(grep -Eo "GETMETA [^ ]*" "$CONFIG" | cut -d " " -f2)
fi

if [ "$CONFREADER" ]; then
	CONFREADER=$(echo $CONFREADER | sed -e 's/getMeta\.pl/readConfig\.pl/')
elif [ ! "$CONFREADER" -a -f "$CONFIG" ]; then
	printf "\nERROR: Environment variable GETMETA was not set.\n\n"
	exit 1
fi
	

if [ "$CONFIG" ]; then
	
	CONFIGPARAMS=$($CONFREADER "${CONFIG}")
	ENV=$(echo $CONFIGPARAMS | cut -d '#' -f2)
	MATRIX=$(echo $CONFIGPARAMS | cut -d '#' -f3)
	MATRIX=$(echo "${MATRIX}" | sed -e 's/_$//' | sed -e 's/_/\\n/g')
	
	CONFIGPARAMS=$(echo $CONFIGPARAMS | cut -d '#' -f1)
	
	
	# Clear $@, take params from config file and append old values of $@.
	# => Command line parameters will have priority over config file parameters
	# => EXCEPT: Environment variables
	# Escaping blanks on the shell works, now
	set -- $CONFIGPARAMS "$@"
fi

# User supplied environment variables in  config file. Default is obtained by printenv --> user set them in shell	
if [ "$ENV" ]; then
	LASTAL=$(echo "$ENV" | grep -Eo "LASTAL=[^ ]*" | cut -d "=" -f2)
	LASTSPLIT=$(echo "$ENV" | grep -Eo "LASTSPLIT=[^ ]*" | cut -d "=" -f2)
	GETMETA=$(echo "$ENV" | grep -Eo "GETMETA=[^ ]*" | cut -d "=" -f2)
	KRONA=$(echo "$ENV" | grep -Eo "KRONA=[^ ]*" | cut -d "=" -f2)
fi

# Process all supplied arguments
readArgs "${@}"
printArgs

# Matrix was not given on CLI, but in config file
if [ -z "$LMAT" -a "$MATRIX" ]; then
	
	# Write matrix
	# $Q can be a directory in case of recursive analysis 
	if [ -d "$Q" ]; then
		MATRIXPATH="$Q"
	else
		MATRIXPATH=$(dirname "$Q")
	fi
	
	printf "${MATRIX}" > "${MATRIXPATH}/temp.matrix" ####output matrix looks strange, but works
	
	LMAT="${MATRIXPATH}/temp.matrix"
	ALIGNINFO="${ALIGNINFO}LAST substitution matrix:\t${LMAT}\n\t"
fi


# User-modified cores or default are printed
ALIGNINFO="${ALIGNINFO}Cores for LAST:\t\t\t${CORES}\n\n"

#==============================================================================================================#
# Check user input in detail
#--------------------------------------------------------------------------------------------------------------#


# Use pre-calculated alignment
# User can give fewer arguments, if --no_align is set
if [ $NOALIGN -eq 1 ]; then
	
	# Perform metagenomic assignment
	if [ $ISFILTER -eq 0 ]; then
		# Check if relevant parameters for workflow were set and are not empty
		checkEssential "$Q" "$LTAX" "$PDBPATH" "$E" "$AC" "$CC" "$M"
		
		INFO="${INFO}${ASSIGNINFO}"
		
		# ...Then by checking arguments that should NOT have been set
		# These won't have any effect using --no_align.
		# Just issue a warning.
		if [ ! -z $LDB ] || [ ! -z $LMAT ] || [ ! -z $LPARAM ] || [ ! -z $LSPLIT ] || [ $CORES != 1 ] || [ $ISREC != 0 ]; then
			INFO="${INFO}\nFollowing parameters have no effect with --no_align and are ignored:\n\n\t${ALIGNINFO}"
		fi	
		
		# Check environment variables for GETMETA and KRONA
		if [ "$GETMETA" = "" -o "$KRONA" = "" ]; then
			printf "\nERROR: Unset essential environment variable. Use -h for help\n\n"
			exit 1
		else
			# Set path to major Krona script
			KRONA="${KRONA}/ImportText.pl"
			
			# Check if scripts exist
			if [ ! -e "$GETMETA" ]; then
				printf "\nERROR: Assignment script not found at ${GETMETA}. Check environmental variable or use -h for help.\n"
				exit 1
			elif [ ! -e "$KRONA" ]; then
				printf "\nERROR: Krona main script not found at ${KRONA}. Check environmental variable or use -h for help.\n"
				exit 1
			fi
		fi
	# Filter reads
	else
		# Check if relevant parameters for workflow were set and are not empty
		checkEssential "$Q"
		
		# Check if scripts exist
		if [ ! -e "$GETMETA" ]; then
			printf "\nERROR: Assignment script not found at ${GETMETA}. Check environmental variable or use -h for help.\n"
			exit 1
		fi
		
		# ...Then by checking arguments that should NOT have been set
		if [ ! -z $LDB ] || [ ! -z $LMAT ] || [ ! -z $LPARAM ] || [ ! -z $LSPLIT ] || [ $CORES != 1 ] || [ $ISREC != 0 ] || [ ! -z $LTAX ] || 
		[ ! -z $PDBPATH ] || [ ! -z $E ] || [ ! -z $AC ] || [ ! -z $CC ] || [ ! -z $M ]; then
			INFO="${INFO}\nFollowing parameters have no effect with --filter and are ignored:\n\n\t${ALIGNINFO}\n\t${ASSIGNINFO}\n"
		fi
	fi

# Don't use pre-calculated alignment
elif [ $NOALIGN -eq 0 ]; then
	
	
	# Perform metagenomic assignment
	if [ $ISFILTER -eq 0 ]; then
		# Check if relevant parameters for workflow were set and are not empty
		checkEssential "$Q" "$LDB" "$LTAX" "$PDBPATH" "$E" "$AC" "$CC" "$M"
		INFO="${INFO}${ASSIGNINFO}"
	else
		# Check if relevant parameters for workflow were set and are not empty
		checkEssential "$Q" "$LDB"
		
		# ...Then by checking arguments that should NOT have been set
		if [ ! -z $LTAX ] || [ ! -z $PDBPATH ] || [ ! -z $E ] || [ ! -z $AC ] || [ ! -z $CC ] || [ ! -z $M ]; then
			INFO="${INFO}\nFollowing parameters have no effect with --filter and are ignored:\n\n\t${ASSIGNINFO}\n"
		fi
	fi	
	
	# Check if environment variables are set
	if [ "$LASTAL" = "" -o "$LASTSPLIT" = "" -o "$GETMETA" = "" ]; then
		printf "\nERROR: Unset essential environment variable(s). Use -h for help\n\n"
		exit 1
	fi
	
	# Check if programs exist
	if [ ! -e "$LASTAL" ]; then
		printf "\nERROR: LASTAL not found at ${LASTAL}. Check environmental variable or use -h for help.\n\n"
		exit 1
		
	# If [-lsplit] is omitted, a LASTSPLIT environment variable is not needed and no checks for LAST-SPLIT are performed	
	elif [ ! -e "$LASTSPLIT" ]; then
		if [ ! "$LSPLIT" = "" ]; then
			printf "\nERROR: LAST-SPLIT not found at ${LASTSPLIT}. Check environmental variable or use -h for help.\n\n"
			exit 1
		fi
	elif [ ! -e "$GETMETA" ]; then
		printf "\nERROR: Assignment script not found at ${GETMETA}. Check environmental variable or use -h for help.\n\n"
		exit 1
	fi

	# I only need Krona outside filter workflow
	if [ $ISFILTER -eq 0 ]; then
		if [ "$KRONA" = "" ]; then
			printf "\nERROR: Unset essential environment variable(s). Use -h for help\n\n"
		else
			# Set path to major Krona script
			KRONA="${KRONA}/ImportText.pl"
		
			if [ ! -e "$KRONA" ]; then
				printf "\nERROR: Krona main script not found at ${KRONA}. Check environmental variable or use -h for help.\n\n"
				exit 1
			fi
		fi		
	fi	
	
	# Check if LASTAL and LAST-SPLIT version is 963 or later
	# If [-lsplit] is omitted, a LASTSPLIT environment variable is not needed and no checks for LAST-SPLIT are performed
	ALVERSION=$($LASTAL --version | grep -i "lastal" | cut -d ' ' -f2)
	if [ ! "$LSPLIT" = "" ]; then
		SPLITVERSION=$($LASTSPLIT --version | grep -i "last-split" | cut -d ' ' -f2)
		if [ $ALVERSION -lt 963 -o $SPLITVERSION -lt 963 ]; then
			printf "\nERROR: LASTAL and LAST-SPLIT version must be 963 or later.\n LASTAL is ${ALVERSION}; LAST-SPLIT is ${SPLITVERSION}.\n"
			printf "Use -h for further information\n\n"
			exit 1
		fi
	else
		if [ $ALVERSION -lt 963 ]; then
			printf "\nERROR: LASTAL version must be 963 or later.\n LASTAL is ${ALVERSION}.\n"
			printf "Use -h for further information\n\n"
			exit 1
		fi	
	fi	
fi
	

# LASTSPLIT was omitted
if [ "$LSPLIT" = "" -a $NOALIGN -eq 0 ]; then
	INFO="${INFO}LAST-SPLIT:\t\t\tomitted\n\t"
fi



# For runtime
STARTTIME=$(date +"%s")



#====================================================================================================================#
# The user already has an alignment and wants to skip LAST alignment with --no_align
#--------------------------------------------------------------------------------------------------------------------#

if [ $NOALIGN -eq 1 ]; then
	
	# Check if alignment file exists
	if [ ! -e "$MAF" ]; then
		printf "\nERROR: Chose to use an existing alignment file, but file not found. Use -h for help\n\n"
		exit 1
	fi
	
	[ -d $Q ] && { printf "\nERROR: Input file must not be a directory!\n\n"; exit 1; }
	[ -f $Q ] || { printf "\nERROR: Input file does not exist or is no file!\n\n"; exit 1; }
	
	DIR=$(dirname "${Q}")
	
	# Optionally transform fastq to fasta. Script must be in same directory as getMeta
	if [ $ISFQ -eq 1 ]; then
		FASTQ2A=$(echo $GETMETA | sed -e 's/getMeta\.pl/fastq2a\.pl/')
		FASTA=$(basename "${Q}" | perl -pne 'chomp($_); if ( $_ =~ m/.fq$/ || $_ =~ m/.fastq$/ ) { $_ =~ s/q$/a/ } else { $_.= ".fasta" }')
		FASTA="${DIR}/calc.${FASTA}"
		
		"$FASTQ2A" "$Q" 1> "$FASTA" 2>/dev/null || { printf "\nERROR: Invalid fastq query file ${Q}\n\n"; rm "$FASTA"; exit 1; }
	
		# Analyze fasta file
		Q="${FASTA}"
	else
		testFasta $Q || { printf "ERROR: Invalid fasta query file ${Q}\n"; exit 1; }
	fi
	
	
	# Perform metagenomic assignment
	if [ $ISFILTER -eq 0 ]; then
		# Print environment variables
		printf "${INFO}Accessing helpers from:\n\n\tGETMETA: ${GETMETA}\n\tKRONA: ${KRONA}\n\n"
		
		# Get taxa
		printf "\n[1/1] Starting taxonomic assignment on given alignment.\n\n"
		getTaxa
		
		# Pack unmatched reads for user
		packUnmatched
	else
		# Print environment variables
		printf "${INFO}Accessing helpers from:\n\n\tGETMETA: ${GETMETA}\n\n"
		
		printf "[1/1] Filtering reads.\n\n"
		filterReads
	fi
	
	# Runtime
	ENDTIME=$(date +"%s")
	CALCTIME=$(($ENDTIME-$STARTTIME))
	
	printf "DONE\n"
	printf 'Finished in: %02dh:%02dm:%02ds\n\n' $(($CALCTIME/3600)) $(($CALCTIME%3600/60)) $(($CALCTIME%60))
	exit 0
fi


#=================================================================================================================#
# The user has no pre-calculated alignment and wants to use LAST. No --no_align flag.
#-----------------------------------------------------------------------------------------------------------------#

if [ "$LSPLIT" = "" ]; then
	LASTSPLIT="omitted"
fi


# Print environment variables. I don't need KRONA in filter workflow.
if [ $ISFILTER -eq 0 ]; then
	printf "${INFO}${ALIGNINFO}Accessing helpers from:\n\n\tLASTAL: ${LASTAL}\n\tLASTSPLIT: ${LASTSPLIT}\n\tGETMETA: ${GETMETA}\n\tKRONA: ${KRONA}\n\n"
else
	printf "${INFO}${ALIGNINFO}Accessing helpers from:\n\n\tLASTAL: ${LASTAL}\n\tLASTSPLIT: ${LASTSPLIT}\n\tGETMETA: ${GETMETA}\n\n"
fi

# User want recursive processing of all query files in a directory
if [ $ISREC -eq 1 ]; then
	
	# Directory is valid?
	[ -d $Q ] || { printf "\nERROR: Input directory invalid\n\n"; exit 1; }
	
	DIR=$( echo $Q | sed -e 's/\/$//' )
	
	TOTAL=$( ls -l $DIR | grep -G "^-" | grep -v "temp.matrix" | wc -l)
	COUNTER=1
	
	for Q in ${DIR}/*
	do
	OUTDIR=$( echo $Q | sed -e 's/\.[^\/]*$//' )
	OUTDIR="${OUTDIR}_${COUNTER}"
		
		# Only process files
		[ ! -f "$Q" -o "$Q" = "temp.matrix" ] && continue
		
		# Skip matrix file
		[ "$(echo $Q | grep -Eo 'temp.matrix')" ] && continue
		
		printf "\n# Processing file ${COUNTER} of ${TOTAL}\n"
		COUNTER=$(( COUNTER + 1 ))
		
		# Optionally transform fastq to fasta. Script must be in same directory as getMeta
		if [ $ISFQ -eq 1 ]; then
			FASTQ2A=$(echo $GETMETA | sed -e 's/getMeta\.pl/fastq2a\.pl/')
			FASTA=$(basename "${Q}" | perl -pne 'chomp($_); if ( $_ =~ m/.fq$/ || $_ =~ m/.fastq$/ ) { $_ =~ s/q$/a/ } else { $_.= ".fasta" }')
			FASTA="${DIR}/calc.${FASTA}"
			FORMERQ="${Q}"
						
			"$FASTQ2A" "$Q" 1> "$FASTA" 2>/dev/null || { printf "WARNING: Invalid fastq query file ${Q}. Continuing with other files.\n"; rm "$FASTA"; continue; }
		
			# Analyze fasta file
			Q="${FASTA}"
		else
			testFasta $Q || { printf "WARNING: Invalid fasta query file ${Q}. Continuing with other files.\n"; continue; }
			FORMERQ="${Q}"
		fi
		
		
		# Create output directory for each query file
		mkdir $OUTDIR || { printf "\nERROR: Could not create output directory\n\n"; exit 1; }
		
		# Calculate alignment
		calcAlign
		
		# Perform metagenomic assignment
		if [ $ISFILTER -eq 0 ]; then	
			printf "[2/2] Starting taxonomic assignment.\n\n"
			getTaxa
		
			# Pack unmatched reads for user
			packUnmatched || exit 1
		else
			printf "[2/2] Filtering reads.\n\n"
			filterReads
		fi
		
		# Move all created files to respective output directory
		mv ${DIR}/calc.* $OUTDIR 
		if [ "$FORMERQ" ]; then mv $FORMERQ $OUTDIR; fi
		
	done
	
else
	[ -f $Q ] || { printf "\nERROR: Input file does not exist!\n\n"; 
	[ -d $Q ] && printf "ERROR: Directory found when file expected: Did you forget to type --recursive?\n\n"; exit 1; }
	
	DIR=$(dirname "${Q}")
	
	# Optionally transform fastq to fasta. Script must be in same directory as getMeta
	if [ $ISFQ -eq 1 ]; then
		FASTQ2A=$(echo $GETMETA | sed -e 's/getMeta\.pl/fastq2a\.pl/')
		FASTA=$(basename "${Q}" | perl -pne 'chomp($_); if ( $_ =~ m/.fq$/ || $_ =~ m/.fastq$/ ) { $_ =~ s/q$/a/ } else { $_.= ".fasta" }')
		FASTA="${DIR}/calc.${FASTA}"
		
		"$FASTQ2A" "$Q" 1> "$FASTA" 2>/dev/null || { printf "\nERROR: Invalid fastq query file ${Q}\n\n"; rm "$FASTA"; exit 1; }
	
		# Analyze fasta file
		Q="${FASTA}"
	else
		testFasta $Q || { printf "ERROR: Invalid fasta query file ${Q}\n"; exit 1; }
	fi
	
	# Calculate alignment
	calcAlign	
		
	# Perform metagenomic assignment
	if [ $ISFILTER -eq 0 ]; then
		printf "[2/2] Starting taxonomic assignment.\n\n"
		getTaxa
		
		# Pack unmatched reads for user
		packUnmatched
	else
		printf "[2/2] Filtering reads.\n\n"
		filterReads
	fi
fi


# Runtime
ENDTIME=$(date +"%s")
CALCTIME=$(($ENDTIME-$STARTTIME))

printf "DONE\n"
printf 'Finished in: %02dh:%02dm:%02ds\n\n' $(($CALCTIME/3600)) $(($CALCTIME%3600/60)) $(($CALCTIME%60))
