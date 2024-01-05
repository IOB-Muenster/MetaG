#!/usr/bin/env bash

# AUTHORS

# Felix Manske, felix.manske@uni-muenster.de
# Marc-Nicolas Bendisch
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

#===================================================================================#
# Train MetaG on a database and sample
#-----------------------------------------------------------------------------------#
# 
# DESCRIPTION
#
# Given a query file, an expected classification and a database, this
# script will find the best LAST alignment and MetaG parameters using the
# Bray-Curtis similarity index or the Rk (a multiclass generalization of the 
# Matthews correlation coefficient). Calculations are performed for a given
# rank and a given parameter range of the MetaG parameters -e, -ac and -cc.
#
#
# USAGE
#
# metag-train.sh -q  QUERY READS -exp EXPECTED CLASSIFICATIONS -rank RANK TO ANALYZE \ 
#	-metric bc | rk -ldb DATABASE -ltax DATABASE TAXONOMY \
#	[-lsplit FILTERING WITH LAST-SPLIT] [-cores NUMBER OF PARALLEL PROCESSES] \
#	-e RANGE FOR E-VALUE CUTOFF -ac RANGE FOR ALIGNMENT SCORE CUTOFF \
#	-cc RANGE FOR CONFIDENCE CUTOFF [--verbose] [--fastq] [-h]
#
# The following environment variables need to be set (export VAR=PATH in Bash):
#
#	LASTAL				Path of LASTAL executable.
#	LASTSPLIT			Path of LASTSPLIT executable.
#	LASTTRAIN			Path of LASTTRAIN executable.
#	GETMETA				Path of getMeta.pl script from MetaG.
#	METAG				Path of metag.sh script from MetaG.
#	KRONA				Path to ImportText.pl script from KronaTools.
#	
# More information on the structure of the EXPECTED CLASSIFICATIONS is provided
# in metag-train.py.
#
#
# OUTPUT
#
# In the directory, of the query file, this script will create an output directory
# named after the values for -e, -ac, -cc and -lsplit. The following output files
# will be created:
# 	calc.QUERY READS	If query reads were in FASTQ format, this file will
# 						contain the FASTA reads. Named after the file with the
# 						query reads.
# 	tmp.patric.tax 		Mock pathogen database.
#	tmp.lasttrain   	Raw output from LAST-TRAIN.
#	calc.bestalign.mat  Best alignment parameters and matrix from tmp.lasttrain.
#	calc.lastalign.maf  Alignment (optionally with LAST-SPLIT) from query to database
#						using calc.bestalign.mat
#	metag_log			Concatenated output of MetaG messages. This is not thread safe.
#	calc.METRIC			Benchmarking results for all parameter combinations. Named after
#						the benchmarking metric.
#	calc.best_*.config  Config file for MetaG with the best parameter combination from
#						the current run. Includes the environment variables used for
#						this run. If a viral analysis was conducted, the line
#						"#metag -vir +" needs to be added.
#	README				Best value for the benchmarking metric.
#	tmp/
#        tmp.run 		Read level classification results by MetaG over all
#        				parameter combinations
#
#
# IMPORTANT
# 
# For optimal results, multiple training routines should be conducted with varying values
# of LAST-SPLIT.
#
#
# DEPENDENCIES
#
# GNU Parallel
# KronaTools 2.7
# metag-train.py
# Perl 5.22 or more recent
# Python 3 or more recent
#
#
# CITE
#
# Internally, this script uses GNU Parallel. If you use this script,
# please cite the GNU Parallel publication:
# 
# 	Tange, O. GNU Parallel 2018. (Ole Tange, 2018).
# 	doi:10.5281/zenodo.1146014.
#
#===================================================================================#

Q=""
EXP=""
RANK=""
METRIC=""
LDB=""
LTAX=""
CORES=1
PDBPATH=""
LSPLIT=""
ERANGE=""
ACRANGE=""
CCRANGE=""

ISFQ=0
ARGCOUNT=0

HELPMSG='
#===================================================================================#
 Train MetaG on a database and sample
#-----------------------------------------------------------------------------------#
 
 DESCRIPTION

 Given a query file, an expected classification and a database, this
 script will find the best LAST alignment and MetaG parameters using the
 Bray-Curtis similarity index or the Rk (a multiclass generalization of the 
 Matthews correlation coefficient). Calculations are performed for a given
 rank and a given parameter range of the MetaG parameters -e, -ac and -cc.


 USAGE

 metag-train.sh -q  QUERY READS -exp EXPECTED CLASSIFICATIONS -rank RANK TO ANALYZE \ 
	-metric bc | rk -ldb DATABASE -ltax DATABASE TAXONOMY \
	[-lsplit FILTERING WITH LAST-SPLIT] [-cores NUMBER OF PARALLEL PROCESSES] \
	-e RANGE FOR E-VALUE CUTOFF -ac RANGE FOR ALIGNMENT SCORE CUTOFF \
	-cc RANGE FOR CONFIDENCE CUTOFF [--verbose] [--fastq] [-h]

 The following environment variables need to be set (export VAR=PATH in Bash):

	LASTAL				Path of LASTAL executable.
	LASTSPLIT			Path of LASTSPLIT executable.
	LASTTRAIN			Path of LASTTRAIN executable.
	GETMETA				Path of getMeta.pl script from MetaG.
	METAG				Path of metag.sh script from MetaG.
	KRONA				Path to ImportText.pl script from KronaTools.
	
 More information on the structure of the EXPECTED CLASSIFICATIONS is provided
 in metag-train.py.


 OUTPUT

 In the directory, of the query file, this script will create an output directory
 named after the values for -e, -ac, -cc and -lsplit. The following output files
 will be created:
 	calc.QUERY READS	If query reads were in FASTQ format, this file will
 						contain the FASTA reads. Named after the file with the
 						query reads.
 	tmp.patric.tax 		Mock pathogen database.
	tmp.lasttrain   	Raw output from LAST-TRAIN.
	calc.bestalign.mat  Best alignment parameters and matrix from tmp.lasttrain.
	calc.lastalign.maf  Alignment (optionally with LAST-SPLIT) from query to database
						using calc.bestalign.mat
	metag_log			Concatenated output of MetaG messages. This is not thread safe.
	calc.METRIC			Benchmarking results for all parameter combinations. Named after
						the benchmarking metric.
	calc.best_*.config  Config file for MetaG with the best parameter combination from
						the current run. Includes the environment variables used for
						this run. If a viral analysis was conducted, the line
						"#metag -vir +" needs to be added.
	README				Best value for the benchmarking metric.
	tmp/
        tmp.run 		Read level classification results by MetaG over all
        				parameter combinations


 IMPORTANT
 
 For optimal results, multiple training routines should be conducted with varying values
 of LAST-SPLIT.


 DEPENDENCIES

 GNU Parallel
 KronaTools 2.7
 metag-train.py
 Perl 5.22 or more recent
 Python 3 or more recent


 CITE

 Internally, this script uses GNU Parallel. If you use this script,
 please cite the GNU Parallel publication:
 
 	Tange, O. GNU Parallel 2018. (Ole Tange, 2018).
 	doi:10.5281/zenodo.1146014.
'


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


# Run MetaG and analyze
trainMetaG () {
	TMP="e$1_ac$2_cc$3"
	mkdir ${TMP}
	cd ${TMP}
	
	ln -s $4 "input.query"
	${METAG} -q input.query -ltax $5 -pdbPath "../../tmp.patric.tax" -e $1 -ac $2 -cc $3 -m "williams" \
		--no_align "../../calc.lastalign.maf" >> "../../metag_log" 2>&1  || \
		{ printf "ERROR: MetaG returned an exception\n"; exit 1; }
	
	# Save the identifications for each read.
	# Use a header line which uniquely identifies
	# Each parameter combination: @@e_ac_cc	
	cat <(printf "@@_$1_$2_$3\n") "calc.LIN.txt" > "tmp.res"
	cat "tmp.res" >&1

	${SCRIPTPATH}/metag-train.py --query input.query --exp $6 --train "tmp.res" --rank $7 --metric $8 --name ${TMP} >&2 || \
        { printf "ERROR: Could not analyze benchmarking output tmp.benchmark\n"; exit 1; }

	cd ..
	rm -r ${TMP}
}


# I don't need the -vir flag from MetaG, since this affects only the pathogen analysis.
# During training, I am not interested in that.
while [ "$1" != "" ]; do
		case $1 in
			-q)		shift
					testValue $1
					countOpt "$Q" "2"
					Q=$(realpath $1)
					;;
			-exp)	shift
					testValue $1
					countOpt "$EXP" "2"
					EXP=$(realpath $1)
					;;
			-rank)	shift
					testValue $1
					countOpt "$RANK" "2"
					RANK=$(printf ${1} | tr '[:upper:]' '[:lower:]')
					;;
			-metric)	shift
					testValue $1
					countOpt "$METRIC" "2"
					METRIC=$(printf ${1} | tr '[:upper:]' '[:lower:]')
					;;
			--fastq)
					countSwitch "$ISFQ" "1"
					ISFQ=1
					;;
			-ldb)	shift
					testValue $1
					countOpt "$LDB" "2"
					LDB=$(realpath $1)
	
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
					LTAX=$(realpath $1)		
					;;
			-cores)		shift
					testValue $1
					countOpt "$CORES" "2"
					CORES=$1
					;;
			-lsplit)	shift
					testValue $1
					countOpt "$LSPLIT" "2"
					LSPLIT=$1
					;;
			-e)		shift
					testValue $1
					countOpt "$ERANGE" "2"
					ERANGE=$1
					;;
			-ac)	shift
					testValue $1
					countOpt "$ACRANGE" "2"
					ACRANGE=$1
					;;
			-cc)	shift
					testValue $1
					countOpt "$CCRANGE" "2"
					CCRANGE=$1
					;;
			-h)		printf "${HELPMSG}" | less
					exit 0
					;;
			*)		printf "\nERROR: Unknown parameter ${1}. Use -h for help.\n\n"
					exit 1
					;;
		esac
		shift
		
done


#------------------------------------------------------------------#
# Get and check environment variables
#------------------------------------------------------------------#
LASTAL=$(printenv LASTAL)
LASTSPLIT=$(printenv LASTSPLIT)
LASTTRAIN=$(printenv LASTTRAIN)
METAG=$(printenv METAG)
GETMETA=$(printenv GETMETA)
KRONA=$(printenv KRONA)

if [ -z ${LASTAL} ]; then
	printf "ERROR: Empty environment variable LASTAL.\n"
	exit 1
elif [ -z "$(${LASTAL} --help)" ]; then
	printf "ERROR: LASTAL does not exist at ${LASTAL} .\n"
	exit 1
fi

if [ -z ${LASTSPLIT} ]; then
	printf "ERROR: Empty environment variable LASTSPLIT.\n"
	exit 1
elif [ -z "$(${LASTSPLIT} --help)" ]; then
	printf "ERROR: LAST-SPLIT does not exist at ${LASTSPLIT} .\n"
	exit
fi

if [ -z ${LASTTRAIN} ]; then
	printf "ERROR: Empty environment variable LASTTRAIN.\n"
	exit 1
elif [ -z "$(${LASTTRAIN} --help)" ]; then
	printf "ERROR: LAST-TRAIN does not exist at ${LASTTRAIN} .\n"
	exit 1
fi

if [ -z ${METAG} ]; then
	printf "ERROR: Empty environment variable METAG.\n"
	exit 1
elif [ -z "$(${METAG} --help)" ]; then
	printf "ERROR: MetaG does not exist at ${METAG} .\n"
	exit 1
fi

if [ -z ${GETMETA} ]; then
	printf "ERROR: Empty environment variable GETMETA.\n"
	exit 1
elif [ ! -s ${GETMETA} ]; then
	printf "ERROR: getMeta.pl does not exist at ${GETMETA} .\n"
	exit 1
fi

if [ -z ${KRONA} ]; then
	printf "ERROR: Empty environment variable KRONA.\n"
	exit 1
elif [ ! -s "${KRONA}/ImportText.pl" ]; then
	printf "ERROR: ImportText.pl from Krona not found in ${KRONA} .\n"
	exit 1
fi	

# I need GNU Parallel in PATH.
# Just looking for parallel might be ambiguous.
if [ -z "$(parallel --version | grep -i 'GNU Parallel')" ]; then
	printf "ERROR: GNU Parallel not found in PATH.\n"
	exit 1
fi

SCRIPTPATH=$(realpath $(dirname $0))
export SCRIPTPATH


#------------------------------------------------------------------#
# Check parameters
#------------------------------------------------------------------#
if [ -z ${Q} ]; then
	printf "ERROR: Query file cannot be empty.\n"
	exit 1
elif [ ! -s ${Q} ]; then
	printf "ERROR: Query file does not exist or is empty.\n"
	exit 1
fi

if [ -z ${EXP} ]; then
	printf "ERROR: File with expected taxonomy cannot be empty.\n"
	exit 1
elif [ ! -s ${EXP} ]; then
	printf "ERROR: File with expected taxonomy does not exist or is empty.\n"
	exit 1
fi

if [ -z ${LTAX} ]; then
	printf "ERROR: Taxonomy file from database cannot be empty.\n"
	exit 1
elif [ ! -s ${LTAX} ]; then
	printf "ERROR: Taxonomy file from database does not exist or is empty.\n"
	exit 1
fi

if [ -z ${RANK} ]; then
	printf "ERROR: Rank cannot be empty.\n"
	exit 1
fi
	
if [ -z ${METRIC} ]; then
	printf "ERROR: Benchmarking metric cannot be empty.\n"
	exit 1
elif [ ${METRIC} != "rk" -a ${METRIC} != "bc" ]; then
	printf "ERROR: Unrecognized metric ${METRIC}.\n"
	exit 1
fi

if [ -z ${LDB} ]; then
	printf "ERROR: LAST database cannot be empty.\n"
	exit 1
fi

if [ -z ${ERANGE} ]; then
	printf "ERROR: E-value range is not set.\n"
	exit 1
elif [ -z "$(printf \"${ERANGE}\" | grep '_')" ]; then
	printf "ERROR: E-value must be a range LOWER_UPPER\n"
	exit 1
fi

if [ -z ${ACRANGE} ]; then
	printf "ERROR: Alignment score range is not set.\n"
	exit 1
elif [ -z "$(printf \"${ACRANGE}\" | grep '_')" ]; then
	printf "ERROR: Alignment score must be a range LOWER_UPPER\n"
	exit 1
fi

if [ -z ${CCRANGE} ]; then
printf "ERROR: Confidence cutoff range is not set.\n"
	exit 1
elif [ -z "$(printf \"${CCRANGE}\" | grep '_')" ]; then
printf "ERROR: Confidence cutoff must be a range LOWER_UPPER\n"
	exit 1
fi


#------------------------------------------------------------------#
# Prepare the run
#------------------------------------------------------------------#
# Name to identify the run later
RUNNAME=""
if [ -z ${LSPLIT} ]; then
	RUNNAME="E_${ERANGE}_AC_${ACRANGE}_CC_${CCRANGE}"
else
	RUNNAME="E_${ERANGE}_AC_${ACRANGE}_CC_${CCRANGE}_lsplit_${LSPLIT}"
fi

OUTPATH="$(dirname ${Q})/train_${RUNNAME}"

if [ -d ${OUTPATH} ]; then
	printf "ERROR: Training directory ${OUTPATH} exists. Won't override.\n"
	exit 1
else
	mkdir ${OUTPATH} || { printf "ERROR: Cannot create training directory ${OUTPATH}\n"; exit 1; }
	
	if [ $ISFQ -eq 1 ]; then
		FASTQ2A=$(echo $GETMETA | sed -e 's/getMeta\.pl/fastq2a\.pl/')
		FASTA=$(basename "${Q}" | perl -pne 'chomp($_); if ( $_ =~ m/.fq$/ || $_ =~ m/.fastq$/ ) { $_ =~ s/q$/a/ } else { $_.= ".fasta" }')
		FASTA="${OUTPATH}/calc.${FASTA}"
		
		"$FASTQ2A" "$Q" 1> "$FASTA" 2>/dev/null || { printf "\nERROR: Invalid fastq query file ${Q}\n\n"; rm "$FASTA"; exit 1; }
	
		# Analyze fasta file
		Q="${FASTA}"
	else
		testFasta $Q || { printf "ERROR: Invalid fasta query file ${Q}\n"; exit 1; }
	fi
	
	cd ${OUTPATH}
	mkdir "tmp"
	
	# Create a mock pathogen file to avoid crash of MetaG. Pathogen detection not part of training.
	printf "#PATRICgenomeID;lineage;#host;#resistance\n" > tmp.patric.tax
	printf ">1234;0;0;0;0;0;0;0;0;0;0;#Human;#0\n" >> tmp.patric.tax
fi


#------------------------------------------------------------------#
# Use LAST-TRAIN to get best alignment parameters
#------------------------------------------------------------------#
printf "INFO: Step [1/2]: Getting optimal alignment parameters.\n"

${LASTTRAIN} -P ${CORES} ${LDB} ${Q} > "tmp.lasttrain" || { printf "ERROR: LAST-TRAIN returned an exception.\n"; exit 1; }

if [ -s "tmp.lasttrain" ]; then
	perl -e 'my $isPrint=0;
		while(<>) {
			if ($isPrint == 1) {
				print $_
			} 
			else {
				if ($_ =~ m/^#last/){
					$isPrint = 1;
					print $_
				}
			}
		}' "tmp.lasttrain" >  "calc.bestalign.mat"
else
	printf "ERROR: LAST-TRAIN did not produce any output.\n"
	exit 1
fi


#------------------------------------------------------------------#
# Test all combinations of -e, -ac and -cc in specified range and
# calculate performance metric.
# Perform alignment outside of MetaG and use --no_align. -> Faster.
#------------------------------------------------------------------#
printf "INFO: Step [2/2]: Getting optimal MetaG parameters.\n"

if [ -z "${LSPLIT}" ]; then
	${LASTAL} -P ${CORES} -p "calc.bestalign.mat" ${LDB} ${Q} > "calc.lastalign.maf" || \
		{ printf "ERROR: LASTAL returned an exception.\n"; exit 1; }
else
	${LASTAL} -P ${CORES} -p "calc.bestalign.mat" ${LDB} ${Q} | ${LASTSPLIT} -m ${LSPLIT} -n -f MAF > "calc.lastalign.maf" || \
		{ printf "ERROR: LASTAL returned an exception.\n"; exit 1; }
fi

cd "tmp"

ELOWER=$(printf -- "${ERANGE}" | cut -d "_" -f1)
EUPPER=$(printf -- "${ERANGE}" | cut -d "_" -f2)
ACLOWER=$(printf -- "${ACRANGE}" | cut -d "_" -f1)
ACUPPER=$(printf -- "${ACRANGE}" | cut -d "_" -f2)
CCLOWER=$(printf -- "${CCRANGE}" | cut -d "_" -f1)
CCUPPER=$(printf -- "${CCRANGE}" | cut -d "_" -f2)

E=$(seq ${ELOWER} 1 ${EUPPER})
AC=$(seq ${ACLOWER} 0.1 ${ACUPPER})
CC=$(seq ${CCLOWER} 0.1 ${CCUPPER})

# Make function visible to GNU Parallel.
# Run jobs in parallel, but preserve output order.
# Crash, if one of the function calls fails.
export -f trainMetaG
parallel -j ${CORES} --keep-order --halt soon,fail=1  trainMetaG ::: "${E[@]}" ::: "${AC[@]}" ::: "${CC[@]}" ::: "${Q}" ::: "${LTAX}" ::: \
	"${EXP}" ::: "${RANK}" ::: "${METRIC}" 1>> "tmp.run" 2>>"tmp.${METRIC}" || { printf "ERROR: MetaG training returned an exception\n"; exit 1; }

sort -k2,2nr -k1,1Vr "tmp.${METRIC}" > "../calc.${METRIC}" && rm "tmp.${METRIC}"


#------------------------------------------------------------------#
# Print config with best results
#------------------------------------------------------------------#
BEST=$(head -n1 ../calc.${METRIC})
BESTMETRIC=$(printf -- "${BEST}" | cut -d $'\t' -f2)
BESTE=$(printf -- "${BEST}" | cut -d $'\t' -f1 | cut -d '_' -f1)
BESTAC=$(printf -- "${BEST}" | cut -d $'\t' -f1 | cut -d '_' -f2)
BESTCC=$(printf -- "${BEST}" | cut -d $'\t' -f1 | cut -d '_' -f3)
CONFIGF="../calc.best_${METRIC}.config"

printf "Best ${METRIC}: ${BESTMETRIC}\n" > "../README"

# Environment variables + best MetaG
printf "
#env LASTAL ${LASTAL}
#env LASTSPLIT ${LASTSPLIT}
#env GETMETA ${GETMETA}
#env KRONA ${KRONA}
#metag -e ${BESTE}
#metag -ac ${BESTAC}
#metag -cc ${BESTCC}
#metag -m williams
" > ${CONFIGF}

# LAST-SPLIT, if used
if [ ! -z ${LSPLIT} ]; then
	printf "#lastsplit -m ${LSPLIT}\n" >> ${CONFIGF}
fi

# Best LAST params and matrix
grep -iv "score matrix" "../calc.bestalign.mat" >> ${CONFIGF}

printf "INFO: Best ${METRIC}: ${BESTMETRIC} for this analysis.\n"
printf "INFO: MetaG config with best values written to $(realpath ${CONFIGF}).\n"
printf "INFO: If you used a viral database, add this line to the config file: #metag -vir +\n"
printf "INFO: Finding the overall optimum, requires several trainings with different values for LAST-SPLIT\n"
printf "\nDONE\n"
