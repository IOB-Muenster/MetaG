#!/bin/sh


#======================================================================================#
# Setup of MetaG.
# Can install or remove MetaG and its dependencies and write a config file.
#--------------------------------------------------------------------------------------#


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


# Functions
testValue () {
	echo "$1" | grep -qG "^-[a-z]\|^--" && { printf "\nERROR: Argument found when value expected. Did you forget to give a value to the argument before [${1}]?\n\n"; exit 1; }
}

checkAbort () {
	
	if [ "$ANSWER" = "q" ] || [ "$ANSWER" = "Q" ]; then
		printf "INFO: User has aborted the setup!\n"
		exit 0
	fi
}

checkDir () {

	if [ -z "$ANSWER" ]; then
		printf "WARNING: Choose a valid directory.\n"
	elif [ ! -d $ANSWER ] && [ $ISREMOVE -eq 0 ]; then
		
		DIR=${ANSWER}
		
		printf "WARNING: Directory ${DIR} does not exist.\n"
		
		while [ 0=0 ]; do
			printf "MSG: Do you wish to create ${DIR}? (yes, no, q/Q to quit) "
			read ANSWER
			checkAbort
			checkYesNo
		done
		
		if [ $ANSWER -eq 1 ]; then
			mkdir $DIR || { printf "\nERROR: Could not create ${DIR}\n";
				printf "INFO: No software changes have been made!\n\n";
				exit 1; }
			break
		else
			printf "WARNING: ${DIR} was not created. Choose new directory.\n";
		fi
		
	elif [ ! -d $ANSWER ] && [ $ISREMOVE -eq 1 ];then
		printf "WARNING: Can't uninstall from nonexistent directory ${ANSWER}\n"
	else
		DIR=${ANSWER}
		
		if [ $ISREMOVE -eq 0 ] && [ $(ls -1qA $DIR | grep -cG ".") -gt 0 ]; then
			
			printf "WARNING: Installation directory ${DIR} is not empty.\n"
			
			FOLDERS=""
			[ -d "$DIR/tools" ] && FOLDERS="${FOLDERS}tools,"
			[ -d "$DIR/calculate" ] && FOLDERS="${FOLDERS}calculate,"
			[ -d "$DIR/krona" ] && FOLDERS="${FOLDERS}krona,"
			[ -d "$DIR/files" ] && FOLDERS="${FOLDERS}files"
			FOLDERS=$( echo $FOLDERS | sed s/,$//)
			
			# User should only be able to remove critical folders
			# calculate, tools, krona, files are used for MetaG
			# will also silently overwrite metag_setup.sh in target folder.
			if [ "$FOLDERS" != "" ]; then
				while [ 0=0 ]; do
					printf "MSG: Delete folder(s) ${FOLDERS} and all files inside? (yes, no, view, q/Q to quit): " # error
					read ANSWER
					checkAbort
					checkYesNoView
				done
				
				if [ $ANSWER -eq 1 ]; then
					[ -d "$DIR/tools" ] && rm -r $DIR/tools
					[ -d "$DIR/calculate" ] && rm -r $DIR/calculate
					[ -d "$DIR/krona" ] && rm -r $DIR/krona
					[ -d "$DIR/files" ] && rm -r $DIR/files
					printf "INFO: Deleted all files in ${FOLDERS} within ${DIR}\n"
					break
				else
					printf "WARNING: Chose non-empty directory ${DIR}\n"
					printf "WARNING: Chose not to delete files in directory.\n"
					printf "WARNING: Must provide new installation directory.\n"
				fi
			
			else
				printf "INFO: ${DIR} not empty. No critical directories, will not delete anything.\n"
				break
			fi
			
		elif [ $ISREMOVE -eq 1 ]; then
			
			# Check if directory contains MetaG files
			if [ -e "${DIR}/calculate/getMeta.pl" ]; then
				printf "INFO: Valid MetaG directory.\n"
				break
			else
				printf "WARNING: MetaG has not been installed in this directory.\n"
				printf "WARNING: Must provide new installation directory.\n"
			fi	
		else
			break
		fi
		
		
	fi
}

checkYesNoView () {
	
	if [ "$ANSWER" = "yes" ] || [ "$ANSWER" = "Yes" ] || [ "$ANSWER" = "y" ] || [ "$ANSWER" = "Y" ]; then
		ANSWER=1
		break
	elif [ "$ANSWER" = "no" ] || [ "$ANSWER" = "No" ] || [ "$ANSWER" = "n" ] || [ "$ANSWER" = "N" ]; then
		ANSWER=0
		break
	elif [ "$ANSWER" = "view" ] || [ "$ANSWER" = "View" ]; then
		if [ -d $DIR ]; then
			CONTENTS=$(ls -lha $DIR)
			printf "\n${CONTENTS}\n\n"
		else
			printf "WARNING: Choose yes or no.\n"
		fi
	else
		printf "WARNING: Choose yes or no.\n"
	fi
}

checkYesNo () {
	
	if [ "$ANSWER" = "yes" ] || [ "$ANSWER" = "Yes" ] || [ "$ANSWER" = "y" ] || [ "$ANSWER" = "Y" ]; then
		ANSWER=1
		break
	elif [ "$ANSWER" = "no" ] || [ "$ANSWER" = "No" ] || [ "$ANSWER" = "n" ] || [ "$ANSWER" = "N" ]; then
		ANSWER=0
		break
	else
		ANSWER=2
		printf "WARNING: You must choose yes or no.\n"
	fi
		
}

checkLASTdep () {
	
	[ "$(which unzip)" = "" ] && { printf "\nERROR: unzip must be installed\n"; exit 1; }
	[ "$(which cc)" = "" ] && { printf "\nERROR: cc or gcc must be installed\n"; exit 1; }
	
	if [ $OS = "FreeBSD" ] || [ $OS = "Darwin" ]; then
		[ "$(which gmake)" = "" ] && { printf "\nERROR: gmake must be installed\n"; exit 1; }
	fi
}

checkOS () {
	
	OS=`uname -s`
	
	if [ "$OS" != "FreeBSD" ] && [ "$OS" != "Linux" ] && [ "$OS" != "Darwin" ]; then
		printf "\nERROR: Installer requires FreeBSD, Linux or MacOS operating system!\n"
		printf "INFO: No software changes have been made!\n\n"
		exit 1
	fi
	
}

# Tries to auto detect path for parameter
findPATH () {
	
	# Check printenv and which with both upper and lower case
	LOWERPARAM=$(echo "$PARAM" | tr '[:upper:]' '[:lower:]')
	
	if [ "$PARAM" = "LASTSPLIT" ]; then
		ANSWER=$(printenv "LAST-SPLIT" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(printenv "last-split" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(which "LAST-SPLIT" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(which "last-split" 2> /dev/null)
	else
		ANSWER=$(printenv "$PARAM" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(printenv "$LOWERPARAM" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(which "$PARAM" 2> /dev/null)
		if [ "$ANSWER" ]; then return; fi
		ANSWER=$(which "$LOWERPARAM" 2> /dev/null)
	fi
	
}

# Get essential KRONA files for generation of interactive plots
getKrona () {
	
	printf "INFO: Installing Krona\n"
	mkdir krona
	cd krona
	
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/LICENSE.txt' && \
	mkdir scripts lib img src && \
	
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/scripts/ImportText.pl' && \
	chmod a+x ImportText.pl && \
	mv ImportText.pl scripts/ && \
	
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/lib/KronaTools.pm' && \
	mv KronaTools.pm lib/ && \
	
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/src/krona-2.0.js' && \
	mv krona-2.0.js src/ && \
	
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/img/hidden.uri' && \
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/img/loading.uri' && \
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/img/favicon.uri' && \
	wget -q 'https://raw.githubusercontent.com/marbl/Krona/master/KronaTools/img/logo-med.uri' && \
	mv *.uri img/ && \
	
	cd .. && \
	printf "INFO: Installed Krona at ${DIR}/krona\n" || \
	{ printf "ERROR: Could not install Krona\n"; exit 1; }
	
}


installAll () {
	
	[ "$(which perl)" = "" ] && { printf "\nERROR: Perl must be installed\n"; exit 1; }
	[ "$(which zip)" = "" ] && { printf "\nERROR: zip must be installed\n"; exit 1; }
	
	if [ "$(which wget)" = "" ]; then
		
		if [ "$OS" != "Darwin" ]; then
			printf "\nERROR: wget must be installed\n"
			exit 1
		else
			printf "\nERROR: wget must be installed. You may want to use Homebrew.\n"
			exit 1
		fi
	else
		if [ "$OS" = "FreeBSD" ]; then
			ISCERTS=0
			ISCERTS=$(pkg list | grep -c "ca_root")
			
			if [ $ISCERTS -eq 0 ]; then
				printf "\nERROR: ca_root_nss must be installed.\n"
				exit 1
			fi
		fi
	fi
	
	printf "INFO: Installing MetaG in: ${DIR}\n"
	
	# Check dependencies
	printf "INFO: Checking dependencies\n"
	
	INSTALL=$(dirname $0)
	
	# Get the path of calculate with standard utilities, even if user is in install folder.
	# 1) Get dirname of this script (install directory)
	# 2) Go to install directory and get full path via pwd; go back
	# 3) Replace "install" with "calculate". 
	CALC=$(echo $(cd $INSTALL && pwd && cd - > /dev/null) | sed s/install$/calculate/)
	
	if [ ! -d $CALC ]; then
		printf "ERROR: Missing directories for MetaG.\n"
		exit 1
	fi
	
	# Check if lastal and last-split are installed
	LASTAL=$( which lastal )
	LASTSPLIT=$( which last-split )
	
	if [ "$LASTAL" != "" ] && [ "$LASTSPLIT" != "" ]; then
		
		# Check versions
		ALVERSION=$(lastal --version | grep -i "lastal" | cut -d ' ' -f2)
		SPLITVERSION=$(last-split --version | grep -i "last-split" | cut -d ' ' -f2)
		
		if [ $ALVERSION -lt 963 -o $SPLITVERSION -lt 963 ]; then
			cat << EOF		
INFO: Found: LASTAL version ${ALVERSION}; LAST-SPLIT version ${SPLITVERSION}.
WARNING: LASTAL and LAST-SPLIT must be version 963 or later.
EOF
		else
			ISLAST=1
			printf "INFO: Correct version of LAST package already installed.\n"
		fi
	
	else
		printf "INFO: LAST was not found on your computer. Will be installed.\n"
	fi
	
	
	# Install LAST, if needed
	if [ $ISLAST -eq 0 ]; then
		
		checkLASTdep
		
		# Check, if LAST can be downloaded
		ping -c 2 -W 5 $GITLAB > /dev/null 2>&1 || { printf "ERROR: Could not reach ${GITLAB}. Check internet connection.\n"; exit 1; }
	
		cd $DIR
	
		printf "INFO: Installing LAST\n"
		
		mkdir "tools"
		cd "tools"
		mkdir "bin"
		
		# Download lastal, last-split
		# Get index.html to get filename
		wget -q --read-timeout 15 -t 3 "$LASTWEB" || { printf "ERROR: Having trouble reaching ${LASTWEB}.\n";
			cd ..;
			rm -r tools;
			exit 1; }
	
		URL=$(cat tags | grep 'gl-button btn btn-sm btn-confirm' | grep -Eo '\/.*\.zip' | head -n1)
		FILE=$(printf ${URL} | grep -Eo '[^/]*\.zip$')
		rm tags
		
		# Get last*zip and extract
		wget -q --read-timeout 15 -t 3 "${GITLAB}${URL}" || { printf "ERROR: Having trouble reaching ${GITLAB}${URL}.\n";
			cd ..;
			rm -r tools;
			exit 1; }
			
		unzip -q $FILE
		rm $FILE
		mv last-* last
		
		
		printf "INFO: Compiling\n"
		
		# Compile and suppress stdout and stderr of compiler
		# Save stderr temporarily. If status code 0, remove error log
		if [ $OS = "FreeBSD" ]; then
			sed -i.tmp -e '{
	s/= g++/= c++/;
	s/= gcc/= cc/
	}'  "last/makefile"
			cd last
			gmake 1> /dev/null 2> ../metag.comp.err || { printf "ERROR: Error during compiling.\n";
				printf "ERROR: Log file: ${DIR}/tools/metag.comp.err\n";
				cd ..;
				rm -r last;
				exit 1; }
			[ -f "../metag.comp.err" ] && rm ../metag.comp.err
			cd ..
		elif [ $OS = "Linux" ]; then
			cd "last"
			make 1> /dev/null 2> ../metag.comp.err || { printf "ERROR: Error during compiling.\n";
				printf "ERROR: Log file: ${DIR}/tools/metag.comp.err\n";
				cd ..;
				rm -r last;
				exit 1; }
			[ -f "../metag.comp.err" ] && rm ../metag.comp.err
			cd ..
		elif [ $OS = "Darwin" ]; then
			cd "last"
			gmake 1> /dev/null 2> ../metag.comp.err || { printf "ERROR: Error during compiling.\n";
				printf "ERROR: Log file: ${DIR}/tools/metag.comp.err\n";
				cd ..;
				rm -r last;
				exit 1; }
			[ -f "../metag.comp.err" ] && rm ../metag.comp.err
			cd ..
		fi
		
		cp -av last/bin/* bin > /dev/null
		
		rm -rf last*
		
		PWD=$(pwd)
		
		LASTAL="${PWD}/bin/lastal"
		LASTSPLIT="${PWD}/bin/last-split"
		
		printf "INFO: Installed LAST at ${PWD}/bin/\n"

		cd $INDIR
	fi
	
	# Move MetaG files
	mv $CALC/ $DIR
	mv $INSTALL/files/ $DIR
	
	
	cp $0 $DIR
	
	
	cd $DIR/
	
	# Install essential Krona files
	getKrona
	
	PWD=$(pwd)
	
	chmod a+x calculate/*
	chmod a+x krona/*
	
	printf "INFO: Installing Metag.\n"
	
	# Write environment variables to config file 
	# Append to existing standard configs which come with MetaG or make new one
	for CONF in $PWD/files/standard.config.*; do
		if [ ! -f "$CONF" ]; then
			printf "INFO: Writing environment variables to config file ${PWD}/files/standard.config.env\n"
			printf "#env LASTAL ${LASTAL}\n#env LASTSPLIT ${LASTSPLIT}\n#env KRONA ${PWD}/krona/scripts\n#env GETMETA ${PWD}/calculate/getMeta.pl" > "${PWD}/files/standard.config.env"
			break
		else
			printf "INFO: Adding environment variables to config file ${CONF}\n"
			printf "\n#env LASTAL ${LASTAL}\n#env LASTSPLIT ${LASTSPLIT}\n#env KRONA ${PWD}/krona/scripts\n#env GETMETA ${PWD}/calculate/getMeta.pl" >> "${CONF}"
		fi
	done
	
	MANINSTALL=$(manpath -q | cut -d ":" -f1)
	if [ -d "${MANINSTALL}/man1" ]; then
		printf "INFO: Attempting to install man pages to ${MANINSTALL}/man1.\n"
		printf "INFO: You need root privileges to continue.\n"
		while [ 0=0 ]; do
			printf "MSG: Continue? [yes, no]: "
			read ANSWER
			checkYesNo
		done
		
		if [ $ANSWER = 0 ]; then
			printf "WARNING: Man pages not installed. Help functions for metag and metag_setup are limited.\n"
			printf "INFO: Man files ${PWD}/files/metag_setup.1.gz and ${PWD}/files/metag.1.gz should be moved to any man1 directory of\n"
			printf "$(manpath -q | sed 's/:/ or /g')\n"
		else
			printf "INFO: Executing sudo find ${PWD}/files/ -name metag*.1.gz -exec mv '{}' ${MANINSTALL}/man1 \;\n"
			
			if [ "$(which sudo)" = "" ]; then
				printf "\nWARNING: sudo must be installed\n"
				printf "WARNING: Man pages could not be installed. Help functions for metag and metag_setup are limited.\n"
			else	
				sudo find "${PWD}/files/" -name "metag*.1.gz" -exec mv '{}' "${MANINSTALL}/man1" \; || { printf "WARNING: Error while copying man pages. Help functions for metag and metag_setup are limited.\n"; }
			fi
			
		fi
	
	else
		printf "WARNING: Man pages could not be installed. Help functions for metag and metag_setup are limited.\n"
		printf "INFO: Man files ${PWD}/files/metag_setup.1.gz and ${PWD}/files/metag.1.gz should be moved to any man1 directory of\n"
		printf "$(manpath| sed 's/:/ /g')\n"
	fi
	
	
	printf "\nSUCCESS\n\n"	
}

# Defaults
GITLAB="gitlab.com"
LASTWEB="${GITLAB}/mcfrith/last/-/tags"
INDIR=$(pwd)
ISREMOVE=0
ISCONF=0
ISLAST=0
DIR=""

HELP=$(cat <<EOF
MetaG setup(1)                                                                                  Documentation                                                                                  MetaG setup(1)

NAME
       MetaG setup - (un)installs MetaG and optionally the LAST package. Assists with writing config files for MetaG.

SYNOPSIS
       For help: metag_setup.sh -h

       metag_setup.sh [-i directory] | [-u directory] | -c [-f]

DESCRIPTION
       The setup of MetaG can be used to (un)install MetaG either interactively or by specifying the appropriate options on the command line.

       The software also assists users with creating the config files for MetaG. This is done in an interactive manner.

       Supplying no options to the program indicates the interactive mode.

OPTIONS
       -h     OPT: Show this help message and exit.

       -i     OPT:  Install MetaG to this directory [PATH]. Target directory must not contain any directories names tools, krona, calculate or data.  If these exist, -f will prevent the program from crash‐
              ing, but it will delete all contents of these folders.

       -u     OPT: Uninstall MetaG from this directory [PATH]. Removes the folders tools, krona, calculate and data and metag_setup.sh in the directory.

       -f     OPT: Continues installation, even if critical directories are present in the installation directory. See -i.

EXAMPLES
       During the installation process, the setup attempts to write the environment variables for MetaG to all config files defining the standard parameters of MetaG.  These  should  have  been  downloaded
       with MetaG. If they are not present, the variables are written to a separate config file.

       During  the  generation  of a config file in the interactive mode variables can be skipped by pressing ENTER.  The script can try auto detect the correct values for the environment parameters. These
       are looked up internally using printenv.  If this fails, the script attempts to use the which command, internally. In case both fail, you have to manually enter the correct value.

ENVIRONMENT
       The each environment variable is only used, if a config file is created and if auto detection of that environment variable is selected.

       LASTAL Location of LAST aligner [PATH].

       LAST-SPLIT
              Location of LAST-SPLIT [PATH].

       GETMETA
              Location of the getMeta.pl script of MetaG [PATH].

       KRONA  Location of Krona [PATH].

FILES
       The setup will create the following files and directories in your installation path:

       calculate
              DIR: Contains the scripts for the calculations including metag.sh.  Modules are located in the modules subdirectory.

       files  DIR: Standard config files for MetaG, named in the format: standard.config.database.sequencing technology
              It may also contain the man files metag.1.gz and metag_setup.1.gz, if those could or were not installed.

       krona  DIR: Contains the img subdirectory and other important files for the graphical representations using Krona.

       tools  DIR: Contains the LAST package in the bin subdirectory.

       metag_setup.sh
              FILE: This script.

       The above folders and files will be deleted during the uninstallation process. Potentially, conflicting files and folders will be overwritten during the installation process. See also: OPTIONS

DEPENDENCIES
       The setup was tested on Ubuntu 18.04.3 LTS, FreeBSD 12.1 and macOS 10.11.6.

       MetaG depends on LAST and LAST-SPLIT versions 963 or later to perform alignments.  If LAST and LAST-SPLIT are not installed or the versions are too old, an internet connection is mandatory to obtain
       the  LAST package. Besides, wget, unzip, GNU make and cc or gcc are needed for the installation. If GNU make is installed, but not available as gmake, you need to provide a link called gmake and add
       it to your search path.  Perl and zip are mandatory for running MetaG. See also: metag

       KronaTools 2.7 or later is mandatory for creating interactive graphs with metag.  For the installation, only wget is needed.

       If you want to install the LAST package and KronaTools manually, you can find them here:

              https://gitlab.com/mcfrith/last

              https://github.com/marbl/Krona/wiki/KronaTools

       Note that wget on FreeBSD will need ca_root_nss to work properly.

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
       ARE  DISCLAIMED.  IN  NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
       GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  TORT  (INCLUDING  NEGLIGENCE  OR
       OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

SEE ALSO
       metag(1), unzip(1), cc(1), gcc(1), gmake(1), make(1), wget(1), zip(1), Perl(1)

Felix Manske                                                                                         2020                                                                                      MetaG setup(1)
EOF
)

cat <<EOF

#================================================================================#
# SETUP for MetaG -the metagenomics pipeline
#--------------------------------------------------------------------------------#

This script will guide you through the (un)installation or configuration of MetaG.
For installation: Please make sure that you have an internet connection,
as some files may be retrieved from the internet.
	
EOF


#==========================================================================#
# Uninstall or install or write config file?
#--------------------------------------------------------------------------#

# User supplied arguments on the command line
if [ "$#" -ne 0 ]; then
	
	if [ "$#" -gt 3 ]; then
		printf "\nERROR: Too many arguments. Use -h for help.\n\n"
		exit 1
	fi
	
	ISFORCE=0
	
	while [ "$1" != "" ]; do
		case $1 in
			-i)		shift
					ISREMOVE=0
					testValue $1
					DIR=$1
					;; 
			-u)		shift
					ISREMOVE=1
					testValue $1
					DIR=$1
					;;
			-f)		ISFORCE=1
					;;
			-h)		man metag_setup 2>/dev/null || printf "${HELP}" | less
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

	checkOS
	
	if [ $ISREMOVE -eq 1 ]; then
		
		if [ ! -d "$DIR" ]; then
			printf "ERROR: Can't uninstall from nonexistent directory ${ANSWER}\n"
			exit 1
		fi
		
		if [ ! -e "${DIR}/calculate/getMeta.pl" ]; then
			printf "ERROR: MetaG has not been installed in this directory.\n"
			exit 1
		fi
		
		printf "\nINFO: Removing MetaG from ${DIR}\n"
		
		if [ -d "${DIR}/calculate" ]; then rm -r "${DIR}/calculate"; fi
		if [ -d "${DIR}/krona" ]; then rm -r "${DIR}/krona"; fi
		if [ -d "${DIR}/tools" ]; then rm -r "${DIR}/tools"; fi
		if [ -d "${DIR}/files" ]; then rm -r "${DIR}/files"; fi
		if [ -f "${DIR}/metag_setup.sh" ]; then rm "${DIR}/metag_setup.sh" || printf "WARNING: Could not remove ${DIR}/metag_setup.sh Please remove manually.\n"; fi
		
		# If available remove man pages
		MANMETAG=$(man -w metag 2>/dev/null)
		MANSETUP=$(man -w metag_setup 2>/dev/null)
		
		
		if [ $MANMETAG ]; then
			printf "INFO: Found metag man page at ${MANMETAG}\n"
			
			printf "INFO: You need root privileges to remove it.\n"
			while [ 0=0 ]; do
				printf "MSG: Continue? [yes, no]: "
				read ANSWER
				checkYesNo
			done
			
			if [ $ANSWER = 0 ]; then
				printf "INFO: Man page for metag was not removed.\n"
				printf "INFO: You can delete it manually from ${MANMETAG}.\n"
			else
				printf "INFO: Executing sudo rm ${MANMETAG}\n"
					
				if [ "$(which sudo)" = "" ]; then
					printf "\nWARNING: sudo must be installed\n"
					printf "WARNING: Man page could not be removed.\n"
					printf "INFO: You can delete it manually from ${MANMETAG}.\n"
				else
					sudo rm "${MANMETAG}" || { printf "WARNING: Man page could not be removed.\n"; }
				fi
			
			fi
			
		fi
		
		if [ $MANSETUP ]; then
			printf "INFO: Found metag_setup man page at ${MANSETUP}\n"
			
			printf "INFO: You need root privileges to remove it.\n"
			while [ 0=0 ]; do
				printf "MSG: Continue? [yes, no]: "
				read ANSWER
				checkYesNo
			done
			
			if [ $ANSWER = 0 ]; then
				printf "INFO: Man page for metag_setup was not removed.\n"
				printf "INFO: You can delete it manually from ${MANSETUP}.\n"
			else
				printf "INFO: Executing sudo rm ${MANSETUP}\n"
					
				if [ "$(which sudo)" = "" ]; then
					printf "\nWARNING: sudo must be installed\n"
					printf "WARNING: Man page could not be removed.\n"
					printf "INFO: You can delete it manually from ${MANSETUP}.\n"
				else
					sudo rm "${MANSETUP}" || { printf "WARNING: Man page could not be removed.\n"; }
				fi
			fi
			
		fi
		
		printf "\nSUCCESS\n\n"
		exit 0
	
	else
		
		if [ -d "$DIR" ] && [ $(ls -1qA "$DIR" | grep -cG ".") -gt 0 ]; then
			
			if [ $ISFORCE -eq 1 ]; then
				printf "INFO: Installation directory ${DIR} is not empty.\n"
				printf "INFO: Chose to overwrite.\n"
				[ -d "$DIR/tools" ] && rm -r $DIR/tools
				[ -d "$DIR/calculate" ] && rm -r $DIR/calculate
				[ -d "$DIR/krona" ] && rm -r $DIR/krona
				[ -d "$DIR/files" ] && rm -r $DIR/files
			else
				FOLDERS=""
				[ -d "$DIR/tools" ] && FOLDERS="${FOLDERS}tools,"
				[ -d "$DIR/calculate" ] && FOLDERS="${FOLDERS}calculate,"
				[ -d "$DIR/krona" ] && FOLDERS="${FOLDERS}krona"
				[ -d "$DIR/files" ] && FOLDERS="${FOLDERS}files"
				FOLDERS=$( echo $FOLDERS | sed s/,$//)
				
				if [ "$FOLDERS" = "" ]; then
					printf "INFO: ${DIR} not empty. No critical directories, will not delete anything.\n"			
				else
					printf "ERROR: Installation directory ${DIR} is not empty. Use [-f] to force installation.\n"
					exit 1
				fi
			fi
			
		elif [ ! -d "$DIR" ]; then
			
			mkdir "$DIR" || { printf "ERROR: Could not create installation directory ${DIR}"; exit 1; }
			printf "INFO: Created installation directory ${DIR}\n"
		fi
		
		installAll
		exit 0
	fi
		
fi



# User provided no arguments --> interactive installation

while [ 0 = 0 ]; do
	
printf "\nMSG: Choose an operation: i for install; u for uninstall; c to create config file and h for help (q/Q to abort): "
	read ANSWER
	checkAbort
	
	if [ "$ANSWER" = "u" ] || [ "$ANSWER" = "U" ]; then
		ISREMOVE=1
		printf "INFO: Uninstalling MetaG\n"
		break
	elif [ "$ANSWER" = "i" ] || [ "$ANSWER" = "I" ]; then
		ISREMOVE=0
		printf "INFO: Installing MetaG\n"
		break
	elif [ "$ANSWER" = "c" ] || [ "$ANSWER" = "C" ]; then
		ISCONF=1
		printf "INFO: Creating config file for MetaG\n"
		break
	elif [ "$ANSWER" = "h" ] || [ "$ANSWER" = "H" ]; then
		man metag_setup 2>/dev/null || printf "${HELP}" | less
		exit 0
	else
		printf "WARNING: Invalid answer!\n"
	fi
	
done

checkOS

#=============================================================================#
# Write a config file
#-----------------------------------------------------------------------------#

if [ "$ISCONF" -eq 1 ]; then
	
	ENVPARAMS="LASTAL LASTSPLIT GETMETA KRONA"
	PATHPARAMS="-q -ldb -ltax -lmat -pdbPath"
	VALUEPARAMS="-lparam -lsplit -cores -e -ac -cc -m"
	RISKYPARAMS="--no_align --fastq --recursive -vir --verbose"
	
	while [ 0=0 ]; do
	
		CONFIG=""
		MATRIX=""
		
		printf "INFO: Where indicated, type auto to try to detect the variable.\n"
		printf "INFO: Press enter to skip a variable.\n"
		
		# Environment variables
		for PARAM in $ENVPARAMS; do
			
			while [ 0=0 ]; do
				printf "MSG: Enter PATH for ${PARAM} (auto for auto detection): "
				read ANSWER
				
				if [ "$ANSWER" ]; then
					
					if [ "$ANSWER" = "auto" ] || [ "$ANSWER" = "Auto" ]; then
						findPATH
						
						if [ "$ANSWER" ]; then 
							CONFIG="${CONFIG}#env ${PARAM} ${ANSWER}\n"
							break
						else
							printf "WARNING: Could not auto detect ${PARAM}. Please fill out manually.\n"
						fi
					else
						CONFIG="${CONFIG}#env ${PARAM} ${ANSWER}\n"
						break
					fi
				else
					break
				fi
			
			done
			
		done
		
		# Variables that need a path
		for PARAM in $PATHPARAMS; do
				
				while [ 0 = 0 ]; do
					
					printf "MSG: Enter PATH for $PARAM: "
					read ANSWER
					
					if [ "$ANSWER" ]; then
						
						if [ "$PARAM" = "-lmat" ]; then
							
							if [ -s "$ANSWER" ]; then
								MATRIX=$(perl -ne 'if ( $_ !~ m/^#/ ) { $_ =~ s/\n$/_/; print $_ }' "$ANSWER")
								MATRIX=$(echo "${MATRIX}" | sed -e 's/_$//' | sed -e 's/_/\\n/g')
								break
							else
							printf "ERROR: Invalid file ${ANSWER}\n"
							fi
								
						else
							CONFIG="${CONFIG}#metag ${PARAM} ${ANSWER}\n"
							break
						fi
					
					else
						break
					fi
					
				done
	
		done
		
		# Variables that need a "regular" value
		for PARAM in $VALUEPARAMS; do
			printf "MSG: Enter Value for $PARAM: "
			read ANSWER
			
			if [ "$ANSWER" ]; then
				
				if [ "$PARAM" = "-lsplit" ]; then
					CONFIG="${CONFIG}#lastsplit -m ${ANSWER}\n"
				elif [ "$PARAM" = "-lparam" ]; then
					# Remove spaces between single LAST parameters and their values: -a 3 -A 2 => -a3 -A2
					# Transform spaces between parameters to \n
					ANSWER=$(echo "${ANSWER}" | sed -e 's/\(\-[a-zA-Z]\) \([0-9a-zA-Z]*\)/\1\2/g' | sed -e 's/ /\\n#last /g')
					ANSWER="#last ${ANSWER}\n"
					CONFIG="${CONFIG}${ANSWER}"					
				else
					CONFIG="${CONFIG}#metag ${PARAM} ${ANSWER}\n"
				fi
				
			fi
		done
		
		
		# Boolean variables that cannot be overwritten with command line values.
		# Check that user wants to continue and knows what he's doing
		printf "\nWARNING: The following analysis-specific parameters cannot be overwritten on the command line.\n\n"
		
		while [ 0=0 ]; do
			
			printf "MSG: Set these parameters anyway? (yes, no): "
			read ANSWER
			checkYesNo
		
		done
		
		if [ "$ANSWER" -eq "1" ]; then
				
			for PARAM in $RISKYPARAMS; do
				
				while [ 0=0 ]; do
					printf "MSG: Set $PARAM (yes, no): "
					read ANSWER
					
					if [ ! "$ANSWER" ]; then ANSWER="0"; break; fi
					checkYesNo
				done
				
				if [ "$ANSWER" -eq "1" ]; then
					# --no_align needs maf path
					if [ "$PARAM" = "--no_align" ];  then
						
						while [ 0=0 ]; do 
							printf "MSG: Enter alignment path. [ENTER] will skip ${PARAM}: "
							read ANSWER
							
							if [ ! "$ANSWER" ]; then
								break
							else
								CONFIG="${CONFIG}#metag ${PARAM} ${ANSWER}\n"	
								break
							fi
							
						done
						
					else
						CONFIG="${CONFIG}#metag ${PARAM} +\n"
					fi
				fi
				
			done
				
		fi	
			
		# Final check by user
		while [ 0=0 ]; do
			
			printf "MSG: Are all answers correct? Press q to quit. (yes, no, q): "
			
			read ANSWER
			checkAbort
			checkYesNo
		
		done
		
		if [ "$ANSWER" -eq "1" ]; then
			break
		fi
		
		printf "INFO: Restarting setup.\n"
		
	done
	
	
	# Write the config file
	while [ 0=0 ]; do
		
		# No values given
		if [ ! "$CONFIG" ] && [ ! "$MATRIX" ]; then
			printf "WARNING: Nothing to write. Exiting.\n"
			exit 0
		fi
		
		printf "MSG: Enter a PATH to write the config file: "
		read ANSWER
		
		if [ -e "$ANSWER" ]; then
			
			printf "ERROR: File exists, won't overwrite it.\n"
		
		else
			
			printf "${CONFIG}${MATRIX}" 2> /dev/null 1> "$ANSWER" 
			
			if [ "$?" -ne "0" ]; then
				printf "ERROR: Cannot write file. Make sure the location exists.\n"
			else
				printf "INFO: Wrote config to ${ANSWER}\n"
				break
			fi
			
		fi
	
	done
	
	printf "\nSUCCESS\n\n"


	exit 0
fi


#===========================================================================#
# Installation directory
#---------------------------------------------------------------------------#

while [ 0 = 0 ]; do
	
	printf "MSG: Enter the installation directory of MetaG (q/Q to abort): "
	read ANSWER
	checkAbort
	checkDir
	
done




#============================================================================#
# Uninstall MetaG LAST
#----------------------------------------------------------------------------#

if [ $ISREMOVE -eq 1 ]; then
	
	while [ 0=0 ]; do
		printf "WARNING: Uninstalling MetaG will lead to loss of all data in calculate, krona, tools and files within ${DIR}\n"
		printf "MSG: Continue? (yes, no, view) "
		read ANSWER
		checkYesNoView
	done
	
	if [ $ANSWER -eq 0 ]; then
		printf "\nINFO: USER aborted uninstallation.\n"
		printf "INFO: No software changes were made.\n\n"
		exit 0
	elif [ $ANSWER -eq 1 ]; then
		if [ -d "${DIR}/calculate" ]; then rm -r "${DIR}/calculate"; fi
		if [ -d "${DIR}/krona" ]; then rm -r "${DIR}/krona"; fi
		if [ -d "${DIR}/tools" ]; then rm -r "${DIR}/tools"; fi
		if [ -d "${DIR}/files" ]; then rm -r "${DIR}/files"; fi
		if [ -f "${DIR}/metag_setup.sh" ]; then rm "${DIR}/metag_setup.sh" || printf "WARNING: Could not remove ${DIR}/metag_setup.sh Please remove manually.\n"; fi
		
		# If available remove man pages
		MANMETAG=$(man -w metag 2>/dev/null)
		MANSETUP=$(man -w metag_setup 2>/dev/null)
		
		
		if [ $MANMETAG ]; then
			printf "INFO: Found metag man page at ${MANMETAG}\n"
			
			printf "INFO: You need root privileges to remove it.\n"
			while [ 0=0 ]; do
				printf "MSG: Continue? [yes, no]: "
				read ANSWER
				checkYesNo
			done
			
			if [ $ANSWER = 0 ]; then
				printf "INFO: Man page for metag was not removed.\n"
				printf "INFO: You can delete it manually from ${MANMETAG}.\n"
			else
				printf "INFO: Executing sudo rm ${MANMETAG}\n"
				sudo rm "${MANMETAG}" || { printf "WARNING: Man page could not be removed.\n"; }
			fi
			
		fi
		
		if [ $MANSETUP ]; then
			printf "INFO: Found metag_setup man page at ${MANSETUP}\n"
			
			printf "INFO: You need root privileges to remove it.\n"
			while [ 0=0 ]; do
				printf "MSG: Continue? [yes, no]: "
				read ANSWER
				checkYesNo
			done
			
			if [ $ANSWER = 0 ]; then
				printf "INFO: Man page for metag_setup was not removed.\n"
				printf "INFO: You can delete it manually from ${MANSETUP}.\n"
			else
				printf "INFO: Executing sudo rm ${MANSETUP}\n"
				sudo rm "${MANSETUP}" || { printf "WARNING: Man page could not be removed.\n"; }
			fi
			
		fi
		
		printf "\nSUCCESS\n\n"
	fi
	
	exit 0
fi


#============================================================================#
# Installation of dependencies and MetaG
#----------------------------------------------------------------------------#

installAll
