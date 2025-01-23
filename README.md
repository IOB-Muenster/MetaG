# MetaG -the metagenomics pipeline
## What is MetaG?
MetaG is a pipeline for the analysis of marker genes (16S, 28S rRNA) and for the
analysis of whole genomes of viruses. It can handle the distinct sequencing
error profiles of short-read and long-read sequencing data.
MetaG has been shown to perform better than state-of-the-art classifiers,
especially at lower taxonomic ranks.

Apart from its outstanding performance, MetaG was designed for easy usage.
This is realized by an installation/configuration via an interactive installer, as
well as a thorough documentation which is also available in man page format.

If you want to learn more about MetaG, read its [publication](https://doi.org/10.1101/2020.03.13.991190).
The program was tested on Ubuntu 18.04.3 LTS, macOS 10.15.3 and FreeBSD 12.1.

Free yourself from hardware restrictions: Use [MetaG online](http://www.bioinformatics.uni-muenster.de/tools/metag)!

**IMPORTANT: We tested MetaG with LAST version <=1256 . In some later versions, the command line syntax of LAST has changed which will lead to unexpected behaviour of MetaG.**

## Major updates
### 14.06.2024 - Description of parameter training
The description, how parameters were trained for nanopore data, is available in the Supplementary Methods 4.1 to 4.2
of this [manuscript](https://doi.org/10.1016/j.csbj.2024.12.031).

### 08.01.2024 - New training
We completely rewrote the training scripts, used new training data, and updated the standards for RDP, MTX, RefSeq, and
the T2T filter database.

### 04.08.2021 - Filter workflow
We implemented a new workflow (--filter) which allows you to remove reads from a specific organism from the sample. Why do you
need this? Let's have a look at two samples: An environmental sample and a patient sample. In the environmental sample from a pond
you may find a lot of reads matching algae, but only few matching your taxon of interest. To make the analysis and display more
straightforward, you may wish to exclude the algae reads first, before analyzing the remaining reads for your taxon of interest.
In case of a virome analysis from a human patient sample, you may wish to exclude human reads first, to avoid any false positive
assignments of human reads.
Please see the help or man page for details.

## Directory and files
### metag_src
Download this folder to install MetaG locally. See [README](../../blob/master/metag_src/install/README) in this directory
to learn how to get started.

### supplemental
Supplemental files and scripts. See [README](../../blob/master/supplemental/README.md) in this directory for further details.

## Authors
* Felix Manske, felix.manske@uni-muenster.de
* Marc-Nicolas Bendisch
* Norbert Grundmann
* Author for correspondence: Wojciech Makalowski, wojmak@uni-muenster.de


## COPYRIGHT
Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO,  PROCUREMENT  OF  SUBSTITUTE GOODS  OR  SERVICES;
LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

