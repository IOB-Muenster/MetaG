# supplemental
## files
### db
Databases used in the MetaG paper.
### query
Samples used in the MetaG paper.

## scripts
These scripts were used to generate the databases for MetaG and to train
the program. You can find information on how to use the scripts in the
[publication](https://doi.org/10.1101/2020.03.13.991190) (especially in the Supplemental Materials).
The training scripts have received a major update, thus the description in the
publication is not up-to-date. Please run metag-train.sh -h for further information.
### db
Scripts used to generate the database files needed for MetaG.
Please note that we manually currated the database files in /files/db .
However, this affects only a minor portion of entries. 
### train
Scripts to obtain the optimal parameters for MetaG either using the
known abundance of taxa or the known origin of each read.
