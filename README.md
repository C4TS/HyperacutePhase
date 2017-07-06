# HyperacutePhase

All raw and normalized data can be found in E-MTAB-5882_Data.zip. This folder includes the following files:

* ChaussabelModules.csv
This file maps Chaussabel modules to Illumina probes and nuIDs.

* Clinical.csv
This is the clinical data file with information on subject, time, batch, group, MODS, infection, ISS, sex, and age for every sample in the data.

* Sample_probe_Profile_5_12_13.txt
* Control_probe_Profile_5_12_13.txt
These files contain the raw expression data and can be turned into a single EListRaw object in limma by a call to read.ilmn().

* NormData.csv
This is a gene level normalized expression matrix.

* ProbeMap.csv
This file maps array address to Illumina IDs and nuIDs.

Analysis scripts are in the Scripts directory, with files covering probe and module level analyses, in addition to the figures from our paper.

The Results directory contains summary output from our module analysis.
