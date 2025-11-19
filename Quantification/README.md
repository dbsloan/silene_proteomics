# Quantification

LC-MS/MS quantification data.

- The abundance.csv files contain ProteomeDiscoverer output (reformatted as csv) that quantifies peak abundance for each sample.

- The [PSMs](PSMs) directory contains ProteomeDiscoverer output with PSM counts by peptide for each sample. The [sum_psms.R](sum_psms.R) script is used to process PSM files in a directory and then write the psms.csv files for each species. *Note: the file name labeling for PSM data is misleading. The "nuclear" files really contain data from total lead samples. They are not specific to the nucleus or nuclear-encoded proteins.*

