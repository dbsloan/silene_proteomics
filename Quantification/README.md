# Quantification

LC-MS/MS quantification data.

- The [ion_intensity](ion_intensity) directory contain ProteomeDiscoverer output that quantifies peak abundance for each sample. The "cleaned" versions were processed with the script: [clean_shared_peptides_and_peaks_found.pl](clean_shared_peptides_and_peaks_found.pl)

- The [PSMs](PSMs) directory contains ProteomeDiscoverer output with PSM counts by peptide for each sample. The [sum_psms.R](sum_psms.R) script is used to process PSM files in a directory and then write the psms.csv files for each species. *Note: the file name labeling for PSM data is misleading. The "nuclear" files really contain data from total leaf samples. They are not specific to the nucleus or nuclear-encoded proteins.*

