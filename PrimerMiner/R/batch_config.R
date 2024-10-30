batch_config <- function(file){

cat(file=file, append=F)


sys <- Sys.info()[['sysname']]

if(sys=="Darwin"){sys <- "MacOSX"}
if(sys=="Windows"){
	stop("WARNING: If you runing this package on a Windows PC, clustering will NOT work! Vsearch is used for OTU clustering and is not available for windows. Please run this package on a Mac or Linuxs based system. Sorry for the inconvinience :( - comand execution was stopped")
}


# config text
text_conf <- c("# Configuration file for batch download",
paste("Version =", packageVersion("PrimerMiner"), "# you might need to regenerate this config file with earlyer or later versions of PrimerMiner!"),
"",
"# General settings",
"Taxon_table = \"path/to/csv/name_of_table.csv\" # csv containing taxa to download",
"Taxon_sep = \",\" # table entries seperated by , ",
"Marker = c(\"COi\", \"CO1\", \"COXi\", \"COX1\") # specify target gene",
"Download = T # if FALSE, no data will be downloaded ",
"Merge_and_Cluster_data = T # if set to FALSE, sequences are not merged / clustered",
"Skip_if_complete = T # if the data for a group was completely downloaded + clustered PrimerMiner will not download / cluster the data in the respecitve folder again. Set \"Skip_if_complete\" and \"Download\" to \"F\" if you like all data to be reclustered but not downloaded again",
"",
"# NCBI download, see ?Download_GB for details",
"download_GB = T",
"merge_GB = T",
"maxlength_GB = 2000",
"custom_query_GB = NULL",
"clipping_left_GB = 0",
"clipping_rigth_GB = 0",
"GB_Subset=NULL # Enter the maximum number of sequences to download to reduce download time",
"",
"# Mitochondria download, see ?Download_mito and ?Mito_GB2fasta for details",
"download_mt = T",
"merge_mt = T",
"minlength_mt = 2001",
"maxlength_mt = 80000",
"custom_query_mt = NULL",
"clipping_left_mt = 0",
"clipping_rigth_mt = 0",
"add_mt = 100",
"rm_dup = T",
"no_marker = T",
"mt_subset = NULL",
"",
"# BOLD download, see ?Download_BOLD for details",
"download_bold = T",
"save_bold_tsv = T # Indicate if bold table should also be saved",
"apikey_bold = \"00000000-0000-0000-0000-000000000000\" # login to your bold account, see user settings",
"marker_bold = \"COI-5P\" # COI marker code, other bold marker codes can be used (see downloaded bold tsv for examples)",
"merge_bold = T",
"clipping_left_bold = 0",
"clipping_rigth_bold = 0",
"subset_bold=NULL # Enter the maximum number of sequences to download to reduce download time. Obtained sequences might be less than specified if e.g. having a different marker code",
"",
"# Clustering sequences, see ?Clustering for details",
paste("operating_system= \"", sys, "\" # autodetected, can be \"MacOSX\" or \"Linux\"", sep=""), 
"vsearchpath = \"integrated\" # uses the Vsearch that comes with the R package. Use \"Vsearch\" for your local Vsearch install, or paste the path to the executable.",
"id = 0.97",
"cmd = NULL",
"threshold =  \"Majority\"",
"",
"# Write summary statistics",
"summstats = T")



cat(text_conf, file=file, append=T, sep="\n")

}







