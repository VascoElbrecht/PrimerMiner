## ---- eval = F-----------------------------------------------------------
#  setwd("your/path/to/PrimerMiner/Sample_Data")
#  
#  evaluate_primer("primer_scoring/01_Plecoptera_subset.fasta", "GGTCAACAAATCATAAAGATATTGG", 1, 25, save="save_evaluation_table_LCO.csv", mm_position ="primer_scoring/Position_v1.csv", adjacent=2, mm_type="primer_scoring/Type_v1.csv")
#  
#  evaluate_primer("primer_scoring/01_Plecoptera_subset.fasta", "ARYATDGTRATDGCHCCDGC", 585, 604, save="save_evaluation_table_BR1.csv", mm_position ="primer_scoring/Position_v1.csv", adjacent=2, mm_type="primer_scoring/Type_v1.csv", forward=F)

## ---- eval = F-----------------------------------------------------------
#  primer_threshold("save_evaluation_table_LCO.csv", "save_evaluation_table_BR1.csv", threshold=120)

