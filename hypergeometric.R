##Probability of getting an overlap of at least a particular size given two lists and the hypergeometric distribution

#### Data input is of the format:

## column_1		column_2		column_3		column_4
## overlap_size 	size_of_list_A		Population_size		size_of_list_B
## overlap_size 	size_of_list_A		Population_size		size_of_list_B
## overlap_size 	size_of_list_A		Population_size		size_of_list_B
## overlap_size 	size_of_list_A		Population_size		size_of_list_B
## overlap_size 	size_of_list_A		Population_size		size_of_list_B
## ....			....			....			.... 

## where: 
## overlap_size is the intersection of list A and list B
## size_of_list_A is the number of proteins in list A (e.g. the number of proteins in the reference ##proteome annotated as having a particular feature)
## Population_size is the entire population size to select from (e.g. the reference proteome size or ##the number of proteins in the 'conservome')
## size_of_list_B is the number of proteins in list B (e.g. the number of proteins in the mRBPome ##annotated as having a particular feature).
##
##
## read in csv file containing data input, in above format
## obtain a raw p-value for obtaining an overlap using hypergeometic distribution
## adjust raw p-value using Benjami-Hochberg false discovery rate correction
## write csv file with the results, including the counts 

data_input <- read.csv("data_input_file.csv")
raw_p <- 1-phyper(data_input[,1]-1,data_input[,2],data_input[,3]-data_input[,2],data_input[,4])
BH_adjusted_p_values <- p.adjust(raw_p,method="BH")
write.csv(cbind(data_input,BH_adjusted_p_values), file="my_output_results.csv",quote=FALSE)
