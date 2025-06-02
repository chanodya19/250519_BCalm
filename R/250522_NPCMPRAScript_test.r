library("BCalm")
library(dplyr)
library(ggplot2)
library(kableExtra)


##PREPROCESSING

#DEFINE FIRST INPUT FILE
BC_Oligo_file_path <- "250219_MPRA_SAGCQR1327/MPRAsnakeflow/results/experiments/npcBbmapSwapped/reporter_experiment.barcode.NPC.fromFile.default.min_oligo_threshold_10.tsv.gz"
file.exists(BC_Oligo_file_path) #check if file path is recognisable
BC_to_oligo_data <- read.csv(gzfile(BC_Oligo_file_path), header = TRUE, sep = "\t") #make dataobject of counts vs CRS_name
nr_reps=3 #establish number of reps

#few quick checks on the data
str(BC_to_oligo_data)  # Check the structure
head(BC_to_oligo_data)  # View the first few rows
dim(BC_to_oligo_data)   # Check number of rows and columns

#make a nice table in html format
kable(head(BC_to_oligo_data), "html") %>% 
    kable_styling("striped") %>% 
    scroll_box(width = "100%") #shown in table 
#Problem - can only see first 6 rows no mtter what i try
#okay figures this out its because I have head LOL KMN

#DEFINE SECOND INPUT FILE
#Makeing the Ref/Alt to ID file
Ref_Alt_to_ID_file_path <- "250219_MPRA_SAGCQR1327/metadata/"

#make vector for refvs alt for each SNP ID - CRS name as cell results
Ref_Alt_to_ID_data <- read.csv(gzfile(Ref_Alt_to_ID_file_path), header = TRUE, sep = "\t")
#few quick checks on the data
str(Ref_Alt_to_ID_data)  # Check the structure
head(Ref_Alt_to_ID_data)  # View the first few rows
dim(Ref_Alt_to_ID_data)   # Check number of rows and column

#make a nice table in html format
kable(head(Ref_Alt_to_ID_data), "html") %>% 
    kable_styling("striped") %>% 
    scroll_box(width = "100%")

##COMBINE INPUTY FILES TOGETHER
var_df <- create_var_df(BC_to_oligo_data,Ref_Alt_to_ID_data)
kable(head(var_df), "html") %>%  #Show the data nicely formatted
    kable_styling("striped") %>% 
    scroll_box(width = "100%")
#Problem : may have issues as the first file the names ae  actuallly IDs in 


#OPTIONAL:DOWNSAMPLING TO REDUCE BARCODE
var_df <- downsample_barcodes(var_df, id_column_name="ID")
var_df

#PROPROCESS FOR BC_TO_OLIGO_DATA FOR MPRASET
##Separate RNA and DNA count for MPRAset. Make both dna_var and rna_var dataframes by using function create_dna_df and createrna_df 
dna_var <- create_dna_df(var_df) #DNA
kable(head(dna_var), "html") %>% 
    kable_styling("striped") %>% 
    scroll_box(width = "100%")
rna_var <- create_rna_df(var_df) #RNA
kable(head(rna_var), "html") %>% 
    kable_styling("striped") %>% 
    scroll_box(width = "100%")

#MAKE THE INPUT TO BCALM - MPRASET 
bc_to_oligo_MPRASet <- MPRASet(DNA = dna_var,
                                RNA = rna_var, 
                                eid = row.names(dna_var), 
                                barcode = NULL)


#ELEMENT TESTING
#load tables with labes for ested/control 
oligo_label_file_path <- "250219_MPRA_SAGCQR1327/metadata/"#define file path
oligo_label <-read.csv(gzfile(oligo_label_file_path), header = TRUE, sep = "\t")

#view file
data(oligo_label)
table(oligo_label)
kable(head(oligo_label), "html") %>% 
    kable_styling("striped") %>% 
    scroll_box(width = "100%") #make table

#downsample label file. Idk why we need to downsample this
elem_df <- downsample_barcodes(oligo_label)

#make unique label files for DNA and rna
dna_elem <- create_dna_df(elem_df, id_column_name="name")
rna_elem <- create_rna_df(elem_df, id_column_name="name")

#Add labels to the MPRASet
bc_label_MPRASet  <- MPRASet(DNA = dna_elem, 
                            RNA = rna_elem, 
                            eid = row.names(dna_elem), 
                            barcode = NULL, 
                            label=oligo_label)

##ANALYSIS
#VARIANT ANALYSIS
bcs_var <- ncol(dna_var) / nr_reps #check block size
design <- data.frame(intcpt = 1, 
                    alt = grepl("alt", colnames(bc_to_oligo_MPRASet))) #design matrix, get columns with alt 
block_vector <- rep(1:nr_reps, each=bcs_var)#make a block vector for correlation groups so assign reps to groups
mpralm_fit_var <- mpralm(object = bc_to_oligo_MPRASet, #run on linear modeling
                        design = design, 
                        aggregate = "none", 
                        normalize = TRUE, #normalisation 
                        model_type = "corr_groups", 
                        plot = FALSE, 
                        block = block_vector) #considers replicates

top_var <- topTable(mpralm_fit_var, 
                    coef = 2, #can adjust this
                    number = Inf) #extract top results with all sig for coeeficient 2
kable(head(rna_var), "html") %>% 
kable_styling("striped") %>% 
scroll_box(width = "100%")

#plot this as volcano plot 
ggplot(top_var, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6)

#ELEMNT ANALYSIS
bcs_elem <- ncol(dna_elem) / nr_reps
block_vector <- rep(1:nr_reps, each=bcs_elem)
mpralm_fit_elem <- fit_elements(object = bc_label_MPRASet, 
                                normalize=TRUE, 
                                block = block_vector, 
                                plot = FALSE)

#VISUALISATION
#results from fir_elements is used here to compare log ratios
plot_groups(mpralm_fit_elem, 0.975, neg_label="control", test_label="tested")
#makes histogram with logration values

#Do T-test
treat <- mpra_treat(mpralm_fit_elem, 0.975, neg_label="control")
result <- topTreat(treat, coef = 1, number = Inf)
head(result)

#attempted to do a plot 
ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = P.Value < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log Fold Change (logFC)", y = "-log10(P-value)") +
  theme_minimal()