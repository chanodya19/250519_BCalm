library(BCalm) #No need" 
library(dplyr)
library(ggplot2)
library(kableExtra)


##PREPROCESSING

#DEFINE FIRST INPUT FILE
data("BcSetExample")
nr_reps=3 #establish number of reps

#CR-optional testing
str(BcSetExample) # Check the structure
head(BcSetExample) #checkfirst few rows
dim(BcSetExample) #dimensions

#make a nice table in html format
kable(head(BcSetExample), "html") %>% 
  kable_styling("striped") %>% 
  scroll_box(width = "100%") #shown in table 


#DEFINE SECOND INPUT FILE
data("MapExample")
kable(head(MapExample), "html") %>%
  kable_styling("striped") %>%
  scroll_box(width = "100%")

#make a nice table in html format
kable(head(MapExample), "html") %>%
  kable_styling("striped") %>% 
  scroll_box(width = "100%")


##COMBINE INPUTY FILES TOGETHER
# create the variant dataframe by adding the variant ID to the DNA and RNA counts
var_df <- create_var_df(BcSetExample, MapExample)
# show the data
kable(head(var_df), "html") %>% 
  kable_styling("striped") %>% 
  scroll_box(width = "100%")


#OPTIONAL:DOWNSAMPLING TO REDUCE BARCODE
var_df <- downsample_barcodes(var_df, id_column_name="variant_id")


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
# create the variant specific MPRAset
BcVariantMPRASetExample <- MPRASet(DNA = dna_var,
                                   RNA = rna_var,
                                   eid = row.names(dna_var), 
                                   barcode = NULL)


#ELEMENT TESTING
#load tables with labes fortested/control 
#view file
data(LabelExample)
table(LabelExample)
kable(head(LabelExample), "html") %>%
  kable_styling("striped") %>% 
  scroll_box(width = "100%")


#downsample label file. Idk why we need to downsample this
elem_df <- downsample_barcodes(BcSetExample)

#make unique label files for DNA and rna
dna_elem <- create_dna_df(elem_df, id_column_name="name")
rna_elem <- create_rna_df(elem_df, id_column_name="name")

#Add labels to the MPRASet
BcLabelMPRASetExample <- MPRASet(DNA = dna_elem, 
                                 RNA = rna_elem, 
                                 eid = row.names(dna_elem), 
                                 barcode = NULL, label=LabelExample)


##ANALYSIS
#VARIANT ANALYSIS
bcs <- ncol(dna_var) / nr_reps #check block size
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(BcVariantMPRASetExample)))#design matrix, get columns with alt 
block_vector <- rep(1:nr_reps, each=bcs)#make a block vector for correlation groups so assign reps to groups
mpralm_fit_var <- mpralm(object = BcVariantMPRASetExample,#run on linear modeling
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

#ELEMeNT ANALYSIS
bcs_elem <- ncol(dna_elem) / nr_reps
block_vector <- rep(1:nr_reps, each=bcs_elem)
mpralm_fit_elem <- fit_elements(object = BcLabelMPRASetExample,
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

#tried to do a plot
ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = P.Value < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log Fold Change (logFC)", y = "-log10(P-value)") +
  theme_minimal()