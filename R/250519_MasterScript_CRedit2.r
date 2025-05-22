library("BCalm")
library(dplyr)
library(ggplot2)
library(kableExtra)

#load dataframe with barcode, CRS name DNAcount 1-3 and RNA count 1-3, Dna_1 and RNA_1 and DNA_2 versa
data("BcSetExample") #for us use zcat reporter_experiment.barcode.NPC.fromFile.default.min_oligo_threshold_10.tsv.gz | head
nr_reps = 3 
# show the data. does not make a dataframe object
kable(head(BcSetExample), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%") #make html file and make it scrollable


#Variant testing
#load the variant map
data("MapExample") 
# show the data
kable(head(MapExample), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")
#BC set is for counts of each barcode in each oligo, mapexample finding variant ID so RS number in column one and then Ref or Alt in column headings for each CRS epi_test_001

# create the variant dataframe by adding the variant ID to the DNA and RNA counts
?create_var()
var_df <- create_var_df(BcSetExample, MapExample) 
#merge the two datafgrames
# show the data
kable(head(var_df), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")


#establish the exact tables required to be asse
#learninhng to doownsample barides to reduce the number of barcdoes - simplification step 
##DOWNSAMPLING IS OPTIONAL
var_df <- downsample_barcodes(var_df, id_column_name="variant_id")
var_df

#prepare data for analysis. Make both dna_var and rna_var dataframes by using function create_dna_df and createrna_df 
dna_var <- create_dna_df(var_df)
kable(head(dna_var), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")
rna_var <- create_rna_df(var_df)
kable(head(rna_var), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")
#OUTPUT COLUMNS SAMPLE_COUNT_1_BC1_REF

#mpraset input is made 
# creatimg the variant specific MPRASet
BcVariantMPRASetExample <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)


# ELEMENT TESTING
# first need to add labels between creference and alt SNPs. Easy for differential group testing
data(LabelExample)
table(LabelExample)
kable(head(LabelExample), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")
#OUTPUT , ROW HEADING IS OLIGO IDENTITY. x AXIS TESTED VERSUS CONTROL

#then again we down sample. Makes a new dtaframe as output call elem_df 
elem_df <- downsample_barcodes(BcSetExample)

#CHANGE THE NAMES OF THE COLUMN idS 
dna_elem <- create_dna_df(elem_df, id_column_name="name")
rna_elem <- create_rna_df(elem_df, id_column_name="name")

#To compare between test and control, we need to add the labels to the MPRASet.
BcLabelMPRASetExample <- MPRASet(DNA = dna_elem, RNA = rna_elem, eid = row.names(dna_elem), barcode = NULL, label=LabelExample)


###ANALYSIS
##Variant analysis
bcs <- ncol(dna_var) / nr_reps #DIVIDE number of columns in dna_var df from number of reps - 3
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(BcVariantMPRASetExample))) #column intcpt made and retrief column anmes from the df. grep for word search alt in column name vetor
block_vector <- rep(1:nr_reps, each=bcs) #store output pf  length in block length
mpralm_fit_var <- mpralm(object = BcVariantMPRASetExample, 
                         design = design, 
                         aggregate = "none", 
                         normalize = TRUE, 
                         model_type = "corr_groups", 
                         plot = FALSE, 
                         block = block_vector) #CHECK LB warning message
#output of mpralm 
top_var <- topTable(mpralm_fit_var, coef = 2, number = Inf)
kable(head(rna_var), "html") %>% kable_styling("striped") %>% scroll_box(width = "100%")

#use data from top_var to make a plot. Log_FC in x axis and p value in . geompoint to add points to ploy
#alpha for transparency level of the dots. this is a volcano plot
ggplot(top_var, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6)

##Element Analysis
#fit_elements function to take the MPRASett object as input and apply stat modeling
#block_Vector is reference on which BC_to_CRS association
#normalise the total count in RNA and dna libraries 
bcs <- ncol(dna_elem) / nr_reps
block_vector <- rep(1:nr_reps, each=bcs)
mpralm_fit_elem <- fit_elements(object = BcLabelMPRASetExample, 
                                normalize=TRUE, block = 
                                  block_vector, plot = FALSE)


#synthesis of plots 
##output probides a bar gra[h woyj yje mehative control and percentiles ]
plot_groups(mpralm_fit_elem, 0.975, #fit elements in the 0.975th percentile of negative controls 
            neg_label="control",
            test_label="tested")

#mpra_treat applies limma's treat() function. Does the important T-test 
treat <- mpra_treat(mpralm_fit_elem, 0.975, neg_label="control")
result <- topTreat(treat, coef = 1, number = Inf)
head(result) #IMPORTANT 

#attempt at a volvcano plot from result
ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = P.Value < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log Fold Change (logFC)", y = "-log10(P-value)") +
  theme_minimal()

