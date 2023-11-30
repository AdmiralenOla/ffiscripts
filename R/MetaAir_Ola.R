################################
#           METAAIR            #
################################


require("dplyr")
require("tidyverse")
require("readr")
library("phyloseq")
library("umap")
library("plotly")

## FUNCTIONS ##

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

## USED TO DROP COLUMNS THAT DOESNT SUM TO AT LEAST X
drop_otus_below_count <- function(df, sum) {
  df_filt <- df %>%
    select_if(is.numeric) %>%
    select(where(~ sum(.x) > sum))
  return(df_filt)
}

## USED TO DROP COLUMNS THAT ARE NON-ZERO IN LESS THAN P% OF SAMPLES
drop_rare_otus <- function(df, prevalence) {
  df_filt <- df %>%
    select_if(is.numeric) %>%
    select(where(~ sum(. != 0) /length(.) >= prevalence))
  return(df_filt)
}

## USED TO DROP COLUMNS THAT DOESNT REPRESENT AT LEAST X% OF A SAMPLE. THESE ARE SET TO 0
## AFTER THIS, USE drop_constant_rows AGAIN TO REMOVE THESE.
set_low_frac_to_0 <- function(df, frac){
  ## All taxa that doesnt represent at least frac of a sample has their count set to 0
  ## After this, use drop_otus_below_count and drop_constant_rows again. NOT drop_rare_otus (what's left is OTUS very prevalent in some samples)
  # EXPLANATION (from chatGPT):
  # - rowwise() groups the data frame by rows (allows us to perform row-wise operations)
  # - mutate() then applies a transformation to each row
  # - across() is used to apply the same transformation to all columns in the row. Notice that all further statements are within across()
  # - everything() is used within across to select all columns
  # - the second argument to across is the conditional ifelse statement
  # - if THIS value (.) is less that frac times the row-sum (sum(c_across())), it is set to 0, otherwise, it is set to THIS value (.)
  # NOTE - Upon testing this function does seem a bit slow - Maybe high runtime
  df_filt <- df %>%
    select_if(is.numeric) %>%
    rowwise %>%
    mutate(across(everything(), ~ifelse(. < frac*sum(c_across()), 0, .)))
  
  return(df_filt)
}


# set_low_frac_to_0_alt <- function(df, frac){
#   df_filt <- df %>%
#     select_if(is.numeric) %>%
#     mutate(TotalCount=rowSums(.)) %>%
#     mutate(across(everything(), ~ifelse(. < frac*TotalCount,0,.))) %>%
#     select(-TotalCount)
#   return(df_filt)
# }

## USED TO DROP ALL COLS WHERE THE VARIANCE IS 0 (USUALLY THIS MEANS ALL POINTS ARE 0)
drop_constant_cols <- function(df) {
  bind_cols(df %>% select(sample),
            df %>% select_if(is.numeric) %>% select(-where(~ var(.) == 0)))
}

# THIS FUNCTION NO LONGER USED
normalize_zscore_OTU_data <- function(df) {
  df %>%
    select_if(is.numeric) %>%
    mutate(Center=rowMeans(across(everything()))) %>%
    rowwise() %>%
    mutate(SD=sd(c_across(-last_col()))) %>%
    ungroup() %>%
    mutate(across(-c(Center,SD), ~ (.-Center)/SD)) %>%
    select(-c(Center,SD))
}

## USED TO NORMALIZE COLS ACCORDING TO TSS - SCALE TO 1M READS PER SAMPLE
normalize_TSS_OTU_data <- function(df) {
  bind_cols(
    df %>% select(sample),
    df %>%
    select_if(is.numeric) %>%
    rowwise %>%
    mutate(TotalSum=sum(c_across(everything()))) %>%
    ungroup %>%
    mutate(across(-TotalSum, ~ floor(./TotalSum * 10000000))) %>% # Normalize to 10M reads
    select(-TotalSum)
  )
}

## USED TO CALCULATE THE DIFFERENCE BETWEEN THE PAIRED SAMPLES - REMEMBER TO NORMALIZE FIRST!
createDifferenceBetweenPairs <- function(df,Pair){
  df %>%
    select_if(is.numeric) %>%
    mutate(Pair=Pair) %>%
    group_by(Pair) %>%
    filter(n() > 1) %>%
    summarise_if(is.numeric, function(x) x[2] - x[1])
}

# FUNCTION FOR FINDING THE INDEX POSITIONS IN A LOWER TRI MATRIX WITH N ROWS AND (N-1) COLUMNS, ROWNR IS I, COLNR IS J, 
find_dist_indexes <- function(dists,j,i){
  N <- attr(dists, "Size")
  index <- (j-1)*(N-j/2) + i - j 
  return(index)
}



###############################
#           CODE              #
###############################

# SOME ACCESSORY DATA
NEG_CTRLS_METAAIR <- c("SL335732","SL335733","SL335734","SL335735","SL335736","SL335737","SL335738","SL335754","SL470202","SL470204","SL470208","SL470276","SL470277","SL470278","SL470279","SL470294","SL470295","SL470297","SL470298","SL470335","SL470337","SL470338")
POS_CTRLS_METAAIR <- c("SL335740","SL335741","SL470280","SL470281","SL470339")
CTRLS_METAAIR <- c(NEG_CTRLS_METAAIR, POS_CTRLS_METAAIR)
MISSING <- c("SL310938","SL310939","SL310940","SL310941","SL310942")

# IMPORT METAAIR DATA - PHYLOSEQ OR OLA METHOD

metaair <- read_tsv("/media/ubuntu/4TB-1/Ola/METAAIR/EXCLUDE_NEGATIVES/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE)
#taxa_names <- unlist(metaair[,1])
#metaair <- metaair[,2:ncol(metaair)]
#rownames(metaair) <- taxa_names
#metaair_otu <- otu_table(object = metaair,taxa_are_rows = TRUE)
metaair <- transpose_df(metaair)
colnames(metaair) <- metaair[1,]
metaair <- metaair[2:nrow(metaair),]
colnames(metaair)[1] <- "sample"

metaair_taxa_names <- read_tsv("/media/ubuntu/4TB-1/Ola/METAAIR/taxid_to_name_2023-11-07.tsv")
translation <- metaair_taxa_names$name[match(colnames(metaair)[2:ncol(metaair)],metaair_taxa_names$taxonomy_id)]
colnames(metaair) <-  c("sample",translation)

metaair_meta <- read_tsv("/media/ubuntu/4TB-1/Ola/METAAIR/META/MetaAir_metadata_final_.tsv",col_names = TRUE,col_types = "cfffffffffdf?dfdf??iididd")
metaair_meta_posremoved <- metaair_meta %>% filter(!(ID %in% POS_CTRLS_METAAIR))
metaair_meta_posremoved_negremoved <- metaair_meta %>% filter(!(ID %in% CTRLS_METAAIR))
# IGNORE WARNINGS - THEY ARE WRONG AND THE DATA IS READ CORRECTLY

# GENERATE SOME MATRIXES THAT ARE USED OFTEN
metaair_norm <- normalize_TSS_OTU_data(metaair)
# TAKING LOG TSS COUNTS INSTEAD OF RAW
metaair_log <- bind_cols(metaair_norm %>% select(sample), metaair_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
metaair_log_var <- drop_constant_cols(metaair_log)


################################
#     ANALYSIS                 #
################################

# BETA-DIVERSITY WORKS GREAT - VERY GOOD SIGNAL FROM CITY
umap.meta <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.layout <- umap.meta[["layout"]] 
umap.meta.layout <- data.frame(umap.meta.layout) 
umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta <- plot_ly(umap.meta.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 

# ALPHA-DIVERSITY - JUST ADD SAMPLE DATA TO METADATA

library(vegan)

metaair_meta_posremoved_negremoved$shannondiv <- diversity(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% select_if(is.numeric))
ggplot(metaair_meta_posremoved_negremoved, aes(x = CITY, y = shannondiv,color=YEAR)) + 
  geom_boxplot() + ylab("Shannon diversity")

# PERFORM TUKEYS HONEST SIGNIFICANT DIFFERENCE TEST
library(agricolae)
anova_result <- aov(shannondiv ~ CITY, metaair_meta_posremoved_negremoved)
tukey_result <- HSD.test(anova_result, "CITY", group = TRUE)
group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggplot(metaair_meta_posremoved_negremoved, aes(x = CITY, y = shannondiv, color=YEAR)) + # ADD COLOR YEAR?
  geom_text(data = data.frame(),
            aes(x = rownames(group_data), y = max(metaair_meta_posremoved_negremoved$shannondiv) + 1, label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot() +
  ggtitle("Shannon diversity by City") +
  xlab("City") +
  ylab("Shannon diversity index")


## STARTING TO WORK ON CORE DEFINITIONS

# IDEA:
# Function must include minimum threshold AVERAGE abundance (e.g. 0.0005), should be possible to adjust
# Function should allow you to not include contaminants
# Function should allow you to cut species where prevalence in negatives is similar to prevalence in samples. (Proportional test for difference?)
# Need a step function for prevalence to calculate core at different prevalences. Prevalence should be on Y axis
# The number of species with prevalence >= than Y is plotted on X axis (OR - reverse axes so that response-outcome is more logical)
# Output should simply be a data frame of x and y coordinates?

filter_core <- function(df,includeContaminants=TRUE,includePankitome=TRUE,minimumAverageAbundance=0.00005,differentialAbundanceFromNegs=TRUE,significance=0.05){
  if(!includeContaminants){
    df <- df %>%
      filter(Contaminants == 0)
  }
  if(!includePankitome){
    df <- df %>%
      filter(Pankitome == 0)
  }
  if(differentialAbundanceFromNegs){
    N_samples <- 787 #rep(787,nrow(df))
    N_negatives <- 22 #rep(22,nrow(df))
    X_samples <- ceiling(df$frequency_perc * N_samples)
    X_negatives <- ceiling(df$`Negative Controls` * N_negatives)
    #p <- prop.test(x=c(,),n=c(N_samples,N_negatives))
    a <- c(1:nrow(df))
    myrowprop <- function(i) {
      suppressWarnings(prop.test(c(X_samples[i],X_negatives[i]),c(787,22))$p.value)
    }
    ListofResults <- unlist(lapply(a,myrowprop))
    df$p.value <- ListofResults
    df <- df %>%
      filter(p.value < significance)
  }
  
  # FILTER MIN ABUNDANCE
  df <- df %>%
    filter(counts_perc_avg >= minimumAverageAbundance)
}

calculate_core <- function(df){
  prev_range <- seq(1.0,0.5,-0.01)
  abundance_range <- c(0.00005,0.003,0.005,0.01)
  combinations <- expand_grid(prev_range, abundance_range)
  
  #num_core <- sapply(X=prev_range,FUN= function(x) length(which(x<df$frequency_perc)))
  num_core <- apply(combinations,MARGIN=1,FUN=num_core_function,df=df)
  return(data.frame(X=num_core,Y=combinations$prev_range,Z=factor(combinations$abundance_range)))
}

num_core_function <- function(combinations,df){
   return(length(which(combinations[1] < df$frequency_perc & combinations[2] < df$counts_perc_avg)))
}

## COMBINED ##
core_mat_combined <- read_tsv(file="/media/ubuntu/4TB-1/Ola/METAAIR/RESULTS/GRIMER_COMBINED_PREVALENCE.tsv",col_names = TRUE)
core_mat_combined_filtered <- filter_core(core_mat_combined,includeContaminants = FALSE)
core_mat_combined_numbers <- calculate_core(core_mat_combined_filtered)
ggplot(data=core_mat_combined_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

## BACTERIA ##
core_mat_bacteria <- read_tsv(file="/media/ubuntu/4TB-1/Ola/METAAIR/RESULTS/GRIMER_BACTERIA_PREVALENCE.tsv",col_names = TRUE)
core_mat_bacteria_filtered <- filter_core(core_mat_bacteria,includeContaminants = FALSE)
core_mat_bacteria_numbers <- calculate_core(core_mat_bacteria_filtered)
ggplot(data=core_mat_bacteria_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

## FUNGI ##
core_mat_fungi <- read_tsv(file="/media/ubuntu/4TB-1/Ola/METAAIR/RESULTS/GRIMER_FUNGI_PREVALENCE.tsv",col_names = TRUE)
core_mat_fungi_filtered <- filter_core(core_mat_fungi,includeContaminants = FALSE)
core_mat_fungi_numbers <- calculate_core(core_mat_fungi_filtered)
ggplot(data=core_mat_fungi_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
# 
