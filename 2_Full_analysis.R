# FULL ANALYSIS (ALL SAMPLES ON 96 WELL - TO SEE HOW THEY CLUSTER)

require("dplyr")
require("tidyverse")
require("readr")


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


## USED TO DROP ALL COLS WHERE THE VARIANCE IS 0 (USUALLY THIS MEANS ALL POINTS ARE 0)
drop_constant_cols <- function(df) {
  col_variances <- df %>%
    select(where(is.numeric)) %>%
    map_dbl(var)
  zero_variance_cols <- names(col_variances[col_variances==0])
  df %>%
    select(-(any_of(zero_variance_cols)))
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
  df %>%
    rowwise %>%
    mutate(TotalSum=sum(c_across(where(is.numeric)))) %>%
    ungroup %>%
    mutate(across(-c(where(is.character),TotalSum), ~ floor(./TotalSum * 1000000))) %>% # Normalize to 1M reads
    select(-TotalSum)
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

##################### CODE ###############################

stok_full <- read_tsv("/media/ubuntu/4TB-1/Ola/Kari_stokastisitet/AGGREGATED_DATA/SPECIES/2023-03-16_Full_data_species.tsv",col_names=TRUE)
stok_full_t <- transpose_df(stok_full)
colnames(stok_full_t) <- stok_full_t[1,]
stok_full_t <- stok_full_t[2:nrow(stok_full_t),]
colnames(stok_full_t)[1] <- "sample"

# FIX NAMES
full_names <- stok_full_t$sample
full_names[83:85] <- c("95_SEA","96_POS","97_NEG")

# FIRST - NON-NORMALIZED
fulltest <- hclust(vegdist(stok_full_t[,2:ncol(stok_full_t)]),method = "ward.D2")
plot(fulltest,labels=full_names,axes=FALSE,ann=FALSE,ylab="",main="Cluster dendrogram of pairwise Bray-Curtis distances",cex=.5)

# THEN - NORMALIZED  - NOTE THAT BRAY-CURTIS SHOULD PERHAPS NOT BE NORMALIZED, SINCE THE DIVISOR IN THE FORMULA WILL BECOME THE SAME (2M) ALWAYS
# OTOH, SAMPLES SUCH AS [1,2,3,4,5] AND [2,4,6,8,10] WILL HAVE A SIGNIFICANT BC DISTANCE DESPITE BEING COMPLETELY EQUAL IN RELATIVE ABUNDANCE
# BRAY-CURTIS MEASURES "SIZE", NOT JUST "SHAPE" OF THE COMPOSITION. (SIZE=ROWSUMS, THE TOTAL OUTPUT PER SAMPLE). CHI-SQUARE IS NEALY PURE SHAPE
# FOR NORMALIZED COUNTS WE CAN USE CHI-SQUARE DISTANCE ON THIS "PROFILE" OF RELATVE ABUNDANCES
stok_full_norm <- normalize_TSS_OTU_data(stok_full_t)
fullnormtest <- hclust(vegdist(stok_full_norm[,2:ncol(stok_full_norm)],method = "chisq"),method = "ward.D2")
plot(fullnormtest,labels=full_names,axes=FALSE,ann=FALSE,ylab="",main="Cluster dendrogram of pairwise Bray-Curtis distances",cex=.5)

