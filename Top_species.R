

require("dplyr")
require("tidyverse")
require("readr")
require("ggplot2")

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

# Get the top N species and choose whether to filter out contaminants
getTopSpecies <- function(df, N, metadata, includeNegsKit=TRUE,includeNegsWater=TRUE,includeContams=TRUE,includeHuman=TRUE, includeDecontam = TRUE, report_percentage=TRUE){
  exclude_columns <- c()
  if (!includeContams){
    rowIndex <- which(df$Sample == "Is_contaminant")
    exclude_columns_contam_ind <- df %>%
      select(where(is.numeric)) %>%
      slice(rowIndex)# %>%
      #unlist %>%
      #as.logical()
    exclude_columns_contam <- colnames(exclude_columns_contam_ind)[exclude_columns_contam_ind==1]
    exclude_columns <- union(exclude_columns, exclude_columns_contam)
  }
  if (!includeHuman){
    rowIndex <- which(df$Sample == "Is_human_related")
    exclude_columns_human_ind <- df %>%
      select(where(is.numeric)) %>%
      slice(rowIndex)
    exclude_columns_human <- colnames(exclude_columns_human_ind)[exclude_columns_human_ind==1]
    exclude_columns <- union(exclude_columns, exclude_columns_human)
  }
  if (!includeNegsKit){
    # FOR KIT, CONSIDER ADDING MINIMUM READ THRESHOLD
    rowIndex <- which(df$Sample == "NEG_KIT")
    exclude_columns_negkit_ind <- df %>%
      select(where(is.numeric)) %>%
      slice(rowIndex)
    exclude_columns_negkit<- colnames(exclude_columns_negkit_ind)[exclude_columns_negkit_ind>0]
    exclude_columns <- union(exclude_columns,exclude_columns_negkit)
  }
  if (!includeNegsWater){
    # FOR WATER, CONSIDER ADDING MINIMUM READ THRESHOLD
    rowIndex <- which(df$Sample == "NEG_WATER")
    exclude_columns_negwater_ind <- df %>%
      select(where(is.numeric)) %>%
      slice(rowIndex)
    exclude_columns_negwater <- colnames(exclude_columns_negwater_ind)[exclude_columns_negwater_ind>0]
    exclude_columns <- union(exclude_columns,exclude_columns_negwater)
  }
  if (!includeDecontam){
    rowIndex <- which(df$Sample == "Decontam")
    exclude_columns_decontam_ind <- df %>%
      select(where(is.numeric)) %>%
      slice(rowIndex)
    exclude_columns_decontam <- colnames(exclude_columns_decontam_ind)[exclude_columns_decontam_ind==1]
    exclude_columns <- union(exclude_columns,exclude_columns_decontam)
  }
  # CALCULATE ROWSUMS
  df <- df %>%
    mutate(rowSums = select(.,-Sample) %>% rowSums(.))

  top_species_group <- df %>%
    select(-one_of(exclude_columns)) %>% # Select only columns to be included - VERIFY that this works!
    filter(Sample %in% metadata) %>%
    summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>% #.names = "{.col}_Percentage"
    select(-rowSums) %>%
    unlist() %>%
    sort(decreasing=TRUE) %>%
    head(N)
  species_to_get_data_for <- names(top_species_group)
  top_species_group <- addOtherSpecies(top_species_group)
  results <- data.frame(Fractions=top_species_group, Species=names(top_species_group),Sample="Group")
  
  for(sample in metadata){
    top_species <- df %>%
      # select(-one_of(exclude_columns)) %>% No longer necessary. Rowsums calculated
      filter(Sample == sample) %>%
      summarise(across(where(is.numeric), ~sum(.) / rowSums)) %>%
      select(all_of(species_to_get_data_for)) %>%
      unlist() %>%
      addOtherSpecies()
    # Now return only those species that are already in the group top
    # Finally, row bind results
    tmp_df <- data.frame(Fractions=top_species, Species=names(top_species),Sample=sample)
    results <- bind_rows(results, tmp_df)
  }
  #results <- data.frame(Fractions=top_species, Species=names(top_species))
  results <- reOrderDataForPlot(results)
  # Add legend colors
  #Legend.color <- '#000000'
  contamIndex <- which(df$Sample == "Is_contaminant")
  contamrow <- df[contamIndex,]
  contamrow$Other <- 0
  Legend.color <- rep('#000000',nrow(results))
  for(i in 1:nrow(results)){
    #Species <- row
    Lookup <- results$Species[i]
    Lookupvalue <- contamrow %>% select(any_of(Lookup))
    Legend.color[i] <- ifelse(as.logical(Lookupvalue),'#FF0000','#000000')
  }
  results$Legend.color <- Legend.color
  return(results)
}

addOtherSpecies <- function(vec){
  missing <- 1.0 - sum(vec)
  vec <- c(vec, missing)
  names(vec)[length(vec)] <- "Other"
  return(vec)
}

reOrderDataForPlot <- function(df){
  #tmp_df <- df %>% filter(Sample == "Group")
  #tmp_species <- reorder(factor(tmp_df$Species),X=tmp_df$Fractions)
  #tmp_species <- relevel(tmp_species, ref=nrow(tmp_df))
  #print(tmp_species)
  #print(levels(tmp_species))
  # Reorder factor levels so that level 1 is the least common and so forth.
  #df$Species <- reorder(factor(tmp_df$Species),X=tmp_df$Fractions)
  df$Species <- reorder(factor(df$Species),X=df$Fractions)
  df$Species <- relevel(df$Species,ref="Other")
  # The level for "Other" should always be 1
  
  return(df)
}

plotStackedBarChart <- function(df, metadata){
  # OUTLINE
  # her vil condition tilsvare species , value være prosent og specie være conditions (f.eks summer, spring etc)
  # Legg til at Contams har rød skrift

  ggplot(df, aes(fill=Species, y=Fractions, x=Sample,legend.text=Legend.color)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = my.colors) +
    theme(legend.position = "bottom", legend.margin = unit(c(1,1,1,1), "cm"), legend.text=element_text(size=5),axis.text.x = element_text(angle=90))
}
########## COLORS ###################################
my.colors <- c('#999999','#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')


########## CODE #####################################

kari_full <- read_tsv("/media/ubuntu/Pandora/KARI21_ANALYSE/AGGREGATED_DATA/BRACKEN_FBAV/SEA/AGGREGATED_0.005_W_CTRLS_CONTAMFLAGGED_DECONTAM.tsv",col_names=TRUE)
kari_full_t <- transpose_df(kari_full)
colnames(kari_full_t) <- kari_full_t[1,]
kari_full_t <- kari_full_t[2:nrow(kari_full_t),]
colnames(kari_full_t)[1] <- "Sample"


kari_s_names <- read_tsv("/media/ubuntu/Pandora/KARI21_ANALYSE/AGGREGATED_DATA/BRACKEN_FBAV/SEA/taxid_to_name.tsv")

# MATCHING NAMES
translation <- kari_s_names$name[match(colnames(kari_full_t)[2:ncol(kari_full_t)],kari_s_names$taxonomy_id)]
# Gonna try to set this as the colnames
colnames(kari_full_t) <-  c("Sample",translation)

kari_long <- pivot_longer(kari_full_t, cols=2:ncol(kari_full_t))

# Reading metadata
kari_meta <- read_tsv("/media/ubuntu/Pandora/KARI21_ANALYSE/AGGREGATED_DATA/BRACKEN_FBAV/SEA/Kari_SEA_metadata_.tsv",col_names = TRUE,col_types="cfffffddDDDf")


# Expand so that table comes in ggplot2 ready form, eg
# SUMMER  SPECIES_1 FRAC_1
# SUMMER  SPECIES_2 FRAC_2
# WINTER  SPECIES_1 FRAC_1

SEA_1_data <- getTopSpecies(kari_full_t,20,"SEA_1",includeContams = FALSE,includeDecontam = FALSE)
SEA_1_ggplot <- data.frame(Species=names(SEA_1_data),Fractions=SEA_1_data)
ggplot(SEA_1_ggplot, aes(fill=Species,y=Fractions,x="SEA_1")) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = my.colors)
# NOTE - Make Sure "Others" is always at bottom or top.
