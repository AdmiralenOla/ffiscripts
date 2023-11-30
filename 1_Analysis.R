## PRE-CONFIG - REQUIRED TO INSTALL SOME PACKAGES ##

#Sys.setenv(CMAKE_BIN="/home/ubuntu/miniconda2/envs/python3.7_environment/bin/cmake")
#install.packages("car")
#install.packages("lme4")
#install.packages("factoextra")
#install.packages("umap")
#install.packages("plotly")
## 

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
  df %>%
    select_if(is.numeric) %>%
    select(-where(~ var(.) == 0))
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
    select_if(is.numeric) %>%
    rowwise %>%
    mutate(TotalSum=sum(c_across(everything()))) %>%
    ungroup %>%
    mutate(across(-TotalSum, ~ floor(./TotalSum * 1000000))) %>% # Normalize to 1M reads
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

# FUNCTION FOR FINDING THE INDEX POSITIONS IN A LOWER TRI MATRIX WITH N ROWS AND (N-1) COLUMNS, ROWNR IS I, COLNR IS J, 
find_dist_indexes <- function(dists,j,i){
  N <- attr(dists, "Size")
  index <- (j-1)*(N-j/2) + i - j 
  return(index)
}

## CODE ##

#stok_s <- read_tsv("/media/ubuntu/4TB-21/Ola/Kari_stokastisitet/AGGREGATED_DATA/SPECIES/Aggregated_data_species.tsv",col_names=TRUE) # OLD DATASET
#stok_s <- read_tsv("/media/ubuntu/4TB-1/Ola/Kari_stokastisitet/AGGREGATED_DATA/2023-10-18_AGGREGATED_FBAV_0.005.tsv",col_names=TRUE) # CURRENT DATASET
stok_s <- read_tsv("/media/ubuntu/4TB-1/Ola/Kari_stokastisitet/AGGREGATED_DATA/2023-10-18_AGGREGATED_FBAV_0.005_NOKITOME.tsv",col_names=TRUE) # NO_KITOME DATASET

stok_s_t <- transpose_df(stok_s)
colnames(stok_s_t) <- stok_s_t[1,]
stok_s_t <- stok_s_t[2:nrow(stok_s_t),]
colnames(stok_s_t)[1] <- "sample"

#stok_s_names <- read_tsv("/media/ubuntu/4TB-21/Ola/Kari_stokastisitet/AGGREGATED_DATA/SPECIES/taxid_to_name_species.tsv") # OLD VERSION
stok_s_names <- read_tsv("/media/ubuntu/4TB-1/Ola/Kari_stokastisitet/AGGREGATED_DATA/taxid_to_name.tsv")

# MATCHING NAMES
translation <- stok_s_names$name[match(colnames(stok_s_t)[2:ncol(stok_s_t)],stok_s_names$taxonomy_id)]
# Gonna try to set this as the colnames
colnames(stok_s_t) <-  c("sample",translation)



# Getting correct IDS - NO LONGER NEEDED OCT 2023
#stok_s_t <-stok_s_t %>%
#  mutate(sample=substring(sample,first=4))

# Reading metadata
stok_meta <- read_tsv("/media/ubuntu/4TB-1/Ola/Kari_stokastisitet/METADATA/Kari_stokastisitet_metadata.tsv",col_names = TRUE,col_types="cffffddDDDDff")

# Drop meta for empty sample 79
stok_meta <- stok_meta %>% filter(id != "SEA-79")

## NOTE - ALSO REMOVE UNPAIRED SAMPLE SEA_76
stok_meta <- stok_meta %>% filter(id != "SEA-76") %>% dplyr::arrange(id)
stok_s_t <- stok_s_t %>% filter(sample != "SEA-76")

# GENERATE SOME MATRIXES THAT ARE USED OFTEN
stok_norm <- normalize_TSS_OTU_data(stok_s_t)
# TAKING LOG TSS COUNTS INSTEAD OF RAW
stok_log <- bind_cols(stok_norm %>% select(sample), stok_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))

###############
#    NOTE     # 
###############
# NOTE THAT BRAY-CURTIS SHOULD PERHAPS NOT BE NORMALIZED, SINCE THE DIVISOR IN THE FORMULA WILL BECOME THE SAME (2M) ALWAYS
# OTOH, SAMPLES SUCH AS [1,2,3,4,5] AND [2,4,6,8,10] WILL HAVE A SIGNIFICANT BC DISTANCE DESPITE BEING COMPLETELY EQUAL IN RELATIVE ABUNDANCE
# BRAY-CURTIS MEASURES "SIZE", NOT JUST "SHAPE" OF THE COMPOSITION. (SIZE=ROWSUMS, THE TOTAL OUTPUT PER SAMPLE). CHI-SQUARE IS NEALY PURE SHAPE
# FOR NORMALIZED COUNTS WE CAN USE CHI-SQUARE DISTANCE ON THIS "PROFILE" OF RELATVE ABUNDANCES
###############

# stok_s_t_var <- drop_constant_cols(stok_s_t) # No longer needed since all remaining taxa have variances due to Bracken filtering


# Add 1 and 2 to a pair - NOTE, this removes control samples
stok_meta_ctrls <- stok_meta %>% filter(Pair %in% c("Reagent","H2O","Standard")) %>% mutate(Pairtype=0)
stok_meta <- stok_meta %>%
  group_by(Pair) %>%
  filter(n() > 1) %>%
  mutate(Pairtype=c(1,2))
stok_meta <- bind_rows(stok_meta,stok_meta_ctrls) # Adding ctrls back again

# Doing some summing of OTUs to investigate how much will be cut off at different thresholds
colsums <- apply(stok_s_t[,2:ncol(stok_s_t)],2,sum)
quantile(colsums,c(0.10,0.50,0.90))
ecdf(colsums)(100) # ecdf(colsums) is now a *function*. Investigate what percentile 100 corresponds to

#rowsums <- apply(stok_s_t[,2:ncol(stok_s_t)],1,sum) # NO longer needed Oct 2023

#prevalence <- function(x){sum((x!=0))/length(x)}
prev <- apply(stok_s_t[,2:ncol(stok_s_t)],2, function(x) sum(x!=0)/length(x))
ecdf(prev)(0.05)
# Original dims 6405x39
# Consider cutoffs for the following
# - OTU count sum less than 100 (approx 17%)
# - OTUs that are present in only a few samples (e.g. 10%? 5%)

## APPLY FURTHER CUTOFFS HERE ##


################################


## INITIAL PCA

#stok_s_t_var <- drop_constant_cols(stok_s_t)
#stok_s_t_var <- drop_constant_cols(stok_norm) # COMMENT ME OUT WHEN DONE
stok_s_t_var <- drop_constant_cols(stok_log)

pca1 <- prcomp(stok_s_t_var %>% select(where(is.numeric)),scale.=TRUE,center=TRUE)
plot(pca1) # Great Scree plot
biplot(pca1) # Pretty decent spread #
library("factoextra")
fviz_eig(pca1) # Awesome Scree plot
fviz_pca_ind(pca1,habillage = stok_meta$Pair,geom.ind = "point") # <- kind of ugly plot. But it shows some degree of stochasticity
fviz_pca_ind(pca1,habillage = stok_meta$Season,geom.ind = "point",ellipses=TRUE)
plot1 <- fviz_pca_ind(pca1) # TRY to add shapes, colors etc manually. Ie color = season, shape=time of day, pairs have lines between them
pca.layout <- cbind(pca1$x[,1:2],stok_meta)

require(plotly)
fig.pca <- plot_ly(pca.layout, x = ~PC1, y = ~PC2, color = ~Season, symbol = ~Time, type = 'scatter', mode = 'markers', size=10) %>% # colors = c('#636EFA','#EF553B','#00CC96')
#fig.pca <- plot_ly(pca.layout, x = ~PC1, y = ~PC2, color = ~Season, type = 'scatter', mode='markers',size=10,legendgroup="Season", legend=FALSE) %>%
#  add_trace(x = ~PC1, y = ~PC2, symbol=~Time, type = 'scatter', mode='markers',size=10,legendgroup="Time of day") #%>%
    layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "PC1 (19%)"),  
    yaxis = list( 
      title = "PC2 (11%)"))

add_segments(fig.pca,
           x = ~PC1[Pairtype == 1], 
           xend= ~ PC1[Pairtype == 2], 
           y = ~PC2[Pairtype == 1], 
           yend = ~PC2[Pairtype == 2],
           inherit=FALSE,
           color=I("black"),
           #line = list(width=1,dash="dot"),
           line = list(width=1),
           #line = list(dash = "dash"),
           #color = ~ pairtype[pairtype == "B"], 
           #symbol = ~ pairtype[pairtype == "B"], 
           #line=list(color="black"), 
           showlegend=FALSE) # %>%
 # add_annotations(text=~id,legendtitle=TRUE, showarrow=FALSE)



fig.pca
### END INITIAL PCA ###

#####################################
### OCT 2023 - NO LONGER RELEVANT ###
#####################################
## PCA with more stringent removal of data. Summed OTU count needs to be 500 and prevalence 20% : ###
stok_s_t_var <- drop_constant_rows(stok_s_t)
#stok_s_t_var <- drop_constant_cols(stok_norm) # COMMENT OUT WHEN DONE
stok_s_t_var_filt <- drop_otus_below_count(stok_s_t_var, 500)
stok_s_t_var_filt_common <- drop_rare_otus(stok_s_t_var_filt, 0.2)

pca2 <- prcomp(stok_s_t_var_filt_common,scale.=TRUE,center=TRUE)

pca2.layout <- cbind(pca2$x[,1:2],stok_meta)

fig.pca2 <- plot_ly(pca2.layout, x = ~PC1, y = ~PC2, color = ~Season, symbol = ~Time, type = 'scatter', mode = 'markers', size=10) %>% # colors = c('#636EFA','#EF553B','#00CC96')
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "PC1 (27%)"),  
    yaxis = list( 
      title = "PC2 (14%)"))

add_segments(fig.pca2,
             x = ~PC1[Pairtype == 1], 
             xend= ~ PC1[Pairtype == 2], 
             y = ~PC2[Pairtype == 1], 
             yend = ~PC2[Pairtype == 2],
             inherit=FALSE,
             color=I("black"),
             #line = list(width=1,dash="dot"),
             line = list(width=1),
             #line = list(dash = "dash"),
             #color = ~ pairtype[pairtype == "B"], 
             #symbol = ~ pairtype[pairtype == "B"], 
             #line=list(color="black"), 
             showlegend=FALSE) # %>%

### END PCA2

### REALLY STRICT PCA ###
#stok_s_t_var <- drop_constant_cols(stok_norm)
stok_s_t_var <- drop_constant_cols(stok_s_t)
stok_s_t_top <- set_low_frac_to_0_alt(stok_s_t_var, 0.001)
stok_s_t_top_var <- drop_constant_rows(stok_s_t_top)
#stok_s_t_var_filt <- drop_otus_below_count(stok_s_t_var, 500)
stok_s_t_top_var_common <- drop_rare_otus(stok_s_t_top_var, 0.05) # This ensures taxa are present in at least 2 samples

pca3 <- prcomp(stok_s_t_top_var_common,scale.=TRUE,center=TRUE)
#pca3 <- prcomp(stok_s_t_top_var_common) # NOTE - If we do NOT normalize, we get a much better PCA. PC1 takes 34% and PC2 22% of variance

pca3.layout <- cbind(pca3$x[,1:2],stok_meta)

fig.pca3 <- plot_ly(pca3.layout, x = ~PC1, y = ~PC2, color = ~Season, symbol = ~Time, type = 'scatter', mode = 'markers', size=10) %>% # colors = c('#636EFA','#EF553B','#00CC96')
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "PC1 (15%)"),  
    yaxis = list( 
      title = "PC2 (11%)"))

add_segments(fig.pca3,
             x = ~PC1[Pairtype == 1], 
             xend= ~ PC1[Pairtype == 2], 
             y = ~PC2[Pairtype == 1], 
             yend = ~PC2[Pairtype == 2],
             inherit=FALSE,
             color=I("black"),
             line = list(width=1),
             showlegend=FALSE) # %>%

#######################################
#######################################


## EXPERIMENTS WITH UMAP INSTEAD
library("umap")
library("plotly")
#stok_s_t_top.2 <- set_low_frac_to_0_alt(stok_s_t_var, 0.001) # DONT DO THIS ANYMORE OCT 2023
# stok_s_t_top_var.2 <- drop_constant_rows(stok_s_t_top.2)
# stok_s_t_top_var_common.2 <- drop_rare_otus(stok_s_t_top_var.2, 0.05) # This ensures taxa are present in at least 2 samples

#umap.1.default <- umap(stok_s_t_top_var_common.2,n_neighbors=8,min_dist=0.1) 

# 2023-10-18 - THIS IS THE ONE IM CURRENTLY USING. I DONT NEED FURTHER FILTERING DUE TO BRACKEN CUTOFFS AND KITOME BLANKING
# NOTE: There is one winter pair that is clustering with summer samples. This is the DAY samples from winter. Might indicate that this is largely driven by human activity? Winter day = summer at any time?
umap.1.default <- umap(stok_s_t_var %>% filter (!(sample %in% c("NEG-KIT","NEG-WATER","POS"))) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.layout <- umap.1.default[["layout"]] 
umap.layout <- data.frame(umap.layout) 
umap.layout <- cbind(umap.layout, stok_meta %>% filter (!(id %in% c("NEG-KIT","NEG_KIT","NEG-WATER","NEG_WATER","POS")))) 

fig.umap <- plot_ly(umap.layout, x = ~X1, y = ~X2, color=~Season, symbol = ~Year, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 

add_segments(fig.umap,
             x = ~X1[Pairtype == 1], 
             xend= ~ X1[Pairtype == 2], 
             y = ~X2[Pairtype == 1], 
             yend = ~X2[Pairtype == 2],
             inherit=FALSE,
             color=I("black"),
             line = list(width=1),
             showlegend=FALSE)



### STRICTEST POSSIBLE UMAP - NO LONGER RELEVANT OCT 2023

stok_s_t_top <- set_low_frac_to_0_alt(stok_s_t_var, 0.01)
stok_s_t_top_var <- drop_constant_rows(stok_s_t_top)
#stok_s_t_var_filt <- drop_otus_below_count(stok_s_t_var, 500)
stok_s_t_top_var_common <- drop_rare_otus(stok_s_t_top_var, 0.05) # This ensures taxa are present in at least 2 samples

umap.1.default <- umap(stok_s_t_top_var_common,n_neighbors=8,min_dist=0.1) # NOTE - Seems extremely robust to changes in n_neighbors and min_dist
umap.layout <- umap.1.default[["layout"]] 
umap.layout <- data.frame(umap.layout) 
umap.layout <- cbind(umap.layout, stok_meta) 

fig.umap <- plot_ly(umap.layout, x = ~X1, y = ~X2, color = ~Season, symbol = ~Time, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 

add_segments(fig.umap,
             x = ~X1[Pairtype == 1], 
             xend= ~ X1[Pairtype == 2], 
             y = ~X2[Pairtype == 1], 
             yend = ~X2[Pairtype == 2],
             inherit=FALSE,
             color=I("black"),
             line = list(width=1),
             showlegend=FALSE)



#####################################
# DIFFERENCE BETWEEN PAIRS ANALYSIS #
#####################################

#stok_diff <- createDifferenceBetweenPairs(stok_s_t,stok_meta$Pair)
## NOTE - NEED TO NORMALIZE FIRST!
stok_norm <- normalize_TSS_OTU_data(stok_s_t)

stok_diff <- createDifferenceBetweenPairs(stok_norm, stok_meta$Pair)

stok_diff_means <- stok_diff %>%
  summarise(across(where(is.numeric),function(x) mean(x,na.rm=TRUE))) %>%
  unlist
plot(density(stok_diff_means),xlim=c(-1000,1000))
# CLEAR MEAN AROUND 0 AT LEAST

#########################
# BRAY-CURTIS DISTANCES #
#########################

require(vegan)
pairwise_BC_dist <- vegdist(stok_s_t_var %>% select(where(is.numeric)) %>% select_if(colSums(.) != 0), method = "chisq")
#pairwise_rho_dist <- get_dist(stok_norm,"spearman")
#pairwise_euclidian_dist <- get_dist(stok_norm,"euclidian")
# Now, fish out within-pair distances and compare them with between-pair
my_labels=str_remove(stok_meta$Pair,"Par")
my_labels=str_c(str_remove(stok_meta$Pair,"Par"),"_",str_remove(stok_meta$id,"SEA_"))
plot(hclust(pairwise_BC_dist,method="ward.D2"),labels=my_labels,axes=FALSE,ann=FALSE,ylab="",main="Cluster dendrogram of pairwise Bray-Curtis distances")
#require(adegenet)
#tmp <- pairDist(pairwise_BC_dist,grp = as.factor(my_labels),within=TRUE) # Not a very beautiful plot, too many comparisons



# Try creating a sideways violinplot in ggplot2, within-pair and between-pair Bray-Curtis distance
require(usedist)
# USEDIST has functions for this. 
# - Specifically dist_subset(distances, c("A","B","C")) gets distances between only A B and C
# - dist_get(distances, c("A","B"), c("D","E")) gets A-D and B-E distances.

# First, from this 820 (=N*(N-1)/2) long double, need to find numbers that correspond to pairs and numbers that do not correspond to pairs
# Index 1 holds distance between 1 and 2, index 2 between 1 and 3 etc. (First N correspond to the distances from 1 to 2..N). So for example,
# to find the distance between 9 and 16 (i >= j -> i=16, j=9), we would need to find object number (j-1)(N-j/2) + i - j = (9-1)(41-9/2)+16-9= 299 (3.58)

# First, find the within-pairs indexes:
stok_meta_index <- stok_meta %>%
  ungroup() %>%
  mutate(index=1:nrow(.)) %>%
  group_by(Pair) %>%
  filter(n() > 1) %>%
  arrange(Pair) %>%
  ungroup()
indexes <- stok_meta_index$index
indexes_from <- indexes[seq(1,length(indexes),2)]
indexes_to <- indexes[seq(2,length(indexes),2)]

pairwise_BC_dist_noctrl <- vegdist(stok_s_t_var %>% filter (!(sample %in% c("NEG-KIT","NEG-WATER","POS"))) %>% select(where(is.numeric)) %>% select_if(colSums(.) != 0) , method = "chisq")


dists_within <- dist_get(pairwise_BC_dist_noctrl,indexes_from,indexes_to)

# ALTERNATIVE, MANUAL WAY OF FINDING
indexes_withinpair <- find_dist_indexes(pairwise_BC_dist_noctrl,indexes_from,indexes_to)
distances_withinpair <- pairwise_BC_dist_noctrl[indexes_withinpair]
distances_betweenpair <- pairwise_BC_dist_noctrl[-indexes_withinpair]
distances_df <- data.frame("Distances"=c(distances_betweenpair,distances_withinpair),"Type"=c(rep(factor("Between"),length(distances_betweenpair)), rep(factor("Within"),length(distances_withinpair))))

p <- ggplot(distances_df, aes(x=Type,y=Distances)) + geom_violin() # can add fill=Type in aes for colors, but I think it looks better without
p + geom_boxplot(width=0.1)

# SOME BETWEEN-DISTS ARE LESS THAN 2. FIND OUT WHICH
indexes_lowdist <- which(pairwise_BC_dist_noctrl < 2)
low_diff_indexes_in_between <- setdiff(indexes_lowdist,indexes_withinpair)
#   1  39 109 208 514 517 533 534 536 551 552 568 584 677
# Corresponds to the following non-pairs:
# SEA-13, SEA-14 (Par01, Par02)   - 1 (j=1, i = 2)
# SEA-14, SEA-16 (Par01, Par02)   - 39  (j=2, i=5)
# SEA-16, SEA-17 (Par01, Par02)   - 109 (j=4, i=5)
# SEA-31, SEA-32 (Par04, Par05)   - 208 (j=7, i=8)
# SEA-43, SEA-44 (Par10, Par11) - 514 (j=19, i=20)
# SEA-43, SEA-48 (Par10, Par11) - 517 (j=19, i=23)
# SEA-44, SEA-45 (Par11, Par12) - 533 (j=20, i=21)
# SEA-44, SEA-47 (Par11, Par10) - 534 (j=20, i=22)
# SEA-44, SEA-48 (Par10, Par11) - 536 (j=20, i=23)
# SEA-45, SEA-47 (Par12, Par10) - 551 (j=21, i=22)
# SEA-45, SEA-48 (Par12, Par11) - 552 (j=21, i=23)
# SEA-47, SEA-48 (Par10, Par11) - 568 (j=22, i=23)
# SEA-48, SEA-49 (Par11, Par12) - 584 (j=23, i=24)
# SEA-68, SEA-70 (Par16, Par18) - 677 (j=31, i=33)







###########################
# NOTHING TO DO FROM HERE #
###########################


# Also need - Heatmap with all pairwise distances AND clustering tree
# NOTE - OCT 2023 - This approach is ugly. Do dendrograms instead
#pairwise_BC_dist_sym <- vegdist(stok_norm,diag=TRUE)
#hclust_rows <- hclust(pairwise_BC_dist_sym,method="ward.D2")
#hclust_rows$order <- 1:38
hclust_rows <- as.dendrogram(hclust(pairwise_BC_dist,method="ward.D2"))  # Calculate hclust dendrograms
hclust_cols <- as.dendrogram(hclust(t(pairwise_BC_dist),method="ward.D2"))
heatmap(as.matrix(pairwise_BC_dist),symm=TRUE,labRow=my_labels,labCol=my_labels,Rowv=hclust_rows,Colv = hclust_cols) #Rowv=as.factor(my_labels),Colv="Rowv"
# LABELS ARE IDS
my_labels <- stok_meta$id

# CREATE FAKE, "IDEAL" LABELS - WHAT WE WOULD EXPECT IF THERE WERE NO SAMPLE MIXUPS
fake_pairs <- stok_meta$Pair
fake_pairs[7:12] <- c("Par04","Par04","Par05","Par05","Par06","Par06")
fake_pairs[13:18] <- c("Par07","Par07","Par08","Par08","Par09","Par09")
fake_pairs[33] <- "Par16"
fake_pairs[34] <- "Par18"

### 2023 OCT NOT NEEDED
stok_diff_filt <- createDifferenceBetweenPairs(stok_s_t_var_filt_common,stok_meta$Pair)
stok_diff_filt_means <- stok_diff_filt %>%
  summarise(across(where(is.numeric),function(x) mean(x,na.rm=TRUE))) %>%
  unlist
plot(density(stok_diff_filt_means),xlim=c(-1000,1000))


stok_diff_veryfilt <- createDifferenceBetweenPairs(stok_s_t_top_var_common,stok_meta$Pair)
stok_diff_veryfilt_means <- stok_diff_veryfilt %>%
  summarise(across(where(is.numeric),function(x) mean(x,na.rm=TRUE))) %>%
  unlist
plot(density(stok_diff_veryfilt_means),xlim=c(-10000,10000))
