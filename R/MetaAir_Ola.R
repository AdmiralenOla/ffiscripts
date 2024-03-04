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
#BOBCAT <- metaair_meta$ID[!is.na(metaair_meta$EXCLUDE)]
#FILTER <- c(BOBCAT, "SL467228") # <- ADD WEIRD SAMPLES TO EXCLUDE, THIS ONE CLUSTERS WITH NEWYORK DESPITE BEING FROM LONDON AND HAS VERY LOW DNA CONC
FILTER <- metaair_meta$ID[metaair_meta$EXCLUDE] # BOBCAT, WEIRD LONDON SAMPLE, STOCKHOLM SEPT, BUT NOT CONTROLS
UNDERGROUND <- metaair_meta %>% filter(GROUNDLEVEL == "Underground") %>% filter(!ID %in% FILTER) %>% pull(ID)
MISSING <- c("SL310938","SL310939","SL310940","SL310941","SL310942")

# IMPORT METAAIR DATA - PHYLOSEQ OR OLA METHOD

metaair <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE)
#taxa_names <- unlist(metaair[,1])
#metaair <- metaair[,2:ncol(metaair)]
#rownames(metaair) <- taxa_names
#metaair_otu <- otu_table(object = metaair,taxa_are_rows = TRUE)
metaair <- transpose_df(metaair)
colnames(metaair) <- metaair[1,]
metaair <- metaair[2:nrow(metaair),]
colnames(metaair)[1] <- "sample"

metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
translation <- metaair_taxa_names$name[match(colnames(metaair)[2:ncol(metaair)],metaair_taxa_names$taxonomy_id)]
colnames(metaair) <-  c("sample",translation)

metaair_meta <- read_tsv("/media/ubuntu/Pandora/METAAIR/META/Metadata_metaair_final.tsv",col_names = TRUE,col_types = "clffffddfffdddiiddddiididid") # ID Exclude Year Type City Station lat long station groundlevel closedopen temp relhumidrawreads cleanreads percremoved coverage percunclassifiedncbi percunclassifiedfbav2 classifiedreads contamreads contamperc contamreadsbac contampercbac contamreadsfungi contampercfungi
metaair_meta_posremoved <- metaair_meta %>% filter(!(ID %in% POS_CTRLS_METAAIR))
metaair_meta_posremoved_negremoved <- metaair_meta %>% filter(!(ID %in% CTRLS_METAAIR))
metaair_meta_posremoved_negremoved_filtered <- metaair_meta_posremoved_negremoved %>% filter(!EXCLUDE)
# NOTE, UNDERGROUND IS ONLY NEEDED FOR KARIS ANALYSIS
metaair_meta_posremoved_negremoved_filtered_underground <- metaair_meta_posremoved_negremoved_filtered %>% filter(GROUNDLEVEL == "Underground")
# IGNORE WARNINGS - THEY ARE WRONG AND THE DATA IS READ CORRECTLY

# GENERATE SOME MATRIXES THAT ARE USED OFTEN
metaair_norm <- normalize_TSS_OTU_data(metaair)
# TAKING LOG TSS COUNTS INSTEAD OF RAW
metaair_log <- bind_cols(metaair_norm %>% select(sample), metaair_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
metaair_log_var <- drop_constant_cols(metaair_log)

# For METAAIR PAPER, include both above and underground. For Karis purposes, use only underground
metaair_final <- metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
names(metaair_final)[1] <- "ID"

################################
#     ANALYSIS                 #
################################


## FIGURES CONCERNING WHAT'S REMOVED
# DO STACKED BARCHARTS WITH PERC CONTAMINATION
# PER YEAR TOTAL + PER YEAR-CITY

# Need to create long-data first
contam_df <- metaair_meta_posremoved_negremoved_filtered
contam_df$CLASSIFIEDPERC <- 100 - contam_df$CONTAMPERC
contam_df_long <- contam_df %>%
  pivot_longer(cols = c(CLASSIFIEDPERC,CONTAMBACPERC,CONTAMFUNGIPERC),names_to="CONTAMTYPE", values_to = "CONTAMVALUEPERC")

contam_df_plot <- contam_df_long %>% summarise(mean(CONTAMVALUEPERC), .by=c(CONTAMTYPE,CITY,YEAR))
names(contam_df_plot)[4] <- "CONTAMVALUEPERC"
contam_df_plot$CONTAMTYPE[contam_df_plot$CONTAMTYPE == "CLASSIFIEDPERC"] <- "Retained"
contam_df_plot$CONTAMTYPE[contam_df_plot$CONTAMTYPE == "CONTAMBACPERC"] <- "Bacterial contam."
contam_df_plot$CONTAMTYPE[contam_df_plot$CONTAMTYPE == "CONTAMFUNGIPERC"] <- "Fungal contam."
contam_df_plot$CITY <- factor(contam_df_plot$CITY, levels=c("Denver","Hong Kong", "London","New York","Oslo","Stockholm"))
#contam_df_long$GROUPS <- interaction(metaair_meta_posremoved_negremoved_filtered$CITY, metaair_meta_posremoved_negremoved_filtered$YEAR)

ggplot(contam_df_plot, aes(fill=CONTAMTYPE, y=CONTAMVALUEPERC,x=interaction(CITY,YEAR,lex.order=TRUE,drop=TRUE, sep=" ") )) + #x=interaction(CITY,YEAR,lex.order=TRUE)
  geom_bar(position="stack", stat="identity") + xlab("CITY") + ylab("Percentage") + labs(fill="Read was") + 
  theme(axis.text.x = element_text(angle=45))





# BETA-DIVERSITY WORKS GREAT - VERY GOOD SIGNAL FROM CITY
umap.meta <- umap(metaair_final %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.layout <- umap.meta[["layout"]] 
umap.meta.layout <- data.frame(umap.meta.layout, "ID"=metaair_final %>% pull(sample)) 
umap.meta.layout <- inner_join(umap.meta.layout, metaair_meta_posremoved_negremoved_filtered)
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta <- plot_ly(umap.meta.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10,
                         text = ~ID, hovertemplate=paste(
                           "<b>%{text}</b><br><br>",
                           "%{x}, %{y}")) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.meta


## Repeat for bacteria and fungi individually
metaair_bac <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__bacteria.tsv",col_names=TRUE)
metaair_fungi <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__fungi.tsv",col_names=TRUE)

metaair_bac <- transpose_df(metaair_bac)
metaair_fungi <- transpose_df(metaair_fungi)
colnames(metaair_bac) <- metaair_bac[1,]
colnames(metaair_fungi) <- metaair_fungi[1,]
metaair_bac <- metaair_bac[2:nrow(metaair_bac),]
metaair_fungi <- metaair_fungi[2:nrow(metaair_fungi),]

metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
translation_bac <- metaair_taxa_names$name[match(colnames(metaair_bac)[2:ncol(metaair_bac)],metaair_taxa_names$taxonomy_id)]
translation_fungi <- metaair_taxa_names$name[match(colnames(metaair_fungi)[2:ncol(metaair_fungi)],metaair_taxa_names$taxonomy_id)]
colnames(metaair_bac) <-  c("sample",translation)
colnames(metaair_fungi) <-  c("sample",translation)

metaair_bac_norm <- normalize_TSS_OTU_data(metaair_bac)
metaair_bac_log <- bind_cols(metaair_bac_norm %>% select(sample), metaair_bac_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
metaair_bac_log_var <- drop_constant_cols(metaair_bac_log)

metaair_fungi_norm <- normalize_TSS_OTU_data(metaair_fungi %>% filter(sample != "SL470335"))
# NOTE - WEIRD BEHAVIOR FROM SL470335 - Happens after normalization, and only for fungi. This sample has 0 fungi reads. It is also one of the weird Stockholm samples
metaair_fungi_log <- bind_cols(metaair_fungi_norm %>% select(sample), metaair_fungi_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
metaair_fungi_log_var <- drop_constant_cols(metaair_fungi_log)

metaair_bac_final <- metaair_bac_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
metaair_fungi_final <- metaair_fungi_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
names(metaair_bac_final)[1] <- "ID"
names(metaair_fungi_final)[1] <- "ID"

umap.bac.meta <- umap(metaair_bac_final %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.bac.meta.layout <- umap.bac.meta[["layout"]] 
umap.bac.meta.layout <- data.frame(umap.bac.meta.layout, "ID"=metaair_bac_final %>% pull(ID)) # ID?
umap.bac.meta.layout <- inner_join(umap.bac.meta.layout, metaair_meta_posremoved_negremoved_filtered)
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.bac.meta <- plot_ly(umap.bac.meta.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10,
                             text = ~ID, hovertemplate=paste(
                               "<b>%{text}</b><br><br>",
                               "%{x}, %{y}")) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.bac.meta

umap.fungi.meta <- umap(metaair_fungi_final %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.fungi.meta.layout <- umap.fungi.meta[["layout"]] 
umap.fungi.meta.layout <- data.frame(umap.fungi.meta.layout, "ID"=metaair_fungi_final %>% pull(ID)) # ID?
umap.fungi.meta.layout <- inner_join(umap.fungi.meta.layout, metaair_meta_posremoved_negremoved_filtered)
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.fungi.meta <- plot_ly(umap.fungi.meta.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10,
                               text = ~ID, hovertemplate=paste(
                                 "<b>%{text}</b><br><br>",
                                 "%{x}, %{y}")) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')),
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.fungi.meta


## TRY THIS ANALYSIS WITHOUT REMOVING CONTAMINANTS
metaair_w_kitome <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_0.005_.tsv",col_names=TRUE)
metaair_w_kitome <- transpose_df(metaair_w_kitome)
colnames(metaair_w_kitome) <- metaair_w_kitome[1,]
metaair_w_kitome <- metaair_w_kitome[2:nrow(metaair_w_kitome),]

metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
translation_w_kitome <- metaair_taxa_names$name[match(colnames(metaair_w_kitome)[2:ncol(metaair_w_kitome)],metaair_taxa_names$taxonomy_id)]
colnames(metaair_w_kitome) <-  c("sample",translation)

metaair_w_kitome_norm <- normalize_TSS_OTU_data(metaair_w_kitome)
metaair_w_kitome_log <- bind_cols(metaair_w_kitome_norm %>% select(sample), metaair_w_kitome_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
metaair_w_kitome_log_var <- drop_constant_cols(metaair_w_kitome_log)

metaair_w_kitome_final <- metaair_w_kitome_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
names(metaair_w_kitome_final)[1] <- "ID"

umap.kitome.meta <- umap(metaair_w_kitome_final %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.kitome.meta.layout <- umap.kitome.meta[["layout"]] 
umap.kitome.meta.layout <- data.frame(umap.kitome.meta.layout, "ID"=metaair_w_kitome_final %>% pull(ID)) # ID?
umap.kitome.meta.layout <- inner_join(umap.kitome.meta.layout, metaair_meta_posremoved_negremoved_filtered)
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.kitome.meta <- plot_ly(umap.kitome.meta.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10,
                             text = ~ID, hovertemplate=paste(
                               "<b>%{text}</b><br><br>",
                               "%{x}, %{y}")) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.kitome.meta

# SAME ANALYSIS, BUT REMOVE FILTERED AND NON-UNDERGROUND
# ONLY NEEDED FOR KARIS PROJECT
# umap.meta.filt.ug <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter((sample %in% UNDERGROUND)) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
# umap.meta.filt.ug.layout <- umap.meta.filt.ug[["layout"]] 
# umap.meta.filt.ug.layout <- data.frame(umap.meta.filt.ug.layout, "ID"=metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter((sample %in% UNDERGROUND)) %>% pull(sample)) 
# umap.meta.filt.ug.layout <- inner_join(umap.meta.filt.ug.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR)) %>% filter(!(ID %in% FILTER)) %>% filter((ID %in% UNDERGROUND)))
# #umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 
# 
# fig.umap.meta.ug.filt <- plot_ly(umap.meta.filt.ug.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
#   layout(
#     plot_bgcolor = "#e5ecf6",
#     legend=list(title=list(text='Sample type')), 
#     xaxis = list( 
#       title = "0"),  
#     yaxis = list( 
#       title = "1")) 
# fig.umap.meta.ug.filt


### ANALYSIS OF ENVIRONMENTAL VARIABLES ###
#    WONDERING IF I SHOULD DO THIS WITH PCA INSTEAD? UMAP USED MOSTLY FOR PLOTTING

umap.meta.layout$TEMPFACTOR = umap.meta.layout$TEMPERATURE >= 25.0
umap.meta.layout$HUMIDFACTOR = umap.meta.layout$HUMIDITY >= 65.0
umap.meta.layout$OCCUPFACTOR = umap.meta.layout$OCCUPANCY > 100


#umap.model <- manova(cbind(X1, X2) ~ YEAR + CITY + GROUNDLEVEL + OCCUPANCY + TEMPERATURE + HUMIDITY, data=umap.meta.layout)
#summary(umap.model)
# EXPECTING ENORMOUS CORRELATION BETWEEN CITY AND ENVIRONMENTAL VARIABLES. 
# BUILDING MODEL FROM GROUND UP. CITY AND YEAR OBVIOUSLY MATTER, KEEP ADDING WITH INTERACTION
# AND CHECK USING F TEST
umap.base.model <- manova(cbind(X1, X2) ~ YEAR + CITY, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL),])
umap.complex.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL, data=umap.meta.layout)
anova(umap.base.model, umap.complex.model)
# CONCLUSION - CITY*GROUNDLEVEL SHOULD BE ADDED TO THE MODEL

umap.base.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL) & !is.na(umap.meta.layout$TEMPFACTOR),])
umap.complex.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL + TEMPFACTOR + CITY:TEMPFACTOR, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL),])
anova(umap.base.model, umap.complex.model)
# CONCLUSION - TEMPFACTOR SHOULD NOT BE ADDED TO MODEL

umap.base.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL) & !is.na(umap.meta.layout$OCCUPFACTOR),])
umap.complex.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL + OCCUPFACTOR + CITY:OCCUPFACTOR, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL),])
anova(umap.base.model, umap.complex.model)
# CONCLUSION - OCCUPFACTOR SHOULD NOT BE ADDED TO MODEL

umap.base.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL) & !is.na(umap.meta.layout$HUMIDFACTOR),])
umap.complex.model <- manova(cbind(X1, X2) ~ YEAR + CITY*GROUNDLEVEL + HUMIDFACTOR + CITY:HUMIDFACTOR, data=umap.meta.layout[!is.na(umap.meta.layout$GROUNDLEVEL),])
anova(umap.base.model, umap.complex.model)
# CONCLUSION - HUMIDFACTOR CAN BE ADDED TO MODEL. CITY:GROUNDLEVEL IS NO LONGER SIGNIFICANT, MEANING HUMIDITY IS PROBABLY MORE IMPORTANT


# PCA INSTEAD, MORE ROBUST THAN UMAP??
PCA.components <- prcomp(drop_constant_cols(metaair_final) %>% select(where(is.numeric)), scale. = TRUE, center = TRUE)
library("factoextra")
fviz_eig(PCA.components) # Awesome Scree plot
pca.layout <- data.frame(PCA.components$x[,1:10],ID=umap.meta.layout$ID,CITY=umap.meta.layout$CITY, YEAR=umap.meta.layout$YEAR, GROUNDLEVEL=umap.meta.layout$GROUNDLEVEL, TEMPFACTOR=umap.meta.layout$TEMPFACTOR, HUMIDFACTOR=umap.meta.layout$HUMIDFACTOR, OCCUPFACTOR=umap.meta.layout$OCCUPFACTOR)

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY, data=pca.layout[!is.na(pca.layout$GROUNDLEVEL),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# NOTE - INTERACTION TERM NOT SIGNIFICANT
# CONCLUSION - GROUNDLEVEL SHOULD BE ADDED

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$TEMPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + TEMPFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - TEMPFACTOR SIGNIFICANCE NOT LOWER THAN 1e-3 (0.025)

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$HUMIDFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + HUMIDFACTOR + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - HUMIDFACTOR NOT SIGNIFICANT ENOUGH (0.021)

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$OCCUPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + OCCUPFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - OCCUPFACTOR NOT SIGNIFICANT ENOUGH (0.036)


# REPEAT ANALYSIS FOR BACTERIA AND FUNGI SEPARATELY
dat = metaair_bac_final
names(dat)[1] <- "sample"
PCA.components <- prcomp(drop_constant_cols(dat) %>% select(where(is.numeric)), scale. = TRUE, center = TRUE)
pca.layout <- data.frame(PCA.components$x[,1:10],ID=umap.meta.layout$ID,CITY=umap.meta.layout$CITY, YEAR=umap.meta.layout$YEAR, GROUNDLEVEL=umap.meta.layout$GROUNDLEVEL, TEMPFACTOR=umap.meta.layout$TEMPFACTOR, HUMIDFACTOR=umap.meta.layout$HUMIDFACTOR, OCCUPFACTOR=umap.meta.layout$OCCUPFACTOR)

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY, data=pca.layout[!is.na(pca.layout$GROUNDLEVEL),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - GROUNDLEVEL SHOULD BE INCLUDED FOR BAC

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$TEMPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + TEMPFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - TEMP SHOULD NOT BE INCLUDED FOR BAC

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$OCCUPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + OCCUPFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - OCCUPANCY SHOULD NOT BE INCLUDED FOR BAC

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$HUMIDFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + HUMIDFACTOR + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - HUMIDFACTOR SHOULD! BE INCLUDED FOR BAC

dat = metaair_fungi_final
names(dat)[1] <- "sample"
PCA.components <- prcomp(drop_constant_cols(dat) %>% select(where(is.numeric)), scale. = TRUE, center = TRUE)
pca.layout <- data.frame(PCA.components$x[,1:10],ID=umap.meta.layout$ID,CITY=umap.meta.layout$CITY, YEAR=umap.meta.layout$YEAR, GROUNDLEVEL=umap.meta.layout$GROUNDLEVEL, TEMPFACTOR=umap.meta.layout$TEMPFACTOR, HUMIDFACTOR=umap.meta.layout$HUMIDFACTOR, OCCUPFACTOR=umap.meta.layout$OCCUPFACTOR)

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY, data=pca.layout[!is.na(pca.layout$GROUNDLEVEL),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - GROUNDLEVEL SHOULD BE INCLUDED FOR FUNGI

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL, data=pca.layout[!is.na(pca.layout$HUMIDFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + HUMIDFACTOR + GROUNDLEVEL, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - HUMIDFACTOR SHOULD! BE INCLUDED FOR FUNGI

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + HUMIDFACTOR, data=pca.layout[!is.na(pca.layout$TEMPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + HUMIDFACTOR + TEMPFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - TEMP SHOULD BE INCLUDED FOR FUNGI

pca.base.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + HUMIDFACTOR + TEMPFACTOR + CITY:GROUNDLEVEL + CITY:HUMIDFACTOR, data=pca.layout[!is.na(pca.layout$OCCUPFACTOR),])
pca.complex.model <- manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ YEAR + CITY + GROUNDLEVEL + HUMIDFACTOR + TEMPFACTOR + OCCUPFACTOR + CITY:GROUNDLEVEL + CITY:HUMIDFACTOR, data=pca.layout)
anova(pca.base.model, pca.complex.model)
# CONCLUSION - OCCUPANCY SHOULD BE INCLUDED FOR FUNGI
# ALL FACTORS RELEVANT FOR FUNGI!



# PLOT OF PCA
fig.pca <- plot_ly(pca.layout, x = ~PC1, y = ~PC2, color = ~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% # colors = c('#636EFA','#EF553B','#00CC96')
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "PC1 (15%)"),  
    yaxis = list( 
      title = "PC2 (4%)"))

# REPEAT PER CITY AND INVESTIGATE ENVIRONMENTAL FACTORS
DENVER <- metaair_meta %>% filter(CITY == "Denver") %>% pull("ID")
HONGKONG <- metaair_meta %>% filter(CITY == "Hong Kong") %>% pull("ID")
LONDON <- metaair_meta %>% filter(CITY == "London") %>% pull("ID")
NEWYORK <- metaair_meta %>% filter(CITY == "New York") %>% pull("ID")
OSLO <- metaair_meta %>% filter(CITY == "Oslo") %>% pull("ID")
STOCKHOLM <- metaair_meta %>% filter(CITY == "Stockholm") %>% pull("ID")

## DENVER - CLOSED/OPEN

umap.meta.filt.ug.denver <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% DENVER) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.filt.ug.denver.layout <- umap.meta.filt.ug.denver[["layout"]] 
umap.meta.filt.ug.denver.layout <- data.frame(umap.meta.filt.ug.denver.layout, "ID"=metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% DENVER) %>% pull(sample)) 
umap.meta.filt.ug.denver.layout <- inner_join(umap.meta.filt.ug.denver.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR)) %>% filter(!(ID %in% FILTER)) %>% filter(ID %in% UNDERGROUND) %>% filter(ID %in% DENVER))
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta.ug.denver.filt <- plot_ly(umap.meta.filt.ug.denver.layout, x = ~X1, y = ~X2, color=~`CLOSED/OPEN`, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.meta.ug.denver.filt

## HONGKONG - OCCUPANCY, TEMP, HUMIDITY
umap.meta.filt.ug.hongkong <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% HONGKONG) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.filt.ug.hongkong.layout <- umap.meta.filt.ug.hongkong[["layout"]] 
umap.meta.filt.ug.hongkong.layout <- data.frame(umap.meta.filt.ug.hongkong.layout, "ID"=metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% HONGKONG) %>% pull(sample)) 
umap.meta.filt.ug.hongkong.layout <- inner_join(umap.meta.filt.ug.hongkong.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR)) %>% filter(!(ID %in% FILTER)) %>% filter(ID %in% UNDERGROUND) %>% filter(ID %in% HONGKONG))
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta.ug.hongkong.filt <- plot_ly(umap.meta.filt.ug.hongkong.layout, x = ~X1, y = ~X2, color=~`OCCUPANCY`, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.meta.ug.hongkong.filt

## NEWYORK - CLOSED/OPEN, TEMP, HUMIDITY
umap.meta.filt.ug.newyork <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% NEWYORK) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.filt.ug.newyork.layout <- umap.meta.filt.ug.newyork[["layout"]] 
umap.meta.filt.ug.newyork.layout <- data.frame(umap.meta.filt.ug.newyork.layout, "ID"=metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% NEWYORK) %>% pull(sample)) 
umap.meta.filt.ug.newyork.layout <- inner_join(umap.meta.filt.ug.newyork.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR)) %>% filter(!(ID %in% FILTER)) %>% filter(ID %in% UNDERGROUND) %>% filter(ID %in% NEWYORK))
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta.ug.newyork.filt <- plot_ly(umap.meta.filt.ug.newyork.layout, x = ~X1, y = ~X2, color=~`Relative humidity`, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.meta.ug.newyork.filt

## OSLO - CLOSED/OPEN, TEMP, HUMIDITY, OCCUPANCY
umap.meta.filt.ug.oslo <- umap(metaair_log_var %>% filter (!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% OSLO) %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
umap.meta.filt.ug.oslo.layout <- umap.meta.filt.ug.oslo[["layout"]] 
umap.meta.filt.ug.oslo.layout <- data.frame(umap.meta.filt.ug.oslo.layout, "ID"=metaair_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER)) %>% filter(sample %in% UNDERGROUND) %>% filter(sample %in% OSLO) %>% pull(sample)) 
umap.meta.filt.ug.oslo.layout <- inner_join(umap.meta.filt.ug.oslo.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR)) %>% filter(!(ID %in% FILTER)) %>% filter(ID %in% UNDERGROUND) %>% filter(ID %in% OSLO))
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

fig.umap.meta.ug.oslo.filt <- plot_ly(umap.meta.filt.ug.oslo.layout, x = ~X1, y = ~X2, color=~`CLOSED/OPEN`, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10) %>% 
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Sample type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig.umap.meta.ug.oslo.filt

# ALPHA-DIVERSITY - JUST ADD SAMPLE DATA TO METADATA

library(vegan)
metaair_final_ <- metaair_final
names(metaair_final_)[1] <- "ID"
metaair_final_meta_and_data <- inner_join(metaair_final_,metaair_meta_posremoved_negremoved_filtered,by="ID")

metaair_final_meta_and_data$shannondiv <- diversity(metaair_final_meta_and_data %>% select(2:3738)) #%>% filter (!(sample %in% CTRLS_METAAIR)) %>% select_if(is.numeric)
ggplot(metaair_final_meta_and_data, aes(x = CITY, y = shannondiv,color=YEAR)) + 
  geom_boxplot() + ylab("Shannon diversity")

# PERFORM TUKEYS HONEST SIGNIFICANT DIFFERENCE TEST
library(agricolae)
level_order_city <- c("Denver", "Hong Kong", "London", "New York", "Oslo", "Stockholm")
anova_result <- aov(shannondiv ~ CITY, metaair_final_meta_and_data)
tukey_result <- HSD.test(anova_result, "CITY", group = TRUE)
group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggplot(metaair_final_meta_and_data, aes(x = factor(CITY,levels=level_order_city), y = shannondiv, color=YEAR)) + # ADD COLOR YEAR?
  geom_text(data = data.frame(),
            aes(x = rownames(group_data), y = max(metaair_final_meta_and_data$shannondiv) + 0.5, label = group_data$groups),
            col = 'black',
            size = 5) +
  geom_boxplot() +
  ggtitle("Shannon diversity by City") +
  xlab("City") +
  ylab("Shannon diversity index")

# SAME ANALYSIS, BACTERIA ONLY
metaair_final_meta_and_data_bac <- inner_join(metaair_bac_final,metaair_meta_posremoved_negremoved_filtered,by="ID")

metaair_final_meta_and_data_bac$shannondiv <- diversity(metaair_final_meta_and_data_bac %>% select(2:2661)) #%>% filter (!(sample %in% CTRLS_METAAIR)) %>% select_if(is.numeric)
anova_result_bac <- aov(shannondiv ~ CITY, metaair_final_meta_and_data_bac)
tukey_result_bac <- HSD.test(anova_result_bac, "CITY", group = TRUE)
group_data_bac <- tukey_result_bac$groups[order(rownames(tukey_result_bac$groups)),]
ggplot(metaair_final_meta_and_data_bac, aes(x = factor(CITY,levels=level_order_city), y = shannondiv,color=YEAR)) + 
  geom_boxplot() + ylab("Shannon diversity index") + xlab("City") +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_bac), y = max(metaair_final_meta_and_data_bac$shannondiv) + 0.5, label = group_data_bac$groups),
            col = 'black',
            size = 5)

# SAME ANALYSIS, FUNGI ONLY
metaair_final_meta_and_data_fungi <- inner_join(metaair_fungi_final,metaair_meta_posremoved_negremoved_filtered,by="ID")

metaair_final_meta_and_data_fungi$shannondiv <- diversity(metaair_final_meta_and_data_fungi %>% select(2:1011)) #%>% filter (!(sample %in% CTRLS_METAAIR)) %>% select_if(is.numeric)
anova_result_fungi <- aov(shannondiv ~ CITY, metaair_final_meta_and_data_fungi)
tukey_result_fungi <- HSD.test(anova_result_fungi, "CITY", group = TRUE)
group_data_fungi <- tukey_result_fungi$groups[order(rownames(tukey_result_fungi$groups)),]
ggplot(metaair_final_meta_and_data_fungi, aes(x = factor(CITY,levels=level_order_city), y = shannondiv,color=YEAR)) + 
  geom_boxplot() + ylab("Shannon diversity index") + xlab("City") +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_bac), y = max(metaair_final_meta_and_data_bac$shannondiv) + 0.5, label = group_data_bac$groups),
            col = 'black',
            size = 5)

# WITH KITOME
metaair_final_meta_and_data_kitome <- inner_join(metaair_w_kitome_final,metaair_meta_posremoved_negremoved_filtered,by="ID")
metaair_final_meta_and_data_kitome$shannondiv <- diversity(metaair_final_meta_and_data_kitome %>% select(2:4036)) #%>% filter (!(sample %in% CTRLS_METAAIR)) %>% select_if(is.numeric)
anova_result_kitome <- aov(shannondiv ~ CITY, metaair_final_meta_and_data_kitome)
tukey_result_kitome <- HSD.test(anova_result_kitome, "CITY", group = TRUE)
group_data_kitome <- tukey_result_kitome$groups[order(rownames(tukey_result_kitome$groups)),]
ggplot(metaair_final_meta_and_data_kitome, aes(x = factor(CITY,levels=level_order_city), y = shannondiv,color=YEAR)) + 
  geom_boxplot() + ylab("Shannon diversity index") + xlab("City") +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_kitome), y = max(metaair_final_meta_and_data_kitome$shannondiv) + 0.5, label = group_data_kitome$groups),
            col = 'black',
            size = 5)

## STARTING TO WORK ON CORE DEFINITIONS

# IDEA:
# Function must include minimum threshold AVERAGE abundance (e.g. 0.0005), should be possible to adjust
# Function should allow you to not include contaminants
# Function should allow you to cut species where prevalence in negatives is similar to prevalence in samples. (Proportional test for difference?)
# Need a step function for prevalence to calculate core at different prevalences. Prevalence should be on Y axis
# The number of species with prevalence >= than Y is plotted on X axis (OR - reverse axes so that response-outcome is more logical)
# Output should simply be a data frame of x and y coordinates?

filter_core <- function(df,includeContaminants=TRUE,includePankitome=TRUE,minimumAverageAbundance=0.00005,differentialAbundanceFromNegs=TRUE,significance=0.05,N_samples=787,N_negatives=22){
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
    #X_samples <- ceiling(df$frequency_perc * N_samples)
    X_samples <- ceiling(df$True_prevalence * N_samples)
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
   return(length(which(combinations[1] < df$True_prevalence & combinations[2] < df$counts_perc_avg)))
}

# NOTE 2024 - NEED NEW GRIMER FILES SINCE SOME SAMPLES HAVE BEEN EXCLUDED.

metaair_metadata_final <- metaair_meta_posremoved_negremoved_filtered %>% filter(ID %in% metaair_final$sample)
# NOTE - GRIMER NEEDS FIRST COLUMN AS ID, NOT EXCLUDE. THIS WONT WORK WITHOUT MANUAL WORK
write_tsv(metaair_metadata_final,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/metadata_final.tsv")

metaair_counts <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE)
metaair_counts <- metaair_counts %>% select(c(taxonomy_id,metaair_metadata_final$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_final.tsv")

# NOTE - NEED TO CREATE COUNT FILES FOR BACTERIA AND FUNGI TOO
metaair_counts_bac <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__bacteria.tsv",col_names=TRUE)
metaair_counts_bac <- metaair_counts_bac %>% select(c(taxonomy_id,metaair_metadata_final$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_bac,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final.tsv")

metaair_counts_fungi <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__fungi.tsv",col_names=TRUE)
metaair_counts_fungi <- metaair_counts_fungi %>% select(c(taxonomy_id,metaair_metadata_final$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_fungi,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final.tsv")



## COMBINED ##
core_mat_combined <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_ALL.tsv",col_names = TRUE)
core_mat_combined_filtered <- filter_core(core_mat_combined,includeContaminants = FALSE)
core_mat_combined_numbers <- calculate_core(core_mat_combined_filtered)
ggplot(data=core_mat_combined_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

## BACTERIA ##
core_mat_bacteria <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC.tsv",col_names = TRUE)
core_mat_bacteria_filtered <- filter_core(core_mat_bacteria,includeContaminants = FALSE)
core_mat_bacteria_numbers <- calculate_core(core_mat_bacteria_filtered)
ggplot(data=core_mat_bacteria_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

## FUNGI ##
core_mat_fungi <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI.tsv",col_names = TRUE)
core_mat_fungi_filtered <- filter_core(core_mat_fungi,includeContaminants = FALSE)
core_mat_fungi_numbers <- calculate_core(core_mat_fungi_filtered)
ggplot(data=core_mat_fungi_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
# 

# TO CONTINUE - MAKING CITY_SPECIFIC CORE FIGURES. FIRST MAKE FILES TO GRIMER, THEN EXPORT TABLE, THEN USE ABOVE FUNCTIONS
# WE CAN ALSO USE THESE GRIMER FILES FOR OUR TOP 20 ANALYSES
# NOTE - NEED TO CREATE COUNT FILES FOR BACTERIA AND FUNGI TOO

# UNCOMMENT TO RE-READ
#metaair_counts_bac <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final.tsv")
#metaair_counts_fungi <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final.tsv")

DENVER <- metaair_metadata_final %>% filter(CITY == "Denver") %>% pull("ID")
HONGKONG <- metaair_metadata_final %>% filter(CITY == "Hong Kong") %>% pull("ID")
LONDON <- metaair_metadata_final %>% filter(CITY == "London") %>% pull("ID")
NEWYORK <- metaair_metadata_final %>% filter(CITY == "New York") %>% pull("ID")
OSLO <- metaair_metadata_final %>% filter(CITY == "Oslo") %>% pull("ID")
STOCKHOLM <- metaair_metadata_final %>% filter(CITY == "Stockholm") %>% pull("ID")

# DENVER
metaair_counts_all_denver <- metaair_counts %>% select(c(taxonomy_id,DENVER,NEG_CTRLS_METAAIR))
metaair_counts_bac_denver <- metaair_counts_bac %>% select(c(taxonomy_id,DENVER,NEG_CTRLS_METAAIR))
metaair_counts_fungi_denver <- metaair_counts_fungi %>% select(c(taxonomy_id,DENVER,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_denver,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_denver.tsv")
write_tsv(metaair_counts_bac_denver,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_denver.tsv")
write_tsv(metaair_counts_fungi_denver,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_denver.tsv")

core_mat_all_denver <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_DENVER.tsv",col_names = TRUE)
core_mat_all_denver_filtered <- filter_core(core_mat_all_denver,includeContaminants = FALSE)
core_mat_all_denver_numbers <- calculate_core(core_mat_all_denver_filtered)
ggplot(data=core_mat_all_denver_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Denver (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

core_mat_bac_denver <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_DENVER.tsv",col_names = TRUE)
core_mat_bac_denver_filtered <- filter_core(core_mat_bac_denver,includeContaminants = FALSE,N_samples = 38,N_negatives = 22)
core_mat_bac_denver_numbers <- calculate_core(core_mat_bac_denver_filtered)
ggplot(data=core_mat_bac_denver_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Denver (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100))

core_mat_fungi_denver <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_DENVER.tsv",col_names = TRUE)
core_mat_fungi_denver_filtered <- filter_core(core_mat_fungi_denver,includeContaminants = FALSE,N_samples = 38,N_negatives = 21)
core_mat_fungi_denver_numbers <- calculate_core(core_mat_fungi_denver_filtered)
ggplot(data=core_mat_fungi_denver_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Denver (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
# 
# 
# 

# HONG KONG
metaair_counts_all_hk <- metaair_counts %>% select(c(taxonomy_id,HONGKONG,NEG_CTRLS_METAAIR))
metaair_counts_bac_hk <- metaair_counts_bac %>% select(c(taxonomy_id,HONGKONG,NEG_CTRLS_METAAIR))
metaair_counts_fungi_hk <- metaair_counts_fungi %>% select(c(taxonomy_id,HONGKONG,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_hk,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_hk.tsv")
write_tsv(metaair_counts_bac_hk,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_hk.tsv")
write_tsv(metaair_counts_fungi_hk,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_hk.tsv")

core_mat_all_hk <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_HK.tsv",col_names = TRUE)
core_mat_all_hk_filtered <- filter_core(core_mat_all_hk,includeContaminants = FALSE)
core_mat_all_hk_numbers <- calculate_core(core_mat_all_hk_filtered)
ggplot(data=core_mat_all_hk_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Hong Kong (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_hk <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_HK.tsv",col_names = TRUE)
core_mat_bac_hk_filtered <- filter_core(core_mat_bac_hk,includeContaminants = FALSE)
core_mat_bac_hk_numbers <- calculate_core(core_mat_bac_hk_filtered)
ggplot(data=core_mat_bac_hk_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Hong Kong (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_hk <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_HK.tsv",col_names = TRUE)
core_mat_fungi_hk_filtered <- filter_core(core_mat_fungi_hk,includeContaminants = FALSE)
core_mat_fungi_hk_numbers <- calculate_core(core_mat_fungi_hk_filtered)
ggplot(data=core_mat_fungi_hk_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Hong Kong (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))


# LONDON
metaair_counts_all_london <- metaair_counts %>% select(c(taxonomy_id,LONDON,NEG_CTRLS_METAAIR))
metaair_counts_bac_london <- metaair_counts_bac %>% select(c(taxonomy_id,LONDON,NEG_CTRLS_METAAIR))
metaair_counts_fungi_london <- metaair_counts_fungi %>% select(c(taxonomy_id,LONDON,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_london,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_london.tsv")
write_tsv(metaair_counts_bac_london,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_london.tsv")
write_tsv(metaair_counts_fungi_london,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_london.tsv")

core_mat_all_london <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_LONDON.tsv",col_names = TRUE)
core_mat_all_london_filtered <- filter_core(core_mat_all_london,includeContaminants = FALSE)
core_mat_all_london_numbers <- calculate_core(core_mat_all_london_filtered)
ggplot(data=core_mat_all_london_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species London (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_london <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_LONDON.tsv",col_names = TRUE)
core_mat_bac_london_filtered <- filter_core(core_mat_bac_london,includeContaminants = FALSE)
core_mat_bac_london_numbers <- calculate_core(core_mat_bac_london_filtered)
ggplot(data=core_mat_bac_london_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species London (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_london <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_LONDON.tsv",col_names = TRUE)
core_mat_fungi_london_filtered <- filter_core(core_mat_fungi_london,includeContaminants = FALSE)
core_mat_fungi_london_numbers <- calculate_core(core_mat_fungi_london_filtered)
ggplot(data=core_mat_fungi_london_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species London (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28))


# NEW YORK
metaair_counts_all_ny <- metaair_counts %>% select(c(taxonomy_id,NEWYORK,NEG_CTRLS_METAAIR))
metaair_counts_bac_ny <- metaair_counts_bac %>% select(c(taxonomy_id,NEWYORK,NEG_CTRLS_METAAIR))
metaair_counts_fungi_ny <- metaair_counts_fungi %>% select(c(taxonomy_id,NEWYORK,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_ny,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_ny.tsv")
write_tsv(metaair_counts_bac_ny,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_ny.tsv")
write_tsv(metaair_counts_fungi_ny,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_ny.tsv")

core_mat_all_ny <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_NY.tsv",col_names = TRUE)
core_mat_all_ny_filtered <- filter_core(core_mat_all_ny,includeContaminants = FALSE)
core_mat_all_ny_numbers <- calculate_core(core_mat_all_ny_filtered)
ggplot(data=core_mat_all_ny_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species New York (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_ny <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_NY.tsv",col_names = TRUE)
core_mat_bac_ny_filtered <- filter_core(core_mat_bac_ny,includeContaminants = FALSE)
core_mat_bac_ny_numbers <- calculate_core(core_mat_bac_ny_filtered)
ggplot(data=core_mat_bac_ny_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species New York (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_ny <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_NY.tsv",col_names = TRUE)
core_mat_fungi_ny_filtered <- filter_core(core_mat_fungi_ny,includeContaminants = FALSE)
core_mat_fungi_ny_numbers <- calculate_core(core_mat_fungi_ny_filtered)
ggplot(data=core_mat_fungi_ny_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species New York (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28))



# OSLO
metaair_counts_all_oslo <- metaair_counts %>% select(c(taxonomy_id,OSLO,NEG_CTRLS_METAAIR))
metaair_counts_bac_oslo <- metaair_counts_bac %>% select(c(taxonomy_id,OSLO,NEG_CTRLS_METAAIR))
metaair_counts_fungi_oslo <- metaair_counts_fungi %>% select(c(taxonomy_id,OSLO,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_oslo,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_oslo.tsv")
write_tsv(metaair_counts_bac_oslo,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_oslo.tsv")
write_tsv(metaair_counts_fungi_oslo,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_oslo.tsv")

core_mat_all_oslo <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_OSLO.tsv",col_names = TRUE)
core_mat_all_oslo_filtered <- filter_core(core_mat_all_oslo,includeContaminants = FALSE)
core_mat_all_oslo_numbers <- calculate_core(core_mat_all_oslo_filtered)
ggplot(data=core_mat_all_oslo_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Oslo (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_oslo <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_OSLO.tsv",col_names = TRUE)
core_mat_bac_oslo_filtered <- filter_core(core_mat_bac_oslo,includeContaminants = FALSE)
core_mat_bac_oslo_numbers <- calculate_core(core_mat_bac_oslo_filtered)
ggplot(data=core_mat_bac_oslo_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Oslo (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_oslo <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_OSLO.tsv",col_names = TRUE)
core_mat_fungi_oslo_filtered <- filter_core(core_mat_fungi_oslo,includeContaminants = FALSE)
core_mat_fungi_oslo_numbers <- calculate_core(core_mat_fungi_oslo_filtered)
ggplot(data=core_mat_fungi_oslo_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Oslo (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28))


# STOCKHOLM
metaair_counts_all_stockholm <- metaair_counts %>% select(c(taxonomy_id,STOCKHOLM,NEG_CTRLS_METAAIR))
metaair_counts_bac_stockholm <- metaair_counts_bac %>% select(c(taxonomy_id,STOCKHOLM,NEG_CTRLS_METAAIR))
metaair_counts_fungi_stockholm <- metaair_counts_fungi %>% select(c(taxonomy_id,STOCKHOLM,NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_all_stockholm,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_all_final_stockholm.tsv")
write_tsv(metaair_counts_bac_stockholm,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_bac_final_stockholm.tsv")
write_tsv(metaair_counts_fungi_stockholm,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/counts_fungi_final_stockholm.tsv")

core_mat_all_stockholm <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-19_GRIMER_COUNTS_ALL_STOCKHOLM.tsv",col_names = TRUE)
core_mat_all_stockholm_filtered <- filter_core(core_mat_all_stockholm,includeContaminants = FALSE)
core_mat_all_stockholm_numbers <- calculate_core(core_mat_all_stockholm_filtered)
ggplot(data=core_mat_all_stockholm_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Stockholm (Bacteria and fungi)") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_stockholm <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_BAC_STOCKHOLM.tsv",col_names = TRUE)
core_mat_bac_stockholm_filtered <- filter_core(core_mat_bac_stockholm,includeContaminants = FALSE)
core_mat_bac_stockholm_numbers <- calculate_core(core_mat_bac_stockholm_filtered)
ggplot(data=core_mat_bac_stockholm_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Stockholm (Bacteria only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_stockholm <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS/2024-01-16_GRIMER_COUNTS_FUNGI_STOCKHOLM.tsv",col_names = TRUE)
core_mat_fungi_stockholm_filtered <- filter_core(core_mat_fungi_stockholm,includeContaminants = FALSE)
core_mat_fungi_stockholm_numbers <- calculate_core(core_mat_fungi_stockholm_filtered)
ggplot(data=core_mat_fungi_stockholm_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species Stockholm (Fungi only)") + 
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28))



###################################
# REPEATING WITH NON_PAREIL > 0.5 #
###################################

metaair_metadata_final_highcov <- metaair_metadata_final %>% filter(COVERAGE >= 0.5)

metaair_counts <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE)
metaair_counts_highcov <- metaair_counts %>% select(c(taxonomy_id,metaair_metadata_final_highcov$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_highcov,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/COUNTS_HIGHCOV/counts_final.tsv")

# NOTE - NEED TO CREATE COUNT FILES FOR BACTERIA AND FUNGI TOO
metaair_counts_bac <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__bacteria.tsv",col_names=TRUE)
metaair_counts_bac_highcov <- metaair_counts_bac %>% select(c(taxonomy_id,metaair_metadata_final_highcov$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_bac_highcov,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/COUNTS_HIGHCOV/counts_bac_final.tsv")

metaair_counts_fungi <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED__fungi.tsv",col_names=TRUE)
metaair_counts_fungi_highcov <- metaair_counts_fungi %>% select(c(taxonomy_id,metaair_metadata_final_highcov$ID, NEG_CTRLS_METAAIR))
write_tsv(metaair_counts_fungi_highcov,file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/COUNTS_HIGHCOV/counts_fungi_final.tsv")

core_mat_all_highcov <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS_HIGHCOV/2024-01-22_GRIMER_COUNTS_HIGHCOV_ALL.tsv",col_names = TRUE)
core_mat_all_highcov_filtered <- filter_core(core_mat_all_highcov,includeContaminants = FALSE)
core_mat_all_highcov_numbers <- calculate_core(core_mat_all_highcov_filtered)
ggplot(data=core_mat_all_highcov_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species, Non-pareil >= 50%") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_bac_highcov <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS_HIGHCOV/2024-01-22_GRIMER_COUNTS_HIGHCOV_BAC.tsv",col_names = TRUE)
core_mat_bac_highcov_filtered <- filter_core(core_mat_bac_highcov,includeContaminants = FALSE)
core_mat_bac_highcov_numbers <- calculate_core(core_mat_bac_highcov_filtered)
ggplot(data=core_mat_bac_highcov_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Bacteria only), Non-pareil >= 50%") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120,140,160,180,200))

core_mat_fungi_highcov <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/GRIMER_FINAL/GRIMER_COUNTS_HIGHCOV/2024-01-22_GRIMER_COUNTS_HIGHCOV_FUNGI.tsv",col_names = TRUE)
core_mat_fungi_highcov_filtered <- filter_core(core_mat_fungi_highcov,includeContaminants = FALSE)
core_mat_fungi_highcov_numbers <- calculate_core(core_mat_fungi_highcov_filtered)
ggplot(data=core_mat_fungi_highcov_numbers, aes(x=X, y=Y*100, group=Z)) + # group=1
  geom_line(aes(color=Z)) +
  geom_point(aes(color=Z)) + xlab("Number of species above this threshold") + ylab("Prevalence (%)") + ggtitle("Core species (Fungi only), Non-pareil >= 50%") +
  labs(color="Abundance") + scale_color_hue(labels = c("0.005%","0.300%", "0.500%", "1.000%")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28))

###########################################
###         TOP SPECIES                 ###
###########################################

########## COLORS ###################################
my.colors <- c('#999999','#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')

########## FUNCTIONS ################################
getTopSpecies <- function(df, N, metadata, includeNegsKit=TRUE,includeNegsWater=TRUE,includeContams=TRUE,includeHuman=TRUE, includeDecontam = TRUE, report_percentage=TRUE){
  exclude_columns <- c()
  # NOTE - I SHOULDNT HAVE TO CARE ABOUT THESE FILTERS ANYMORE, USE FILES THAT ARE ALREADY FILTERED
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
    mutate(rowSums = select(.,-ID) %>% rowSums(.))
  
  top_species_group <- df %>%
    select(-all_of(exclude_columns)) %>% # Select only columns to be included - VERIFY that this works!
    filter(ID %in% metadata) %>%
    summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>% #.names = "{.col}_Percentage"
    select(-rowSums) %>%
    unlist() %>%
    sort(decreasing=TRUE) %>%
    head(N)
  species_to_get_data_for <- names(top_species_group)
  top_species_group <- addOtherSpecies(top_species_group) # <- Drop this if we don't want the "other" category
  results <- data.frame(Fractions=top_species_group, Species=names(top_species_group),Sample="Group")
  
  for(sample in metadata$ID){
    top_species <- df %>%
      # select(-one_of(exclude_columns)) %>% No longer necessary. Rowsums calculated
      filter(ID == sample) %>%
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
  
  # NOTE - The following code is only relevant if we have contaminants included, which we dont.
  #contamIndex <- which(df$Sample == "Is_contaminant")
  #contamrow <- df[contamIndex,]
  #contamrow$Other <- 0
  #Legend.color <- rep('#000000',nrow(results))
  #for(i in 1:nrow(results)){
    #Species <- row
  #  Lookup <- results$Species[i]
  #  Lookupvalue <- contamrow %>% select(any_of(Lookup))
  #  Legend.color[i] <- ifelse(as.logical(Lookupvalue),'#FF0000','#000000')
  #}
  #results$Legend.color <- Legend.color
  
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

TOP20_BAC <- getTopSpecies(metaair_bac_final,20,metaair_metadata_final)


#################
# CROSS-KINGDOM #
#################

# First, make files:
data_crosskingdom <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE) %>% select(c(taxonomy_id,metaair_metadata_final$ID))
data_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017]))
data_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018]))
data_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019]))
data_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID))

write_tsv(data_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/data_2017.tsv")
write_tsv(data_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/data_2018.tsv")
write_tsv(data_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/data_2019.tsv")
write_tsv(data_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/data_all.tsv")

# I have the data stores, but find it easier to just manually type in data here

ck_year <- c(rep("Total",4),rep("2017",4),rep("2018",4),rep("2019",4))
ck_kingdom <- rep(c("Bacteria","Fungi","Archaea","Virus"),4)
ck_perc <- c(74.4,25.3,0.3,0.0, 86.6,13.1,0.3,0.0, 70.7,28.9,0.4,0.0, 65.6,34.1,0.3,0.0)
ck_df <- data.frame(year=ck_year,kingdom=ck_kingdom,perc=ck_perc)

ggplot(ck_df, aes(fill=ck_kingdom, y=ck_perc, x=ck_year)) + 
  geom_bar(position="stack", stat="identity") + xlab("Year") + ylab("Percentage") + labs(fill="Kingdom")



# REPEAT PER CITY
denver_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "Denver"]))
denver_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "Denver"]))
denver_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "Denver"]))
denver_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "Denver"]))
write_tsv(denver_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/denver_2017.tsv")
write_tsv(denver_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/denver_2018.tsv")
write_tsv(denver_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/denver_2019.tsv")
write_tsv(denver_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/denver_all.tsv")

hk_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "Hong Kong"]))
hk_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "Hong Kong"]))
hk_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "Hong Kong"]))
hk_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "Hong Kong"]))
write_tsv(hk_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/hk_2017.tsv")
write_tsv(hk_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/hk_2018.tsv")
write_tsv(hk_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/hk_2019.tsv")
write_tsv(hk_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/hk_all.tsv")

london_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "London"]))
london_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "London"]))
london_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "London"]))
london_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "London"]))
write_tsv(london_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/london_2017.tsv")
write_tsv(london_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/london_2018.tsv")
write_tsv(london_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/london_2019.tsv")
write_tsv(london_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/london_all.tsv")

ny_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "New York"]))
ny_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "New York"]))
ny_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "New York"]))
ny_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "New York"]))
write_tsv(ny_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/ny_2017.tsv")
write_tsv(ny_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/ny_2018.tsv")
write_tsv(ny_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/ny_2019.tsv")
write_tsv(ny_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/ny_all.tsv")

oslo_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "Oslo"]))
oslo_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "Oslo"]))
oslo_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "Oslo"]))
oslo_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "Oslo"]))
write_tsv(oslo_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/oslo_2017.tsv")
write_tsv(oslo_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/oslo_2018.tsv")
write_tsv(oslo_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/oslo_2019.tsv")
write_tsv(oslo_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/oslo_all.tsv")

stockholm_crosskingdom_2017 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2017 & metaair_metadata_final$CITY == "Stockholm"]))
stockholm_crosskingdom_2018 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2018 & metaair_metadata_final$CITY == "Stockholm"]))
stockholm_crosskingdom_2019 <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$YEAR == 2019 & metaair_metadata_final$CITY == "Stockholm"]))
stockholm_crosskingdom_all <- data_crosskingdom %>% select(c(taxonomy_id,metaair_metadata_final$ID[metaair_metadata_final$CITY == "Stockholm"]))
write_tsv(stockholm_crosskingdom_2017,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/stockholm_2017.tsv")
write_tsv(stockholm_crosskingdom_2018,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/stockholm_2018.tsv")
write_tsv(stockholm_crosskingdom_2019,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/stockholm_2019.tsv")
write_tsv(stockholm_crosskingdom_all,file="/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM/stockholm_all.tsv")

# PLOT
ck_city <- c(rep("Total",4),rep("Denver",4),rep("Hong Kong",4),rep("London",4), rep("New York",4),rep("Oslo",4),rep("Stockholm",4))
ck_kingdom <- rep(c("Bacteria","Fungi","Archaea","Virus"),7)
ck_perc <- c(74.4,25.3,0.3,0.0, 63.4,36.6,0.0,0.0, 84.0,15.9,0.1,0.0, 76.3,23.5,0.2,0.0, 77.5,21.9,0.6,0.0, 61.7,38.0,0.3,0.0, 72.9,25.3,1.8,0.0)
ck_df <- data.frame(year=ck_city,kingdom=ck_kingdom,perc=ck_perc)

ggplot(ck_df, aes(fill=ck_kingdom, y=ck_perc, x=ck_city)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Kingdom")


#######################
# FIG 1 CROSS-KINGDOM #
#     (NCBInr db)     #
#######################

# NOTE - Total and 2019 needs to be redone, Elements messing around.
# Order: Total, 2017, 2018, 2019
# NB - These numbers were not normalized! Redone below
# ck_unclass_fig1 <- c(23.25,19.53,22.44,24.53)
# ck_bact_fig1 <- c(60.67,70.10,62.29,57.90)
# ck_archaea_fig1 <- c(0.47, 0.35,0.51,0.44)
# ck_fungi_fig1 <- c(9.30,5.03,8.69,10.42)
# ck_metazoa_fig1 <- c(2.63,2.22,2.58,2.74)
# ck_viruses_fig1 <- c(0.21,0.20,0.21,0.22)
# ck_plantae_fig1 <- c(2.18,1.71,2.05,2.38)
# ck_othereuk_fig1 <- c(1.28, 0.86,1.23,1.38)
ck_unclass_fig1 <- c(21.80,20.08,21.70,23.71)
ck_bact_fig1 <- c(64.23,69.31,63.71,59.46)
ck_archaea_fig1 <- c(0.43, 0.36,0.49,0.43)
ck_fungi_fig1 <- c(7.68,5.08,8.17,9.89)
ck_metazoa_fig1 <- c(2.52,2.30,2.56,2.71)
ck_viruses_fig1 <- c(0.22,0.20,0.21,0.24)
ck_plantae_fig1 <- c(1.98,1.78,1.96,2.22)
ck_othereuk_fig1 <- c(1.14, 0.88,1.20,1.33)

ck_category_fig1 <- c("Total","2017","2018","2019")

ck_palette_8 <- c('#bdbdbd','#6b6ecf', '#d6616b', '#de9ed6', '#637939', '#8ca252', '#b5cf6b','#cedb9c')

fig1_df <- data.frame(Unclassified=ck_unclass_fig1, 
                      Bacteria=ck_bact_fig1, 
                      Archaea=ck_archaea_fig1,
                      Fungi=ck_fungi_fig1,
                      Metazoa=ck_metazoa_fig1, 
                      Viruses=ck_viruses_fig1, 
                      Plants=ck_plantae_fig1,
                      Other_Eukaryotes=ck_othereuk_fig1,
                      Year=ck_category_fig1)
fig1_df_long <- pivot_longer(fig1_df,cols = c(Unclassified,Bacteria,Archaea,Fungi,Metazoa,Viruses,Plants,Other_Eukaryotes))
fig1_df_long$name[fig1_df_long$name=="Other_Eukaryotes"] <- "Other Eukaryotes"
fig1_df_long$name <- factor(fig1_df_long$name,levels = c("Unclassified","Bacteria","Archaea","Viruses","Fungi","Metazoa","Plants","Other Eukaryotes"))
ggplot(fig1_df_long, aes(x=Year,y=value,fill=name)) + 
  geom_bar(position="stack",stat="identity") +xlab("Year") +ylab("Percentage") + labs(fill="Classification") +
  scale_fill_manual(values=ck_palette_8)

## TOP SPECIES ##

trim_topspecies_data <- function(master,subset){
  # Only get data that are in master's `...1` column
  sub <- subset %>% filter(`...1` %in% master$`...1`)
  return(sub)
}

add_remaining <- function(df,city){
  # Sum Percentage_of_total
  remainder <- 100.0 - sum(df$Percentage_of_total)
  newrow <- data.frame(`...1` = "Other", Percentage_of_total=remainder, Prevalence=100.0)
  newdf <- bind_rows(df,newrow)
  newdf$City <- city
  return(newdf)
}

top20.bac.colors <- c('#999999','#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')

top20_all <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/data_all_RESULTS.tsv")
top20_all_denver <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/denver_all_RESULTS.tsv")
top20_all_hk <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/hk_all_RESULTS.tsv")
top20_all_london <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/london_all_RESULTS.tsv")
top20_all_ny <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/ny_all_RESULTS.tsv")
top20_all_oslo <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/oslo_all_RESULTS.tsv")
top20_all_stockholm <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/stockholm_all_RESULTS.tsv")

top20_all_full <- add_remaining(top20_all,"Total")
top20_all_denver_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_denver),"Denver")
top20_all_hk_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_hk),"Hong Kong")
top20_all_london_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_london),"London")
top20_all_ny_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_ny),"New York")
top20_all_oslo_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_oslo),"Oslo")
top20_all_stockholm_full <- add_remaining(trim_topspecies_data(top20_all,top20_all_stockholm),"Stockholm")

top20_all_df <- bind_rows(top20_all_full,top20_all_denver_full,top20_all_hk_full,top20_all_london_full,top20_all_ny_full,top20_all_oslo_full,top20_all_stockholm_full)

ggplot(top20_all_df, aes(fill=relevel(factor(`...1`),ref="Other"), y=Percentage_of_total, x=City)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors)



# Read master data Bacteria and define which top species get colors
top20_bac <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/data_all_bacteria_RESULTS.tsv")

# Read in all other city bacteria datasets
# Note, I have solved this by printing the top 4000/1000 species for all data sets except master in OVERALL
top20_bac_denver <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/denver_all_bacteria_RESULTS.tsv")
top20_bac_hk <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/hk_all_bacteria_RESULTS.tsv")
top20_bac_london <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/london_all_bacteria_RESULTS.tsv")
top20_bac_ny <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/ny_all_bacteria_RESULTS.tsv")
top20_bac_oslo <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/oslo_all_bacteria_RESULTS.tsv")
top20_bac_stockholm <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/stockholm_all_bacteria_RESULTS.tsv")

top20_bac_full <- add_remaining(top20_bac,"Total")
top20_bac_denver_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_denver),"Denver")
top20_bac_hk_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_hk),"Hong Kong")
top20_bac_london_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_london),"London")
top20_bac_ny_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_ny),"New York")
top20_bac_oslo_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_oslo),"Oslo")
top20_bac_stockholm_full <- add_remaining(trim_topspecies_data(top20_bac,top20_bac_stockholm),"Stockholm")

top20_bac_df <- bind_rows(top20_bac_full,top20_bac_denver_full,top20_bac_hk_full,top20_bac_london_full,top20_bac_ny_full,top20_bac_oslo_full,top20_bac_stockholm_full)

ggplot(top20_bac_df, aes(fill=relevel(factor(`...1`),ref="Other"), y=Percentage_of_total, x=City)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors)


## FUNGI
top20_fungi <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/data_all_fungi_RESULTS.tsv")

# Read in all other city bacteria datasets
# Note, I have solved this by printing the top 4000/1000 species for all data sets except master in OVERALL
top20_fungi_denver <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/denver_all_fungi_RESULTS.tsv")
top20_fungi_hk <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/hk_all_fungi_RESULTS.tsv")
top20_fungi_london <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/london_all_fungi_RESULTS.tsv")
top20_fungi_ny <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/ny_all_fungi_RESULTS.tsv")
top20_fungi_oslo <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/oslo_all_fungi_RESULTS.tsv")
top20_fungi_stockholm <- read_tsv("/media/ubuntu/Pandora/METAAIR/CROSSKINGDOM_FIG4/OVERALL/stockholm_all_fungi_RESULTS.tsv")

top20_fungi_full <- add_remaining(top20_fungi,"Total")
top20_fungi_denver_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_denver),"Denver")
top20_fungi_hk_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_hk),"Hong Kong")
top20_fungi_london_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_london),"London")
top20_fungi_ny_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_ny),"New York")
top20_fungi_oslo_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_oslo),"Oslo")
top20_fungi_stockholm_full <- add_remaining(trim_topspecies_data(top20_fungi,top20_fungi_stockholm),"Stockholm")

top20_fungi_df <- bind_rows(top20_fungi_full,top20_fungi_denver_full,top20_fungi_hk_full,top20_fungi_london_full,top20_fungi_ny_full,top20_fungi_oslo_full,top20_fungi_stockholm_full)

ggplot(top20_fungi_df, aes(fill=relevel(factor(`...1`),ref="Other"), y=Percentage_of_total, x=City)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors)

########################################
# TOP 20 FIGURES, BUT FOR CONTAMINANTS #
########################################

contam_bac <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_CONTAMINANTSONLY_BAC.tsv",col_names=TRUE)
contam_fungi <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_CONTAMINANTSONLY_FUNGI.tsv",col_names=TRUE)

contam_bac <- transpose_df(contam_bac)
colnames(contam_bac) <- contam_bac[1,]
contam_bac <- contam_bac[2:nrow(contam_bac),]
colnames(contam_bac)[1] <- "sample"
contam_fungi <- transpose_df(contam_fungi)
colnames(contam_fungi) <- contam_fungi[1,]
contam_fungi <- contam_fungi[2:nrow(contam_fungi),]
colnames(contam_fungi)[1] <- "sample"

#metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
translation_bac <- metaair_taxa_names$name[match(colnames(contam_bac)[2:ncol(contam_bac)],metaair_taxa_names$taxonomy_id)]
translation_fungi <- metaair_taxa_names$name[match(colnames(contam_fungi)[2:ncol(contam_fungi)],metaair_taxa_names$taxonomy_id)]
colnames(contam_bac) <-  c("sample",translation_bac)
colnames(contam_fungi) <- c("sample", translation_fungi)

contam_bac_norm <- normalize_TSS_OTU_data(contam_bac)
contam_fungi_norm <- normalize_TSS_OTU_data(contam_fungi)
# TAKING LOG TSS COUNTS INSTEAD OF RAW
contam_bac_log <- bind_cols(contam_bac_norm %>% select(sample), contam_bac_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
contam_fungi_log <- bind_cols(contam_fungi_norm %>% select(sample), contam_fungi_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
contam_bac_log_var <- drop_constant_cols(contam_bac_log)
contam_fungi_log_var <- drop_constant_cols(contam_fungi_log)

contam_bac_final <- contam_bac_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
names(contam_bac_final)[1] <- "ID"
contam_fungi_final <- contam_fungi_log_var %>% filter(!(sample %in% CTRLS_METAAIR)) %>% filter(!(sample %in% FILTER))
names(contam_fungi_final)[1] <- "ID"

# Need:
# Top 20 species total
# Data frame with NAME, PERCENTAGE_OF_TOTAL, CITY AND POSSIBLY YEAR

getTopSpeciesSimple <- function(df, N=20, Species=NULL,City, Year){
  # CALCULATE ROWSUMS
  df <- df %>%
    mutate(rowSums = select(.,-ID) %>% rowSums(.))
  Prevalence <- df %>%
    reframe(Prevalence = colSums(select(.,-ID,-rowSums) != 0) / (nrow(.))) %>% unlist()
  #print(names(Prevalence))
  #print(Prevalence)
  #print(names(df %>% select(.,-ID,-rowSums)))
  names(Prevalence) <- names(df %>% select(.,-ID,-rowSums))
  if (is.null(Species)){
    top_species_group <- df %>%
      summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>% #.names = "{.col}_Percentage"
      select(-rowSums) %>%
      unlist() %>%
      sort(decreasing=TRUE) %>%
      head(N)
  }
  else {
    top_species_group <- df %>%
    summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>%
    select(-rowSums) %>%
    select(all_of(Species[1:20])) %>%
    unlist()
  }
  # ADD PREVALENCE
  Prevalence_to_add <- as.vector(Prevalence[names(top_species_group)])
  Prevalence_to_add <- c(Prevalence_to_add,1.0)
  names(Prevalence_to_add)[length(Prevalence_to_add)] <- "Other"
  #species_to_get_data_for <- names(top_species_group)
  top_species_group <- addOtherSpecies(top_species_group) # <- Drop this if we don't want the "other" category)
  results <- data.frame(Percentage_of_total=top_species_group, Species=names(top_species_group),City=City,Year=Year,Prevalence=Prevalence_to_add)
  
  #results <- reOrderDataForPlot(results)

  return(as_tibble(results))
}

contam_bac_final_fig <- getTopSpeciesSimple(contam_bac_final,N=20,City="Combined",Year="All Years")
master_species_bac <- contam_bac_final_fig$Species

for (city in c("Denver","Hong Kong","London","New York", "Oslo","Stockholm")){
  new <- getTopSpeciesSimple(contam_bac_final %>% filter(ID %in% (metaair_meta %>% filter(CITY == city)  %>% pull(ID))),City=city,Year="All years",Species=master_species_bac)
  contam_bac_final_fig <- bind_rows(contam_bac_final_fig,new)
}


ggplot(contam_bac_final_fig, aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=City)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))




contam_fungi_final_fig <- getTopSpeciesSimple(contam_fungi_final,N=20,City="Combined",Year="All Years")
master_species_fungi <- contam_fungi_final_fig$Species

for (city in c("Denver","Hong Kong","London","New York", "Oslo","Stockholm")){
  new <- getTopSpeciesSimple(contam_fungi_final %>% filter(ID %in% (metaair_meta %>% filter(CITY == city)  %>% pull(ID))),City=city,Year="All years",Species=master_species_fungi)
  contam_fungi_final_fig <- bind_rows(contam_fungi_final_fig,new)
}

ggplot(contam_fungi_final_fig, aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=City)) + 
  geom_bar(position="stack", stat="identity") + xlab("City") + ylab("Percentage") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))


#############################################
# ABUNDANCE PLOTS CONTROLS BEFORE AND AFTER #
#############################################

NEG_CTRLS_METAAIR <- c("SL335732","SL335733","SL335734","SL335735","SL335736","SL335737","SL335738","SL335754","SL470202","SL470204","SL470208","SL470276","SL470277","SL470278","SL470279","SL470294","SL470295","SL470297","SL470298","SL470335","SL470337","SL470338")
POS_CTRLS_METAAIR <- c("SL335740","SL335741","SL470280","SL470281","SL470339")

getTopSpeciesCtrls <- function(df, N=20, Species=NULL,Samplename="",absolute=TRUE){
  # CALCULATE ROWSUMS
  df <- df %>%
    mutate(rowSums = select(.,-ID) %>% rowSums(.))
  
  if (is.null(Species)){
    top_species_group <- df %>%
      summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>% #.names = "{.col}_Percentage"
      select(-rowSums) %>%
      unlist() %>%
      sort(decreasing=TRUE) %>%
      head(N)
    top_species_group <- addOtherSpecies(top_species_group)
  }
  else {
    if (absolute){
      top_species_group <- df %>%
        summarise(across(where(is.numeric), ~sum(.) )) %>%
        select(-rowSums) %>%
        select(all_of(Species[1:20])) %>%
        unlist()
      top_species_group <- c(top_species_group,df$rowSums)
      names(top_species_group)[length(top_species_group)] <- "Other"
      }
    else{
      top_species_group <- df %>%
        summarise(across(where(is.numeric), ~sum(.) / sum(rowSums) )) %>%
        select(-rowSums) %>%
        select(all_of(Species[1:20])) %>%
        unlist()
      top_species_group <- addOtherSpecies(top_species_group)
    }
  }
  

  results <- data.frame(Percentage_of_total=top_species_group, Species=names(top_species_group),Name=Samplename)
  return(as_tibble(results))
}

metaair_with_kitome <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/AGGREGATED_0.005_.tsv",col_names = TRUE)
metaair_with_kitome <- transpose_df(metaair_with_kitome)
colnames(metaair_with_kitome) <- metaair_with_kitome[1,]
metaair_with_kitome <- metaair_with_kitome[2:nrow(metaair_with_kitome),]
# metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
# translation <- metaair_taxa_names$name[match(colnames(metaair)[2:ncol(metaair)],metaair_taxa_names$taxonomy_id)]
colnames(metaair_with_kitome) <-  c("sample",translation)

metaair_without_kitome <- read_tsv(file="/media/ubuntu/Pandora/METAAIR/AGGREGATED_0.005_NEGSREMOVED_.tsv",col_names = TRUE)
metaair_without_kitome <- transpose_df(metaair_without_kitome)
colnames(metaair_without_kitome) <- metaair_without_kitome[1,]
metaair_without_kitome <- metaair_without_kitome[2:nrow(metaair_without_kitome),]
# metaair_taxa_names <- read_tsv("/media/ubuntu/Pandora/METAAIR/taxid_to_name_2023-11-07.tsv")
# translation <- metaair_taxa_names$name[match(colnames(metaair)[2:ncol(metaair)],metaair_taxa_names$taxonomy_id)]
colnames(metaair_without_kitome) <-  c("sample",translation)

# NOTE - CONSIDER DOING THIS WITHOUT NORMALIZATION - SHOW ACTUAL ABUNDANCE PRIOR TO NORM

#metaair_with_kitome_norm <- normalize_TSS_OTU_data(metaair_with_kitome)
#metaair_with_kitome_log <- bind_cols(metaair_with_kitome_norm %>% select(sample), metaair_with_kitome_norm %>% select(where(is.numeric)) %>% mutate(log(.+1)))
#metaair_with_kitome_log_var <- drop_constant_cols(metaair_with_kitome_log)

metaair_negctrls_kitome <- metaair_with_kitome %>% filter(sample %in% NEG_CTRLS_METAAIR)
metaair_posctrls_kitome <- metaair_with_kitome %>% filter(sample %in% POS_CTRLS_METAAIR)

metaair_negctrls_no_kitome <- metaair_without_kitome %>% filter(sample %in% NEG_CTRLS_METAAIR)
metaair_posctrls_no_kitome <- metaair_without_kitome %>% filter(sample %in% POS_CTRLS_METAAIR)
#names(metaair_final)[1] <- "ID"
names(metaair_negctrls_kitome)[1] <- "ID"
names(metaair_posctrls_kitome)[1] <- "ID"
names(metaair_negctrls_no_kitome)[1] <- "ID"
names(metaair_posctrls_no_kitome)[1] <- "ID"

#### NEGATIVES - KTIOME NOT EXCLUDED
metaair_negctrls_kitome_fig <- getTopSpeciesCtrls(metaair_negctrls_kitome,N=20,Samplename = "COMBINED",absolute = FALSE)
master_species_negctrls_kitome <- metaair_negctrls_kitome_fig$Species

#NOTE - In absolute count mode, the "COMBINED" level should not really be included.

for (s in NEG_CTRLS_METAAIR){
  new <- getTopSpeciesCtrls(metaair_negctrls_kitome %>% filter(ID == s),Samplename=s,Species=master_species_negctrls_kitome,absolute=TRUE)
  metaair_negctrls_kitome_fig <- bind_rows(metaair_negctrls_kitome_fig,new)
}

ggplot(metaair_negctrls_kitome_fig %>% filter(Name != "COMBINED"), aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=Name)) + 
  geom_bar(position="stack", stat="identity") + xlab("Sample") + ylab("Read counts") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))



#### NEGATIVES _ KITOME EXCLUDED
# USE THE SAME SPECIES AS WITH KITOME
metaair_negctrls_kitome_fig <- getTopSpeciesCtrls(metaair_negctrls_kitome,N=20,Samplename = "COMBINED",absolute = FALSE)
master_species_negctrls_kitome <- metaair_negctrls_kitome_fig$Species

for (s in NEG_CTRLS_METAAIR){
  new <- getTopSpeciesCtrls(metaair_negctrls_no_kitome %>% filter(ID == s),Samplename=s,Species=master_species_negctrls_kitome,absolute=TRUE)
  metaair_negctrls_kitome_fig <- bind_rows(metaair_negctrls_kitome_fig,new)
}

ggplot(metaair_negctrls_kitome_fig %>% filter(Name != "COMBINED"), aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=Name)) + 
  geom_bar(position="stack", stat="identity") + xlab("Sample") + ylab("Read counts") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))


## POSITIVES - WITH KITOME
metaair_posctrls_kitome_fig <- getTopSpeciesCtrls(metaair_posctrls_kitome,N=20,Samplename = "COMBINED",absolute = FALSE)
master_species_posctrls_kitome <- metaair_posctrls_kitome_fig$Species

for (s in POS_CTRLS_METAAIR){
  new <- getTopSpeciesCtrls(metaair_posctrls_kitome %>% filter(ID == s),Samplename=s,Species=master_species_posctrls_kitome,absolute=TRUE)
  metaair_posctrls_kitome_fig <- bind_rows(metaair_posctrls_kitome_fig,new)
}

ggplot(metaair_posctrls_kitome_fig %>% filter(Name != "COMBINED"), aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=Name)) + 
  geom_bar(position="stack", stat="identity") + xlab("Sample") + ylab("Read counts") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))

## POSITIVES - WITHOUT KITOME
# USE SAME SPECIES AS WITH KITOME
metaair_posctrls_kitome_fig <- getTopSpeciesCtrls(metaair_posctrls_kitome,N=20,Samplename = "COMBINED",absolute = FALSE)
master_species_posctrls_kitome <- metaair_posctrls_kitome_fig$Species

for (s in POS_CTRLS_METAAIR){
  new <- getTopSpeciesCtrls(metaair_posctrls_no_kitome %>% filter(ID == s),Samplename=s,Species=master_species_posctrls_kitome,absolute=TRUE)
  metaair_posctrls_kitome_fig <- bind_rows(metaair_posctrls_kitome_fig,new)
}

ggplot(metaair_posctrls_kitome_fig %>% filter(Name != "COMBINED"), aes(fill=relevel(factor(Species),ref="Other"), y=Percentage_of_total, x=Name)) + 
  geom_bar(position="stack", stat="identity") + xlab("Sample") + ylab("Read counts") + labs(fill="Species") +
  scale_fill_manual(values=top20.bac.colors) + theme(axis.text.x = element_text(angle=45))



# SOME OTHER USEFUL PLOTS

# Nonpareil coverage and number of reads, structured by year
ggplot(data=metaair_meta_posremoved_negremoved_filtered,aes(x=CLEANPEREADS,y=COVERAGE,color=YEAR)) + geom_point()


##################################################
###                 KARI                       ###
##################################################

# UNCOMMENT TO REDO KARI FILES
# metaair_meta_kari <- metaair_meta_posremoved_negremoved_filtered_underground
# metaair_kari_countdata <- read_tsv("/media/ubuntu/Pandora/METAAIR/AGGREGATED_USE/AGGREGATED_0.005_NEGSREMOVED_POSCTRLSREMOVED_.tsv",col_names=TRUE)
# metaair_kari <- metaair_kari_countdata %>% select(taxonomy_id, metaair_meta_kari$ID)
# 
# denver <- metaair_meta_kari %>% filter(CITY == "Denver") %>% pull(ID)
# denver_17 <- metaair_meta_kari %>% filter(CITY == "Denver" & YEAR == 2017) %>% pull(ID)
# denver_18 <- metaair_meta_kari %>% filter(CITY == "Denver" & YEAR == 2017) %>% pull(ID)
# denver_19 <- metaair_meta_kari %>% filter(CITY == "Denver" & YEAR == 2017) %>% pull(ID)
# hk <- metaair_meta_kari %>% filter(CITY == "Hong Kong") %>% pull(ID)
# hk_17 <- metaair_meta_kari %>% filter(CITY == "Hong Kong" & YEAR == 2017) %>% pull(ID)
# hk_18 <- metaair_meta_kari %>% filter(CITY == "Hong Kong" & YEAR == 2018) %>% pull(ID)
# hk_19 <- metaair_meta_kari %>% filter(CITY == "Hong Kong" & YEAR == 2019) %>% pull(ID)
# london <- metaair_meta_kari %>% filter(CITY == "London") %>% pull(ID)
# london_17 <- metaair_meta_kari %>% filter(CITY == "London" & YEAR == 2017) %>% pull(ID)
# london_18 <- metaair_meta_kari %>% filter(CITY == "London" & YEAR == 2018) %>% pull(ID)
# london_19 <- metaair_meta_kari %>% filter(CITY == "London" & YEAR == 2019) %>% pull(ID)
# ny <- metaair_meta_kari %>% filter(CITY == "New York") %>% pull(ID)
# ny_17 <- metaair_meta_kari %>% filter(CITY == "New York" & YEAR == 2017) %>% pull(ID)
# ny_18 <- metaair_meta_kari %>% filter(CITY == "New York" & YEAR == 2018) %>% pull(ID)
# ny_19 <- metaair_meta_kari %>% filter(CITY == "New York" & YEAR == 2019) %>% pull(ID)
# oslo <- metaair_meta_kari %>% filter(CITY == "Oslo") %>% pull(ID)
# oslo_17 <- metaair_meta_kari %>% filter(CITY == "Oslo" & YEAR == 2017) %>% pull(ID)
# oslo_18 <- metaair_meta_kari %>% filter(CITY == "Oslo" & YEAR == 2018) %>% pull(ID)
# oslo_19 <- metaair_meta_kari %>% filter(CITY == "Oslo" & YEAR == 2019) %>% pull(ID)
# stockholm <- metaair_meta_kari %>% filter(CITY == "Stockholm") %>% pull(ID)
# stockholm_17 <- metaair_meta_kari %>% filter(CITY == "Stockholm" & YEAR == 2017) %>% pull(ID)
# stockholm_18 <- metaair_meta_kari %>% filter(CITY == "Stockholm" & YEAR == 2018) %>% pull(ID)
# year_17 <- metaair_meta_kari %>% filter(YEAR == 2017) %>% pull(ID)
# year_18 <- metaair_meta_kari %>% filter(YEAR == 2018) %>% pull(ID)
# year_19 <- metaair_meta_kari %>% filter(YEAR == 2019) %>% pull(ID)
# 
# metaair_kari_denver <- metaair_kari %>% select(taxonomy_id, denver)
# metaair_kari_denver_17 <- metaair_kari %>% select(taxonomy_id, denver_17)
# metaair_kari_denver_18 <- metaair_kari %>% select(taxonomy_id, denver_18)
# metaair_kari_denver_19 <- metaair_kari %>% select(taxonomy_id,denver_19)
# metaair_kari_hongkong <- metaair_kari %>% select(taxonomy_id, hk)
# metaair_kari_hongkong_17 <- metaair_kari %>% select(taxonomy_id, hk_17)
# metaair_kari_hongkong_18 <- metaair_kari %>% select(taxonomy_id, hk_18)
# metaair_kari_hongkong_19 <- metaair_kari %>% select(taxonomy_id, hk_19)
# metaair_kari_london <- metaair_kari %>% select(taxonomy_id, london)
# metaair_kari_london_17 <- metaair_kari %>% select(taxonomy_id, london_17)
# metaair_kari_london_18 <- metaair_kari %>% select(taxonomy_id, london_18)
# metaair_kari_london_19 <- metaair_kari %>% select(taxonomy_id, london_19)
# metaair_kari_newyork <- metaair_kari %>% select(taxonomy_id, ny)
# metaair_kari_newyork_17 <- metaair_kari %>% select(taxonomy_id, ny_17)
# metaair_kari_newyork_18 <- metaair_kari %>% select(taxonomy_id, ny_18)
# metaair_kari_newyork_19 <- metaair_kari %>% select(taxonomy_id, ny_19)
# metaair_kari_oslo <- metaair_kari %>% select(taxonomy_id, oslo)
# metaair_kari_oslo_17 <- metaair_kari %>% select(taxonomy_id, oslo_17)
# metaair_kari_oslo_18 <- metaair_kari %>% select(taxonomy_id, oslo_18)
# metaair_kari_oslo_19 <- metaair_kari %>% select(taxonomy_id, oslo_19)
# metaair_kari_stockholm <- metaair_kari %>% select(taxonomy_id, stockholm)
# metaair_kari_stockholm_17 <- metaair_kari %>% select(taxonomy_id, stockholm_17)
# metaair_kari_stockholm_18 <- metaair_kari %>% select(taxonomy_id, stockholm_18)
# metaair_kari_17 <- metaair_kari %>% select(taxonomy_id, year_17)
# metaair_kari_18 <- metaair_kari %>% select(taxonomy_id, year_18)
# metaair_kari_19 <- metaair_kari %>% select(taxonomy_id, year_19)
# 
# write_tsv(x=metaair_kari_denver,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/denver.tsv") 
# write_tsv(x=metaair_kari_denver_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/denver_17.tsv") 
# write_tsv(x=metaair_kari_denver_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/denver_18.tsv") 
# write_tsv(x=metaair_kari_denver_19,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/denver_19.tsv") 
# write_tsv(x=metaair_kari_hongkong,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/hongkong.tsv") 
# write_tsv(x=metaair_kari_hongkong_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/hongkong_17.tsv") 
# write_tsv(x=metaair_kari_hongkong_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/hongkong_18.tsv") 
# write_tsv(x=metaair_kari_hongkong_19,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/hongkong_19.tsv") 
# write_tsv(x=metaair_kari_london,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/london.tsv") 
# write_tsv(x=metaair_kari_london_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/london_17.tsv") 
# write_tsv(x=metaair_kari_london_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/london_18.tsv") 
# write_tsv(x=metaair_kari_london_19,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/london_19.tsv") 
# write_tsv(x=metaair_kari_newyork,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/newyork.tsv") 
# write_tsv(x=metaair_kari_newyork_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/newyork_17.tsv") 
# write_tsv(x=metaair_kari_newyork_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/newyork_18.tsv") 
# write_tsv(x=metaair_kari_newyork_19,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/newyork_19.tsv") 
# write_tsv(x=metaair_kari_oslo,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/oslo.tsv") 
# write_tsv(x=metaair_kari_oslo_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/oslo_17.tsv") 
# write_tsv(x=metaair_kari_oslo_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/oslo_18.tsv") 
# write_tsv(x=metaair_kari_oslo_19,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/oslo_19.tsv") 
# write_tsv(x=metaair_kari_stockholm,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/stockholm.tsv") 
# write_tsv(x=metaair_kari_stockholm_17,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/stockholm_17.tsv") 
# write_tsv(x=metaair_kari_stockholm_18,file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/stockholm_18.tsv") 
# write_tsv(x=metaair_kari, file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/total.tsv")
# write_tsv(x=metaair_kari_17, file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/total_17.tsv")
# write_tsv(x=metaair_kari_18, file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/total_18.tsv")
# write_tsv(x=metaair_kari_19, file="/media/ubuntu/Pandora/KARI_ALL/METAAIR_UNDERGROUND/total_19.tsv")

# metaair_kari_notdenver <- metaair_final_ %>% filter(ID %in% (metaair_meta_kari %>% filter(!CITY == "Denver") %>% pull(ID)))
# NOTE - UMAP performed on matrix where rows are samples and columns are species, not other way around

# umap.kari.notdenver.meta <- umap(metaair_kari_notdenver %>% select(where(is.numeric)),n_neighbors=8,min_dist=0.1)
# umap.meta.kari.notdenver.layout <- umap.kari.notdenver.meta[["layout"]] 
# umap.meta.kari.notdenver.layout <- data.frame(umap.meta.kari.notdenver.layout, "ID"=metaair_kari_notdenver$ID) 
# umap.meta.kari.notdenver.layout <- inner_join(umap.meta.kari.notdenver.layout, metaair_meta_posremoved_negremoved_filtered)
#umap.meta.layout <- cbind(umap.meta.layout, metaair_meta %>% filter (!(ID %in% CTRLS_METAAIR))) 

# fig.umap.meta.notdenver <- plot_ly(umap.meta.kari.notdenver.layout, x = ~X1, y = ~X2, color=~CITY, symbol = ~YEAR, type = 'scatter', mode = 'markers', size=10,
#                         text = ~ID, hovertemplate=paste(
#                           "<b>%{text}</b><br><br>",
#                           "%{x}, %{y}")) %>%
#  layout(
#    plot_bgcolor = "#e5ecf6",
#    legend=list(title=list(text='Sample type')), 
#    xaxis = list( 
#      title = "0"),  
#    yaxis = list( 
#      title = "1")) 
#fig.umap.meta.notdenver
