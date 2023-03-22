library(ggplot2)
library(magrittr)
library(dplyr)
library(dummies)
library(ggtree)
library(data.table)

# # how is jaccard-distance implemented? 
# # test:
# vegan::vegdist(x = data.frame(row.names = c("A", "B", "C"),
#                               trait1 = c(1,0,1),
#                               trait2 = c(0,1,1),
#                               trait3 = c(0,2,0)),
#                method = "jaccard",
#                upper = F)
# #result: jaccard-distance is implemented as a distance and not a similiarity measure
# # jaccard distance is wrong!!! it does not respect the order of traits!!!


test_df <- 
data.frame(row.names = c("A", "B", "C"),
           trait1 = c(1,0,1),
           trait2 = c(0,1,1),
           trait3 = c(0,2,0))

# hamming distance
library(e1071)
e1071::hamming.distance(as.matrix(test_df))/ncol(test_df)

# gower distance
vegan::vegdist(x = test_df,
               method = "gower",
               upper = F)

########################################################
# DATA
########################################################

# Meta data
transect_i_df <- readr::read_csv(file.path("1_data", 
                                           "transect_superregion_macroregion.csv"))

taxgroups <- unique(readr::read_csv(file.path("1_data", 
                                              "taxgroups_v9_13Mar23.csv")))

# set output path
output_path <- file.path("3_output", 
                         "42_distance")
dir.create(output_path)


# load color hexcodes
source(file.path("2_scripts", 
                 "color_palette.R"))
higher_order_hexcodes

higher_order_hexcodes_manual <- setNames(higher_order_hexcodes$color_hex,
                                         higher_order_hexcodes$higher_order)

taxgroups_fac_hexcodes <- 
  dplyr::left_join(taxgroups,
                   higher_order_hexcodes) 
taxgroups_fac_hexcodes$higher_order  <-  factor(taxgroups_fac_hexcodes$higher_order,
                                                   levels = forcats::fct_reorder(unique(taxgroups_fac_hexcodes$higher_order),
                                                                                 unique(taxgroups_fac_hexcodes$chrono_order)))


###########
# TOOLS
MASTERTABLE_v7_all_revised_tools_with_TSrange_uniqueID <- 
  readr::read_delim(file.path("1_data", 
                              "1511NAC_Database", 
                              "MASTERTABLE_v7_all_revised_tools_with_TSrange_uniqueID.csv"),
                    delim = ",", 
                    escape_double = T, 
                    trim_ws = TRUE)

tiplabs_tools_df <- dplyr::left_join(MASTERTABLE_v7_all_revised_tools_with_TSrange_uniqueID[,c("TaxUnit_unique", 
                                                                                               "TaxUnit_unique_TS")],
                                     taxgroups_fac_hexcodes) %>% 
  unique()

# tiplabs_tools_df <- dplyr::left_join(tiplabs_tools_df,
#                                      higher_order_hexcodes)

tools <- MASTERTABLE_v7_all_revised_tools_with_TSrange_uniqueID
tools_matrix <- as.matrix(tools[,c(which(names(tools) == "A_p"):which(names(tools) == "D_adze_axe"))]) # columns 7:59 = "A_p":"D_adze_axe"
tools_matrix[which(tools_matrix == "na")] <- 2 # "not applicable" becomes third character state

tools_df <- as.data.frame(apply(tools_matrix, 2, as.numeric))
# tools_df_dummy <- dummies::dummy.data.frame(data = tools_df, dummy.classes = "ALL")
rownames(tools_df) <- tools$TaxUnit_unique_TS 

dist_tools_df_dummy <- vegan::vegdist(x = tools_df,
                                 method = "gower",
                                 upper = F)


#############
# TECHNOLOGY
MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID <- readr::read_delim(file.path("1_data", 
                                                                                           "1511NAC_Database", 
                                                                                           "MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID.csv"),
                                                                                 delim = ",", 
                                                                                 escape_double = T, 
                                                                                 trim_ws = TRUE)

tiplabs_technology_df <- dplyr::left_join(MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID[,c("TaxUnit_unique", 
                                                                                                         "TaxUnit_unique_TS")],
                                          taxgroups_fac_hexcodes)

technology <- MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID
technology_matrix <- as.matrix(technology[,c(which(names(technology) == "LP_reduction_strat_1"):which(names(technology) == "Microburin"))]) # columns 7:59 = "LP_reduction_strat_1":"Microburin"
technology_matrix[which(technology_matrix == "na")] <- 2 # not applicable becomes third character state

technology_df <- as.data.frame(apply(technology_matrix, 2, as.numeric))
# technology_df_dummy <- dummies::dummy.data.frame(data = technology_df, dummy.classes = "ALL")
rownames(technology_df) <- technology$TaxUnit_unique_TS 

dist_technology_df_dummy <- vegan::vegdist(x = technology_df,
                                      method = "gower",
                                      upper = F)


###########
# OUTLINES
outlines_AR_centered_w_metric_measures_scaled_PCA <- readRDS(file.path("3_output", 
                                                                       "20_outline_analysis_PCA",
                                                                       "outlines_AR_centered_w_metric_measures_scaled_PCA.RDS"))
outlines <- outlines_AR_centered_w_metric_measures_scaled_PCA$fac
outlines_dist_matrix <- as.matrix(dist(outlines_AR_centered_w_metric_measures_scaled_PCA$x,
                                       upper = F))
diag(outlines_dist_matrix) <- NA

tiplabs_outlines_df <- 
  dplyr::left_join(unique(do.call(rbind.data.frame, lapply(outlines_AR_centered_w_metric_measures_scaled_PCA$fac$ARTEFACTNAME,
                                                             function(a){
                                                               data.frame(tiplab = a,
                                                                          TaxUnit_unique = paste0(strsplit(a,split = "_")[[1]][[2]],
                                                                                                  "_",
                                                                                                  strsplit(a,split = "_")[[1]][[3]]))}))),
                   taxgroups %>% select(TaxUnit_unique, higher_order) %>% unique(),
                   by = "TaxUnit_unique")

tiplabs_outlines_df <- dplyr::left_join(tiplabs_outlines_df,
                                        unique(MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID[,c("TaxUnit_unique")]),
                                        by = "TaxUnit_unique")
tiplabs_outlines_df <- dplyr::left_join(tiplabs_outlines_df,
                                          higher_order_hexcodes)

# tiplabs_outlines_df <- tiplabs_outlines_df[-which(duplicated(tiplabs_outlines_df)),]


#######################################
# # Standardised Effect Size Test (SES)
#######################################

no_sampling_it_outlines <- 100000
no_sampling_it <- 10000

unique_semantic_groups <- as.vector(na.omit(unique(taxgroups$higher_order)))

cophenet_TaxUnit_df_list_datatype <- list()
pairwise_meandistances_within_df_list_datatype <- list()
pairwise_meandistances_NULL_df_list_datatype <- list()

for(current_class_index in c("Outlines", "Tools", "Technology")){ #  "Outlines", 
  
  if(current_class_index == "Outlines") {
    current_distance_matrix <- outlines_dist_matrix
    tiplabs_currentClass_df_raw <- tiplabs_outlines_df
    
  } else if(current_class_index == "Tools") {
    current_distance_matrix <- as.matrix(dist_tools_df_dummy)
    tiplabs_currentClass_df_raw <- tiplabs_tools_df
    
  } else if(current_class_index == "Technology") {
    current_distance_matrix <- as.matrix(dist_technology_df_dummy)
    tiplabs_currentClass_df_raw <- tiplabs_technology_df
  }
  
  diag(current_distance_matrix) <- NA
  current_distance_matrix[upper.tri(current_distance_matrix)] <- NA
  
  if(current_class_index == "Outlines"){
    current_label <- "tiplab"
  } else {
    current_label <- "TaxUnit_unique_TS"
  }
  
  
  cophenet_TaxUnit_df_list <- list()
  pairwise_meandistances_within_df_list <- list()
  pairwise_meandistances_NULL_df_list <- list()
  
  for(current_semantic_group in unique_semantic_groups){
    
    # WITHIN-GROUP OBSERVATIONS
    
    ## select the labels of the within-group samples
    within_group_sampleNames <- subset(tiplabs_currentClass_df_raw, higher_order == current_semantic_group)[[current_label]]
    
    currentClass_cophen_by_TaxUnit <-
      current_distance_matrix[which(rownames(current_distance_matrix) %in% within_group_sampleNames), # from rownames observed
                              which(colnames(current_distance_matrix) %in% within_group_sampleNames)] # to colnames observed
    
    ## save the distances as vector
    pairwise_distances_within <- as.vector(na.omit(as.vector(currentClass_cophen_by_TaxUnit)))
    
    ## calculate the mean within-group distance
    meandist_within <- mean(pairwise_distances_within)
    
    # calculate the within-group variance
    # meandist_within <- var(pairwise_distances_within)
    
    pairwise_meandistances_within_df_list[[current_semantic_group]] <- 
      data.frame(mean_distance = meandist_within,
                 group = "Within-group",
                 semantic_group = current_semantic_group,
                 class = current_class_index)
    
    
    
    # NULL-DISTRIBUTION
    
    Null_distance_matrix <- current_distance_matrix
    
    ## select all distances that are the within-group pairwise distances and set them to NA
    Null_distance_matrix[which(rownames(Null_distance_matrix) %in% within_group_sampleNames), 
                         which(colnames(Null_distance_matrix) %in% within_group_sampleNames)] <- NA
    ## save all distances as vector
    pairwise_distances_NULL <- as.vector(na.omit(as.vector(Null_distance_matrix)))
    
    ## permutation
    if(current_class_index == "Outlines"){
      no_sampling_it <- no_sampling_it_outlines # more samples for outlines, as it's a larger data set
    } else {
      no_sampling_it <- no_sampling_it
    }
    
    future::plan(future::multisession(workers = parallel::detectCores())) ## Run in parallel on local computer
    tictoc::tic()
    resampled_values <-
      future.apply::future_lapply(1:no_sampling_it,
                                  future.seed = TRUE,
                                  FUN = function(x){
                                    sample <- sample(x = pairwise_distances_NULL, # sample from the distance matrix. 
                                                     size = length(pairwise_distances_within), # the same number of observations as are included in the within-group vector
                                                     replace = FALSE)
                                    mean(sample) # For each iteration, calculate the mean distance.
                                  }
      ) %>%
      unlist()
    tictoc::toc()
    
    pairwise_meandistances_NULL_df_list[[current_semantic_group]] <- 
      data.frame(mean_distance = resampled_values,
                 group = "Null-distribution",
                 semantic_group = current_semantic_group,
                 class = current_class_index)
    
    
    # calculate SES
    SES <- (meandist_within - mean(resampled_values))/sd(resampled_values)
    print(paste(current_class_index, current_semantic_group, "SES:", SES))
    
    
    cophenet_TaxUnit_df_list[[current_semantic_group]] <-
      data.frame(semantic_group = current_semantic_group,
                 meanMeasure_within = meandist_within,
                 meanMeasure_NULL = mean(resampled_values),
                 effect_size = meandist_within - mean(resampled_values),
                 std_effect_size = round(SES, digits = 2),
                 nPairwiseDistances_within = length(pairwise_distances_within),
                 nPairwiseDistances_NULL = length(pairwise_distances_NULL),
                 nWithinGroupObs = length(within_group_sampleNames),
                 nAllObsAvail = nrow(tiplabs_currentClass_df_raw),
                 class = current_class_index)
    
  }
  
  cophenet_TaxUnit_df_list_datatype[[current_class_index]] <- do.call(rbind.data.frame, cophenet_TaxUnit_df_list)
  pairwise_meandistances_within_df_list_datatype[[current_class_index]] <- do.call(rbind.data.frame, pairwise_meandistances_within_df_list) 
  pairwise_meandistances_NULL_df_list_datatype[[current_class_index]] <- do.call(rbind.data.frame, pairwise_meandistances_NULL_df_list)
  
}

cophenet_TaxUnit_df_list_datatype_df <- do.call(rbind.data.frame, cophenet_TaxUnit_df_list_datatype)
cophenet_TaxUnit_df_list_datatype_df

pairwise_meandistances_within_df_list_datatype_df <- do.call(rbind.data.frame, pairwise_meandistances_within_df_list_datatype) %>% 
  dplyr::left_join(., taxgroups_fac_hexcodes,
                   by = c("semantic_group" = "higher_order"))
pairwise_meandistances_NULL_df_list_datatype_df <- do.call(rbind.data.frame, pairwise_meandistances_NULL_df_list_datatype) %>% 
  dplyr::left_join(., taxgroups_fac_hexcodes,
                   by = c("semantic_group" = "higher_order"))

cophenet_TaxUnit_df_list_datatype_df$SES <- paste0("SES = ", cophenet_TaxUnit_df_list_datatype_df$std_effect_size)

library(ggplot2)
ggplot() +
  geom_histogram(data = pairwise_meandistances_NULL_df_list_datatype_df,
                 aes(x = mean_distance)) +
  geom_vline(data = pairwise_meandistances_within_df_list_datatype_df,
             aes(xintercept = mean_distance,
                 color = semantic_group),
             size = 2) +
  geom_label(data = cophenet_TaxUnit_df_list_datatype_df,
             mapping = aes(x = Inf, 
                           y = Inf, 
                           label = SES),
             hjust   = 1.3,
             vjust   = 1) +
  facet_wrap(class~semantic_group,
             scales = "free", 
             ncol = 7) +
  theme_bw() +
  xlab("Mean distance") +
  ylab("N") +
  scale_color_manual(values = higher_order_hexcodes_manual) +
  # guides(color = guide_legend("Semantic group")) +
  theme(legend.position="none")



readr::write_csv(x = cophenet_TaxUnit_df_list_datatype_df,
                 file = file.path(output_path, "42_distance_SES.csv"))

saveRDS(list(Null = pairwise_meandistances_NULL_df_list_datatype_df,
             Within = pairwise_meandistances_within_df_list_datatype_df),
        file = file.path(output_path, "42_distance_SES.RDS"))



##################################### 
# bootstrapped dendrograms
##################################### 


# tools

# plot(hclust(dist_tools_df_dummy))
# ape::plot.phylo(ape::as.phylo(hclust(dist_tools_df_dummy)))
# dist_tools_df_dummy_plot <-
# ggtree(ape::ladderize(ape::as.phylo(hclust(dist_tools_df_dummy))))+
#   geom_tiplab(align = F)

# bootstrapping
dist_fun <- function(x) {
  dista <- vegan::vegdist(x = x,
                          method = "gower")
  dista[is.nan(dista)] <- 0 # for gower distance in vegdist: if two taxa have the same traits, vegdist returns NaN
  return(dista)
}
f <- function(x) {
  # phangorn::upgma(dist_fun(x))
  ape::as.phylo(hclust((dist_fun(x)), method = "ward.D"))
}
tools_starting_tree <- f(tools_df)
tools_X <- ape::boot.phylo(phy = tools_starting_tree, 
                           x = tools_df, 
                           FUN = f, 
                           trees = TRUE,
                           B = 1000) 
tools_tree <- phangorn::plotBS(tools_starting_tree, tools_X$trees)
tools_tree2 <- phangorn::pruneTree(tools_tree, 50)#mayority consensus tree

dist_tools_df_dummy_plot <-
  ggtree(tools_tree2)+
  geom_tiplab(align = T) 

q <- dist_tools_df_dummy_plot
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)

plopp <- q + geom_text(data=d, aes(label=label))

tools_ggtree <- unique(as.data.frame(tiplabs_tools_df[tiplabs_tools_df$TaxUnit_unique_TS %in%
                                                        tools_tree2$tip.label, ]) %>%
                         dplyr::select(TaxUnit_unique_TS, higher_order, Expert_editor, chrono_order, color_hex) %>%
                         dplyr::mutate(higher_order = forcats::fct_reorder(higher_order, chrono_order))) %>% 
  unique()

rownames(tools_ggtree) <- tools_ggtree$TaxUnit_unique_TS

tools_ggtree <- 
tools_ggtree %>% 
  dplyr::select(higher_order)

tools_ggtree <- dummies::dummy.data.frame(tools_ggtree, fun = as.logical)

colnames(tools_ggtree) <- gsub("higher_order", "", colnames(tools_ggtree))

for (i in colnames(tools_ggtree)) {
  
  tools_ggtree[[i]] <- gsub(TRUE, paste0(i), tools_ggtree[[i]])
    
}



tools_bootstrap_color_plot <- 
gheatmap(plopp,
         tools_ggtree,
         offset = 0.1,
         width = 0.3,
         colnames_angle = 45,
         # colnames_offset_y = -0.3,
         colnames_position = "bottom",
         hjust = 1) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0,10))+
  scale_fill_manual(values = c("NA" = "#999999", # NA
                               "Azilian" = "#E69F00", # Azilian
                               "FBBT/LBI" = "#56B4E9", # FBBT
                               "Mesolithic" = "#009E73", # Meso
                               "ABP" = "#F0E442", # ABP/FMG
                               "LTP" = "#0072B2", # LTP/TPC
                               "Magdalenian" = "#D55E00", # Magdalenian
                               "Epigravettian" = "#CC79A7", # Epigravettian,
                               "FALSE" = "grey96")) +
  guides(fill="none") #+
  # ggtitle("Tools")#

tools_bootstrap_color_plot

ggsave(plot = tools_bootstrap_color_plot,
       filename = file.path(output_path, "tools_bootstrap_color_plot.png"), 
       width = 35, height = 40, units = "cm")
# ggsave(plot = tools_bootstrap_color_plot,
#        filename = file.path(output_path, "tools_bootstrap_color_plot.svg"), 
#        width = 35, height = 40, units = "cm")


# technology

# # plot(hclust(dist_technology_df_dummy))
# # ape::plot.phylo(ape::as.phylo(hclust(dist_technology_df_dummy)))
# dist_technology_df_dummy_plot <-
#   ggtree(ape::ladderize(ape::as.phylo(hclust(dist_technology_df_dummy))))+
#   geom_tiplab(align = F)

technology_starting_tree <- f(technology_df)
technology_X <- ape::boot.phylo(phy = technology_starting_tree, 
                                x = technology_df, 
                                FUN = f, 
                                trees = TRUE,
                                B = 1000) 
# phytools::writeNexus(technology_X$trees,
#                      file = file.path(output_path, "technology_trees.trees"))

technology_tree <- phangorn::plotBS(technology_starting_tree, technology_X$trees)
technology_tree2 <- phangorn::pruneTree(technology_tree, 50.01)#mayority consensus tree

dist_technology_df_dummy_plot <-
  ggtree(technology_tree2)+
  geom_tiplab(align = T)

q <- dist_technology_df_dummy_plot
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)

blopp <- q + geom_text(data=d, aes(label=label))

technology_ggtree <- unique(as.data.frame(tiplabs_technology_df[tiplabs_technology_df$TaxUnit_unique_TS %in%
                                                                  technology_tree2$tip.label, ]) %>%
                              select(TaxUnit_unique_TS, higher_order, Expert_editor) %>%
                              mutate(higher_order = as.factor(higher_order)))

rownames(technology_ggtree) <- technology_ggtree$TaxUnit_unique_TS

technology_ggtree <- technology_ggtree %>% 
  select(higher_order)

technology_ggtree <- dummies::dummy.data.frame(technology_ggtree, fun = as.logical)
colnames(technology_ggtree) <- gsub("higher_order", "", colnames(technology_ggtree))

for (i in colnames(technology_ggtree)) {
  
  technology_ggtree[[i]] <- gsub(TRUE, paste0(i), technology_ggtree[[i]])
  
}


technology_bootstrap_color_plot <- 
gheatmap(blopp,
         technology_ggtree,
         offset = 0.1,
         width = 0.3,
         colnames_angle = 45,
         # colnames_offset_y = -0.3,
         colnames_position = "bottom",
         hjust = 1) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0,10))+  
  scale_fill_manual(values = c("NA" = "#999999", # NA
                               "Azilian" = "#E69F00", # Azilian
                               "FBBT/LBI" = "#56B4E9", # FBBT
                               "Mesolithic" = "#009E73", # Meso
                               "ABP" = "#F0E442", # ABP/FMG
                               "LTP" = "#0072B2", # LTP/TPC
                               "Magdalenian" = "#D55E00", # Magdalenian
                               "Epigravettian" = "#CC79A7", # Epigravettian,
                               "FALSE" = "grey96")) +
  guides(fill="none") 

technology_bootstrap_color_plot

ggsave(plot = technology_bootstrap_color_plot,
       filename = file.path(output_path, "technology_bootstrap_color_plot.png"), 
       width = 35, height = 40, units = "cm")

# outlines 


# plot(hclust(dist_tools_df_dummy))
# ape::plot.phylo(ape::as.phylo(hclust(dist_tools_df_dummy)))
hclust_outlines <- hclust(dist(outlines_AR_centered_w_metric_measures_scaled_PCA$x[,c(1:scree_min)],
                               upper = F))

dist_outlines_df_dummy_plot <-
  ggtree(hclust_outlines)


set.seed(123)
stratified_subset_outlines_taxunit <- splitstackshape::stratified(tiplabs_outlines_df,
                                                                  group = "TaxUnit_unique",
                                                                  size = 5)

outlines_subset <- outlines_AR_centered_w_metric_measures_scaled_PCA$x[which(rownames(outlines_AR_centered_w_metric_measures_scaled_PCA$x) %in% stratified_subset_outlines_taxunit$tiplab),c(1:scree_min)]

f_o <- function(x) ape::as.phylo(hclust(dist(x), method = "ward.D2"))
outlines_starting_tree <- f_o(outlines_subset)
outlines_X <- ape::boot.phylo(phy = outlines_starting_tree, 
                              x = outlines_subset, 
                              FUN = f_o, 
                              trees = TRUE,
                              B = 1000) 
# phytools::writeNexus(outlines_X$trees,
#                      file = file.path(output_path, "outlines_trees.trees"))

outlines_tree <- phangorn::plotBS(outlines_starting_tree, outlines_X$trees)
outlines_tree2 <- phangorn::pruneTree(outlines_tree, 50.01)#mayority consensus tree

dist_outlines_df_dummy_plot <-
  ggtree(outlines_tree2) #+
# geom_tiplab(align = T)

q <- dist_outlines_df_dummy_plot
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 50,]

blubb <- q + geom_text(data=d, aes(label=label))

outlines_ggtree <- unique(as.data.frame(tiplabs_outlines_df[tiplabs_outlines_df$tiplab %in%
                                                              stratified_subset_outlines_taxunit$tiplab, ]) %>%
                            select(tiplab, higher_order) %>%
                            mutate(higher_order = as.factor(higher_order)))
rownames(outlines_ggtree) <- outlines_ggtree$tiplab
outlines_ggtree <- outlines_ggtree %>% select(higher_order)
outlines_ggtree <- dummies::dummy.data.frame(outlines_ggtree, fun = as.logical)
colnames(outlines_ggtree) <- gsub("higher_order", "", colnames(outlines_ggtree))


for (i in colnames(outlines_ggtree)) {
  
  outlines_ggtree[[i]] <- gsub(TRUE, paste0(i), outlines_ggtree[[i]])
  
}


outlines_bootstrap_color_plot <- 
  gheatmap(blubb,
           outlines_ggtree,
           offset = 0.1,
           width = 0.3,
           colnames_angle = 45,
           # colnames_offset_y = -0.3,
           colnames_position = "bottom",
           hjust = 1) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0,40))+  
  scale_fill_manual(values = c("NA" = "#999999", # NA
                               "Azilian" = "#E69F00", # Azilian
                               "FBBT/LBI" = "#56B4E9", # FBBT
                               "Mesolithic" = "#009E73", # Meso
                               "ABP" = "#F0E442", # ABP/FMG
                               "LTP" = "#0072B2", # LTP/TPC
                               "Magdalenian" = "#D55E00", # Magdalenian
                               "Epigravettian" = "#CC79A7", # Epigravettian,
                               "FALSE" = "grey96")) +
  guides(fill="none") 


outlines_bootstrap_color_plot


ggsave(plot = outlines_bootstrap_color_plot,
       filename = file.path(output_path, "outlines_bootstrap_color_plot.png"), 
       width = 35, height = 40, units = "cm")



##################################### 
# tanglegram tools-technology 
##################################### 

##### tanglegram
library(dendextend)
# dndlist <- dendextend::dendlist(phytools::force.ultrametric(tools_tree2), phytools::force.ultrametric(technology_tree2))
dndlist <- dendextend::dendlist(as.dendrogram(tools_starting_tree), as.dendrogram(technology_starting_tree))

dendextend::tanglegram(dndlist,
                       highlight_distinct_edges = F,
                       common_subtrees_color_lines = T,
                       common_subtrees_color_lines_default_single_leaf_color = "grey96",
                       common_subtrees_color_branches = F,
                       margin_inner = 10,
                       sort = T,
                       main_left = "Tools",
                       main_right = "Technology")


ggsave(plot = last_plot(),
       filename = file.path(output_path, "tanglegram_tools_tech.png"), 
       width = 35, height = 40, units = "cm")

# tools_to_prune <- labels(as.dendrogram(hclust(dist_tools_df_dummy)))[which(!(labels(as.dendrogram(hclust(dist_tools_df_dummy))) %in%
#                                                                                labels(as.dendrogram(hclust(dist_technology_df_dummy)))))]
# 
# tech_to_prune <- labels(as.dendrogram(hclust(dist_technology_df_dummy)))[which(!(labels(as.dendrogram(hclust(dist_technology_df_dummy))) %in%
#                                                                                    labels(as.dendrogram(hclust(dist_tools_df_dummy)))))]
# 
# tech_dendro <- as.dendrogram(hclust(dist_technology_df_dummy))
# tools_dendro <- as.dendrogram(hclust(dist_tools_df_dummy))
# 
# dndlist2 <- dendextend::dendlist(dendextend::prune(tools_dendro, tools_to_prune),
#                                  dendextend::prune(tech_dendro, tech_to_prune))
# dndlist2 %>% entanglement()
# dndlist2 %>% untangle(method = "DendSer") %>%
#   tanglegram(highlight_distinct_edges = F,
#              common_subtrees_color_lines = T,
#              common_subtrees_color_lines_default_single_leaf_color = "grey",
#              common_subtrees_color_branches = FALSE,
#              margin_inner = 8,
#              sort = T,
#              main_left = "Tools",
#              main_right = "Technology")
# 
# # while (TRUE) tools_dendro <- click_rotate(tools_dendro)
# 
# tools_multi <- ape::di2multi(ape::as.phylo(hclust(dist_tools_df_dummy)), 0.055)
# plot(tools_multi)
# technology_multi <- ape::di2multi(ape::as.phylo(hclust(dist_technology_df_dummy)), 0.0475)
# plot(technology_multi)
# 
# dendextend::dendlist(as.dendrogram(phytools::force.ultrametric(ape::ladderize(tools_multi))),
#                      as.dendrogram(phytools::force.ultrametric(ape::ladderize(technology_multi)))) %>%
#   untangle(method = "DendSer") %>%
#   tanglegram(highlight_distinct_edges = F,
#              common_subtrees_color_lines = T,
#              common_subtrees_color_lines_default_single_leaf_color = "grey",
#              common_subtrees_color_branches = FALSE,
#              margin_inner = 8,
#              sort = T,
#              main_left = "Tools",
#              main_right = "Technology")

# obj2 <- phytools::cophylo(tr1= phytools::force.ultrametric(ape::ladderize(ape::as.phylo((dendextend::prune(tools_multi, tools_to_prune))))),
#                           tr2 = phytools::force.ultrametric(ape::ladderize(ape::as.phylo(dendextend::prune(technology_multi, tech_to_prune)))),
#                           # assoc = matrix(c(labels(as.dendrogram(hclust(dist_technology_df_dummy)))[which((labels(as.dendrogram(hclust(dist_technology_df_dummy))) %in%
#                           #                                                                                   labels(as.dendrogram(hclust(dist_tools_df_dummy)))))],
#                           #                  labels(as.dendrogram(hclust(dist_technology_df_dummy)))[which((labels(as.dendrogram(hclust(dist_technology_df_dummy))) %in%
#                           #                                                                                   labels(as.dendrogram(hclust(dist_tools_df_dummy)))))]),
#                           #                ncol = 2),
#                           use.edge.length= T)
# library(phytools)
# plot(obj2, link.type="curved",link.lwd=4,
#      link.lty="solid",link.col=make.transparent("red",
#                                                 0.25))
# obj <- phytools::cophylo(tr1= ape::as.phylo(hclust(dist_tools_df_dummy)),
#                          tr2 = ape::as.phylo(hclust(dist_technology_df_dummy)),
#                          use.edge.length= F)
# plot(obj, link.type="curved",link.lwd=4,
#      link.lty="solid",link.col=make.transparent("red",
#                                                 0.25))
# 
# 
# obj <- phytools::cophylo(tr1= tools_tree2,
#                          tr2 = technology_tree2,
#                          use.edge.length= F,
#                          rotate = T)
# plot(obj, link.type="curved",link.lwd=1,
#      link.lty="solid",link.col=make.transparent("red",
#                                                 0.25))





##################################### 
# Mantel Test
##################################### 

## Geographic vs. phylogenetic difference, i.e. is ‘just’ similarity-by-distance that structures cultural diversity or is transmission (and with it migration?) causal?
library(vegan)
library(geosphere)
library(readr)


Macroregion_mean_coords_df <- readr::read_csv(file = file.path("1_data",
                                                               "Macroregion_mean_coords.csv"))
Macroregion_mean_coords_df$Macro_region_code <- factor(Macroregion_mean_coords_df$Macro_region_code,
                                                       levels = Macroregion_mean_coords_df$Macro_region_code[order(Macroregion_mean_coords_df$Macro_region_code_mean_Lat, decreasing = T)])
Macroregion_mean_coords_df$Macro_region <- factor(Macroregion_mean_coords_df$Macro_region,
                                                  levels = Macroregion_mean_coords_df$Macro_region[order(Macroregion_mean_coords_df$Macro_region_code_mean_Lat, decreasing = T)])

# technology_TaxUnit_by_timeslice_df <- readr::read_delim("MASTERTABLE_v7_all_revised_technology_with_TSrange_uniqueID.csv",
#                                                         ",", escape_double = T, trim_ws = TRUE)
# 
# tech_discreteTS <- readr::read_delim("datasheet_discreteTS_list_technology.csv",
#                                    ",", escape_double = T, trim_ws = TRUE)
# time_dist <- dist(tech_discreteTS[,c("TaxUnit_unique_TS", "Timeslice_discrete")])

# # View(technology_TaxUnit_by_timeslice_df)
technology_coord <- dplyr::left_join(technology,
                               Macroregion_mean_coords_df[,c("Macro_region_code", 
                                                             "Macro_region_code_n", 
                                                             "Macro_region_code_mean_Long", 
                                                             "Macro_region_code_mean_Lat")],
                               by = "Macro_region_code")
tools_coord <- dplyr::left_join(tools,
                                     Macroregion_mean_coords_df[,c("Macro_region_code", 
                                                                   "Macro_region_code_n", 
                                                                   "Macro_region_code_mean_Long", 
                                                                   "Macro_region_code_mean_Lat")],
                                     by = "Macro_region_code")


#geographic data frame - haversine distance
# "technology"
technology_coord_matrix <- as.data.frame(technology_coord[,c("Macro_region_code_mean_Long", 
                                                             "Macro_region_code_mean_Lat", 
                                                             "TaxUnit_unique_TS")])
rownames(technology_coord_matrix) <- technology_coord_matrix$TaxUnit_unique_TS
if(duplicated(technology_coord_matrix$TaxUnit_unique_TS)==T){
  technology_coord_matrix <- technology_coord_matrix[-which(duplicated(technology_coord_matrix$TaxUnit_unique_TS)),]
}
technology_coord_matrix <- technology_coord_matrix[,-which(colnames(technology_coord_matrix) == "TaxUnit_unique_TS")]

tech_d.geo = geosphere::distm(technology_coord_matrix, 
                         fun = distHaversine) # 'Haversine' great circle distance. The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), according to the 'haversine method'. This method assumes a spherical earth, ignoring ellipsoidal effects.
rownames(tech_d.geo) <- rownames(technology_coord_matrix)
colnames(tech_d.geo) <- rownames(technology_coord_matrix)
tech_dist.geo = as.dist(tech_d.geo)

# tools
tools_coord_matrix <- tools_coord[,c("Macro_region_code_mean_Long", 
                                     "Macro_region_code_mean_Lat", 
                                     "TaxUnit_unique_TS")]
rownames(tools_coord_matrix) <- tools_coord_matrix$TaxUnit_unique_TS

if(duplicated(tools_coord_matrix$TaxUnit_unique_TS)==T){
  tools_coord_matrix <- tools_coord_matrix[-which(duplicated(tools_coord_matrix$TaxUnit_unique_TS)),]
}
tools_coord_matrix <- tools_coord_matrix[,-which(colnames(tools_coord_matrix) == "TaxUnit_unique_TS")]

tools_d.geo = geosphere::distm(tools_coord_matrix, 
                              fun = distHaversine) # 'Haversine' great circle distance. The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), according to the 'haversine method'. This method assumes a spherical earth, ignoring ellipsoidal effects.
rownames(tools_d.geo) <- rownames(tools_coord_matrix)
colnames(tools_d.geo) <- rownames(tools_coord_matrix)
tools_dist.geo = as.dist(tools_d.geo)

# outlines
outlines_coord_matrix <- outlines[,c("Long", "Lat", "ARTEFACTNAME")]
outlines_coord_matrix <- na.omit(outlines_coord_matrix)

outlines_coord_matrix <- subset(outlines_coord_matrix, ARTEFACTNAME %in% rownames(outlines_dist_matrix))
outlines_dist_matrix <- outlines_dist_matrix[which(rownames(outlines_dist_matrix) %in% outlines_coord_matrix$ARTEFACTNAME),
                                             which(colnames(outlines_dist_matrix) %in% outlines_coord_matrix$ARTEFACTNAME)]

if(any(duplicated(outlines_coord_matrix$ARTEFACTNAME)==T)){
  outlines_coord_matrix <- outlines_coord_matrix[-which(duplicated(outlines_coord_matrix$ARTEFACTNAME)),]
}
outlines_coord_matrix_names <- outlines_coord_matrix[,which(colnames(outlines_coord_matrix) == "ARTEFACTNAME")]
outlines_coord_matrix <- outlines_coord_matrix[,-which(colnames(outlines_coord_matrix) == "ARTEFACTNAME")]

outlines_d.geo = geosphere::distm(outlines_coord_matrix, 
                               fun = distHaversine) # 'Haversine' great circle distance. The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), according to the 'haversine method'. This method assumes a spherical earth, ignoring ellipsoidal effects.
rownames(outlines_d.geo) <- outlines_coord_matrix_names$ARTEFACTNAME
colnames(outlines_d.geo) <- outlines_coord_matrix_names$ARTEFACTNAME
outlines_dist.geo = as.dist(outlines_d.geo)


# cophenetic vs geographic distance

# technology
mantel_techno <- 
  ecodist::mantel(dist_technology_df_dummy ~ tech_dist.geo,
                  nperm = 1000,
                  nboot = 1000,
                  mrank = T)
mantel_techno
mantel_techno[1]^2

tech_correlog <- vegan::mantel.correlog(D.eco = dist_technology_df_dummy, 
                                        D.geo = tech_dist.geo/1000,
                                        # XY=technology_coord_matrix, 
                                        cutoff=FALSE, 
                                        r.type="spearman", 
                                        nperm=1000,
                                        mult = "hochberg")
tech_correlog   
plot(tech_correlog)
# outdist <- as.vector(as.dist(dist_technology_df_dummy))
# geodist <- as.vector(tech_dist.geo)/1000
# blubb <- data.frame(outdist = outdist,
#                     geodist = geodist)
# gam <- mgcv::gam(outdist ~ geodist, data = blubb)
# summary(gam)

# tools
mantel_tools <- 
ecodist::mantel(dist_tools_df_dummy ~ tools_dist.geo,
                nperm = 1000,
                nboot = 1000,
                mrank = T)
mantel_tools
mantel_tools[1]^2

tools_correlog <- vegan::mantel.correlog(D.eco = dist_tools_df_dummy, 
                                        D.geo = tools_dist.geo/1000,
                                        # XY=tools_coord_matrix, 
                                        cutoff=FALSE, 
                                        r.type="spearman", 
                                        nperm=1000,
                                        mult = "hochberg")
tools_correlog   
plot(tools_correlog)
# outdist <- as.vector(as.dist(dist_tools_df_dummy))
# geodist <- as.vector(tools_dist.geo)/1000
# blubb <- data.frame(outdist = outdist,
#                     geodist = geodist)
# gam <- mgcv::gam(outdist ~ geodist, data = blubb)
# summary(gam)

##########
# tools vs tech dist
vegan::mantel(dist_tools_df_dummy, dist_technology_df_dummy)

vegan::mantel.partial(dist_tools_df_dummy, 
                      dist_technology_df_dummy,
                      tools_dist.geo/1000,
                      permutations=1000)

# expert editor
exped <- data.frame(Expert_editor = tools_coord$Expert_editor,
                    row.names = tools_coord$TaxUnit_unique_TS)
exped_dummy <- dummies::dummy.data.frame(exped)
exped_dummy_jacc <- vegan::vegdist(exped_dummy,
                                   method = "gower")
plot(hclust(exped_dummy_jacc))

vegan::mantel(dist_tools_df_dummy , exped_dummy_jacc)

vegan::mantel.partial(dist_tools_df_dummy, 
                      dist_technology_df_dummy,
                      exped_dummy_jacc,
                      permutations=1000)

vegan::mantel(dist_technology_df_dummy , exped_dummy_jacc)
vegan::mantel.partial(dist_technology_df_dummy,
                      dist_tools_df_dummy, 
                      exped_dummy_jacc,
                      permutations=1000)


##########
# outlines
mantel_outlines <- 
ecodist::mantel(as.dist(outlines_dist_matrix) ~ outlines_dist.geo,
                nperm = 1000,
                nboot = 1000)
mantel_outlines
mantel_outlines[1]^2

outlines_correlog <- vegan::mantel.correlog(D.eco = outlines_dist_matrix, 
                                         D.geo = outlines_dist.geo/1000,
                                         # XY=outlines_coord_matrix, 
                                         cutoff=FALSE, 
                                           r.type="spearman", 
                                         nperm=1000,
                                         mult = "hochberg")
outlines_correlog   
plot(outlines_correlog)


plot(x = as.vector(outlines_dist.geo)/1000,
     y = as.vector(as.dist(outlines_dist_matrix)))

# outdist <- as.vector(as.dist(outlines_dist_matrix))
# geodist <- as.vector(outlines_dist.geo)/1000
# blubb <- data.frame(outdist = outdist,
#                     geodist = geodist)
# gam <- mgcv::gam(outdist ~ geodist, data = blubb)
# summary(gam)
# plot.gam(gam,
#      pages=1,
#      residuals=TRUE,
#      all.terms=TRUE,
#      shade=TRUE,
#      shade.col=2)

#####

tech_mantel_res_df <- as.data.frame(tech_correlog$mantel.res)
tech_mantel_res_df$data <- "Technology"
tools_mantel_res_df <- as.data.frame(tools_correlog$mantel.res)
tools_mantel_res_df$data <- "Tools"
outlines_mantel_res_df <- as.data.frame(outlines_correlog$mantel.res)
outlines_mantel_res_df$data <- "Outlines"

mantel_res_df <- rbind.data.frame(tech_mantel_res_df, 
                                  tools_mantel_res_df,
                                  outlines_mantel_res_df)

mantel_res_df$Significance <- sapply(mantel_res_df$`Pr(corrected)`, 
                                     function(x){
                                       if(x <= 0.05){
                                         "p≤0.05"
                                       } else { "n.s."}
                                     })
mantel_res_df$Significance <- factor(mantel_res_df$Significance,
                                     levels = c("p≤0.05", "n.s."))



mantel_plot <- 
  ggplot(data = mantel_res_df,
         aes(x = class.index,
             y = Mantel.cor)) +
  geom_point(size = 3, aes(shape = Significance)) +
  scale_shape_manual(values=c(15,0)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "Geographic Distance (km)",
       y = "Mantel Correlation",
       fill = "Significance") + 
  scale_x_continuous(breaks = c(500,1000,1500,2000,2500)#labels = scales::comma
  ) +
  facet_wrap(~data,
             ncol = 1)
  # theme(legend.position =  c(0.9,0.9),
  #       legend.background = element_rect(fill="white",
  #                                        size=0.2, linetype="solid", 
  #                                        colour ="black"))
mantel_plot
ggsave(plot = mantel_plot,
       filename = file.path(output_path, "mantel_plot.png"))



################################################


# mantel test time
# technology; keep only the observations that span a single time-slice
only_singleTS_technology <- 
  technology %>% 
  filter(., Timeslice %in% c(1, 2, 3, 4)) %>% 
  select(TaxUnit_unique_TS, Timeslice_range_from) %>% 
  mutate(Timeslice_range_from = as.integer(Timeslice_range_from))

only_singleTS_technology_matrix <- data.frame(Timeslice= only_singleTS_technology$Timeslice_range_from,
                                            row.names = only_singleTS_technology$TaxUnit_unique_TS)
technology_dist.time <- dist(only_singleTS_technology_matrix,
                           upper = F)

technology_dist_matrix_time <- vegan::vegdist(x = technology_df[which(rownames(technology_df)%in%only_singleTS_technology$TaxUnit_unique_TS),],
                                              method = "gower",
                                              upper = F)
mantel_technology_time <- 
  ecodist::mantel(technology_dist_matrix_time ~ technology_dist.time,
                  nperm = 1000,
                  nboot = 1000)
mantel_technology_time
mantel_technology_time[1]^2

technology_correlog_time <- vegan::mantel.correlog(D.eco = technology_dist_matrix_time, 
                                                 D.geo = technology_dist.time,
                                                 # XY=technology_coord_matrix, 
                                                 cutoff=FALSE, 
                                                 r.type="spearman", 
                                                 nperm=1000,
                                                 mult = "hochberg")
# tools; keep only the observations that span a single time-slice
only_singleTS_tools <- 
  tools %>% 
  filter(., Timeslice %in% c(1, 2, 3, 4)) %>% 
  select(TaxUnit_unique_TS, Timeslice_range_from) %>% 
  mutate(Timeslice_range_from = as.integer(Timeslice_range_from))

only_singleTS_tools_matrix <- data.frame(Timeslice= only_singleTS_tools$Timeslice_range_from,
                                              row.names = only_singleTS_tools$TaxUnit_unique_TS)
tools_dist.time <- dist(only_singleTS_tools_matrix,
                             upper = F)

tools_dist_matrix_time <- vegan::vegdist(x = tools_df[which(rownames(tools_df)%in%only_singleTS_tools$TaxUnit_unique_TS),],
                                              method = "gower",
                                              upper = F)
mantel_tools_time <- 
  ecodist::mantel(tools_dist_matrix_time ~ tools_dist.time,
                  nperm = 1000,
                  nboot = 1000)
mantel_tools_time
mantel_tools_time[1]^2

tools_correlog_time <- vegan::mantel.correlog(D.eco = tools_dist_matrix_time, 
                                                   D.geo = tools_dist.time,
                                                   # XY=tools_coord_matrix, 
                                                   cutoff=FALSE, 
                                                   r.type="spearman", 
                                                   nperm=1000,
                                                   mult = "hochberg")
# outlines; keep only the observations that span a single time-slice
only_singleTS_outlines <- 
  outlines %>% 
  filter(., Timeslice %in% c("TS1", "TS2", "TS3", "TS4")) %>% 
  select(ARTEFACTNAME, Timeslice_range_from) %>% 
  mutate(Timeslice_range_from = as.integer(Timeslice_range_from))
only_singleTS_outlines_matrix <- data.frame(Timeslice= only_singleTS_outlines$Timeslice_range_from,
                                            row.names = only_singleTS_outlines$ARTEFACTNAME)
outlines_dist.time <- dist(only_singleTS_outlines_matrix,
                           upper = F)
outlines_dist_matrix_time <- dist(outlines_AR_centered_w_metric_measures_scaled_PCA$x[which(rownames(outlines_AR_centered_w_metric_measures_scaled_PCA$x)%in%only_singleTS_outlines$ARTEFACTNAME),
                                                           c(1:scree_min)],
       upper = F)
mantel_outlines_time <- 
  ecodist::mantel(outlines_dist_matrix_time ~ outlines_dist.time,
                  nperm = 1000,
                  nboot = 1000)
mantel_outlines_time
mantel_outlines_time[1]^2

outlines_correlog_time <- vegan::mantel.correlog(D.eco = outlines_dist_matrix_time, 
                                            D.geo = outlines_dist.time,
                                            # XY=outlines_coord_matrix, 
                                            cutoff=FALSE, 
                                            r.type="spearman", 
                                            nperm=1000,
                                            mult = "hochberg")



tech_mantel_time_res_df <- as.data.frame(technology_correlog_time$mantel.res)
tech_mantel_time_res_df$data <- "Technology"
tools_mantel_time_res_df <- as.data.frame(tools_correlog_time$mantel.res)
tools_mantel_time_res_df$data <- "Tools"
outlines_mantel_time_res_df <- as.data.frame(outlines_correlog_time$mantel.res)
outlines_mantel_time_res_df$data <- "Outlines"

mantel_time_res_df <- na.omit(rbind.data.frame(tech_mantel_time_res_df[,c(1:4,6)], 
                                               tools_mantel_time_res_df[,c(1:4,6)],
                                               outlines_mantel_time_res_df[,c(1:4,6)]))


mantel_time_res_df$Significance <- sapply(mantel_time_res_df$`Pr(Mantel)`, 
                                     function(x){
                                       if(x <= 0.001){
                                         "p≤0.001"
                                       } else if(x <= 0.05){
                                         "p≤0.05"
                                       } else if(x >= 0.05){ 
                                         "n.s."
                                         }
                                     })
mantel_time_res_df$Significance <- factor(mantel_time_res_df$Significance,
                                     levels = c("p≤0.05", "n.s."))


library(ggplot2)
mantel_time_plot <-
  ggplot(data = mantel_time_res_df,
         aes(x = class.index,
             y = Mantel.cor)) +
  geom_point(size = 3, aes(shape = Significance)) +
  scale_shape_manual(values=c(15,0)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "Time-slice distance",
       y = "Mantel Correlation",
       fill = "Significance") +
    facet_wrap(~data,
               ncol = 1)

mantel_time_plot

ggsave(plot = mantel_time_plot,
       filename = file.path(output_path, "mantel_time_plot.png"))














fm1 <- lm(as.dist(outlines_dist_matrix_time) ~ outlines_dist.time)
fm1
summary(fm1)
nlme::Variogram(resid(fm1), (outlines_dist.time))

library(reshape2)
N <- nrow(outlines_dist.time)
melt(as.matrix(outlines_dist.time))

df_dist <- data.frame(x=as.vector(outlines_dist_matrix_time),
                      y=as.vector(outlines_dist.time))

ggplot2::ggplot(data = df_dist) +
  ggplot2::geom_point(ggplot2::aes(x=x,
                                 y=y))
                
     
     
     
     
     
