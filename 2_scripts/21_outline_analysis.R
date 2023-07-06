# set input folder path for PCA data
input_folder <- file.path("3_output", 
                          "20_outline_analysis_PCA")

# load PCA data
outlines_AR_centered_w_metric_measures_scaled_PCA <- readRDS(file.path(input_folder,
                                                                       "outlines_AR_centered_w_metric_measures_scaled_PCA.RDS"))

# discrete time-slices
datasheet_discreteTS_list_outlines <- readr::read_delim(file.path(input_folder, 
                                                                  "datasheet_discreteTS_list_outlines_AR_wMetrics.csv"),
                                                        delim = ",", 
                                                        escape_double = T, 
                                                        trim_ws = TRUE)
# # stratified subset
# set.seed(123)
# subset_TaxUnit_Region_TSdisc <- 
# as.data.frame(splitstackshape::stratified(datasheet_discreteTS_list_outlines,
#                                           group = c("TaxUnit", "Region", "Timeslice_discrete"),
#                                           size = 10,
#                                           replace = F))

output_path <- file.path("3_output", "21_outline_analysis")
dir.create(output_path)

######################################################
# DISPARITY
######################################################
# custom bins based on period
rownames_DATASETS <- list()
for(i in unique(datasheet_discreteTS_list_outlines$Timeslice_discrete)){
  rownames_DATASETS[[i]] <- as.character(subset(datasheet_discreteTS_list_outlines, 
                                                Timeslice_discrete == i)$ARTEFACTNAME)
}

library(dispRity)
TS_subsets <- custom.subsets(outlines_AR_centered_w_metric_measures_scaled_PCA$x, 
                             group = rownames_DATASETS)

TS_boot <- boot.matrix(TS_subsets, bootstraps = 1000)
TS_disp <- dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

readr::write_csv(summary(TS_disp),
                 file = file.path(output_path, "outlines_AR_disparity_results.csv"))


# Wilcox.test
test.dispRity(TS_disp,
              test = wilcox.test,
              comparisons = "pairwise",
              correction = "bonferroni")

ttest_res <- do.call(cbind.data.frame,
                     test.dispRity(TS_disp,
                                   test = wilcox.test,
                                   comparisons = "pairwise",
                                   correction = "bonferroni"))
ttest_res$comparison <- rownames(ttest_res)

readr::write_csv(ttest_res[,c(3,1,2)],
                 file = file.path(output_path, "outlines_AR_disparity_results_pairwiseTtest_bonf.csv"))

TS_names <- names(TS_disp$disparity)
disparity_df_list <- list()
for(i in TS_names){
  disparity_df_list[[i]] <- data.frame(Period = paste0("TS", i, 
                                                       "\n(n=",nrow(TS_disp$subsets[[i]]$elements),")"),
                                       disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                       nelements = nrow(TS_disp$subsets[[i]]$elements),
                                       TS = i)
}
disparity_df_TSdiscrete_outlines_AR_perTShard <- do.call(rbind.data.frame, disparity_df_list)

disparity_df_TSdiscrete_outlines_AR_perTShard$Period <- gsub("TS1", "I", disparity_df_TSdiscrete_outlines_AR_perTShard$Period)
disparity_df_TSdiscrete_outlines_AR_perTShard$Period <- gsub("TS2", "II", disparity_df_TSdiscrete_outlines_AR_perTShard$Period)
disparity_df_TSdiscrete_outlines_AR_perTShard$Period <- gsub("TS3", "III", disparity_df_TSdiscrete_outlines_AR_perTShard$Period)
disparity_df_TSdiscrete_outlines_AR_perTShard$Period <- gsub("TS4", "IV", disparity_df_TSdiscrete_outlines_AR_perTShard$Period)

text_size <- 18

disparity_TSdiscrete_outlines_AR_ggplot_perTShard <- 
  ggplot(data = disparity_df_TSdiscrete_outlines_AR_perTShard, aes(x = Period, y = disparity)) +
  geom_violin(aes(fill = TS)) + 
  geom_boxplot(notch = T, width = 0.1, fill = "white", color = "black") +
  theme_bw() +
  ggtitle(NULL) +
  xlab("Time-slice") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = text_size, face = "bold"),
        axis.text=element_text(size=text_size), #,face="bold"
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size=text_size), # angle = 45, 
        axis.title.y = element_text(vjust = 0),
        axis.title=element_text(size=text_size)) +
  scale_fill_grey(start = 0.6,
                  end = 0.9) +
  guides(color = FALSE, fill = FALSE)

disparity_TSdiscrete_outlines_AR_ggplot_perTShard

ggsave(filename = file.path(output_path, "Fig_9_outlines_AR_disparity_per_TS.png"),
       plot = disparity_TSdiscrete_outlines_AR_ggplot_perTShard,
       width = 7,
       height = 7)

svg(filename=file.path(output_path , "Fig_9_outlines_AR_disparity_per_TS.svg"), width = 7, height = 7)
disparity_TSdiscrete_outlines_AR_ggplot_perTShard
dev.off()

###########################################################################################################




















