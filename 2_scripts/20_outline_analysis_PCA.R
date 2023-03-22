# outline analysis

## set paths
outlines_input_path <- file.path("1_data", 
                                 "1511NAC_Database", 
                                 "outlines")

output_folder <- file.path("3_output", 
                           "20_outline_analysis_PCA")
dir.create(output_folder)

## load outlines from 1511NAC_Database
outlines_AR_centered_w_metric_measures <- readRDS(file = file.path(outlines_input_path,
                                                                   "outlines_AR_centered_w_metric_measures.RDS"))
  
outlines_BR_centered_w_metric_measures <- readRDS(file = file.path(outlines_input_path,
                                                                   "outlines_BR_centered_w_metric_measures.RDS")) 

outlines_ES_centered_w_metric_measures <- readRDS(file = file.path(outlines_input_path,
                                                                   "outlines_ES_centered_w_metric_measures.RDS")) 

## combine AR, BR, and ES outlines in a single list
outlines_centered_w_metric_measures_list <- list(outlines_AR_centered_w_metric_measures = outlines_AR_centered_w_metric_measures,
                                                 outlines_BR_centered_w_metric_measures = outlines_BR_centered_w_metric_measures,
                                                 outlines_ES_centered_w_metric_measures = outlines_ES_centered_w_metric_measures)


# run PCA on all outlines
outlines_centered_w_metric_measures_scaled_PCA_list <- list()
for(current_outlines_centered_w_metric_measures_index in names(outlines_centered_w_metric_measures_list)) {
  
  current_outlines_centered_w_metric_measures <- outlines_centered_w_metric_measures_list[[current_outlines_centered_w_metric_measures_index]]
  
  current_outlines_centered_w_metric_measures_scaled <- Momocs::coo_scale(current_outlines_centered_w_metric_measures)
  
  # harmonic calibration. Estimates the number of harmonics required for the Fourier methods implemented in Momocs. 
  ## This is the only step in this section that produces the data we need in the subsequent step.
  current_outlines_centered_w_metric_measures_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(current_outlines_centered_w_metric_measures_scaled, 
                                                                                                           plot = F)
  # current_outlines_centered_w_metric_measures_scaled_harmonics$minh
  
  # efourier
  current_outlines_centered_w_metric_measures_scaled_efourier <- Momocs::efourier(current_outlines_centered_w_metric_measures_scaled,
                                                                             nb.h = as.matrix(current_outlines_centered_w_metric_measures_scaled_harmonics[["minh"]])[[4,1]], # choses number of harmonics for 99.9%
                                                                             norm = F) 
  # PCA
  outlines_centered_w_metric_measures_scaled_PCA_list[[current_outlines_centered_w_metric_measures_index]] <- 
    Momocs::PCA(current_outlines_centered_w_metric_measures_scaled_efourier) # PCA on Coe objects, using prcomp.
  
  # save PCA object of each outline type 
  saveRDS(outlines_centered_w_metric_measures_scaled_PCA_list[[current_outlines_centered_w_metric_measures_index]],
          file = file.path(output_folder, paste0(current_outlines_centered_w_metric_measures_index, "_scaled_PCA.RDS")))
}

outlines_AR_centered_w_metric_measures_scaled_PCA <- outlines_centered_w_metric_measures_scaled_PCA_list[["outlines_AR_centered_w_metric_measures"]]
# outlines_AR_centered_w_metric_measures_scaled_PCA <- readRDS(file.path(output_folder, "outlines_AR_centered_w_metric_measures_scaled_PCA.RDS"))


#####################

# discretize AR outlines by timeslice for later per-discrete time slide analysis.

current_datasheet <- outlines_AR_centered_w_metric_measures_scaled_PCA$fac
current_datasheet$Timeslice <- gsub(" ", "", current_datasheet$Timeslice) # remove empty spaces 
current_datasheet$Timeslice <- gsub("TS", "", current_datasheet$Timeslice) # remove empty spaces 
current_datasheet$Timeslice_discrete <- 9999 # create new column with place holder values

current_datasheet_discretizedTS_list <- list()

for (current_sheet_row in 1:nrow(current_datasheet)) {
  
  current_row <- current_datasheet[current_sheet_row,]
  
  if(nchar(current_row$Timeslice) > 1){ # if there is more than one timeslice
    
    current_row_all_timeslices <- strsplit(current_row$Timeslice, 
                                           split = ",")[[1]]
    if(nchar(current_row_all_timeslices[1])>1){
      current_row_all_timeslices <- strsplit(current_row$Timeslice, 
                                             split = "")[[1]]
    }
    
    new_row_with_discrete_ts_list <- list()
    for(current_timeslice in 1:length(current_row_all_timeslices)) {
      current_row$Timeslice_discrete <- as.numeric(current_row_all_timeslices[current_timeslice])
      
      new_row_with_discrete_ts_list[[current_timeslice]] <- current_row
    }
    
    current_datasheet_discretizedTS_list[[current_sheet_row]] <- do.call(rbind.data.frame, new_row_with_discrete_ts_list)
    
    
  } else {
    current_row$Timeslice_discrete <- as.numeric(current_row$Timeslice)
    current_datasheet_discretizedTS_list[[current_sheet_row]] <- current_row
  }
  
  
}
current_df <- do.call(rbind.data.frame, current_datasheet_discretizedTS_list)

readr::write_csv(x = current_df,
                 path = file.path(output_folder, "datasheet_discreteTS_list_outlines_AR_wMetrics.csv"))
#####################





#####################
# plot PCA
library(ggplot2)

## PC contribution
minimum_no_of_pcs_outlines_AR <- Momocs::scree_min(outlines_AR_centered_w_metric_measures_scaled_PCA,
                                              prop = 0.99) 
minimum_no_of_pcs_outlines_AR

outlines_AR_screeplot <- Momocs::scree_plot(outlines_AR_centered_w_metric_measures_scaled_PCA,
                                       nax = 1:minimum_no_of_pcs_outlines_AR)
outlines_AR_screeplot
ggsave(outlines_AR_screeplot,
       filename = file.path(output_folder, "outlines_AR_screeplot.png"),
       width = 10,
       height = 5)

## PC shape variation
gg <- Momocs::PCcontrib(outlines_AR_centered_w_metric_measures_scaled_PCA,
                        nax = 1:5,
                        sd.r = c(-2,-1,0,1,2)) 

outlines_AR_pccontrib <- gg$gg + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

outlines_AR_pccontrib
ggsave(outlines_AR_pccontrib,
       filename = file.path(output_folder, "outlines_AR_pccontrib.png"),
       width = 5,
       height = 5)

PCA_df <- as.data.frame(outlines_AR_centered_w_metric_measures_scaled_PCA$x)
PCA_df$ARTEFACTNAME <- rownames(PCA_df)

ggplot(data = PCA_df) +
  geom_point(aes(x = PC1,
                 y = PC2)) +
  theme_bw()

ggplot(data = PCA_df) +
  geom_point(aes(x = PC2,
                 y = PC3)) +
  theme_bw()

ggplot(data = PCA_df) +
  geom_point(aes(x = PC1,
                 y = PC3)) +
  theme_bw()




