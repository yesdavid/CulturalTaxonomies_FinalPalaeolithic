#####################################################################################
########################### Part 2: create output and plots #########################

output_path <- file.path("3_output", "10_NACmaps")

library(ggplot2)


####################################### MAPS ######################################## 

# set manual shapes + colours for 
# region = number
# taxunit = shape 

datasheet_discreteTS_list <- readRDS(file = file.path("1_data", "datasheet_discreteTS_list.RDS"))

datasheet_discreteTS_list$sites$TaxUnit_shape <- factor(rep(0, times = nrow(datasheet_discreteTS_list$sites)), levels = c(0:15))

for(current_macroregion in unique(datasheet_discreteTS_list$sites$Macro_region_code)) {
  current_macroregion_subset <- subset(datasheet_discreteTS_list$sites, Macro_region_code == current_macroregion)
  current_unique_TaxUnits <- unique(current_macroregion_subset$TaxUnit)
  
  for (current_unique_TaxUnits_index in 1:length(current_unique_TaxUnits)) {
    a <- current_unique_TaxUnits[current_unique_TaxUnits_index]
    datasheet_discreteTS_list$sites[which(datasheet_discreteTS_list$sites[,"Macro_region_code"] == current_macroregion & datasheet_discreteTS_list$sites[,"TaxUnit"] == a),"TaxUnit_shape"] <- 
      factor(current_unique_TaxUnits_index, levels = c(0:15))
  }
}


# get base map data for extent
world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
clipper_europe <- as(raster::extent(-9.969482, 26.97, 35.74455, 57.5), "SpatialPolygons")
sp::proj4string(clipper_europe) <- sp::CRS(sp::proj4string(world))
world_clip <- raster::intersect(world, clipper_europe)
world_clip_f <- ggplot2::fortify(world_clip)

base_map <- 
  ggplot() +
  geom_polygon(data = world_clip_f,
               aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey") +
  coord_map() +
  theme_classic() +  
  xlab("Longitude") +
  ylab("Latitude") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


##################### map with sites in macroregion color
map_all_sites_RegionCol <-
  base_map + 
  labs(color = "Macroregion") +
  geom_text(data = datasheet_discreteTS_list$sites,
            aes(x = Long, y = Lat,
                label=Macro_region_number),
            color = "black",
            # alpha = 0.7,
            size = 2)

map_all_sites_RegionCol

# ggsave(plot = map_all_sites_RegionCol,
#        filename = file.path(output_path, "map_all_sites.svg"), 
#        width = 30, height = 16, units = "cm")
ggsave(plot = map_all_sites_RegionCol,
       filename = file.path(output_path, "map_all_sites.png"), 
       width = 30, height = 16, units = "cm")

#####################
# sort macroregions by latitude

Macroregion_mean_coords_df <- readr::read_csv(file = file.path("1_data", "Macroregion_mean_coords.csv"))

Macroregion_mean_coords_df$Macro_region_number_text <- 
  factor(paste(Macroregion_mean_coords_df$Macro_region_number,
               " ", 
               Macroregion_mean_coords_df$Macro_region),
         levels = paste(Macroregion_mean_coords_df$Macro_region_number,
                        " ", 
                        Macroregion_mean_coords_df$Macro_region)[order(Macroregion_mean_coords_df$Macro_region_code_mean_Lat, 
                                                                       decreasing = T)])

map_all_sites_plus_center <-
  map_all_sites_RegionCol + 
  geom_label(data = Macroregion_mean_coords_df,
             aes(x = Macro_region_code_mean_Long, 
                 y = Macro_region_code_mean_Lat,
                 label = Macro_region_number), 
             size = 5) + 
  geom_point(data = Macroregion_mean_coords_df,
             aes(x = Macro_region_code_mean_Long, 
                 y = Macro_region_code_mean_Lat,
                 fill = Macro_region_number_text),
             alpha = 0) +
  labs(fill = "Macroregion")

map_all_sites_plus_center

# ggsave(plot = map_all_sites_plus_center,
#        filename = file.path(output_path, "map_all_sites_plus_center_numbers.svg"), 
#        width = 30, height = 16, units = "cm")
ggsave(plot = map_all_sites_plus_center,
       filename = file.path(output_path, "map_all_sites_plus_center_numbers.png"), 
       width = 30, height = 16, units = "cm")


##################### map sites with NACs

taxgroups <- readr::read_csv(file.path("1_data", "taxgroups_v6_31Jan23.csv"))

data_w_taxgroupsNAC <-
  dplyr::left_join(datasheet_discreteTS_list$sites,
                   taxgroups, 
                   by = c("TaxUnit_unique", "TaxUnit"))

data_w_taxgroupsNAC$higher_order <- forcats::fct_reorder(data_w_taxgroupsNAC$higher_order, data_w_taxgroupsNAC$chrono_order)

# change TS to roman numerals
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(1, "I", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(2, "II", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(3, "III", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(4, "IV", data_w_taxgroupsNAC$Timeslice_discrete)

data_w_taxgroupsNAC$Timeslice_discrete <- factor(data_w_taxgroupsNAC$Timeslice_discrete,
                                                 levels = c("I", "II", "III", "IV"))


map_with_NACs_per_TS <-
  base_map +
  geom_point(data = data_w_taxgroupsNAC,
             aes(x = Long, y = Lat,
                 fill = higher_order),
             alpha = 0.9,
             size = 3,
             shape = 21) +
  facet_wrap(~Timeslice_discrete)  +
  scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
                                 "Epigravettian" = "#CC79A7",
                                 "Azilian" = "#E69F00", # Azilian
                                 "ABP" = "#F0E442", # ABP/FMG
                                 "LTP" = "#0072B2", # LTP/TPC
                                 "FBBT/LBI" = "#56B4E9", # FBBT
                                 "Mesolithic" = "#009E73", # Meso
                                 "NA" = "#999999" # NA
  )) +
  labs(color = "NAC")

map_with_NACs_per_TS

# ggsave(plot = map_with_NACs_per_TS,
#        filename = file.path(output_path, "map_with_NACs_per_TS.svg"),
#        width = 30, height = 24, units = "cm")
ggsave(plot = map_with_NACs_per_TS,
       filename = file.path(output_path, "map_with_NACs_per_TS.png"),
       width = 30, height = 24, units = "cm")


map_with_NACs <-
  base_map +
  geom_point(data = data_w_taxgroupsNAC,
             aes(x = Long, y = Lat,
                 fill = higher_order),
             alpha = 0.9,
             size = 3,
             shape = 21) +
  scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
                                 "Epigravettian" = "#CC79A7",
                                 "Azilian" = "#E69F00", # Azilian
                                 "ABP" = "#F0E442", # ABP/FMG
                                 "LTP" = "#0072B2", # LTP/TPC
                                 "FBBT/LBI" = "#56B4E9", # FBBT
                                 "Mesolithic" = "#009E73", # Meso
                                 "NA" = "#999999" # NA
  )) +
  labs(color = "NAC")

map_with_NACs

# ggsave(plot = map_with_NACs,
#        filename = file.path(output_path, "map_with_NACs.svg"),
#        width = 30, height = 24, units = "cm")
ggsave(plot = map_with_NACs,
       filename = file.path(output_path, "map_with_NACs.png"),
       width = 30, height = 24, units = "cm")

map_with_NACs_wrapNAC <-
  base_map +
  geom_point(data = data_w_taxgroupsNAC,
             aes(x = Long, y = Lat,
                 fill = higher_order),
             alpha = 0.9,
             size = 3,
             shape = 21) +
  facet_wrap(~higher_order,
             ncol = 2) +
  scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
                                 "Epigravettian" = "#CC79A7",
                                 "Azilian" = "#E69F00", # Azilian
                                 "ABP" = "#F0E442", # ABP/FMG
                                 "LTP" = "#0072B2", # LTP/TPC
                                 "FBBT/LBI" = "#56B4E9", # FBBT
                                 "Mesolithic" = "#009E73", # Meso
                                 "NA" = "#999999" # NA
  )) +
  # labs(color = "NAC") +
  theme(legend.position = "none")

map_with_NACs_wrapNAC

# ggsave(plot = map_with_NACs_wrapNAC,
#        filename = file.path(output_path, "map_with_NACs_wrapNAC.svg"),
#        width = 30, height = 24, units = "cm")
ggsave(plot = map_with_NACs_wrapNAC,
       filename = file.path(output_path, "map_with_NACs_wrapNAC.png"),
       width = 30, height = 24, units = "cm")

###################
# what are the NAs?
base_map +
  geom_point(data = subset(data_w_taxgroupsNAC, is.na(higher_order)),
             aes(x = Long, y = Lat,
                 fill = TaxUnit),
             alpha = 0.9,
             size = 3,
             shape = 21)+
  ggtitle("NA") +
  theme(plot.title = element_text(hjust = 0.5))

library(dplyr)
subset(data_w_taxgroupsNAC, is.na(higher_order)) %>% 
  dplyr::select(TaxUnit, TaxUnit_unique, higher_order) %>% 
  unique() #%>% 
  # readr::write_csv(file.path("1_data", "taxgroups_v5_26Aug22_missing.csv"))
###################

# MAP without NAs

map_with_NACs_wrapNAC_noNA <-
  base_map +
  geom_point(data = subset(data_w_taxgroupsNAC, !is.na(higher_order)),
             aes(x = Long, y = Lat,
                 fill = higher_order),
             alpha = 0.9,
             size = 3,
             shape = 21) +
  facet_wrap(~higher_order,
             ncol = 2) +
  scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
                                 "Epigravettian" = "#CC79A7",
                                 "Azilian" = "#E69F00", # Azilian
                                 "ABP" = "#F0E442", # ABP/FMG
                                 "LTP" = "#0072B2", # LTP/TPC
                                 "FBBT/LBI" = "#56B4E9", # FBBT
                                 "Mesolithic" = "#009E73", # Meso
                                 "NA" = "#999999" # NA
  )) +
  # labs(color = "NAC") +
  theme(legend.position = "none")

map_with_NACs_wrapNAC_noNA

# ggsave(plot = map_with_NACs_wrapNAC,
#        filename = file.path(output_path, "map_with_NACs_wrapNAC.svg"),
#        width = 30, height = 24, units = "cm")
ggsave(plot = map_with_NACs_wrapNAC_noNA,
       filename = file.path(output_path, "map_with_NACs_wrapNAC_noNA.png"),
       width = 30, height = 24, units = "cm")

########
map_with_NACs_wrapNAC_noNA_timeGrid <- 
  base_map +
  geom_point(data = subset(data_w_taxgroupsNAC %>% 
                             mutate(higher_order = forcats::fct_rev(higher_order)), 
                           !is.na(higher_order)),
             aes(x = Long, y = Lat,
                 fill = higher_order),
             alpha = 0.9,
             size = 2,
             shape = 21) +
  facet_grid(higher_order~Timeslice_discrete) +
  scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
                                 "Epigravettian" = "#CC79A7",
                                 "Azilian" = "#E69F00", # Azilian
                                 "ABP" = "#F0E442", # ABP/FMG
                                 "LTP" = "#0072B2", # LTP/TPC
                                 "FBBT/LBI" = "#56B4E9", # FBBT
                                 "Mesolithic" = "#009E73", # Meso
                                 "NA" = "#999999" # NA
  )) +
  # labs(color = "NAC") +
  theme(legend.position = "none")

ggsave(plot = map_with_NACs_wrapNAC_noNA_timeGrid,
       filename = file.path(output_path, "map_with_NACs_wrapNAC_noNA_timeGrid.png"),
       width = 23, height = 23, units = "cm")

# data_w_taxgroupsNAC$Timeslice_discrete_rev <- factor(data_w_taxgroupsNAC$Timeslice_discrete,
#                                                     levels = rev(c("I", "II", "III", "IV")))
# base_map +
#   geom_point(data = subset(data_w_taxgroupsNAC, !is.na(higher_order)),
#              aes(x = Long, y = Lat,
#                  fill = higher_order),
#              alpha = 0.9,
#              size = 3,
#              shape = 21) +
#   facet_grid(Timeslice_discrete_rev~higher_order) +
#   scale_colour_manual(values = c("Magdalenian" = "#D55E00", # Magdalenian
#                                  "Epigravettian" = "#CC79A7",
#                                  "Azilian" = "#E69F00", # Azilian
#                                  "ABP" = "#F0E442", # ABP/FMG
#                                  "LTP" = "#0072B2", # LTP/TPC
#                                  "FBBT/LBI" = "#56B4E9", # FBBT
#                                  "Mesolithic" = "#009E73", # Meso
#                                  "NA" = "#999999" # NA
#   )) +
#   # labs(color = "NAC") +
#   theme(legend.position = "none")
