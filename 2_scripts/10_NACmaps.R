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
             alpha = 0)

map_all_sites_plus_center 

ggsave(plot = map_all_sites_plus_center +
         theme(legend.position = "none"),
       filename = file.path(output_path, "Fig_3_map_all_sites_plus_center_numbers_NO_legend.png"), 
       width = 20, height = 16, units = "cm")


ggsave(plot = map_all_sites_plus_center +
         labs(fill = "Macroregion"),
       filename = file.path(output_path, "Fig_3_map_all_sites_plus_center_numbers_WITH_legend.png"), 
       width = 30, height = 16, units = "cm")


##################### map sites with NACs
taxgroups <- readr::read_csv(file.path("1_data", "taxgroups_v9_13Mar23.csv"))
########################################################
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
                                                levels = rev(forcats::fct_reorder(unique(taxgroups_fac_hexcodes$higher_order),
                                                                                  unique(taxgroups_fac_hexcodes$chrono_order))))

########################################################
data_w_taxgroupsNAC <-
  dplyr::left_join(datasheet_discreteTS_list$sites,
                   taxgroups)

data_w_taxgroupsNAC$higher_order <- forcats::fct_reorder(data_w_taxgroupsNAC$higher_order, data_w_taxgroupsNAC$chrono_order)

# change TS to Roman numerals
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(1, "I", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(2, "II", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(3, "III", data_w_taxgroupsNAC$Timeslice_discrete)
data_w_taxgroupsNAC$Timeslice_discrete <- gsub(4, "IV", data_w_taxgroupsNAC$Timeslice_discrete)

data_w_taxgroupsNAC$Timeslice_discrete <- factor(data_w_taxgroupsNAC$Timeslice_discrete,
                                                 levels = c("I", "II", "III", "IV"))

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
  scale_colour_manual(values = higher_order_hexcodes_manual) +
  # labs(color = "NAC") +
  theme(legend.position = "none")

map_with_NACs_wrapNAC_noNA_timeGrid

ggsave(plot = map_with_NACs_wrapNAC_noNA_timeGrid,
       filename = file.path(output_path, "SI_map_with_NACs_wrapNAC_noNA_timeGrid.png"),
       width = 23, height = 23, units = "cm")

