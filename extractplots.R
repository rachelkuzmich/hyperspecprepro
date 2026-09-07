library(dplyr)
library(terra)
library(data.table)
library(readxl)
library(purrr)

#the goal is to crop all plots
#extracted the 400m2 raster (and a buffered 800m2 version)
#also as individual site_plot tables, and a master table

#set path to path
folder_path <- "R:/Users/rjkuz/transect_mosaic"

#list all tif files in the directory
file_list <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE)

#import and name each file
for (file in file_list) {
  #extract the filename without the extension
  var_name <- tools::file_path_sans_ext(basename(file))
  
  #read the hyperspectral image
  assign(var_name, terra::rast(file), envir = .GlobalEnv)
}

#centre coordinates####

#dynamically parse plot center coordinates from inventory spreadsheet
#ultimately, these hyperspectral plots will be aligned with the inventory

inventory <- "C:/Users/rakuz6966/OneDrive - Norwegian University of Life Sciences/Desktop/0_forests4society_forestinventory.xlsx"

selected_sheets <- c("plot1(ref)", "plot2(dry)",                  
                     "plot3(ref)", "plot4(dry)", "plot5(ref)",                  
                     "plot6(dry)", "plot7(ref)", "plot8(dry) ",                 
                     "plot9(ref)", "plot10(dry)", "plot11(ref)",                 
                     "plot12(dry)", "plot13(ref)", "pot14(dry)",                  
                     "plot15(ref)","plot16(dry)","plot17(ref)",
                     "plot18(dry)", "plot19(ref)","plot20(dry)")

#read metadata and loop tree inventory sheets
meta_sheet <- read_excel(inventory, sheet = "coordinates and instal dates", col_types = "text") %>%
  janitor::clean_names()

combined_inventory <- map_df(
  selected_sheets,
  function(s) {
    read_excel(inventory, sheet = s, col_types = "text")
  },
  .id = "sheet_name"
)

#tidy up metadata
#clean and parse metadata
meta_clean <- meta_sheet %>%
  filter(!grepl("Plot|number|EUREF|UTM", x1, ignore.case = TRUE)) %>%
  filter(!is.na(x1)) %>%
  mutate(
    Plot_number = stringr::str_trim(as.character(x1)),
    site        = stringr::str_to_lower(stringr::str_trim(x2)),
    
    #clean text strings before matching to capture typos
    raw_treatment = stringr::str_to_lower(stringr::str_trim(x3)),
    
    #map raw text values cleanly
    treatment   = case_when(
      raw_treatment %in% c("reference", "ref", "plot1(ref)", "typical") ~ "typical",
      raw_treatment %in% c("dry", "limited")                            ~ "limited",
      TRUE                                                              ~ raw_treatment
    ),
    plot_E      = as.numeric(from_gnss),
    plot_N      = as.numeric(x5) 
  ) %>%
  select(Plot_number, site, treatment, plot_E, plot_N)

#clean the naming tags to match conventions
meta_spatial_prep <- meta_clean %>%
  filter(!is.na(plot_E) & !is.na(plot_N)) %>%
  mutate(
    clean_site = case_when(
      site == "finstad north" ~ "finsted",
      site == "elvål south"   ~ "elval",
      site == "rømskog"       ~ "romskog",
      site == "våler"         ~ "valer",
      site == "åsnes"         ~ "asnes",
      TRUE                    ~ site
    ),
    #generates a clean key strictly ending in _typical or _limited
    plot_name_key = paste(clean_site, treatment, sep = "_")
  )

#build the named list of terra vectors
plots <- list()

for (i in 1:nrow(meta_spatial_prep)) {
  row_data <- meta_spatial_prep[i, ]
  
  #bind the UTM coordinates
  plot_vector <- terra::vect(
    cbind(row_data$plot_E, row_data$plot_N), 
    type = "points",
    crs  = "EPSG:25832" # ETRS89 / UTM zone 32N
  )
  
  #assign to list matching exact lower loop variable structure
  plots[[row_data$plot_name_key]] <- plot_vector
}

cat("Successfully generated named list for", length(plots), "plots with valid EPSG structures.\n")

#calculate the specific radii needed for the target areas
r_400 <- sqrt(400 / pi) # ~11.284 meters
r_800 <- sqrt(800 / pi) # ~15.958 meters

#process each plot
#disaggregate pixels and process each plot with spatial smoothing
for (plot_name in names(plots)) {
  
  #determine which raster corresponds to which plot based on the prefix name
  site_prefix <- sub("_(typical|limited)$", "", plot_name)
  raster_name <- paste0(site_prefix, "_hyperspec_mosaic")
  
  #double check if the target raster actually exists
  if (exists(raster_name, envir = .GlobalEnv)) {
    raster_obj <- get(raster_name, envir = .GlobalEnv)
    point_obj  <- plots[[plot_name]]
    
    #reproject centre point to match the specific raster's CRS
    point_projected <- terra::project(point_obj, terra::crs(raster_obj))
    
    #create the circular buffers in meters
    buf_400 <- terra::buffer(point_projected, width = r_400)
    buf_800 <- terra::buffer(point_projected, width = r_800)
    
    #crop to the bounding box of the 800m2 plot 
    # this prevents R from wasting memory smoothing the entire massive mosaic
    raster_local <- terra::crop(raster_obj, buf_800)
    
    #subdivide 1.25m pixels into 0.25m pixels using bilinear interpolation
    #this just makes for smoother pixels in the crop to plot
    raster_smooth <- terra::disagg(raster_local, fact = 5, method = "bilinear")
    
    #crop and mask the newly smoothed micro-raster to the precise circular shapes
    cropped_400 <- terra::crop(raster_smooth, buf_400, mask = TRUE)
    cropped_800 <- terra::crop(raster_smooth, buf_800, mask = TRUE)
    
    #save the output clipped rasters back into the global environment
    assign(paste0(plot_name, "_400m2"), cropped_400, envir = .GlobalEnv)
    assign(paste0(plot_name, "_800m2"), cropped_800, envir = .GlobalEnv)
    
    message("Successfully extracted smoothed plots for: ", plot_name)
  } else {
    warning("Could not find matching raster object named: ", raster_name)
  }
}

#check resolution and dimensions to confirm success
cat("\n--- Final Validation Metrics ---\n")
print(res(rena_typical_400m2))
print(dim(rena_typical_400m2))

#visual check####
#define the 10 core site prefixes
sites <- c("aremark", "asnes", "eidskog", "elval", "elverum", 
           "finsted", "koppang", "rena", "romskog", "valer")

#loop through each site and generate a 4-panel comparison plot
for (site in sites) {
  
  if (exists(paste0(site, "_typical_400m2"), envir = .GlobalEnv)) {
    
    dev.new(width = 10, height = 10)
    
    #set up a 2x2 layout grid
    par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
    
    #typical 400m2
    obj_name <- paste0(site, "_typical_400m2")
    if (exists(obj_name)) {
      terra::plotRGB(get(obj_name), r=3, g=2, b=1, stretch="lin", 
                     main = paste(site, "Typical - 400m2"))
    }
    
    #typical 800m2
    obj_name <- paste0(site, "_typical_800m2")
    if (exists(obj_name)) {
      terra::plotRGB(get(obj_name), r=3, g=2, b=1, stretch="lin", 
                     main = paste(site, "Typical - 800m2"))
    }
    
    #limited 400m2
    obj_name <- paste0(site, "_limited_400m2")
    if (exists(obj_name)) {
      terra::plotRGB(get(obj_name), r=3, g=2, b=1, stretch="lin", 
                     main = paste(site, "Limited - 400m2"))
    }
    
    #limited 800m2
    obj_name <- paste0(site, "_limited_800m2")
    if (exists(obj_name)) {
      terra::plotRGB(get(obj_name), r=3, g=2, b=1, stretch="lin", 
                     main = paste(site, "Limited - 800m2"))
    }
    
  } else {
    warning("No processed objects found for site: ", site)
  }
}

#save as rasters and tables####
#define output directories
dir_400 <- "R:/Users/rjkuz/transect_mosaic/plots_400m2"
dir_800 <- "R:/Users/rjkuz/transect_mosaic/plots_800m2"

#create the folders if they do not exist
#if (!dir.exists(dir_400)) dir.create(dir_400, recursive = TRUE)
#if (!dir.exists(dir_800)) dir.create(dir_800, recursive = TRUE)

#get names of all objects currently in your global environment
all_objects <- ls(envir = .GlobalEnv)

#filter for SpatRasters and export them
# for (obj_name in all_objects) {
# 
#   #fetch the object from environment
#   obj <- get(obj_name, envir = .GlobalEnv)
# 
#   #check if it's a SpatRaster object before trying to save it, in case there are centre points in there too
#   if (inherits(obj, "SpatRaster")) {
# 
#     #do the 400m2 objects
#     if (grepl("_400m2$", obj_name)) {
#       file_path <- file.path(dir_400, paste0(obj_name, ".tif"))
#       terra::writeRaster(obj, filename = file_path, overwrite = TRUE)
#       message("Saved full SpatRaster data matrix to: ", file_path)
#     }
# 
#     #and now do the 800m2 objects
#     if (grepl("_800m2$", obj_name)) {
#       file_path <- file.path(dir_800, paste0(obj_name, ".tif"))
#       terra::writeRaster(obj, filename = file_path, overwrite = TRUE)
#       message("Saved full SpatRaster data matrix to: ", file_path)
#     }
#   }
# }


#save in tables###
#this makes a table per site_plot
#define output directory for tables
dir_tables <- "R:/Users/rjkuz/transect_mosaic/table_plots_400m2"

#create directory if it does not exist
#if (!dir.exists(dir_tables)) dir.create(dir_tables, recursive = TRUE)

#list all objects in the environment
all_objects <- ls(envir = .GlobalEnv)

#filter, extract data matrices, and save
# for (obj_name in all_objects) {
# 
#   #target only the 400m2 objects, these are the actual plots
#   if (grepl("_400m2$", obj_name)) {
#     obj <- get(obj_name, envir = .GlobalEnv)
# 
#     #double check that it is a SpatRaster spatial data object
#     if (inherits(obj, "SpatRaster")) {
# 
#       #convert raster pixels to a data frame (xy = TRUE keeps coordinates,
#       #na.rm = TRUE drops the masked NA cells outside the circular plot boundary)
#       pixel_table <- terra::as.data.frame(obj, xy = TRUE, na.rm = TRUE)
# 
#       #define output text file path
#       csv_path <- file.path(dir_tables, paste0(obj_name, "_matrix.csv"))
# 
#       #save out to disk
#       write.csv(pixel_table, file = csv_path, row.names = FALSE)
#       message("Exported spectral matrix table: ", csv_path)
#     }
#   }
# }

#make a master table too
#define input and output structures
dir_tables <- "R:/Users/rjkuz/transect_mosaic/table_plots_400m2"
all_objects <- ls(envir = .GlobalEnv)

#initialize an empty list to store individual data frames efficiently
master_list <- list()

#extract tables from environment and append the source tracking column
for (obj_name in all_objects) {
  if (grepl("_400m2$", obj_name)) {
    obj <- get(obj_name, envir = .GlobalEnv)
    
    if (inherits(obj, "SpatRaster")) {
      #extract pixel rows, omitting the background NA cells outside the circle
      pixel_table <- terra::as.data.frame(obj, xy = TRUE, na.rm = TRUE)
      
      #clean up the object name to get the plot identifier (e.g., "aremark_typical")
      plot_id <- sub("_400m2$", "", obj_name)
      
      #add the tracking column at the front of the data matrix
      pixel_table <- cbind(plot_id = plot_id, pixel_table)
      
      #append to our list storage
      master_list[[plot_id]] <- pixel_table
    }
  }
}

#bind all tables together and export
# if (length(master_list) > 0) {
#   # rbindlist handles mismatched column names or variations efficiently
#   master_table <- data.table::rbindlist(master_list, fill = TRUE)
# 
#   output_path <- file.path(dir_tables, "master_plots_400m2_matrix.csv")
#   write.csv(master_table, file = output_path, row.names = FALSE)
# 
#   message("Successfully compiled master matrix table at: ", output_path)
# } else {
#   warning("No 400m2 SpatRaster objects were found in the workspace.")
# }

#STOP here####
#STOP here####
#STOP here####

#this is the previous approach to centre coordinates
#I had originally tried the coordinates this way
#but there may have been errors
# aremark_typical <- terra::vect(cbind(c(11.6690),
#                                      c(59.2290)), 
#                                crs="+proj=longlat")
# 
# aremark_limited <- terra::vect(cbind(c(11.6520),
#                                      c(59.2330)), 
#                                crs="+proj=longlat")
# 
# 
# asnes_typical <- terra::vect(cbind(c(12.357),
#                                    c(60.624)),
#                              crs="+proj=longlat")
# 
# asnes_limited <- terra::vect(cbind(c(12.355),
#                                    c(60.624)),
#                              crs="+proj=longlat")
# 
# 
# eidskog_typical <- terra::vect(cbind(c(11.8280),
#                                      c(59.9240)), 
#                                crs="+proj=longlat")
# 
# eidskog_limited <- terra::vect(cbind(c(11.8320),
#                                      c(59.9230)), 
#                                crs="+proj=longlat")
# 
# 
# 
# elval_typical <- terra::vect(cbind(c(10.8890),
#                                    c(61.9200)), 
#                              crs="+proj=longlat")
# 
# elval_limited <- terra::vect(cbind(c(10.8900),
#                                    c(61.9200)), 
#                              crs="+proj=longlat")
# 
# 
# elverum_typical <- terra::vect(cbind(c(11.4700),
#                                      c(60.9480)), 
#                                crs="+proj=longlat")
# 
# elverum_limited <- terra::vect(cbind(c(11.4670),
#                                      c(60.9470)), 
#                                crs="+proj=longlat")
# 
# 
# finsted_typical <- terra::vect(cbind(c(11.0060),
#                                      c(62.1230)), 
#                                crs="+proj=longlat")
# 
# finsted_limited <- terra::vect(cbind(c(11.0050),
#                                      c(62.1240)), 
#                                crs="+proj=longlat")
# 
# 
# koppang_typical <- terra::vect(cbind(c(11.1630),
#                                      c(61.5930)), 
#                                crs="+proj=longlat")
# 
# koppang_limited <- terra::vect(cbind(c(11.1610),
#                                      c(61.5930)), 
#                                crs="+proj=longlat")
# 
# 
# rena_typical <- terra::vect(cbind(c(11.2220),
#                                   c(61.2740)), 
#                             crs="+proj=longlat")
# 
# rena_limited <- terra::vect(cbind(c(11.2170),
#                                   c(61.2720)), 
#                             crs="+proj=longlat")
# 
# 
# romskog_typical <- terra::vect(cbind(c(11.9170),
#                                      c(59.7140)), 
#                                crs="+proj=longlat")
# 
# romskog_limited <- terra::vect(cbind(c(11.9210),
#                                      c(59.7160)), 
#                                crs="+proj=longlat")
# 
# 
# valer_typical <- terra::vect(cbind(c(10.8740),
#                                    c(59.5160)), 
#                              crs="+proj=longlat")
# 
# valer_limited <- terra::vect(cbind(c(10.8730),
#                                    c(59.5100)), 
#                              crs="+proj=longlat")
# 
# 
# #combine point vectors into a single named list
# plots <- list(
#   aremark_typical = aremark_typical, aremark_limited = aremark_limited,
#   asnes_typical   = asnes_typical,   asnes_limited   = asnes_limited,
#   eidskog_typical = eidskog_typical, eidskog_limited = eidskog_limited,
#   elval_typical   = elval_typical,   elval_limited   = elval_limited,
#   elverum_typical = elverum_typical, elverum_limited = elverum_limited,
#   finsted_typical = finsted_typical, finsted_limited = finsted_limited,
#   koppang_typical = koppang_typical, koppang_limited = koppang_limited,
#   rena_typical    = rena_typical,    rena_limited    = rena_limited,
#   romskog_typical = romskog_typical, romskog_limited = romskog_limited,
#   valer_typical   = valer_typical,   valer_limited   = valer_limited
# )
# 
# #calculate the specific radii needed for the target areas
# r_400 <- sqrt(400 / pi) # ~11.284 meters
# r_800 <- sqrt(800 / pi) # ~15.958 meters
# 
# #process each plot
# #disaggregate pixels and process each plot with spatial smoothing
# for (plot_name in names(plots)) {
#   
#   #determine which raster corresponds to which plot based on the prefix name
#   site_prefix <- sub("_(typical|limited)$", "", plot_name)
#   raster_name <- paste0(site_prefix, "_hyperspec_mosaic")
#   
#   #double check if the target raster actually exists
#   if (exists(raster_name, envir = .GlobalEnv)) {
#     raster_obj <- get(raster_name, envir = .GlobalEnv)
#     point_obj  <- plots[[plot_name]]
#     
#     #reproject centre point to match the specific raster's CRS
#     point_projected <- terra::project(point_obj, terra::crs(raster_obj))
#     
#     #create the circular buffers in meters
#     buf_400 <- terra::buffer(point_projected, width = r_400)
#     buf_800 <- terra::buffer(point_projected, width = r_800)
#     
#     # crop to the bounding box of the 800m2 plot 
#     #this prevents R from wasting memory smoothing the entire massive mosaic
#     raster_local <- terra::crop(raster_obj, buf_800)
#     
#     #subdivide 1.25m pixels into 0.25m pixels using bilinear interpolation
#     #this just makes for smoother pixels in the crop to plot
#     raster_smooth <- terra::disagg(raster_local, fact = 5, method = "bilinear")
#     
#     #crop and mask the newly smoothed micro-raster to the precise circular shapes
#     cropped_400 <- terra::crop(raster_smooth, buf_400, mask = TRUE)
#     cropped_800 <- terra::crop(raster_smooth, buf_800, mask = TRUE)
#     
#     #save the output clipped rasters back into your global environment
#     assign(paste0(plot_name, "_400m2"), cropped_400, envir = .GlobalEnv)
#     assign(paste0(plot_name, "_800m2"), cropped_800, envir = .GlobalEnv)
#     
#     message("Successfully extracted smoothed plots for: ", plot_name)
#   } else {
#     warning("Could not find matching raster object named: ", raster_name)
#   }
# }
# 
# #check resolution and dimensions
# res(rena_typical_400m2)
# dim(rena_typical_400m2)

#above, pixels were made smalled
#but, to keep the original pixel size - use this code 

# for (plot_name in names(plots)) {
#   
#   #determine which raster corresponds to this plot based on the prefix name
#   # (e.g., "aremark_typical" looks for an object named "aremark_hyperspec_mosaic")
#   site_prefix <- sub("_(typical|limited)$", "", plot_name)
#   raster_name <- paste0(site_prefix, "_hyperspec_mosaic")
#   
#   #check if the target raster actually exists in your R environment
#   if (exists(raster_name, envir = .GlobalEnv)) {
#     raster_obj <- get(raster_name, envir = .GlobalEnv)
#     point_obj  <- plots[[plot_name]]
#     
#     #reproject the longlat point to match the specific raster's CRS
#     point_projected <- terra::project(point_obj, terra::crs(raster_obj))
#     
#     #create the circular buffers in meters
#     buf_400 <- terra::buffer(point_projected, width = r_400)
#     buf_800 <- terra::buffer(point_projected, width = r_800)
#     
#     #crop and mask the rasters to the circles (mask = TRUE makes the outside NA)
#     cropped_400 <- terra::crop(raster_obj, buf_400, mask = TRUE)
#     cropped_800 <- terra::crop(raster_obj, buf_800, mask = TRUE)
#     
#     #save the output clipped rasters back into your global environment
#     assign(paste0(plot_name, "_400m2"), cropped_400, envir = .GlobalEnv)
#     assign(paste0(plot_name, "_800m2"), cropped_800, envir = .GlobalEnv)
#     
#     message("Successfully extracted plots for: ", plot_name)
#   } else {
#     warning("Could not find matching raster object named: ", raster_name)
#   }
# }
# #check resolution and dimensions
# res(rena_typical_400m2)
# dim(rena_typical_400m2)

