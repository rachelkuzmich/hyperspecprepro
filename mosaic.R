library(terra)
library(prospectr)
library(signal)
library(RStoolbox)
library(ggplot2)
library(tidyr)
library(dplyr)

#Aremark####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/aremark/CS")

#list reflectance files
aremark_files <- list.files(pattern="_atm\\.bsq$")
print(aremark_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(aremark_files)){
  f <- aremark_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])
# plot(flightlines[[5]][[50]])
# plot(flightlines[[6]][[50]])
# plot(flightlines[[7]][[50]])
# plot(flightlines[[8]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])
# plot(flightlines_sg[[5]][[20]])
# plot(flightlines_sg[[6]][[20]])
# plot(flightlines_sg[[7]][[20]])
# plot(flightlines_sg[[8]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])
# plot(flightlines_normalized[[5]][[20]])
# plot(flightlines_normalized[[6]][[20]])
# plot(flightlines_normalized[[7]][[20]])
# plot(flightlines_normalized[[8]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic, 
#   filename = "R:/Users/rjkuz/transect_mosaic/aremark_hyperspec_mosaic.tif", 
#   gdal = c("COMPRESS=LZW"), 
#   overwrite = TRUE
# )
#done

#Asnes####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/asnes/CS")

#list reflectance files
asnes_files <- list.files(pattern="_atm\\.bsq$")
print(asnes_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(asnes_files)){
  f <- asnes_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/asnes_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done

#Eidskog####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/eidskog/CS")

#list reflectance files
eidskog_files <- list.files(pattern="_atm\\.bsq$")
print(eidskog_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(eidskog_files)){
  f <- eidskog_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/eidskog_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done

#Elval####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/elval/CS")

#list reflectance files
elval_files <- list.files(pattern="_atm\\.bsq$")
print(elval_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(elval_files)){
  f <- elval_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/elval_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done

#Elverum####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/elverum/CS")

#list reflectance files
elverum_files <- list.files(pattern="_atm\\.bsq$")
print(elverum_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(elverum_files)){
  f <- elverum_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/elverum_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done


#Finsted####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/finsted/CS")

#list reflectance files
finsted_files <- list.files(pattern="_atm\\.bsq$")
print(finsted_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(finsted_files)){
  f <- finsted_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/finsted_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done


#Koppang####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/koppang/CS")

#list reflectance files
koppang_files <- list.files(pattern="_atm\\.bsq$")
print(koppang_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(koppang_files)){
  f <- koppang_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/koppang_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done


#Rena####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/rena/CS")

#list reflectance files
rena_files <- list.files(pattern="_atm\\.bsq$")
print(rena_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(rena_files)){
  f <- rena_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/rena_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done

#Romskog####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/romskog/CS")

#list reflectance files
romskog_files <- list.files(pattern="_atm\\.bsq$")
print(romskog_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(romskog_files)){
  f <- romskog_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/romskog_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done

#Valer####

#prep reflectance####
#set working directory
setwd("R:/Data/Forest4Society/2023/3_RemoteSensingData/valer/CS")

#list reflectance files
valer_files <- list.files(pattern="_atm\\.bsq$")
print(valer_files) #check

#extract wavelengths
get_wavelengths <- function(hdr_file){
  
  txt <- readLines(hdr_file)
  
  start <- grep("^wavelength", txt, ignore.case=TRUE)
  end <- start + grep("}", txt[start:length(txt)])[1] - 1
  
  wl_lines <- txt[start:end]
  
  wl_txt <- paste(wl_lines, collapse=" ")
  wl_txt <- sub(".*\\{","", wl_txt)
  wl_txt <- sub("\\}.*","", wl_txt)
  
  as.numeric(trimws(strsplit(wl_txt,",")[[1]]))
}

#load flightlignes
flightlines <- list()
orig_bandnames <- list()

for(i in seq_along(valer_files)){
  f <- valer_files[i]
  r <- rast(f)
  r <- r / 10000
  hdr <- sub("\\.bsq$", ".hdr", f)
  wl <- get_wavelengths(hdr)
  orig_bandnames[[i]] <- names(r)
  names(r) <- paste0("wl_", round(wl,2))
  flightlines[[i]] <- r
}

length(flightlines) #check
plot(flightlines[[1]][[50]])
# plot(flightlines[[2]][[50]])
# plot(flightlines[[3]][[50]])
# plot(flightlines[[4]][[50]])

#Czech Globe report notes that there are flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  #original names
  orig_names <- orig_bandnames[[i]] 
  
  #remove interpolated bands
  interp <- startsWith(orig_names, "band*") 
  r_clean <- r[[!interp]]
  
  #extract wavelengths
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) 
  
  #remove spectral tails
  keep <- wl_clean >= 450 & wl_clean <= 2300 
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#Savitzky-Golay####
#signal package: https://cran.r-project.org/web/packages/signal/index.html
#prospectr package: https://antoinestevens.github.io/prospectr/
#note different parameterization in this analysis https://doi.org/10.1016/j.rse.2025.114907 
#see SavitzkyGolay_testing R script for details on package and parameter testing and implementation

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#define segments
#VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#fully vectorized SG function
#no row loops, pixel-by-pixel was taking too long
spectral_smooth_vectorized <- function(r, groups) {
  #extract data matrix (Rows = pixels, Cols = bands)
  v_in <- terra::values(r, mat = TRUE)
  v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #parameters
  w_size <- 5
  p_order <- 2
  
  #loop through the five distinct spectral segments
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    n_bands <- length(band_idx)
    
    if (n_bands >= w_size) {
      #extract the standard (5x5) SG filter weights
      sg_filter <- signal::sgolay(p = p_order, n = w_size)
      
      #build a full-segment transformation matrix (n_bands x n_bands)
      #this maps the filter window perfectly across internal bands and edge cases
      W <- matrix(0, nrow = n_bands, ncol = n_bands)
      
      #fill the edge case weights (first and last (w-1)/2 rows)
      k <- (w_size - 1) / 2
      for (i in 1:k) {
        W[i, 1:w_size] <- sg_filter[i, ]
        W[n_bands - k + i, (n_bands - w_size + 1):n_bands] <- sg_filter[w_size - k + i, ]
      }
      
      #fill the internal moving window weights
      for (i in (k + 1):(n_bands - k)) {
        W[i, (i - k):(i + k)] <- sg_filter[k + 1, ]
      }
      
      #compute the entire segment for millions of rows
      #multiply the pixel data by the transpose of our transformation matrix
      v_out[, band_idx] <- segment_data %*% t(W)
      
    } else {
      v_out[, band_idx] <- segment_data
    }
  }
  
  #rebuild the SpatRaster
  r_out <- r
  terra::values(r_out) <- v_out
  return(r_out)
}

#loop with progress tracker
total_images <- length(flightlines_filter)
start_total_time <- Sys.time()

flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
  cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
  
  start_single_time <- Sys.time()
  smoothed_img <- spectral_smooth_vectorized(flightlines_filter[[i]], groups = group_ids)
  end_single_time <- Sys.time()
  
  run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
  cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
  
  return(smoothed_img)
})

end_total_time <- Sys.time()
total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

#SG check####
plot(flightlines_sg[[1]][[20]])
# plot(flightlines_sg[[2]][[20]])
# plot(flightlines_sg[[3]][[20]])
# plot(flightlines_sg[[4]][[20]])

#quick check to compare raw vs SG for ONE pixel
#extract wavelengths and segment groups
my_wavelengths <- c(
  seq(455, 498, length.out = 4),   # Blue
  seq(512, 598, length.out = 7),   # Green
  seq(612, 669, length.out = 5),   # Red
  seq(683, 740, length.out = 5),   # Red Edge
  seq(755, 883, length.out = 10),  # NIR
  seq(1018, 1118, length.out = 7), # SWIR 1
  seq(1212, 1318, length.out = 8), # SWIR 2
  seq(1498, 1722, length.out = 16),# SWIR 3
  seq(2068, 2292, length.out = 16) # SWIR 4
)
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#extract raw and smoothed pixel values
img_raw <- flightlines_filter[[1]]
img_sm  <- flightlines_sg[[1]]

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2)) #can change the pixel here

raw_pixel <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel  <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel
)

#diagnostic plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  #raw data shown as discrete points to highlight sensor chatter
  geom_point(aes(y = Raw, color = "Raw Data"), size = 1.8, alpha = 0.6) +
  #smooth tracking lines spanning edge-to-edge within segments
  geom_line(aes(y = Smoothed, color = "Vectorized SG (w=5, p=2)"), size = 1.0) +
  scale_color_manual(values = c("Raw Data" = "gray40", "Vectorized SG (w=5, p=2)" = "#009E73")) +
  theme_minimal() +
  labs(
    title = "Spectral Trajectory: Raw vs. Vectorized Savitzky-Golay",
    subtitle = "Verifying edge corrections and noise suppression across data gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#brightness normalization####
#again, this is vectorized for speed
brightness_normalize_raster <- function(r) {
  #extract the smoothed data matrix (Rows = pixels, Cols = 78 bands)
  v_in <- terra::values(r, mat = TRUE)
  
  #calculate the brightness norm for each individual pixel (row-wise)
  #this computes the square root of the sum of squared values across all bands
  pixel_norms <- sqrt(rowSums(v_in^2, na.rm = TRUE))
  
  #avoid division by zero for completely black or masked pixels
  pixel_norms[pixel_norms == 0] <- 1
  
  #divide each band by its pixel's brightness norm
  v_norm <- v_in / pixel_norms
  
  #rebuild the SpatRaster with the normalized values
  r_out <- r
  terra::values(r_out) <- v_norm
  return(r_out)
}

#normalize the entire list of smoothed flightlines
total_images <- length(flightlines_sg)
start_time <- Sys.time()

flightlines_normalized <- lapply(seq_along(flightlines_sg), function(i) {
  cat(paste0("\n[", Sys.time(), "] Brightness Normalizing Flightline ", i, " of ", total_images, "...\n"))
  
  norm_img <- brightness_normalize_raster(flightlines_sg[[i]])
  return(norm_img)
})

end_time <- Sys.time()
cat(paste0("\n=== Normalization complete in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds! ===\n"))

#norm check####
plot(flightlines_normalized[[1]][[20]])
# plot(flightlines_normalized[[2]][[20]])
# plot(flightlines_normalized[[3]][[20]])
# plot(flightlines_normalized[[4]][[20]])

#diagnostic plot again
#extract values from all three stages of the pipeline
img_raw  <- flightlines_filter[[1]]
img_sm   <- flightlines_sg[[1]]
img_norm <- flightlines_normalized[[1]] # Extract the newly normalized raster

mid_cell <- cellFromRowCol(img_raw, round(nrow(img_raw)/2), round(ncol(img_raw)/2))

raw_pixel  <- as.vector(terra::values(img_raw, mat = TRUE)[mid_cell, ])
sm_pixel   <- as.vector(terra::values(img_sm, mat = TRUE)[mid_cell, ])
norm_pixel <- as.vector(terra::values(img_norm, mat = TRUE)[mid_cell, ])

#structure data frame for plotting
plot_df <- data.frame(
  Wavelength = my_wavelengths,
  Group      = group_ids,
  Raw        = raw_pixel,
  Smoothed   = sm_pixel,
  Normalized = norm_pixel
)

#generate the 3-Way Diagnostic Plot (Using Dual Y-Axes)
#brightness values are scaled by a factor of 4 to match the primary reflectance axis visually
scale_factor <- 4

ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40",
    "2. Vectorized SG (w=5, p=2)" = "#009E73",
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(
    name = "Reflectance Value (Raw & Smoothed)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Normalized Unit Scale (Vector Length = 1)")
  ) +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel from raw signal to brightness-leveled output",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(color = "#D55E00"),
    axis.text.y.right = element_text(color = "#D55E00")
  )

#replot the 3-Way Diagnostic Plot on a single shared Y-axis
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG (w=5, p=2)" = "#009E73", 
    "3. Brightness Normalized" = "#D55E00"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(
    title = "Complete Spectral Pipeline: Smoothing & Normalization",
    subtitle = "Tracking a single forest canopy pixel plotted directly on a single, shared axis",
    x = "Wavelength (nm)",
    color = "Processing Stage"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 10)
  )

#mosaicking####
#robust threshold converter to Catch near-zero background values
clean_null_backgrounds <- function(r) {
  #extract values matrix
  v <- terra::values(r, mat = TRUE)
  
  #find rows where the maximum value across all 78 bands is extremely low
  #this isolates the background mask from real forest data
  is_background <- rowSums(v > 0.001, na.rm = TRUE) == 0
  
  #force those entire background rows to be explicit NAs
  v[is_background, ] <- NA
  
  #rebuild raster
  r_clean <- r
  terra::values(r_clean) <- v
  return(r_clean)
}

#run the background cleaner
cat(paste0("[", Sys.time(), "] Stripping near-zero backgrounds dynamically...\n"))
flightlines_na_corrected <- lapply(flightlines_normalized, clean_null_backgrounds)

#recreate the collection and restitch the mosaic
flightline_collection <- terra::sprc(flightlines_na_corrected)

cat(paste0("[", Sys.time(), "] Re-stitching master mosaic...\n"))
forest_mosaic <- terra::mosaic(flightline_collection, fun = "mean")

#check
terra::plotRGB(
  x = forest_mosaic,
  r = 22, g = 14, b = 7, # CIR Composite
  stretch = "lin",
  main = "Fixed CASI-SASI Forest Mosaic"
)

terra::plotRGB(
  x = forest_mosaic,
  r = 13, g = 8, b = 3, # RGB Composite
  stretch = "hist",
  main = "Fixed CASI-SASI Forest Mosaic"
)

#export the mosaic###
#use LZW compression to reduce the hyperspectral file size
# terra::writeRaster(
#   x = forest_mosaic,
#   filename = "R:/Users/rjkuz/transect_mosaic/valer_hyperspec_mosaic.tif",
#   gdal = c("COMPRESS=LZW"),
#   overwrite = TRUE
# )
#done