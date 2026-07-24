library(terra)
library(prospectr)
library(signal)
library(RStoolbox)
library(ggplot2)
library(tidyr)
library(dplyr)

#using Aremark to test Savitzky–Golay for spectral smoothing packages and parameters

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

#CG report notes flightlines to be removed from CASI only
#this analysis is using their fused CS dataset
#I assume that they fused the data with only the appropriate flightlines

#remove tails and interpolated bands####
flightlines_filter <- vector("list", length(flightlines))
for(i in seq_along(flightlines)){
  r <- flightlines[[i]]
  
  orig_names <- orig_bandnames[[i]] #original names
  
  interp <- startsWith(orig_names, "band*") #remove interpolated bands
  r_clean <- r[[!interp]]
  
  wl_clean <- as.numeric(sub("wl_", "", names(r_clean))) #extract wavelengths
  
  keep <- wl_clean >= 450 & wl_clean <= 2300 #remove spectral tails
  
  flightlines_filter[[i]] <- r_clean[[keep]]
}

#check
flightlines_filter[[1]]
nlyr(flightlines_filter[[1]]) #78 bands, looks good!

#SG with prospectr package: https://antoinestevens.github.io/prospectr/####

#checks before running
r <-  flightlines_filter[[1]]
v <- terra::values(r, mat = TRUE)
dim(v)
ncell(r) #should match nrow(v) from dim(v)
nlyr(r) #should match ncol(v) from dim

#note parameterization here https://doi.org/10.1016/j.rse.2025.114907 

#testing all together, not considering segments
#extract a single real spectrum from your first flightline, from the center of the raster to avoid edges
img <- flightlines_filter[[1]]
mid_row <- round(nrow(img) / 2) #changed the /N value to test different pixels
mid_col <- round(ncol(img) / 2)#changed the /N value to test different pixels
real_pixel <- as.vector(terra::values(img, mat = TRUE)[cellFromRowCol(img, mid_row, mid_col), ])
#create a base data frame for plotting
wavelength_axis <- 1:length(real_pixel) 
results_df <- data.frame(Wavelength = wavelength_axis, Raw = real_pixel)

#grid search loop over all valid combinations (w = 3,5,7,9,11; p = 2,3,4)
window_sizes <- c(3, 5, 7, 9, 11)
poly_orders  <- c(2, 3, 4)
for (w in window_sizes) {
  for (p in poly_orders) {
    #skip mathematically invalid combinations (polynomial must be smaller than window)
    if (p >= w) next
    #run Savitzky-Golay on our single pixel matrix
    sm_pixel <- prospectr::savitzkyGolay(X = matrix(real_pixel, nrow = 1), m = 0, p = p, w = w)
    #calculate pad size to realign the trimmed tails
    k <- (length(real_pixel) - length(sm_pixel)) / 2
    padded_pixel <- c(rep(NA, k), as.vector(sm_pixel), rep(NA, k))
    #save column with a clear name (e.g., "w5_p2")
    col_name <- paste0("w", w, "_p", p)
    results_df[[col_name]] <- padded_pixel
  }
}

#pivot data for ggplot
df_plot <- results_df %>%
  pivot_longer(cols = -Wavelength, names_to = "Setting", values_to = "Reflectance")

#diagnostic plot (as facet_wrap)
ggplot(df_plot, aes(x = Wavelength, y = Reflectance)) +
  geom_point(data = filter(df_plot, Setting == "Raw"), color = "gray60", alpha = 0.6, size = 1.2) +
  geom_line(data = filter(df_plot, Setting != "Raw"), aes(color = Setting), size = 0.8) +
  facet_wrap(~Setting, ncol = 3) + # Organizes all 13 combinations into a clean grid
  theme_minimal() +
  labs(
    title = "Real Data Grid Search: Savitzky-Golay Parameters",
    subtitle = "Reviewing valid combinations of window size (w) and polynomial order (p)",
    x = "Band Number / Wavelength",
    y = "Reflectance Value"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

#diagnostic plot (all together)
ggplot(df_plot, aes(x = Wavelength, y = Reflectance)) +
  #raw data points plotted as a background anchor
  geom_point(data = filter(df_plot, Setting == "Raw"), color = "gray40", alpha = 0.5, size = 1.5) +
  #smooth lines overlaid together
  geom_line(data = filter(df_plot, Setting != "Raw"), aes(color = Setting), size = 0.8) +
  theme_minimal() +
  labs(
    title = "Overlaid Savitzky-Golay Parameter Comparison",
    subtitle = "Comparing all 13 valid smoothing combinations against raw forest data",
    x = "Band Number / Wavelength",
    y = "Reflectance Value",
    color = "Parameters"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

#testing with segments because SG smooths over gaps (where interpolated bands were removed)
#build wavelengths vector
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

#assign group IDs based on the five continuous segments
group_ids <- c(
  rep(1, 31), # VNIR continuous block (4 + 7 + 5 + 5 + 10 bands)
  rep(2, 7),  # SWIR 1
  rep(3, 8),  # SWIR 2
  rep(4, 16), # SWIR 3
  rep(5, 16)  # SWIR 4
)

#extract  pixel data
img <- flightlines_filter[[1]]
mid_cell <- cellFromRowCol(img, round(nrow(img)/2), round(ncol(img)/2))
real_pixel <- as.vector(terra::values(img, mat = TRUE)[mid_cell, ])

results_df <- data.frame(Wavelength = my_wavelengths, Group = group_ids, Raw = real_pixel)

#loop function using strict check (segment length must be strictly greater than w)
run_segmented_sg <- function(pixel_vec, groups, w, p) {
  out_vec <- rep(NA, length(pixel_vec))
  
  for (g in unique(groups)) {
    idx <- which(groups == g)
    seg <- pixel_vec[idx]
    
    # CRITICAL: segment length must be strictly greater than window size (length > w)
    if (length(seg) > w && p < w) {
      sm <- as.vector(savitzkyGolay(X = matrix(seg, nrow = 1), m = 0, p = p, w = w))
      k <- (length(seg) - length(sm)) / 2
      out_vec[idx[(k + 1):(length(seg) - k)]] <- sm
    }
  }
  return(out_vec)
}

#grid search over configurations
configs <- list(
  c(3,2), c(5,2), c(7,2), c(9,2), c(11,2), 
  c(5,3), c(7,3), c(9,3), c(11,3), 
  c(5,4), c(7,4), c(9,4), c(11,4)
)

for (cfg in configs) {
  w <- cfg[1]; p <- cfg[2]
  col_name <- paste0("w", w, "_p", p)
  results_df[[col_name]] <- run_segmented_sg(real_pixel, group_ids, w, p)
}

#pivot and plot diagnostics
df_plot <- results_df %>% 
  pivot_longer(cols = -c(Wavelength, Group), names_to = "Setting", values_to = "Reflectance")

ggplot(df_plot, aes(x = Wavelength, y = Reflectance, group = interaction(Setting, Group))) +
  geom_point(data = filter(df_plot, Setting == "Raw"), color = "gray60", alpha = 0.4) +
  geom_line(data = filter(df_plot, Setting != "Raw"), aes(color = Setting), size = 0.8) +
  theme_minimal() +
  labs(
    title = "Corrected Segmented Savitzky-Golay Diagnostics",
    subtitle = "Visible-to-NIR grouped continuously; shorter SWIR windows auto-skipped when invalid",
    x = "Wavelength (nm)",
    y = "Reflectance"
  )

#second version of plot as facet_wrap
#create a separate data frame for the raw points so they appear as a background on every facet
df_raw_bg <- results_df %>% 
  select(Wavelength, Group, Reflectance = Raw)

#filter out the "Raw" category from the line dataset for cleaner faceting
df_lines_only <- df_plot %>% 
  filter(Setting != "Raw")

ggplot() +
  #background raw points on every facet
  geom_point(data = df_raw_bg, aes(x = Wavelength, y = Reflectance, group = Group), 
             color = "gray75", alpha = 0.6, size = 1.2) +
  #segmented smooth lines
  geom_line(data = df_lines_only, aes(x = Wavelength, y = Reflectance, color = Setting, group = Group), 
            size = 0.8) +
  #split into a clean grid by configuration
  facet_wrap(~Setting, ncol = 3) + 
  theme_minimal() +
  labs(
    title = "Faceted Segmented Savitzky-Golay Diagnostics",
    subtitle = "Comparing each valid configuration independently against raw forest data",
    x = "Wavelength (nm)",
    y = "Reflectance Value"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(1, "lines")
  )

#SG with signal package: https://cran.r-project.org/web/packages/signal/index.html####
#setup exact wavelengths and continuous segment groups 
#same as above, not re-run if the above was ran
# my_wavelengths <- c(
#   seq(455, 498, length.out = 4),   # Blue
#   seq(512, 598, length.out = 7),   # Green
#   seq(612, 669, length.out = 5),   # Red
#   seq(683, 740, length.out = 5),   # Red Edge
#   seq(755, 883, length.out = 10),  # NIR
#   seq(1018, 1118, length.out = 7), # SWIR 1
#   seq(1212, 1318, length.out = 8), # SWIR 2
#   seq(1498, 1722, length.out = 16),# SWIR 3
#   seq(2068, 2292, length.out = 16) # SWIR 4
# )
# group_ids <- c(
#   rep(1, 31), # VNIR continuous block (4 + 7 + 5 + 5 + 10 bands)
#   rep(2, 7),  # SWIR 1
#   rep(3, 8),  # SWIR 2
#   rep(4, 16), # SWIR 3
#   rep(5, 16)  # SWIR 4
# )
# 
# #extract the diagnostic pixel spectrum
# img <- flightlines_filter[]
# mid_cell <- cellFromRowCol(img, round(nrow(img)/2), round(ncol(img)/2))
# real_pixel <- as.vector(terra::values(img, mat = TRUE)[mid_cell, ])
# 
# results_df <- data.frame(Wavelength = my_wavelengths, Group = group_ids, Raw = real_pixel)

#segmented loop function with signal::sgolayfilt (safely covers edge cases)
#signal automatically handles edges by fitting progressively smaller, asymmetric polynomials as it approaches segment boundaries
run_segmented_signal <- function(pixel_vec, groups, w, p) {
  out_vec <- rep(NA, length(pixel_vec))
  
  for (g in unique(groups)) {
    idx <- which(groups == g)
    seg <- pixel_vec[idx]
    
    # signal requires window length 'n' to be strictly greater than polynomial order 'p'
    # It also handles edge points automatically over the full segment length
    if (length(seg) >= w && p < w) {
      out_vec[idx] <- signal::sgolayfilt(seg, p = p, n = w)
    }
  }
  return(out_vec)
}

#grid search configuration matrix
#same as above, not re-run if the above was ran
# configs <- list(
#   c(3,2), c(5,2), c(7,2), c(9,2), c(11,2), 
#   c(5,3), c(7,3), c(9,3), c(11,3), 
#   c(5,4), c(7,4), c(9,4), c(11,4)
# )

for (cfg in configs) {
  w <- cfg[1]; p <- cfg[2]
  col_name <- paste0("w", w, "_p", p)
  results_df[[col_name]] <- run_segmented_signal(real_pixel, group_ids, w, p)
}

#reshape data for plotting
df_plot <- results_df %>% 
  pivot_longer(cols = -c(Wavelength, Group), names_to = "Setting", values_to = "Reflectance")

df_raw_bg <- results_df %>% 
  select(Wavelength, Group, Reflectance = Raw)

df_lines_only <- df_plot %>% 
  filter(Setting != "Raw")

#generate the edge-corrected facet plot
ggplot() +
  geom_point(data = df_raw_bg, aes(x = Wavelength, y = Reflectance, group = Group), 
             color = "gray75", alpha = 0.6, size = 1.2) +
  geom_line(data = df_lines_only, aes(x = Wavelength, y = Reflectance, color = Setting, group = Group), 
            size = 0.8) +
  facet_wrap(~Setting, ncol = 3) + 
  theme_minimal() +
  labs(
    title = "Faceted Segmented 'Signal' Savitzky-Golay Diagnostics",
    subtitle = "Edge-corrected curves running across full segment lengths without gaps",
    x = "Wavelength (nm)",
    y = "Reflectance Value"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(1, "lines")
  )

#w5_p2 is best, there are no artifacts on edges
#lines are stable right up to the edges, smooths entire segments without introducing weird flares

#implementation tests####
# #this takes too long because it is running pixel-by pixel
# #define the continuous segments
# # VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
# group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))
# 
# #optimized smoothing function (based on earlier diagnostics)
# spectral_smooth_signal <- function(r, groups) {
#   #extract data matrix (Rows = pixels, Cols = bands)
#   v_in <- terra::values(r, mat = TRUE)
#   v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
#   
#   #loop and smooth each continuous segment independently
#   for (g in unique(groups)) {
#     band_idx <- which(groups == g)
#     segment_data <- v_in[, band_idx, drop = FALSE]
#     n_bands <- length(band_idx)
#     
#     #parameters
#     w_size <- 5
#     p_order <- 2
#     #m = 0 for smoothing is the default
#     
#     #fallback adjust if an individual segment ever drops below window size
#     #it should not happen with these data, but just in case the data change in some future analysis
#     if (n_bands <= w_size) {
#       w_size <- if (n_bands %% 2 == 0) n_bands - 1 else n_bands
#       w_size <- max(3, w_size)
#     }
#     
#     #fast row-wise application across all raster pixels
#     if (n_bands >= w_size) {
#       sm_segment <- t(apply(segment_data, 1, function(pixel_row) {
#         signal::sgolayfilt(pixel_row, p = p_order, n = w_size)
#       }))
#       
#       #populate the master output matrix
#       v_out[, band_idx] <- sm_segment
#     } else {
#       #if segment is too small for any filter, preserve raw values
#       v_out[, band_idx] <- segment_data
#     }
#   }
#   
#   #rebuild the SpatRaster object with smoothed data
#   r_out <- r
#   terra::values(r_out) <- v_out
#   return(r_out)
# }

#apply across flightlines list
# flightlines_sg <- lapply(flightlines_filter, function(img) {
#   spectral_smooth_signal(img, groups = group_ids)
# })

#alternate version with progress tracker
#also takes too long
#first, get total number of images to process
# total_images <- length(flightlines_filter)
# start_total_time <- Sys.time()
# 
# flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
#   #print progress update to console
#   cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
#   
#   start_single_time <- Sys.time()
#   
#   #run the smoothing function on the current image
#   smoothed_img <- spectral_smooth_signal(flightlines_filter[[i]], groups = group_ids)
#   
#   #print completion stats for this single image
#   end_single_time <- Sys.time()
#   run_duration <- round(difftime(end_single_time, start_single_time, units = "mins"), 2)
#   cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " minutes.\n"))
#   
#   return(smoothed_img)
# })
# 
# end_total_time <- Sys.time()
# total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
# cat(paste0("\n=== All processing complete! Total elapsed time: ", total_duration, " minutes. ===\n"))

#trying a faster version
#instead of pixel-by-pixel, rewrite the function to perform vectorized matrix math
#continuous segment grouping for your 78 bands
#group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

#optimized matrix-vectorized smoothing function
# spectral_smooth_signal_fast <- function(r, groups) {
#   # Extract full data matrix (Rows = pixels, Cols = bands)
#   v_in <- terra::values(r, mat = TRUE)
#   v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
#   
#   #loop through segments (only 5 iterations total, instead of millions)
#   for (g in unique(groups)) {
#     band_idx <- which(groups == g)
#     segment_data <- v_in[, band_idx, drop = FALSE]
#     n_bands <- length(band_idx)
#     
#     w_size <- 5
#     p_order <- 2
#     
#     #mathematical fallback check
#     if (n_bands <= w_size) {
#       w_size <- if (n_bands %% 2 == 0) n_bands - 1 else n_bands
#       w_size <- max(3, w_size)
#     }
#     
#     if (n_bands >= w_size) {
#       #extract the exact filter coefficient matrix from signal package (speed fix)
#       #this matrix handles the asymmetric edge-case weights
#       sg_filter_matrix <- signal::sgolay(p = p_order, n = w_size)
#       
#       #perform ultra-fast matrix multiplication across all millions of pixels instantly
#       # %*% handles the internal smoothing, but transpose (t) to keep the matrix shapes aligned
#       sm_segment <- t(sg_filter_matrix %*% t(segment_data))
#       
#       #populate output matrix
#       v_out[, band_idx] <- sm_segment
#     } else {
#       v_out[, band_idx] <- segment_data
#     }
#   }
#   
#   #rebuild the SpatRaster object
#   r_out <- r
#   terra::values(r_out) <- v_out
#   return(r_out)
# }
# 
# #loop with progress tracker
# total_images <- length(flightlines_filter)
# start_total_time <- Sys.time()
# 
# flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
#   cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
#   
#   start_single_time <- Sys.time()
#   smoothed_img <- spectral_smooth_signal_fast(flightlines_filter[[i]], groups = group_ids)
#   end_single_time <- Sys.time()
#   
#   run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
#   cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
#   
#   return(smoothed_img)
# })
# 
# end_total_time <- Sys.time()
# total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
# cat(paste0("\n=== Complete! Total time: ", total_duration, " minutes. ===\n"))

#still too slow
#30 mins and the first flightline hasn't even run
# #there are gaps due to interpolated bands which were removed
# #need to pass each segment with prospectr::savitzkyGolay
# #then use signal::sgolayfilt on the segment edges to fill in the missing tail pieces
# #best of both worlds: faster and smoothed edges without boundaries
# 
# #define continuous segments for the 78 bands
# # VNIR (31 bands), SWIR 1 (7 bands), SWIR 2 (8 bands), SWIR 3 (16 bands), SWIR 4 (16 bands)
# group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))
# 
# #matrix-Vectorized smoothing function with edge corrections
# spectral_smooth_hybrid_fast <- function(r, groups) {
#   #extract full data matrix (Rows = pixels, Cols = bands)
#   v_in <- terra::values(r, mat = TRUE)
#   v_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
#   
#   #set parameters
#   w_size <- 5
#   p_order <- 2
#   
#   #loop through segments (there are 5 segments/iterations total)
#   for (g in unique(groups)) {
#     band_idx <- which(groups == g)
#     segment_data <- v_in[, band_idx, drop = FALSE]
#     n_bands <- length(band_idx)
#     
#     # run prospectr's fast vectorized filter for the segment core
#     sm_core <- prospectr::savitzkyGolay(X = segment_data, m = 0, p = p_order, w = w_size)
#     
#     #calculate how many bands are trimmed at the tails
#     k <- (n_bands - ncol(sm_core)) / 2
#     smoothed_cols <- (k + 1):(n_bands - k)
#     
#     #insert the fast-smoothed center core into the output matrix
#     v_out[, band_idx[smoothed_cols]] <- sm_core
#     
#     #clean up the edges using signal package (only loop over the 2+2 edge bands)
#     edge_cols <- setdiff(1:n_bands, smoothed_cols)
#     
#     if (length(edge_cols) > 0) {
#       #apply signal row-by-row but ONLY to the k edge positions to maximize speed
#       edge_data <- t(apply(segment_data, 1, function(pixel_row) {
#         full_smooth <- signal::sgolayfilt(pixel_row, p = p_order, n = w_size)
#         return(full_smooth[edge_cols])
#       }))
#       
#       v_out[, band_idx[edge_cols]] <- edge_data
#     }
#   }
#   
#   #rebuild the SpatRaster object
#   r_out <- r
#   terra::values(r_out) <- v_out
#   return(r_out)
# }
# 
# #loop with progress tracker
# total_images <- length(flightlines_filter)
# start_total_time <- Sys.time()
# 
# flightlines_sg <- lapply(seq_along(flightlines_filter), function(i) {
#   cat(paste0("\n[", Sys.time(), "] Processing Flightline ", i, " of ", total_images, "...\n"))
#   
#   start_single_time <- Sys.time()
#   smoothed_img <- spectral_smooth_hybrid_fast(flightlines_filter[[i]], groups = group_ids)
#   end_single_time <- Sys.time()
#   
#   run_duration <- round(difftime(end_single_time, start_single_time, units = "secs"), 2)
#   cat(paste0("--> Finished Flightline ", i, " in ", run_duration, " seconds.\n"))
#   
#   return(smoothed_img)
# })
# 
# end_total_time <- Sys.time()
# total_duration <- round(difftime(end_total_time, start_total_time, units = "mins"), 2)
# cat(paste0("\n=== Complete! Total processing time: ", total_duration, " minutes. ===\n"))

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
  #raw data points
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  
  #smoothed SG line
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  
  #brightness Normalized line (scaled for the primary axis)
  geom_line(aes(y = Normalized * scale_factor, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  
  #map custom colors
  scale_color_manual(values = c(
    "Raw Data" = "gray40", 
    "Vectorized SG (w=5, p=2)" = "#009E73", 
    "Brightness Normalized" = "#D55E00"
  )) +
  
  #configure dual Y-axes to handle the different value scales cleanly
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
  #raw data points
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.5) +
  
  #smoothed SG line
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG (w=5, p=2)"), size = 1.0) +
  
  #normalized brightness  line
  geom_line(aes(y = Normalized, color = "3. Brightness Normalized"), 
            size = 1.0, linetype = "dashed") +
  
  #map custom colors
  scale_color_manual(values = c(
    "Raw Data" = "gray40", 
    "Vectorized SG (w=5, p=2)" = "#009E73", 
    "Brightness Normalized" = "#D55E00"
  )) +
  
  #clean single Y-axis definition
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

#segment-wise brightness normalization function test (vectorized to make it faster)
brightness_normalize_segmented <- function(r, groups) {
  v_in <- terra::values(r, mat = TRUE)
  v_norm_out <- matrix(NA, nrow = nrow(v_in), ncol = ncol(v_in))
  
  #normalize each continuous segment independently
  for (g in unique(groups)) {
    band_idx <- which(groups == g)
    segment_data <- v_in[, band_idx, drop = FALSE]
    
    #calculate the localized vector norm for this segment only
    segment_norms <- sqrt(rowSums(segment_data^2, na.rm = TRUE))
    segment_norms[segment_norms == 0] <- 1
    
    #scale only the bands belonging to this segment
    v_norm_out[, band_idx] <- segment_data / segment_norms
  }
  
  r_out <- r
  terra::values(r_out) <- v_norm_out
  return(r_out)
}

#normalize the entire list of smoothed flightlines
group_ids <- c(rep(1, 31), rep(2, 7), rep(3, 8), rep(4, 16), rep(5, 16))

flightlines_seg_norm <- lapply(seq_along(flightlines_sg), function(i) {
  brightness_normalize_segmented(flightlines_sg[[i]], groups = group_ids)
})

#diagnostic plot
#extract the new segment-normalized pixel value
img_seg_norm <- flightlines_seg_norm[[1]]
seg_norm_pixel <- as.vector(terra::values(img_seg_norm, mat = TRUE)[mid_cell, ])

#add to your existing plot data frame
plot_df$Seg_Normalized <- seg_norm_pixel

#generate the 4-way plot
ggplot(plot_df, aes(x = Wavelength, group = Group)) +
  geom_point(aes(y = Raw, color = "1. Raw Data"), size = 1.8, alpha = 0.4) +
  geom_line(aes(y = Smoothed, color = "2. Vectorized SG"), size = 1.0) +
  geom_line(aes(y = Normalized, color = "3. Global Brightness Norm"), size = 1.0, linetype = "dashed") +
  geom_line(aes(y = Seg_Normalized, color = "4. Segment-Wise Norm"), size = 1.0, linetype = "dotdash") +
  scale_color_manual(values = c(
    "1. Raw Data" = "gray40", 
    "2. Vectorized SG" = "#009E73", 
    "3. Global Brightness Norm" = "#D55E00",
    "4. Segment-Wise Norm" = "#0072B2"
  )) +
  scale_y_continuous(name = "Value Scale") +
  theme_minimal() +
  labs(title = "Global vs. Segment-Wise Normalization", x = "Wavelength (nm)", color = "Processing Stage") +
  theme(legend.position = "bottom")
#not using this approach: nearer to raw scale, but  it does not preserve relative shape
#the goal brightness normalization is not to keep the values near the raw scale
#it is to eliminate shading and cross-track illumination variations across flightlines 
#while protecting the genuine spectral shape. The global approach does exactly this
