install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Cimba1 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Cimba1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Cimba1.csv", header=TRUE, stringsAsFactors=FALSE)
Cimba1=Cimba1[c(1:5956),] # Oligocene-Miocene Unconformity
head(Cimba1)
plot(Cimba1, type="l", xlim = c(150, 1400), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Cimba1_interpolated <- linterp(Cimba1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Cmean = Gmean(Cimba1_interpolated$GR)
Cstd = Gsd(Cimba1_interpolated$GR)
Cimba1_scaled = (Cimba1_interpolated$GR - Cmean)/Cstd
Cimba1_rescaled = data.frame(Cimba1_interpolated$DEPT, Cimba1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Cimba1_scaled = mwStats(Cimba1_rescaled, cols = 2, win=3, ends = T)
Cimba1_standardized = data.frame(Cimba1_scaled$Center_win, Cimba1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Cimba1_standardized, type="l", xlim = c(150, 1400), ylim = c(-20, 20), xlab = "Cimba1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_c1_p1_ap1 <- dtw(Cimba1_standardized$Cimba1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_c1_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Cimba1_on_Picard1_depth = tune(Cimba1_standardized, cbind(Cimba1_standardized$Cimba1_scaled.Center_win[al_c1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_c1_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(Cimba1_on_Picard1_depth, col = "red")

# DTW Distance
al_c1_p1_ap1$normalizedDistance
al_c1_p1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Cimba1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Cimba1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 260))
base_1_y <- which.min(abs(Cimba1_standardized[,1] - 270))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 370))
base_2_y <- which.min(abs(Cimba1_standardized[,1] - 370))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 430))
base_3_y <- which.min(abs(Cimba1_standardized[,1] - 475))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 545))
base_4_y <- which.min(abs(Cimba1_standardized[,1] - 560))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 1010))
base_5_y <- which.min(abs(Cimba1_standardized[,1] - 980))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Cimba1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Cimba1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+200):nrow(Cimba1_standardized),1:(base_2_x-200)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+200):nrow(Cimba1_standardized),1:(base_3_x-200)] <- 0
compare.window[1:(base_3_y-200),(base_3_x+200):ncol(compare.window)] <- 0

compare.window[(base_4_y+200):nrow(Cimba1_standardized),1:(base_4_x-200)] <- 0
compare.window[1:(base_4_y-200),(base_4_x+200):ncol(compare.window)] <- 0

compare.window[(base_5_y+400):nrow(Cimba1_standardized),1:(base_5_x-400)] <- 0
compare.window[1:(base_5_y-400),(base_5_x+400):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Cimba1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Cimba1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_c1_p1_ap1 <- dtw(Cimba1_standardized$Cimba1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_c1_p1_ap1, type = "threeway")

# DTW Distance measure
al_c1_p1_ap1$normalizedDistance
al_c1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Cimba1_standardized[,1], z = compare.window, useRaster = T)
lines(Cimba1_standardized$Cimba1_scaled.Center_win[al_c1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_c1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Cimba1_on_Picard1_depth = tune(Cimba1_standardized, cbind(Cimba1_standardized$Cimba1_scaled.Center_win[al_c1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_c1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Cimba-1)")
lines(Cimba1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Cimba1_originalGR_on_Picard1_depth = data.frame(Cimba1_on_Picard1_depth$X1, Cimba1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Cimba-1)")
lines(Cimba1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Cimba1 

U1463Age_on_Cimba1_depth = tune(Cimba1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Cimba1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Cimba1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Cimba1_depth) <- new_column_names
write.csv(U1463Age_on_Cimba1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Cimba 1.csv", row.names = FALSE)
