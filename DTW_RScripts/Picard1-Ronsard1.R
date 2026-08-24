install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Ronsard1 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Ronsard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Ronsard1.csv", header=TRUE, stringsAsFactors=FALSE)
Ronsard1=Ronsard1[c(1:8728),]# Oligocene-Miocene
head(Ronsard1)
plot(Ronsard1, type="l", xlim = c(150, 1500), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Ronsard1_interpolated <- linterp(Ronsard1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Rmean = Gmean(Ronsard1_interpolated$GR)
Rstd = Gsd(Ronsard1_interpolated$GR)
Ronsard1_scaled = (Ronsard1_interpolated$GR - Rmean)/Rstd
Ronsard1_rescaled = data.frame(Ronsard1_interpolated$DEPT, Ronsard1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Ronsard1_scaled = mwStats(Ronsard1_rescaled, cols = 2, win=3, ends = T)
Ronsard1_standardized = data.frame(Ronsard1_scaled$Center_win, Ronsard1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Ronsard1_standardized, type="l", xlim = c(150, 1700), ylim = c(-20, 20), xlab = "Ronsard1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_r1_p1_ap1 <- dtw(Ronsard1_standardized$Ronsard1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_r1_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Ronsard1_on_Picard1_depth = tune(Ronsard1_standardized, cbind(Ronsard1_standardized$Ronsard1_scaled.Center_win[al_r1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_r1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(Ronsard1_on_Picard1_depth, col = "red")

# DTW Distance 
al_r1_p1_ap1$normalizedDistance
al_r1_p1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Ronsard1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Ronsard1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 240))
base_1_y <- which.min(abs(Ronsard1_standardized[,1] - 340))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 370))
base_2_y <- which.min(abs(Ronsard1_standardized[,1] - 420))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 430))
base_3_y <- which.min(abs(Ronsard1_standardized[,1] - 495))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 545))
base_4_y <- which.min(abs(Ronsard1_standardized[,1] - 660))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 1010))
base_5_y <- which.min(abs(Ronsard1_standardized[,1] - 1225))

# Depth values for sixth datum
base_6_x <- which.min(abs(Picard1_standardized[,1] - 1237))
base_6_y <- which.min(abs(Ronsard1_standardized[,1] - 1427))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Ronsard1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Ronsard1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+200):nrow(Ronsard1_standardized),1:(base_2_x-200)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+200):nrow(Ronsard1_standardized),1:(base_3_x-200)] <- 0
compare.window[1:(base_3_y-200),(base_3_x+200):ncol(compare.window)] <- 0

compare.window[(base_4_y+200):nrow(Ronsard1_standardized),1:(base_4_x-200)] <- 0
compare.window[1:(base_4_y-200),(base_4_x+200):ncol(compare.window)] <- 0

compare.window[(base_5_y+200):nrow(Ronsard1_standardized),1:(base_5_x-200)] <- 0
compare.window[1:(base_5_y-200),(base_5_x+200):ncol(compare.window)] <- 0

compare.window[(base_6_y+100):nrow(Ronsard1_standardized),1:(base_6_x-100)] <- 0
compare.window[1:(base_6_y-100),(base_6_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Ronsard1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Ronsard1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_r1_p1_ap1 <- dtw(Ronsard1_standardized$Ronsard1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_r1_p1_ap1, type = "threeway")

# DTW Distance measure
al_r1_p1_ap1$normalizedDistance
al_r1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Ronsard1_standardized[,1], z = compare.window, useRaster = T)
lines(Ronsard1_standardized$Ronsard1_scaled.Center_win[al_r1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_r1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Ronsard1_on_Picard1_depth = tune(Ronsard1_standardized, cbind(Ronsard1_standardized$Ronsard1_scaled.Center_win[al_r1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_r1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Ronsard-1)")
lines(Ronsard1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Ronsard1_originalGR_on_Picard1_depth = data.frame(Ronsard1_on_Picard1_depth$X1, Ronsard1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Ronsard-1)")
lines(Ronsard1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Ronsard1 

U1463Age_on_Ronsard1_depth = tune(Ronsard1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Ronsard1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Ronsard1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Ronsard1_depth) <- new_column_names
write.csv(U1463Age_on_Ronsard1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Ronsard 1.csv", row.names = FALSE)
