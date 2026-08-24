install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Whitetail1 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Whitetail1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Whitetail1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
W1 = Gmean(Whitetail1[c(1:320),2])
W2 = Gmean(Whitetail1[c(300:500),2])
SD1 = Gsd(Whitetail1[c(1:320),2])
SD2 = Gsd(Whitetail1[c(300:500),2])
Whitetail1[c(1:320),2]=(Whitetail1[c(1:320),2]+(W2-W1))*(SD1/SD2)

Whitetail1=Whitetail1[c(1:5430),] # Oligocene-Miocene
head(Whitetail1)
plot(Whitetail1, type="l", xlim = c(1000, 2100), ylim = c(0, 100))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Whitetail1_interpolated <- linterp(Whitetail1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Wmean = Gmean(Whitetail1_interpolated$GR)
Wstd = Gsd(Whitetail1_interpolated$GR)
Whitetail1_scaled = (Whitetail1_interpolated$GR - Wmean)/Wstd
Whitetail1_rescaled = data.frame(Whitetail1_interpolated$DEPT, Whitetail1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Whitetail1_scaled = mwStats(Whitetail1_rescaled, cols = 2, win=3, ends = T)
Whitetail1_standardized = data.frame(Whitetail1_scaled$Center_win, Whitetail1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Whitetail1_standardized, type="l", xlim = c(1000, 1900), ylim = c(-20, 20), xlab = "Whitetail1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_w1_p1_ap1 <- dtw(Whitetail1_standardized$Whitetail1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_w1_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Whitetail1_on_Picard1_depth = tune(Whitetail1_standardized, cbind(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(Whitetail1_on_Picard1_depth, col = "red")

# DTW Distance
al_w1_p1_ap1$normalizedDistance
al_w1_p1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Whitetail1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 260))
base_1_y <- which.min(abs(Whitetail1_standardized[,1] - 1100))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 370))
base_2_y <- which.min(abs(Whitetail1_standardized[,1] - 1190))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 430))
base_3_y <- which.min(abs(Whitetail1_standardized[,1] - 1230))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 545))
base_4_y <- which.min(abs(Whitetail1_standardized[,1] - 1350))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 1010))
base_5_y <- which.min(abs(Whitetail1_standardized[,1] - 1750))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Whitetail1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Whitetail1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Whitetail1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(Whitetail1_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(Whitetail1_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-400),(base_4_x+400):ncol(compare.window)] <- 0

compare.window[(base_5_y+500):nrow(Whitetail1_standardized),1:(base_5_x-500)] <- 0
compare.window[1:(base_5_y-500),(base_5_x+500):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Whitetail1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_w1_p1_ap1 <- dtw(Whitetail1_standardized$Whitetail1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_w1_p1_ap1, type = "threeway")

# DTW Distance measure
al_w1_p1_ap1$normalizedDistance
al_w1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Whitetail1_standardized[,1], z = compare.window, useRaster = T)
lines(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Whitetail1_on_Picard1_depth = tune(Whitetail1_standardized, cbind(Whitetail1_standardized$Whitetail1_scaled.Center_win[al_w1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_w1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Whitetail1)")
lines(Whitetail1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Whitetail1_originalGR_on_Picard1_depth = data.frame(Whitetail1_on_Picard1_depth$X1, Whitetail1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Whitetail-1)")
lines(Whitetail1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Whitetail1 

U1463Age_on_Whitetail1_depth = tune(Whitetail1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Whitetail1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Whitetail-1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Whitetail1_depth) <- new_column_names
write.csv(U1463Age_on_Whitetail1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Whitetail 1.csv", row.names = FALSE)
