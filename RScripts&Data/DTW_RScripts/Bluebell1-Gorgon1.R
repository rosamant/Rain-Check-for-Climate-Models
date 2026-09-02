install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Bluebell1 and Gorgon1 datasets

Bluebell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Bluebell1.csv", header=TRUE, stringsAsFactors=FALSE)
Bluebell1=Bluebell1[c(1:11479),] # Oligocene-Miocene
head(Bluebell1)
plot(Bluebell1, type="l", xlim = c(200, 2000), ylim = c(0, 60))

Gorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgon1=Gorgon1[c(1:9280),] # Oligocene-Miocene
head(Gorgon1)
plot(Gorgon1, type="l", xlim = c(300, 1700), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Bluebell1_interpolated <- linterp(Bluebell1, dt = 0.2, genplot = F)
Gorgon1_interpolated <- linterp(Gorgon1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Bluebell1_interpolated$GR)
Bstd = Gsd(Bluebell1_interpolated$GR)
Bluebell1_scaled = (Bluebell1_interpolated$GR - Bmean)/Bstd
Bluebell1_rescaled = data.frame(Bluebell1_interpolated$DEPT, Bluebell1_scaled)

Gmean = Gmean(Gorgon1_interpolated$GR)
Gstd = Gsd(Gorgon1_interpolated$GR)
Gorgon1_scaled = (Gorgon1_interpolated$GR - Gmean)/Gstd
Gorgon1_rescaled = data.frame(Gorgon1_interpolated$DEPT, Gorgon1_scaled)

# Resampling the data using moving window statistics
Bluebell1_scaled = mwStats(Bluebell1_rescaled, cols = 2, win=3, ends = T)
Bluebell1_standardized = data.frame(Bluebell1_scaled$Center_win, Bluebell1_scaled$Average)

Gorgon1_scaled = mwStats(Gorgon1_rescaled, cols = 2, win=3)
Gorgon1_standardized = data.frame(Gorgon1_scaled$Center_win, Gorgon1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Bluebell1_standardized, type="l", xlim = c(200, 2000), ylim = c(-20, 20), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR")
plot(Gorgon1_standardized, type="l", xlim = c(300, 1700), ylim = c(-20, 20), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")

# DTW with asymmetricP01 and open-close ends

system.time(al_g1_b1_ap1 <- dtw(Gorgon1_standardized$Gorgon1_scaled.Average, Bluebell1_standardized$Bluebell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_g1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Gorgon1_on_Bluebell1_depth = tune(Gorgon1_standardized, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Bluebell1_standardized, type = "l", ylim = c(-20, 20), xlim = c(200, 2000), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR")
lines(Gorgon1_on_Bluebell1_depth, col = "red")

# DTW Distance and RMSE

al_g1_b1_ap1$normalizedDistance
al_g1_b1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Gorgon1_standardized),ncol=nrow(Bluebell1_standardized))
image(x=Bluebell1_standardized[,1],y=Gorgon1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Bluebell1_standardized[,1] - 310))
base_1_y <- which.min(abs(Gorgon1_standardized[,1] - 380))

# Depth values for second datum
base_2_x <- which.min(abs(Bluebell1_standardized[,1] - 450))
base_2_y <- which.min(abs(Gorgon1_standardized[,1] - 525))

# Depth values for third datum
base_3_x <- which.min(abs(Bluebell1_standardized[,1] - 550))
base_3_y <- which.min(abs(Gorgon1_standardized[,1] - 625))

# Depth values for fourth datum
base_4_x <- which.min(abs(Bluebell1_standardized[,1] - 900))
base_4_y <- which.min(abs(Gorgon1_standardized[,1] - 870))

# Depth values for fifth datum
base_5_x <- which.min(abs(Bluebell1_standardized[,1] - 1625))
base_5_y <- which.min(abs(Gorgon1_standardized[,1] - 1550))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Gorgon1_standardized), ncol = nrow(Bluebell1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(Gorgon1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Gorgon1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-500),(base_2_x+500):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(Gorgon1_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(Gorgon1_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-300),(base_4_x+300):ncol(compare.window)] <- 0

compare.window[(base_5_y+400):nrow(Gorgon1_standardized),1:(base_5_x-400)] <- 0
compare.window[1:(base_5_y-300),(base_5_x+300):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Bluebell1_standardized[,1],y=Gorgon1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Bluebell1_standardized[,1],y=Gorgon1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_g1_b1_ap1 <- dtw(Gorgon1_standardized$Gorgon1_scaled.Average, Bluebell1_standardized$Bluebell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))

# DTW Distance measure
al_g1_b1_ap1$normalizedDistance
al_g1_b1_ap1$distance

plot(al_g1_b1_ap1, type = "threeway")

image(y = Bluebell1_standardized[,1], x = Gorgon1_standardized[,1], z = compare.window, useRaster = T)
lines(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Gorgon1_on_Bluebell1_depth = tune(Gorgon1_standardized, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Bluebell1_standardized, type = "l", ylim = c(-20, 20), xlim = c(200, 2000), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR (Gorgon-1)")
lines(Gorgon1_on_Bluebell1_depth, col = "red")

# Tuning Gorgon1 data on Fisher1 depth
Gorgon1_on_Fisher1_depth = tune(Gorgon1_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Gorgon1 data on Finucane1 depth
Gorgon1_on_Finucane1_depth = tune(Gorgon1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Gorgon1 data on Picard1 depth
Gorgon1_on_Picard1_depth = tune(Gorgon1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Gorgon-1)")
lines(Gorgon1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Gorgon1_originalGR_on_Picard1_depth = data.frame(Gorgon1_on_Picard1_depth$X1, Gorgon1_interpolated[1:7054,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Gorgon-1)")
lines(Gorgon1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Gorgon1 

U1463Age_on_Gorgon1_depth = tune(Gorgon1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Gorgon1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Gorgon1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Gorgon1_depth) <- new_column_names
write.csv(U1463Age_on_Gorgon1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Gorgon 1.csv", row.names = FALSE)
