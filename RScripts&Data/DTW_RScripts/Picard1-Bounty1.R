install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Bounty1 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Bounty1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Bounty1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
B1 = Gmean(Bounty1[c(1:1961),2])
B2 = Gmean(Bounty1[c(1961:4000),2])
SD1 = Gsd(Bounty1[c(1:1961),2])
SD2 = Gsd(Bounty1[c(1961:4000),2])
Bounty1[c(1:1961),2]=(Bounty1[c(1:1961),2]+(B2-B1))*(SD1/SD2)

Bounty1=Bounty1[c(1:10146),] # Oligocene-Miocene
head(Bounty1)
plot(Bounty1, type="l", xlim = c(150, 1800), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Bounty1_interpolated <- linterp(Bounty1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Bmean = Gmean(Bounty1_interpolated$GR)
Bstd = Gsd(Bounty1_interpolated$GR)
Bounty1_scaled = (Bounty1_interpolated$GR - Bmean)/Bstd
Bounty1_rescaled = data.frame(Bounty1_interpolated$DEPT, Bounty1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Bounty1_scaled = mwStats(Bounty1_rescaled, cols = 2, win=3, ends = T)
Bounty1_standardized = data.frame(Bounty1_scaled$Center_win, Bounty1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Bounty1_standardized, type="l", xlim = c(150, 1700), ylim = c(-20, 20), xlab = "Bounty1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_b1_p1_ap1 <- dtw(Bounty1_standardized$Bounty1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_b1_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Bounty1_on_Picard1_depth = tune(Bounty1_standardized, cbind(Bounty1_standardized$Bounty1_scaled.Center_win[al_b1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_b1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(Bounty1_on_Picard1_depth, col = "red")

# DTW Distance 
al_b1_p1_ap1$normalizedDistance
al_b1_p1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Bounty1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Bounty1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 260))
base_1_y <- which.min(abs(Bounty1_standardized[,1] - 320))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 370))
base_2_y <- which.min(abs(Bounty1_standardized[,1] - 375))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 430))
base_3_y <- which.min(abs(Bounty1_standardized[,1] - 460))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 545))
base_4_y <- which.min(abs(Bounty1_standardized[,1] - 630))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 1010))
base_5_y <- which.min(abs(Bounty1_standardized[,1] - 1000))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Bounty1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Bounty1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Bounty1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(Bounty1_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+500):nrow(Bounty1_standardized),1:(base_4_x-500)] <- 0
compare.window[1:(base_4_y-400),(base_4_x+400):ncol(compare.window)] <- 0

compare.window[(base_5_y+500):nrow(Bounty1_standardized),1:(base_5_x-500)] <- 0
compare.window[1:(base_5_y-500),(base_5_x+500):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Bounty1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Bounty1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_b1_p1_ap1 <- dtw(Bounty1_standardized$Bounty1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_b1_p1_ap1, type = "threeway")

# DTW Distance measure
al_b1_p1_ap1$normalizedDistance
al_b1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Bounty1_standardized[,1], z = compare.window, useRaster = T)
lines(Bounty1_standardized$Bounty1_scaled.Center_win[al_b1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_b1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Bounty1_on_Picard1_depth = tune(Bounty1_standardized, cbind(Bounty1_standardized$Bounty1_scaled.Center_win[al_b1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_b1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bounty-1)")
lines(Bounty1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Bounty1_originalGR_on_Picard1_depth = data.frame(Bounty1_on_Picard1_depth$X1, Bounty1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bounty-1)")
lines(Bounty1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Bounty1 

U1463Age_on_Bounty1_depth = tune(Bounty1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Bounty1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Bounty-1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Bounty1_depth) <- new_column_names
write.csv(U1463Age_on_Bounty1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Bounty 1.csv", row.names = FALSE)
