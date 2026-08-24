install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Gorgon1 and Griffin7 datasets

Gorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgon1=Gorgon1[c(1:9280),] # Data required till Eocene-Miocene Unconformity
head(Gorgon1)
plot(Gorgon1, type="l", xlim = c(300, 1700), ylim = c(0, 50))

Griffin7 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Griffin7.csv", header=TRUE, stringsAsFactors=FALSE)
Griffin7 = Griffin7[c(1:6857),] # Data required till Eocene-Miocene Unconformity
head(Griffin7)
plot(Griffin7, type="l", xlim = c(200, 1550), ylim = c(0, 80))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgon1_interpolated <- linterp(Gorgon1, dt = 0.2, genplot = F)
Griffin7_interpolated <- linterp(Griffin7, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgon1_interpolated$GR)
Gstd = Gsd(Gorgon1_interpolated$GR)
Gorgon1_scaled = (Gorgon1_interpolated$GR - Gmean)/Gstd
Gorgon1_rescaled = data.frame(Gorgon1_interpolated$DEPT, Gorgon1_scaled)

Grmean = Gmean(Griffin7_interpolated$GR)
Grstd = Gsd(Griffin7_interpolated$GR)
Griffin7_scaled = (Griffin7_interpolated$GR - Grmean)/Grstd
Griffin7_rescaled = data.frame(Griffin7_interpolated$DEPT, Griffin7_scaled)

# Resampling the data using moving window statistics
Gorgon1_scaled = mwStats(Gorgon1_rescaled, cols = 2, win = 3, ends = T)
Gorgon1_standardized = data.frame(Gorgon1_scaled$Center_win, Gorgon1_scaled$Average)

Griffin7_scaled = mwStats(Griffin7_rescaled, cols = 2, win = 3, ends = T)
Griffin7_standardized = data.frame(Griffin7_scaled$Center_win, Griffin7_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgon1_standardized, type="l", xlim = c(300, 1700), ylim = c(-20, 20), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
plot(Griffin7_standardized, type="l", xlim = c(200, 1550), ylim = c(-20, 20), xlab = "Griffin7 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_gr1_g1_ap1 <- dtw(Griffin7_standardized$Griffin7_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_gr1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Griffin7_on_Gorgon1_depth = tune(Griffin7_standardized, cbind(Griffin7_standardized$Griffin7_scaled.Center_win[al_gr1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_gr1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
lines(Griffin7_on_Gorgon1_depth, col = "red")

# DTW Distance measure
al_gr1_g1_ap1$normalizedDistance
al_gr1_g1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Griffin7_standardized),ncol=nrow(Gorgon1_standardized))
image(x=Gorgon1_standardized[,1],y=Griffin7_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Gorgon1_standardized[,1] - 370))
base_1_y <- which.min(abs(Griffin7_standardized[,1] - 320))

# Depth values for second datum
base_2_x <- which.min(abs(Gorgon1_standardized[,1] - 430))
base_2_y <- which.min(abs(Griffin7_standardized[,1] - 365))

# Depth values for third datum
base_3_x <- which.min(abs(Gorgon1_standardized[,1] - 545))
base_3_y <- which.min(abs(Griffin7_standardized[,1] - 480))

# Depth values for fourth datum
base_4_x <- which.min(abs(Gorgon1_standardized[,1] - 1010))
base_4_y <- which.min(abs(Griffin7_standardized[,1] - 810))

# Depth values for fifth datum
base_5_x <- which.min(abs(Gorgon1_standardized[,1] - 1190))
base_5_y <- which.min(abs(York1_standardized[,1] - 940))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Griffin7_standardized), ncol = nrow(Gorgon1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(Griffin7_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-300),(base_1_x+300):ncol(compare.window)] <- 0

compare.window[(base_2_y+400):nrow(Griffin7_standardized),1:(base_2_x-400)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(Griffin7_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(Griffin7_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-300),(base_4_x+300):ncol(compare.window)] <- 0

compare.window[(base_5_y+500):nrow(Griffin7_standardized),1:(base_5_x-500)] <- 0
compare.window[1:(base_5_y-400),(base_5_x+400):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Gorgon1_standardized[,1],y=Griffin7_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Gorgon1_standardized[,1],y=Griffin7_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_gr1_g1_ap1 <- dtw(Griffin7_standardized$Griffin7_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_gr1_g1_ap1, type = "threeway")

# DTW Distance measure
al_gr1_g1_ap1$normalizedDistance
al_gr1_g1_ap1$distance

image(y = Gorgon1_standardized[,1], x = Griffin7_standardized[,1], z = compare.window, useRaster = T)
lines(Griffin7_standardized$Griffin7_scaled.Center_win[al_gr1_g1_ap1$index1], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_gr1_g1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Griffin7_on_Gorgon1_depth = tune(Griffin7_standardized, cbind(Griffin7_standardized$Griffin7_scaled.Center_win[al_gr1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_gr1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR (Griffin-7)")
lines(Griffin7_on_Gorgon1_depth, col = "red")

# Tuning Griffin7 data on Bluebell1 depth
Griffin7_on_Bluebell1_depth = tune(Griffin7_on_Gorgon1_depth, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

# Tuning Griffin7 data on Fisher1 depth
Griffin7_on_Fisher1_depth = tune(Griffin7_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Griffin7 data on Finucane1 depth
Griffin7_on_Finucane1_depth = tune(Griffin7_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Griffin7 data on Picard1 depth
Griffin7_on_Picard1_depth = tune(Griffin7_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Griffin-7)")
lines(Griffin7_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Griffin7_originalGR_on_Picard1_depth = data.frame(Griffin7_on_Picard1_depth$X1, Griffin7_interpolated[1:6820,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 70), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Griffin-7)")
lines(Griffin7_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Griffin7 

U1463Age_on_Griffin7_depth = tune(Griffin7_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Griffin7_depth, type = "l", ylim = c(0, 70), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Griffin7")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Griffin7_depth) <- new_column_names
write.csv(U1463Age_on_Griffin7_depth, file = "RScripts&Data/Sites Data_Age-NGR/Griffin 7.csv", row.names = FALSE)
