install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Angel2 and Aurora1 datasets

Angel2 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Angel2.csv", header=TRUE, stringsAsFactors=FALSE)
Angel2=Angel2[c(1:8419),] # Oligocene-Miocene
head(Angel2)
plot(Angel2, type="l", xlim = c(100, 1600), ylim = c(0, 50))

Aurora1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Aurora1.csv", header=TRUE, stringsAsFactors=FALSE)
Aurora1=Aurora1[c(1:5614),] # Eocene-Miocene
head(Aurora1)
plot(Aurora1, type="l", xlim = c(150, 1000), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Angel2_interpolated <- linterp(Angel2, dt = 0.2, genplot = F)
Aurora1_interpolated <- linterp(Aurora1, dt = 0.2, genplot = F)

# Scaling the data
Amean = Gmean(Angel2_interpolated$GR)
Astd = Gsd(Angel2_interpolated$GR)
Angel2_scaled = (Angel2_interpolated$GR - Amean)/Astd
Angel2_rescaled = data.frame(Angel2_interpolated$DEPT, Angel2_scaled)

Aumean = Gmean(Aurora1_interpolated$GR)
Austd = Gsd(Aurora1_interpolated$GR)
Aurora1_scaled = (Aurora1_interpolated$GR - Aumean)/Austd
Aurora1_rescaled = data.frame(Aurora1_interpolated$DEPT, Aurora1_scaled)

# Resampling the data using moving window statistics
Angel2_scaled = mwStats(Angel2_rescaled, cols = 2, win=3, ends = T)
Angel2_standardized = data.frame(Angel2_scaled$Center_win, Angel2_scaled$Average)

Aurora1_scaled = mwStats(Aurora1_rescaled, cols = 2, win=3, ends = T)
Aurora1_standardized = data.frame(Aurora1_scaled$Center_win, Aurora1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Angel2_standardized, type="l", xlim = c(100, 1400), ylim = c(-20, 20), xlab = "Angel2 Resampled Depth", ylab = "Normalized GR")
plot(Aurora1_standardized, type="l", xlim = c(150, 1000), ylim = c(-20, 20), xlab = "Aurora1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_au1_a1_ap1 <- dtw(Aurora1_standardized$Aurora1_scaled.Average, Angel2_standardized$Angel2_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_au1_a1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Aurora1_on_Angel2_depth = tune(Aurora1_standardized, cbind(Aurora1_standardized$Aurora1_scaled.Center_win[al_au1_a1_ap1$index1s], Angel2_standardized$Angel2_scaled.Center_win[al_au1_a1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Angel2_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1400), xlab = "Angel2 Resampled Depth", ylab = "Normalized GR")
lines(Aurora1_on_Angel2_depth, col = "red")

# DTW Distance 
al_au1_a1_ap1$normalizedDistance
al_au1_a1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Aurora1_standardized),ncol=nrow(Angel2_standardized))
image(x=Angel2_standardized[,1],y=Aurora1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Angel2_standardized[,1] - 225))
base_1_y <- which.min(abs(Aurora1_standardized[,1] - 270))

# Depth values for second datum
base_2_x <- which.min(abs(Angel2_standardized[,1] - 350))
base_2_y <- which.min(abs(Aurora1_standardized[,1] - 385))

# Depth values for third datum
base_3_x <- which.min(abs(Angel2_standardized[,1] - 430))
base_3_y <- which.min(abs(Aurora1_standardized[,1] - 450))

# Depth values for fourth datum
base_4_x <- which.min(abs(Angel2_standardized[,1] - 630))
base_4_y <- which.min(abs(Aurora1_standardized[,1] - 620))

# Depth values for fifth datum
base_5_x <- which.min(abs(Angel2_standardized[,1] - 1275))
base_5_y <- which.min(abs(Aurora1_standardized[,1] - 875))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Aurora1_standardized), ncol = nrow(Angel2_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(Aurora1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-300),(base_1_x+300):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Aurora1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(Aurora1_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+300):nrow(Aurora1_standardized),1:(base_4_x-300)] <- 0
compare.window[1:(base_4_y-300),(base_4_x+300):ncol(compare.window)] <- 0

compare.window[(base_5_y+350):nrow(Aurora1_standardized),1:(base_5_x-350)] <- 0
compare.window[1:(base_5_y-200),(base_5_x+200):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Angel2_standardized[,1],y=Aurora1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Angel2_standardized[,1],y=Aurora1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_au1_a1_ap1 <- dtw(Aurora1_standardized$Aurora1_scaled.Average, Angel2_standardized$Angel2_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_au1_a1_ap1, type = "threeway")

# DTW Distance measure
al_au1_a1_ap1$normalizedDistance
al_au1_a1_ap1$distance

image(y = Angel2_standardized[,1], x = Aurora1_standardized[,1], z = compare.window, useRaster = T)
lines(Aurora1_standardized$Aurora1_scaled.Center_win[al_au1_a1_ap1$index1], Angel2_standardized$Angel2_scaled.Center_win[al_au1_a1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Aurora1_on_Angel2_depth = tune(Aurora1_standardized, cbind(Aurora1_standardized$Aurora1_scaled.Center_win[al_au1_a1_ap1$index1s], Angel2_standardized$Angel2_scaled.Center_win[al_au1_a1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Angel2_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1400), xlab = "Angel2 Resampled Depth", ylab = "Normalized GR (Aurora-1)")
lines(Aurora1_on_Angel2_depth, col = "red")

# Tuning Aurora1 data on Picard1 depth
Aurora1_on_Picard1_depth = tune(Aurora1_on_Angel2_depth, cbind(Angel2_standardized$Angel2_scaled.Center_win[al_a1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_a1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Aurora1 Resampled Depth", ylab = "Normalized GR (Aurora-1)")
lines(Aurora1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Aurora1_originalGR_on_Picard1_depth = data.frame(Aurora1_on_Picard1_depth$X1, Aurora1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Aurora-1)")
lines(Aurora1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Aurora1 

U1463Age_on_Aurora1_depth = tune(Aurora1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Aurora1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Aurora1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Aurora1_depth) <- new_column_names
write.csv(U1463Age_on_Aurora1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Aurora 1.csv", row.names = FALSE)
