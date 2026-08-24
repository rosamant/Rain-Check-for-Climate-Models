install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Venture1 and MaitlandNorth1 datasets

Venture1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Venture1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
V1 = Gmean(Venture1[c(1:1766),2])
V2 = Gmean(Venture1[c(800:1500),2])
SD1 = Gsd(Venture1[c(1:1766),2])
SD2 = Gsd(Venture1[c(800:1500),2])
Venture1[c(1:1766),2]=(Venture1[c(1:1766),2]+(V2-V1))*(SD1/SD2)

Venture1=Venture1[c(1:7223),] # Oligocene-Miocene
head(Venture1)
plot(Venture1, type="l", xlim = c(100, 1550), ylim = c(0, 50))

MaitlandNorth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/MaitlandNorth1.csv", header=TRUE, stringsAsFactors=FALSE)
MaitlandNorth1=MaitlandNorth1[c(1:5855),] # Oligocene-Miocene
head(MaitlandNorth1)
plot(MaitlandNorth1, type="l", xlim = c(100, 1000), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Venture1_interpolated <- linterp(Venture1, dt = 0.2, genplot = F)
MaitlandNorth1_interpolated <- linterp(MaitlandNorth1, dt = 0.2, genplot = F)

# Scaling the data
Vmean = Gmean(Venture1_interpolated$GR)
Vstd = Gsd(Venture1_interpolated$GR)
Venture1_scaled = (Venture1_interpolated$GR - Vmean)/Vstd
Venture1_rescaled = data.frame(Venture1_interpolated$DEPT, Venture1_scaled)

MNmean = Gmean(MaitlandNorth1_interpolated$GR)
MNstd = Gsd(MaitlandNorth1_interpolated$GR)
MaitlandNorth1_scaled = (MaitlandNorth1_interpolated$GR - MNmean)/MNstd
MaitlandNorth1_rescaled = data.frame(MaitlandNorth1_interpolated$DEPT, MaitlandNorth1_scaled)

# Resampling the data using moving window statistics
Venture1_scaled = mwStats(Venture1_rescaled, cols = 2, win=3, ends = T)
Venture1_standardized = data.frame(Venture1_scaled$Center_win, Venture1_scaled$Average)

MaitlandNorth1_scaled = mwStats(MaitlandNorth1_rescaled, cols = 2, win=3, ends = T)
MaitlandNorth1_standardized = data.frame(MaitlandNorth1_scaled$Center_win, MaitlandNorth1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Venture1_standardized, type="l", xlim = c(100, 1550), ylim = c(-20, 20), xlab = "Venture1 Resampled Depth", ylab = "Normalized GR")
plot(MaitlandNorth1_standardized, type="l", xlim = c(100, 1000), ylim = c(-20, 20), xlab = "MaitlandNorth1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_mn1_v1_ap1 <- dtw(MaitlandNorth1_standardized$MaitlandNorth1_scaled.Average, Venture1_standardized$Venture1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_mn1_v1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
MaitlandNorth1_on_Venture1_depth = tune(MaitlandNorth1_standardized, cbind(MaitlandNorth1_standardized$MaitlandNorth1_scaled.Center_win[al_mn1_v1_ap1$index1s], Venture1_standardized$Venture1_scaled.Center_win[al_mn1_v1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Venture1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1550), xlab = "Venture1 Resampled Depth", ylab = "Normalized GR")
lines(MaitlandNorth1_on_Venture1_depth, col = "red")

# DTW Distance 
al_mn1_v1_ap1$normalizedDistance
al_mn1_v1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(MaitlandNorth1_standardized),ncol=nrow(Venture1_standardized))
image(x=Venture1_standardized[,1],y=MaitlandNorth1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Venture1_standardized[,1] - 175))
base_1_y <- which.min(abs(MaitlandNorth1_standardized[,1] - 140))

# Depth values for second datum
base_2_x <- which.min(abs(Venture1_standardized[,1] - 250))
base_2_y <- which.min(abs(MaitlandNorth1_standardized[,1] - 200))

# Depth values for third datum
base_3_x <- which.min(abs(Venture1_standardized[,1] - 320))
base_3_y <- which.min(abs(MaitlandNorth1_standardized[,1] - 250))

# Depth values for fourth datum
base_4_x <- which.min(abs(Venture1_standardized[,1] - 570))
base_4_y <- which.min(abs(MaitlandNorth1_standardized[,1] - 500))

# Depth values for fifth datum
base_5_x <- which.min(abs(Venture1_standardized[,1] - 1500))
base_5_y <- which.min(abs(MaitlandNorth1_standardized[,1] - 900))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(MaitlandNorth1_standardized), ncol = nrow(Venture1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(MaitlandNorth1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-100),(base_1_x+100):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(MaitlandNorth1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(MaitlandNorth1_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(MaitlandNorth1_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-600),(base_4_x+600):ncol(compare.window)] <- 0

compare.window[(base_5_y+400):nrow(MaitlandNorth1_standardized),1:(base_5_x-400)] <- 0
compare.window[1:(base_5_y-200),(base_5_x+200):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Venture1_standardized[,1],y=MaitlandNorth1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Venture1_standardized[,1],y=MaitlandNorth1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_mn1_v1_ap1 <- dtw(MaitlandNorth1_standardized$MaitlandNorth1_scaled.Average, Venture1_standardized$Venture1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_mn1_v1_ap1, type = "threeway")

# DTW Distance measure
al_mn1_v1_ap1$normalizedDistance
al_mn1_v1_ap1$distance

image(y = Venture1_standardized[,1], x = MaitlandNorth1_standardized[,1], z = compare.window, useRaster = T)
lines(MaitlandNorth1_standardized$MaitlandNorth1_scaled.Center_win[al_mn1_v1_ap1$index1], Venture1_standardized$Venture1_scaled.Center_win[al_mn1_v1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
MaitlandNorth1_on_Venture1_depth = tune(MaitlandNorth1_standardized, cbind(MaitlandNorth1_standardized$MaitlandNorth1_scaled.Center_win[al_mn1_v1_ap1$index1s], Venture1_standardized$Venture1_scaled.Center_win[al_mn1_v1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Venture1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1550), xlab = "Venture1 Resampled Depth", ylab = "Normalized GR (Maitland North-1)")
lines(MaitlandNorth1_on_Venture1_depth, col = "red")

# Tuning MaitlandNorth1 data on Fisher1 depth
MaitlandNorth1_on_Fisher1_depth = tune(MaitlandNorth1_on_Venture1_depth, cbind(Venture1_standardized$Venture1_scaled.Center_win[al_v1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_v1_fi1_ap1$index2s]), extrapolate = F)

# Tuning MaitlandNorth1 data on Finucane1 depth
MaitlandNorth1_on_Finucane1_depth = tune(MaitlandNorth1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning MaitlandNorth1 data on Picard1 depth
MaitlandNorth1_on_Picard1_depth = tune(MaitlandNorth1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Maitland North-1)")
lines(MaitlandNorth1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
MaitlandNorth1_originalGR_on_Picard1_depth = data.frame(MaitlandNorth1_on_Picard1_depth$X1, MaitlandNorth1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Spar-1)")
lines(MaitlandNorth1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to MaitlandNorth1 

U1463Age_on_MaitlandNorth1_depth = tune(MaitlandNorth1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_MaitlandNorth1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "MaitlandNorth1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_MaitlandNorth1_depth) <- new_column_names
write.csv(U1463Age_on_MaitlandNorth1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Maitland North 1.csv", row.names = FALSE)
