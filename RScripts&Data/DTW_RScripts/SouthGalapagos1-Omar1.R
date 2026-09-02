install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import South Galapagos 1 and Omar 1 datasets

SouthGalapagos1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1.csv", header=TRUE, stringsAsFactors=FALSE)
SouthGalapagos1=SouthGalapagos1[c(30:3710),] # Oligocene-Miocene
head(SouthGalapagos1)
plot(SouthGalapagos1, type="l", xlim = c(500, 1200), ylim = c(0, 70))

Omar1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Omar1.csv", header=TRUE, stringsAsFactors=FALSE)
Omar1=Omar1[c(260:6481),] # Oligocene-Miocene
head(Omar1)
plot(Omar1, type="l", xlim = c(450, 1800), ylim = c(0, 100))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
SouthGalapagos1_interpolated <- linterp(SouthGalapagos1, dt = 0.2, genplot = F)
Omar1_interpolated <- linterp(Omar1, dt = 0.2, genplot = F)

# Scaling the data
Smean = Gmean(SouthGalapagos1_interpolated$GR)
Sstd = Gsd(SouthGalapagos1_interpolated$GR)
SouthGalapagos1_scaled = (SouthGalapagos1_interpolated$GR - Smean)/Sstd
SouthGalapagos1_rescaled = data.frame(SouthGalapagos1_interpolated$DEPT, SouthGalapagos1_scaled)

Omean = Gmean(Omar1_interpolated$GR)
Ostd = Gsd(Omar1_interpolated$GR)
Omar1_scaled = (Omar1_interpolated$GR - Omean)/Ostd
Omar1_rescaled = data.frame(Omar1_interpolated$DEPT, Omar1_scaled)

# Resampling the data using moving window statistics
SouthGalapagos1_scaled = mwStats(SouthGalapagos1_rescaled, cols = 2, win=3, ends = T)
SouthGalapagos1_standardized = data.frame(SouthGalapagos1_scaled$Center_win, SouthGalapagos1_scaled$Average)

Omar1_scaled = mwStats(Omar1_rescaled, cols = 2, win=3, ends = T)
Omar1_standardized = data.frame(Omar1_scaled$Center_win, Omar1_scaled$Average)

# Plotting the rescaled and resampled data
plot(SouthGalapagos1_standardized, type="l", xlim = c(500, 1200), ylim = c(-20, 20), xlab = "South Galapagos1 Resampled Depth", ylab = "Normalized GR")
plot(Omar1_standardized, type="l", xlim = c(450, 1700), ylim = c(-20, 30), xlab = "Omar1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_o1_sg1_ap1 <- dtw(Omar1_standardized$Omar1_scaled.Average, SouthGalapagos1_standardized$SouthGalapagos1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_o1_sg1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Omar1_on_SouthGalapagos1_depth = tune(Omar1_standardized, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap1$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "Omar1 Resampled Depth", ylab = "Normalized GR")
lines(Omar1_on_SouthGalapagos1_depth, col = "red")

# DTW Distance
al_o1_sg1_ap1$normalizedDistance
al_o1_sg1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Omar1_standardized),ncol=nrow(SouthGalapagos1_standardized))
image(x=SouthGalapagos1_standardized[,1],y=Omar1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 679))
base_1_y <- which.min(abs(Omar1_standardized[,1] - 773))

# Depth values for second datum
base_2_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 731))
base_2_y <- which.min(abs(Omar1_standardized[,1] - 841))

# Depth values for third datum
base_3_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 824))
base_3_y <- which.min(abs(Omar1_standardized[,1] - 1197))

# Depth values for fourth datum
base_4_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 912))
base_4_y <- which.min(abs(Omar1_standardized[,1] - 1239))

# Depth values for fifth datum
base_5_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 1021))
base_5_y <- which.min(abs(Omar1_standardized[,1] - 1442))

# Depth values for sixth datum
base_6_x <- which.min(abs(SouthGalapagos1_standardized[,1] - 1108))
base_6_y <- which.min(abs(Omar1_standardized[,1] - 1642))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Omar1_standardized), ncol = nrow(SouthGalapagos1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+50):nrow(Omar1_standardized),1:(base_1_x-50)] <- 0
compare.window[1:(base_1_y-50),(base_1_x+50):ncol(compare.window)] <- 0

compare.window[(base_2_y+100):nrow(Omar1_standardized),1:(base_2_x-100)] <- 0
compare.window[1:(base_2_y-100),(base_2_x+100):ncol(compare.window)] <- 0

compare.window[(base_3_y+200):nrow(Omar1_standardized),1:(base_3_x-200)] <- 0
compare.window[1:(base_3_y-200),(base_3_x+200):ncol(compare.window)] <- 0

compare.window[(base_4_y+40):nrow(Omar1_standardized),1:(base_4_x-40)] <- 0
compare.window[1:(base_4_y-40),(base_4_x+40):ncol(compare.window)] <- 0

compare.window[(base_5_y+30):nrow(Omar1_standardized),1:(base_5_x-30)] <- 0
compare.window[1:(base_5_y-30),(base_5_x+30):ncol(compare.window)] <- 0

compare.window[(base_6_y+50):nrow(Omar1_standardized),1:(base_6_x-50)] <- 0
compare.window[1:(base_6_y-50),(base_6_x+50):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=SouthGalapagos1_standardized[,1],y=Omar1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=SouthGalapagos1_standardized[,1],y=Omar1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_o1_sg1_ap2 <- dtw(Omar1_standardized$Omar1_scaled.Average, SouthGalapagos1_standardized$SouthGalapagos1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_o1_sg1_ap2, type = "threeway")

# DTW Distance measure
al_o1_sg1_ap2$normalizedDistance
al_o1_sg1_ap2$distance

image(y = SouthGalapagos1_standardized[,1], x = Omar1_standardized[,1], z = compare.window, useRaster = T)
lines(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Omar1_on_SouthGalapagos1_depth1 = tune(Omar1_standardized, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "South Galapagos-1 Resampled Depth", ylab = "Normalized GR (Omar-1)")
lines(Omar1_on_SouthGalapagos1_depth1, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Omar1_originalGR_on_SouthGalapagos1_depth = data.frame(Omar1_on_SouthGalapagos1_depth1$X1, Omar1_interpolated$GR)

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "SouthGalapagos1 Depth", ylab = "GR (Omar1)")
lines(Omar1_originalGR_on_SouthGalapagos1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Omar1 

SG1Age_on_Omar1_depth = tune(Omar1_originalGR_on_SouthGalapagos1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Omar1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Omar1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Omar1_depth) <- new_column_names
SG1Age_on_Omar1_depth[,1] = SG1Age_on_Omar1_depth[,1] * 1000
write.csv(SG1Age_on_Omar1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Omar 1.csv", row.names = FALSE)
