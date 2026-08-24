install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Omar1 and Calliance2 datasets

Omar1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Omar1.csv", header=TRUE, stringsAsFactors=FALSE)
Omar1=Omar1[c(260:6481),] # Oligocene-Miocene
head(Omar1)
plot(Omar1, type="l", xlim = c(450, 1800), ylim = c(0, 100))

Calliance2 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Calliance 2.csv", header=TRUE, stringsAsFactors=FALSE)
Calliance2=Calliance2[c(1:14202),] # Oligocene-Miocene
head(Calliance2)
plot(Calliance2, type="l", xlim = c(600, 2000), ylim = c(0, 80))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Omar1_interpolated <- linterp(Omar1, dt = 0.2, genplot = F)
Calliance2_interpolated <- linterp(Calliance2, dt = 0.2, genplot = F)

# Scaling the data

Omean = Gmean(Omar1_interpolated$GR)
Ostd = Gsd(Omar1_interpolated$GR)
Omar1_scaled = (Omar1_interpolated$GR - Omean)/Ostd
Omar1_rescaled = data.frame(Omar1_interpolated$DEPT, Omar1_scaled)

C2mean = Gmean(Calliance2_interpolated$GR)
C2std = Gsd(Calliance2_interpolated$GR)
Calliance2_scaled = (Calliance2_interpolated$GR - C2mean)/C2std
Calliance2_rescaled = data.frame(Calliance2_interpolated$DEPT, Calliance2_scaled)

# Resampling the data using moving window statistics

Omar1_scaled = mwStats(Omar1_rescaled, cols = 2, win=3, ends = T)
Omar1_standardized = data.frame(Omar1_scaled$Center_win, Omar1_scaled$Average)

Calliance2_scaled = mwStats(Calliance2_rescaled, cols = 2, win=3, ends = T)
Calliance2_standardized = data.frame(Calliance2_scaled$Center_win, Calliance2_scaled$Average)

# Plotting the rescaled and resampled data
plot(Omar1_standardized, type="l", xlim = c(450, 1800), ylim = c(-20, 30), xlab = "Omar1 Resampled Depth", ylab = "Normalized GR")
plot(Calliance2_standardized, type="l", xlim = c(600, 2000), ylim = c(-20, 30), xlab = "Calliance2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_c1_o1_ap1 <- dtw(Calliance2_standardized$Calliance2_scaled.Average, Omar1_standardized$Omar1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_c1_o1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Calliance2_on_Omar1_depth = tune(Calliance2_standardized, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c1_o1_ap1$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c1_o1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Omar1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1800), xlab = "Calliance2 Resampled Depth", ylab = "Normalized GR")
lines(Calliance2_on_Omar1_depth, col = "red")

# DTW Distance
al_c1_o1_ap1$normalizedDistance
al_c1_o1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Calliance2_standardized),ncol=nrow(Omar1_standardized))
image(x=Omar1_standardized[,1],y=Calliance2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Omar1_standardized[,1] - 750))
base_1_y <- which.min(abs(Calliance2_standardized[,1] - 760))

# Depth values for second datum
base_2_x <- which.min(abs(Omar1_standardized[,1] - 841))
base_2_y <- which.min(abs(Calliance2_standardized[,1] - 884))

# Depth values for third datum
base_3_x <- which.min(abs(Omar1_standardized[,1] - 1197))
base_3_y <- which.min(abs(Calliance2_standardized[,1] - 1305))

# Depth values for fourth datum
base_4_x <- which.min(abs(Omar1_standardized[,1] - 1239))
base_4_y <- which.min(abs(Calliance2_standardized[,1] - 1431))

# Depth values for fifth datum
base_5_x <- which.min(abs(Omar1_standardized[,1] - 1442))
base_5_y <- which.min(abs(Calliance2_standardized[,1] - 1592))

# Depth values for sixth datum
base_6_x <- which.min(abs(Omar1_standardized[,1] - 1642))
base_6_y <- which.min(abs(Calliance2_standardized[,1] - 1812))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Calliance2_standardized), ncol = nrow(Omar1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+20):nrow(Calliance2_standardized),1:(base_1_x-20)] <- 0
compare.window[1:(base_1_y-20),(base_1_x+20):ncol(compare.window)] <- 0

compare.window[(base_2_y+20):nrow(Calliance2_standardized),1:(base_2_x-20)] <- 0
compare.window[1:(base_2_y-20),(base_2_x+20):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(Calliance2_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+100):nrow(Calliance2_standardized),1:(base_4_x-100)] <- 0
compare.window[1:(base_4_y-100),(base_4_x+100):ncol(compare.window)] <- 0

compare.window[(base_5_y+100):nrow(Calliance2_standardized),1:(base_5_x-100)] <- 0
compare.window[1:(base_5_y-100),(base_5_x+100):ncol(compare.window)] <- 0

compare.window[(base_6_y+100):nrow(Calliance2_standardized),1:(base_6_x-100)] <- 0
compare.window[1:(base_6_y-100),(base_6_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Omar1_standardized[,1],y=Calliance2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Omar1_standardized[,1],y=Calliance2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_c2_o1_ap2 <- dtw(Calliance2_standardized$Calliance2_scaled.Average, Omar1_standardized$Omar1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_c2_o1_ap2, type = "threeway")

# DTW Distance measure
al_c2_o1_ap2$normalizedDistance
al_c2_o1_ap2$distance

image(y = Omar1_standardized[,1], x = Calliance2_standardized[,1], z = compare.window, useRaster = T)
lines(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Calliance2_on_Omar1_depth1 = tune(Calliance2_standardized, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Omar1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1800), xlab = "Omar-1 Resampled Depth", ylab = "Normalized GR (Calliance-2)")
lines(Calliance2_on_Omar1_depth1, col = "red")

# Tuning Calliance2 data on SG1 depth
Calliance2_on_SG1_depth = tune(Calliance2_on_Omar1_depth1, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Calliance-2)")
lines(Calliance2_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Calliance2_originalGR_on_SG1_depth = data.frame(Calliance2_on_SG1_depth$X1, Calliance2_interpolated$GR)

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Calliance-2)")
lines(Calliance2_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Calliance2 

SG1Age_on_Calliance2_depth = tune(Calliance2_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Calliance2_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Calliance2")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Calliance2_depth) <- new_column_names
SG1Age_on_Calliance2_depth[,1] = SG1Age_on_Calliance2_depth[,1] * 1000
write.csv(SG1Age_on_Calliance2_depth, file = "RScripts&Data/Sites Data_Age-NGR/Calliance 2.csv", row.names = FALSE)
