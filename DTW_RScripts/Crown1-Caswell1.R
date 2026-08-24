install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Caswell1 and Crown1 datasets

Caswell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Caswell1.csv", header=TRUE, stringsAsFactors=FALSE)
Caswell1=Caswell1[c(1:6434),] # Oligocene-Miocene
head(Caswell1)
plot(Caswell1, type="l", xlim = c(350, 1750), ylim = c(0, 50))

Crown1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Crown 1.csv", header=TRUE, stringsAsFactors=FALSE)
Crown1=Crown1[c(1:8673),] # Oligocene-Miocene
head(Crown1)
plot(Crown1, type="l", xlim = c(450, 2250), ylim = c(0, 90))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Caswell1_interpolated <- linterp(Caswell1, dt = 0.2, genplot = F)
Crown1_interpolated <- linterp(Crown1, dt = 0.2, genplot = F)

# Scaling the data
Cmean = Gmean(Caswell1_interpolated$GR)
Cstd = Gsd(Caswell1_interpolated$GR)
Caswell1_scaled = (Caswell1_interpolated$GR - Cmean)/Cstd
Caswell1_rescaled = data.frame(Caswell1_interpolated$DEPT, Caswell1_scaled)

Crmean = Gmean(Crown1_interpolated$GR)
Crstd = Gsd(Crown1_interpolated$GR)
Crown1_scaled = (Crown1_interpolated$GR - Crmean)/Crstd
Crown1_rescaled = data.frame(Crown1_interpolated$DEPT, Crown1_scaled)

# Resampling the data using moving window statistics
Caswell1_scaled = mwStats(Caswell1_rescaled, cols = 2, win=3, ends = T)
Caswell1_standardized = data.frame(Caswell1_scaled$Center_win, Caswell1_scaled$Average)

Crown1_scaled = mwStats(Crown1_rescaled, cols = 2, win=3, ends = T)
Crown1_standardized = data.frame(Crown1_scaled$Center_win, Crown1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Caswell1_standardized, type="l", xlim = c(350, 1700), ylim = c(-20, 20), xlab = "Caswell1 Resampled Depth", ylab = "Normalized GR")
plot(Crown1_standardized, type="l", xlim = c(450, 2250), ylim = c(-20, 40), xlab = "Crown1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_cr1_c1_ap1 <- dtw(Crown1_standardized$Crown1_scaled.Average, Caswell1_standardized$Caswell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_cr1_c1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Crown1_on_Caswell1_depth = tune(Crown1_standardized, cbind(Crown1_standardized$Crown1_scaled.Center_win[al_cr1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_cr1_c1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Caswell1_standardized, type = "l", ylim = c(-20, 40), xlim = c(350, 1700), xlab = "Caswell1 Resampled Depth", ylab = "Normalized GR")
lines(Crown1_on_Caswell1_depth, col = "red")

# DTW Distance

al_cr1_c1_ap1$normalizedDistance
al_cr1_c1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Crown1_standardized),ncol=nrow(Caswell1_standardized))
image(x=Caswell1_standardized[,1],y=Crown1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Omar1_standardized[,1] - 706))
base_1_y <- which.min(abs(Calliance2_standardized[,1] - 893))

# Depth values for second datum
base_2_x <- which.min(abs(Caswell1_standardized[,1] - 870))
base_2_y <- which.min(abs(Crown1_standardized[,1] - 981))

# Depth values for third datum
base_3_x <- which.min(abs(Caswell1_standardized[,1] - 1098))
base_3_y <- which.min(abs(Crown1_standardized[,1] - 1383))

# Depth values for fourth datum
base_4_x <- which.min(abs(Caswell1_standardized[,1] - 1228))
base_4_y <- which.min(abs(Crown1_standardized[,1] - 1618))

# Depth values for fifth datum
base_5_x <- which.min(abs(Caswell1_standardized[,1] - 1352))
base_5_y <- which.min(abs(Crown1_standardized[,1] - 1793))

# Depth values for sixth datum
base_6_x <- which.min(abs(Caswell1_standardized[,1] - 1511))
base_6_y <- which.min(abs(Crown1_standardized[,1] - 1991))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Crown1_standardized), ncol = nrow(Caswell1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(Crown1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-300),(base_1_x+300):ncol(compare.window)] <- 0

compare.window[(base_2_y+20):nrow(Crown1_standardized),1:(base_2_x-20)] <- 0
compare.window[1:(base_2_y-20),(base_2_x+20):ncol(compare.window)] <- 0

compare.window[(base_3_y+20):nrow(Crown1_standardized),1:(base_3_x-20)] <- 0
compare.window[1:(base_3_y-20),(base_3_x+20):ncol(compare.window)] <- 0

compare.window[(base_4_y+20):nrow(Crown1_standardized),1:(base_4_x-20)] <- 0
compare.window[1:(base_4_y-20),(base_4_x+20):ncol(compare.window)] <- 0

compare.window[(base_5_y+20):nrow(Crown1_standardized),1:(base_5_x-20)] <- 0
compare.window[1:(base_5_y-20),(base_5_x+20):ncol(compare.window)] <- 0

compare.window[(base_6_y+100):nrow(Crown1_standardized),1:(base_6_x-100)] <- 0
compare.window[1:(base_6_y-100),(base_6_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Caswell1_standardized[,1],y=Crown1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Caswell1_standardized[,1],y=Crown1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_cr1_c1_ap2 <- dtw(Crown1_standardized$Crown1_scaled.Average, Caswell1_standardized$Caswell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_cr1_c1_ap2, type = "threeway")

# DTW Distance measure
al_cr1_c1_ap2$normalizedDistance
al_cr1_c1_ap2$distance

image(y = Caswell1_standardized[,1], x = Crown1_standardized[,1], z = compare.window, useRaster = T)
lines(Crown1_standardized$Crown1_scaled.Center_win[al_cr1_c1_ap2$index1], Caswell1_standardized$Caswell1_scaled.Center_win[al_cr1_c1_ap2$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Crown1_on_Caswell1_depth1 = tune(Crown1_standardized, cbind(Crown1_standardized$Crown1_scaled.Center_win[al_cr1_c1_ap2$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_cr1_c1_ap2$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Caswell1_standardized, type = "l", ylim = c(-20, 40), xlim = c(350, 1700), xlab = "Caswell-1 Resampled Depth", ylab = "Normalized GR (Crown-1)")
lines(Crown1_on_Caswell1_depth1, col = "red")

# Tuning Crown1 data on Calliance2 depth
Crown1_on_Calliance2_depth = tune(Crown1_on_Caswell1_depth1, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Crown1 data on Omar1 depth
Crown1_on_Omar1_depth = tune(Crown1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Crown1 data on SG1 depth
Crown1_on_SG1_depth = tune(Crown1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Crown-1)")
lines(Crown1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Crown1_originalGR_on_SG1_depth = data.frame(Crown1_on_SG1_depth$X1, Crown1_interpolated[1312:8494,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Crown-1)")
lines(Crown1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Crown1 

SG1Age_on_Crown1_depth = tune(Crown1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Crown1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Crown1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Crown1_depth) <- new_column_names
SG1Age_on_Crown1_depth[,1] = SG1Age_on_Crown1_depth[,1] * 1000
write.csv(SG1Age_on_Crown1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Crown 1.csv", row.names = FALSE)

Crown1_age_depth <- approx(x = Crown1_interpolated$GR,
                             y = Crown1_interpolated$DEPT,
                             xout = SG1Age_on_Crown1_depth$GR,
                             rule = 1)$y

Crown1_agemodel = data.frame(Crown1_age_depth, SG1Age_on_Crown1_depth$AGE)
Crown1_agemodel <- na.omit(Crown1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Crown1_agemodel) <- new_column_names1
plot(Crown1_agemodel)
write.csv(Crown1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Crown1_DepthAge.csv", row.names = FALSE)
