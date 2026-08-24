install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Gorgon1 and York1 datasets

Gorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgon1=Gorgon1[c(1:9280),] # Data required till Eocene-Miocene Unconformity
head(Gorgon1)
plot(Gorgon1, type="l", xlim = c(300, 1700), ylim = c(0, 50))

York1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/York1.csv", header=TRUE, stringsAsFactors=FALSE)
York1 = York1[c(1:5654),] # Data required till Eocene-Miocene Unconformity
head(York1)
plot(York1, type="l", xlim = c(400, 1500), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgon1_interpolated <- linterp(Gorgon1, dt = 0.2, genplot = F)
York1_interpolated <- linterp(York1, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgon1_interpolated$GR)
Gstd = Gsd(Gorgon1_interpolated$GR)
Gorgon1_scaled = (Gorgon1_interpolated$GR - Gmean)/Gstd
Gorgon1_rescaled = data.frame(Gorgon1_interpolated$DEPT, Gorgon1_scaled)

Ymean = Gmean(York1_interpolated$GR)
Ystd = Gsd(York1_interpolated$GR)
York1_scaled = (York1_interpolated$GR - Ymean)/Ystd
York1_rescaled = data.frame(York1_interpolated$DEPT, York1_scaled)

# Resampling the data using moving window statistics
Gorgon1_scaled = mwStats(Gorgon1_rescaled, cols = 2, win = 3, ends = T)
Gorgon1_standardized = data.frame(Gorgon1_scaled$Center_win, Gorgon1_scaled$Average)

York1_scaled = mwStats(York1_rescaled, cols = 2, win = 3, ends = T)
York1_standardized = data.frame(York1_scaled$Center_win, York1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgon1_standardized, type="l", xlim = c(300, 1700), ylim = c(-20, 20), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
plot(York1_standardized, type="l", xlim = c(400, 1500), ylim = c(-20, 20), xlab = "York1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_y1_g1_ap1 <- dtw(York1_standardized$York1_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_y1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
York1_on_Gorgon1_depth = tune(York1_standardized, cbind(York1_standardized$York1_scaled.Center_win[al_y1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_y1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
lines(York1_on_Gorgon1_depth, col = "red")

# DTW Distance measure
al_y1_g1_ap1$normalizedDistance
al_y1_g1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(York1_standardized),ncol=nrow(Gorgon1_standardized))
image(x=Gorgon1_standardized[,1],y=York1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Gorgon1_standardized[,1] - 370))
base_1_y <- which.min(abs(York1_standardized[,1] - 450))

# Depth values for second datum
base_2_x <- which.min(abs(Gorgon1_standardized[,1] - 430))
base_2_y <- which.min(abs(York1_standardized[,1] - 500))

# Depth values for third datum
base_3_x <- which.min(abs(Gorgon1_standardized[,1] - 545))
base_3_y <- which.min(abs(York1_standardized[,1] - 600))

# Depth values for fourth datum
base_4_x <- which.min(abs(Gorgon1_standardized[,1] - 1010))
base_4_y <- which.min(abs(York1_standardized[,1] - 1000))

# Depth values for fifth datum
base_5_x <- which.min(abs(Gorgon1_standardized[,1] - 1190))
base_5_y <- which.min(abs(York1_standardized[,1] - 1170))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(York1_standardized), ncol = nrow(Gorgon1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(York1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-300),(base_1_x+300):ncol(compare.window)] <- 0

compare.window[(base_2_y+400):nrow(York1_standardized),1:(base_2_x-400)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+400):nrow(York1_standardized),1:(base_3_x-400)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(York1_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-300),(base_4_x+300):ncol(compare.window)] <- 0

compare.window[(base_5_y+500):nrow(York1_standardized),1:(base_5_x-500)] <- 0
compare.window[1:(base_5_y-400),(base_5_x+400):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Gorgon1_standardized[,1],y=York1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Gorgon1_standardized[,1],y=York1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_y1_g1_ap1 <- dtw(York1_standardized$York1_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_y1_g1_ap1, type = "threeway")

# DTW Distance measure
al_y1_g1_ap1$normalizedDistance
al_y1_g1_ap1$distance

image(y = Gorgon1_standardized[,1], x = York1_standardized[,1], z = compare.window, useRaster = T)
lines(York1_standardized$York1_scaled.Center_win[al_y1_g1_ap1$index1], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_y1_g1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
York1_on_Gorgon1_depth = tune(York1_standardized, cbind(York1_standardized$York1_scaled.Center_win[al_y1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_y1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR (York-1)")
lines(York1_on_Gorgon1_depth, col = "red")

# Tuning York1 data on Bluebell1 depth
York1_on_Bluebell1_depth = tune(York1_on_Gorgon1_depth, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

# Tuning York1 data on Fisher1 depth
York1_on_Fisher1_depth = tune(York1_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning York1 data on Finucane1 depth
York1_on_Finucane1_depth = tune(York1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning York1 data on Picard1 depth
York1_on_Picard1_depth = tune(York1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (York-1)")
lines(York1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
York1_originalGR_on_Picard1_depth = data.frame(York1_on_Picard1_depth$X1, York1_interpolated[1:5662,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (York-1)")
lines(York1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to York1 

U1463Age_on_York1_depth = tune(York1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_York1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "York1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_York1_depth) <- new_column_names
write.csv(U1463Age_on_York1_depth, file = "RScripts&Data/Sites Data_Age-NGR/York 1.csv", row.names = FALSE)

York1_age_depth <- approx(x = York1_interpolated$GR,
                                  y = York1_interpolated$DEPT,
                                  xout = U1463Age_on_York1_depth$GR,
                                  rule = 1)$y

York1_agemodel = data.frame(York1_age_depth, U1463Age_on_York1_depth$AGE)
York1_agemodel <- na.omit(York1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(York1_agemodel) <- new_column_names1
plot(York1_agemodel)
write.csv(York1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/York1_DepthAge.csv", row.names = FALSE)
