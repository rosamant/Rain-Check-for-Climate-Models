install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and U1464 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

U1464 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/U1464-HSGR.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
M1 = Gmean(U1464[c(1:551),2])
M2 = Gmean(U1464[c(535:900),2])
SD1 = Gsd(U1464[c(1:551),2])
SD2 = Gsd(U1464[c(535:900),2])
U1464[c(1:553),2]=(U1464[c(1:553),2]+(M2-M1))*(SD1/SD2)

head(U1464)
plot(U1464, type="l", xlim = c(0, 800), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
U1464_interpolated <- linterp(U1464, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Umean = Gmean(U1464_interpolated$HSGR)
Ustd = Gsd(U1464_interpolated$HSGR)
U1464_scaled = (U1464_interpolated$HSGR - Umean)/Ustd
U1464_rescaled = data.frame(U1464_interpolated$DEPTH_WMSF, U1464_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

U1464_scaled = mwStats(U1464_rescaled, cols = 2, win=3, ends = T)
U1464_standardized = data.frame(U1464_scaled$Center_win, U1464_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(U1464_standardized, type="l", xlim = c(0, 800), ylim = c(-20, 20), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_U1464_p1_ap1 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_U1464_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth = tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth, col = "red")

# DTW Distance
al_U1464_p1_ap1$normalizedDistance
al_U1464_p1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(U1464_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 190))
base_1_y <- which.min(abs(U1464_standardized[,1] - 50))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 270))
base_2_y <- which.min(abs(U1464_standardized[,1] - 120))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 310))
base_3_y <- which.min(abs(U1464_standardized[,1] - 184))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 390))
base_4_y <- which.min(abs(U1464_standardized[,1] - 275))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 1010))
base_5_y <- which.min(abs(U1464_standardized[,1] - 700))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(U1464_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(U1464_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+400):nrow(U1464_standardized),1:(base_2_x-400)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(U1464_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-400),(base_3_x+400):ncol(compare.window)] <- 0

compare.window[(base_4_y+300):nrow(U1464_standardized),1:(base_4_x-300)] <- 0
compare.window[1:(base_4_y-500),(base_4_x+500):ncol(compare.window)] <- 0

compare.window[(base_5_y+100):nrow(U1464_standardized),1:(base_5_x-100)] <- 0
compare.window[1:(base_5_y-100),(base_5_x+100):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_U1464_p1_ap1 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_U1464_p1_ap1, type = "threeway")

# DTW Distance measure
al_U1464_p1_ap1$normalizedDistance
al_U1464_p1_ap1$distance

image(y = Picard1_standardized[,1], x = U1464_standardized[,1], z = compare.window, useRaster = T)
lines(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth = tune(U1464_standardized, cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (U1464)")
lines(U1464_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
U1464_originalGR_on_Picard1_depth = data.frame(U1464_on_Picard1_depth$X1, U1464_interpolated$HSGR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (U1464)")
lines(U1464_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to U1464 

U1463Age_on_U1464_depth = tune(U1464_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_U1464_depth, type = "l", ylim = c(0, 60), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "U1464")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_U1464_depth) <- new_column_names
write.csv(U1463Age_on_U1464_depth, file = "C:/Users/Rohit/OneDrive - Universität Münster/Maps/Base Map/U1463_Age-Site_GR/U1464.csv", row.names = FALSE)

U1464_age_depth <- approx(x = U1464_interpolated$HSGR,
                            y = U1464_interpolated$DEPTH_WMSF,
                            xout = U1463Age_on_U1464_depth$GR,
                            rule = 1)$y

U1464_agemodel = data.frame(U1464_age_depth, U1463Age_on_U1464_depth$AGE)
U1464_agemodel <- na.omit(U1464_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(U1464_agemodel) <- new_column_names1
plot(U1464_agemodel)
write.csv(U1464_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/U1464_DepthAge.csv", row.names = FALSE)
