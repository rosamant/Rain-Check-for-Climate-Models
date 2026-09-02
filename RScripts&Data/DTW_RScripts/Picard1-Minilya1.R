install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and Minilya1 datasets

Picard1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PICARD 1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7750),] # Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1, type="l", xlim = c(150, 1300), ylim = c(0, 50))

Minilya1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Minilya_1.csv", header=TRUE, stringsAsFactors=FALSE)
Minilya1=Minilya1[c(1:5613),] # Eocene-Miocene Unconformity
head(Minilya1)
plot(Minilya1, type="l", xlim = c(150, 1100), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Picard1_interpolated <- linterp(Picard1, dt = 0.2, genplot = F)
Minilya1_interpolated <- linterp(Minilya1, dt = 0.2, genplot = F)

# Scaling the data
Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)
Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled)

Mmean = Gmean(Minilya1_interpolated$GR)
Mstd = Gsd(Minilya1_interpolated$GR)
Minilya1_scaled = (Minilya1_interpolated$GR - Mmean)/Mstd
Minilya1_rescaled = data.frame(Minilya1_interpolated$DEPT, Minilya1_scaled)

# Resampling the data using moving window statistics
Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win = 3, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average)

Minilya1_scaled = mwStats(Minilya1_rescaled, cols = 2, win = 3, ends = T)
Minilya1_standardized = data.frame(Minilya1_scaled$Center_win, Minilya1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Picard1_standardized, type="l", xlim = c(150, 1300), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
plot(Minilya1_standardized, type="l", xlim = c(150, 1100), ylim = c(-20, 20), xlab = "Minilya1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_m1_p1_ap1 <- dtw(Minilya1_standardized$Minilya1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_m1_p1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Minilya1_on_Picard1_depth = tune(Minilya1_standardized, cbind(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(Minilya1_on_Picard1_depth, col = "red")

# DTW Distance measure
al_m1_p1_ap1$normalizedDistance
al_m1_p1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Minilya1_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Picard1_standardized[,1] - 260))
base_1_y <- which.min(abs(Minilya1_standardized[,1] - 300))

# Depth values for second datum
base_2_x <- which.min(abs(Picard1_standardized[,1] - 370))
base_2_y <- which.min(abs(Minilya1_standardized[,1] - 420))

# Depth values for third datum
base_3_x <- which.min(abs(Picard1_standardized[,1] - 430))
base_3_y <- which.min(abs(Minilya1_standardized[,1] - 510))

# Depth values for fourth datum
base_4_x <- which.min(abs(Picard1_standardized[,1] - 545))
base_4_y <- which.min(abs(Minilya1_standardized[,1] - 620))

# Depth values for fifth datum
base_5_x <- which.min(abs(Picard1_standardized[,1] - 993))
base_5_y <- which.min(abs(Minilya1_standardized[,1] - 790))

# Depth values for sixth datum
base_6_x <- which.min(abs(Picard1_standardized[,1] - 1190))
base_6_y <- which.min(abs(Minilya1_standardized[,1] - 940))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Minilya1_standardized), ncol = nrow(Picard1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+100):nrow(Minilya1_standardized),1:(base_1_x-100)] <- 0
compare.window[1:(base_1_y-100),(base_1_x+100):ncol(compare.window)] <- 0

compare.window[(base_2_y+100):nrow(Minilya1_standardized),1:(base_2_x-100)] <- 0
compare.window[1:(base_2_y-100),(base_2_x+100):ncol(compare.window)] <- 0

compare.window[(base_3_y+100):nrow(Minilya1_standardized),1:(base_3_x-100)] <- 0
compare.window[1:(base_3_y-100),(base_3_x+100):ncol(compare.window)] <- 0

compare.window[(base_4_y+100):nrow(Minilya1_standardized),1:(base_4_x-100)] <- 0
compare.window[1:(base_4_y-100),(base_4_x+100):ncol(compare.window)] <- 0

compare.window[(base_5_y+200):nrow(Minilya1_standardized),1:(base_5_x-200)] <- 0
compare.window[1:(base_5_y-200),(base_5_x+200):ncol(compare.window)] <- 0

compare.window[(base_6_y+450):nrow(Minilya1_standardized),1:(base_6_x-450)] <- 0
compare.window[1:(base_6_y-200),(base_6_x+200):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Picard1_standardized[,1],y=Minilya1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_m1_p1_ap1 <- dtw(Minilya1_standardized$Minilya1_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = F))
plot(al_m1_p1_ap1, type = "threeway")

# DTW Distance measure
al_m1_p1_ap1$normalizedDistance
al_m1_p1_ap1$distance

image(y = Picard1_standardized[,1], x = Minilya1_standardized[,1], z = compare.window, useRaster = T, xlab = "Minilya-1", ylab = "Picard-1", cex.lab = 1.25, cex.axis = 1.25)
lines(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Minilya1_on_Picard1_depth = tune(Minilya1_standardized, cbind(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Minilya-1)")
lines(Minilya1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Minilya1_originalGR_on_Picard1_depth = data.frame(Minilya1_on_Picard1_depth$X1, Minilya1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Minilya-1)")
lines(Minilya1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Picard1 

U1463Age_on_Picard1_depth = tune(Picard1_originalGR, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Picard1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Picard-1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

# Tuning the age model data to Minilya1 

U1463Age_on_Minilya1_depth = tune(Minilya1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Minilya1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Minilya-1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Picard1_depth) <- new_column_names
write.csv(U1463Age_on_Picard1_depth, file = "C:/Users/Rohit/OneDrive - Universität Münster/Maps/Base Map/U1463_Age-Site_GR/Picard 1.csv", row.names = FALSE)

colnames(U1463Age_on_Minilya1_depth) <- new_column_names
write.csv(U1463Age_on_Minilya1_depth, file = "C:/Users/Rohit/OneDrive - Universität Münster/Maps/Base Map/U1463_Age-Site_GR/Minilya 1.csv", row.names = FALSE)

Picard1_age_depth <- approx(x = Picard1_interpolated$GR,
                              y = Picard1_interpolated$DEPT,
                              xout = U1463Age_on_Picard1_depth$GR,
                              rule = 1)$y

Picard1_agemodel = data.frame(Picard1_age_depth, U1463Age_on_Picard1_depth$AGE)
Picard1_agemodel <- na.omit(Picard1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Picard1_agemodel) <- new_column_names1
plot(Picard1_agemodel)
write.csv(Picard1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Picard1_DepthAge.csv", row.names = FALSE)
