install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Minilya1 and Phoenix1 datasets

Minilya1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Minilya_1.csv", header=TRUE, stringsAsFactors=FALSE)
Minilya1=Minilya1[c(1:5613),] # Eocene-Miocene Unconformity
head(Minilya1)
plot(Minilya1, type="l", xlim = c(150, 1100), ylim = c(0, 50))

Phoenix1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/PHOENIX_1.csv", header=TRUE, stringsAsFactors=FALSE)
Phoenix1=Phoenix1[c(1:4140),] # Eocene-Miocene Unconformity
head(Phoenix1)
plot(Phoenix1, type="l", xlim = c(150, 800), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Minilya1_interpolated <- linterp(Minilya1, dt = 0.2, genplot = F)
Phoenix1_interpolated <- linterp(Phoenix1, dt = 0.2, genplot = F)

# Scaling the data

Mmean = Gmean(Minilya1_interpolated$GR)
Mstd = Gsd(Minilya1_interpolated$GR)
Minilya1_scaled = (Minilya1_interpolated$GR - Mmean)/Mstd
Minilya1_rescaled = data.frame(Minilya1_interpolated$DEPT, Minilya1_scaled)

Phmean = Gmean(Phoenix1_interpolated$GR)
Phstd = Gsd(Phoenix1_interpolated$GR)
Phoenix1_scaled = (Phoenix1_interpolated$GR - Phmean)/Phstd
Phoenix1_rescaled = data.frame(Phoenix1_interpolated$DEPT, Phoenix1_scaled)

# Resampling the data using moving window statistics

Minilya1_scaled = mwStats(Minilya1_rescaled, cols = 2, win = 3, ends = T)
Minilya1_standardized = data.frame(Minilya1_scaled$Center_win, Minilya1_scaled$Average)

Phoenix1_scaled = mwStats(Phoenix1_rescaled, cols = 2, win = 3, ends = T)
Phoenix1_standardized = data.frame(Phoenix1_scaled$Center_win, Phoenix1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Minilya1_standardized, type="l", xlim = c(150, 1100), ylim = c(-20, 20), xlab = "Minilya1 Resampled Depth", ylab = "Normalized GR")
plot(Phoenix1_standardized, type="l", xlim = c(150, 800), ylim = c(-20, 20), xlab = "Phoenix1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_ph1_m1_ap1 <- dtw(Phoenix1_standardized$Phoenix1_scaled.Average, Minilya1_standardized$Minilya1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_ph1_m1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Phoenix1_on_Minilya1_depth = tune(Phoenix1_standardized, cbind(Phoenix1_standardized$Phoenix1_scaled.Center_win[al_ph1_m1_ap1$index1s], Minilya1_standardized$Minilya1_scaled.Center_win[al_ph1_m1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Minilya1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1100), xlab = "Minilya1 Resampled Depth", ylab = "Normalized GR")
lines(Phoenix1_on_Minilya1_depth, col = "red")

# DTW Distance measure
al_ph1_m1_ap1$normalizedDistance
al_ph1_m1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Phoenix1_standardized),ncol=nrow(Minilya1_standardized))
image(x=Minilya1_standardized[,1],y=Phoenix1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Minilya1_standardized[,1] - 300))
base_1_y <- which.min(abs(Phoenix1_standardized[,1] - 320))

# Depth values for second datum
base_2_x <- which.min(abs(Minilya1_standardized[,1] - 420))
base_2_y <- which.min(abs(Phoenix1_standardized[,1] - 430))

# Depth values for third datum
base_3_x <- which.min(abs(Minilya1_standardized[,1] - 510))
base_3_y <- which.min(abs(Phoenix1_standardized[,1] - 480))

# Depth values for fourth datum
base_4_x <- which.min(abs(Minilya1_standardized[,1] - 680))
base_4_y <- which.min(abs(Phoenix1_standardized[,1] - 590))

# Depth values for fifth datum
base_5_x <- which.min(abs(Minilya1_standardized[,1] - 810))
base_5_y <- which.min(abs(Phoenix1_standardized[,1] - 660))

# Depth values for sixth datum
base_6_x <- which.min(abs(Minilya1_standardized[,1] - 990))
base_6_y <- which.min(abs(Phoenix1_standardized[,1] - 780))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Phoenix1_standardized), ncol = nrow(Minilya1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+100):nrow(Phoenix1_standardized),1:(base_1_x-100)] <- 0
compare.window[1:(base_1_y-100),(base_1_x+100):ncol(compare.window)] <- 0

compare.window[(base_2_y+100):nrow(Phoenix1_standardized),1:(base_2_x-100)] <- 0
compare.window[1:(base_2_y-100),(base_2_x+100):ncol(compare.window)] <- 0

compare.window[(base_3_y+100):nrow(Phoenix1_standardized),1:(base_3_x-100)] <- 0
compare.window[1:(base_3_y-100),(base_3_x+100):ncol(compare.window)] <- 0

compare.window[(base_4_y+50):nrow(Phoenix1_standardized),1:(base_4_x-50)] <- 0
compare.window[1:(base_4_y-50),(base_4_x+50):ncol(compare.window)] <- 0

compare.window[(base_5_y+50):nrow(Phoenix1_standardized),1:(base_5_x-50)] <- 0
compare.window[1:(base_5_y-50),(base_5_x+50):ncol(compare.window)] <- 0

compare.window[(base_6_y+50):nrow(Phoenix1_standardized),1:(base_6_x-50)] <- 0
compare.window[1:(base_6_y-50),(base_6_x+50):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Minilya1_standardized[,1],y=Phoenix1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Minilya1_standardized[,1],y=Phoenix1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_ph1_m1_ap1 <- dtw(Phoenix1_standardized$Phoenix1_scaled.Average, Minilya1_standardized$Minilya1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_ph1_m1_ap1, type = "threeway")

# DTW Distance measure
al_ph1_m1_ap1$normalizedDistance
al_ph1_m1_ap1$distance

image(y = Minilya1_standardized[,1], x = Phoenix1_standardized[,1], z = compare.window, useRaster = T, xlab = "Phoenix-1", ylab = "Minilya-1", cex.lab = 1.25, cex.axis = 1.25)
lines(Phoenix1_standardized$Phoenix1_scaled.Center_win[al_ph1_m1_ap1$index1], Minilya1_standardized$Minilya1_scaled.Center_win[al_ph1_m1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Phoenix1_on_Minilya1_depth = tune(Phoenix1_standardized, cbind(Phoenix1_standardized$Phoenix1_scaled.Center_win[al_ph1_m1_ap1$index1s], Minilya1_standardized$Minilya1_scaled.Center_win[al_ph1_m1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Minilya1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1100), xlab = "Minilya1 Resampled Depth", ylab = "Normalized GR (Phoenix-1)")
lines(Phoenix1_on_Minilya1_depth, col = "red")

# Tuning Phoenix1 data on Picard1 depth
Phoenix1_on_Picard1_depth = tune(Phoenix1_on_Minilya1_depth, cbind(Minilya1_standardized$Minilya1_scaled.Center_win[al_m1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_m1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Fisher-1)")
lines(Phoenix1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Phoenix1_originalGR_on_Picard1_depth = data.frame(Phoenix1_on_Picard1_depth$X1, Phoenix1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 60), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Fisher-1)")
lines(Phoenix1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Phoenix1 

U1463Age_on_Phoenix1_depth = tune(Phoenix1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Phoenix1_depth, type = "l", ylim = c(0, 60), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Phoenix1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Phoenix1_depth) <- new_column_names
write.csv(U1463Age_on_Phoenix1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Phoenix 1.csv", row.names = FALSE)
