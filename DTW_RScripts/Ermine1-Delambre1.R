install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Ermine1 and Delambre1 datasets

Ermine1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Ermine1.csv", header=TRUE, stringsAsFactors=FALSE)
Ermine1=Ermine1[c(559:14259),] # Oligocene-Miocene
head(Ermine1)
plot(Ermine1, type="l", xlim = c(600, 2000), ylim = c(0, 70))

Delambre1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Delambre1.csv", header=TRUE, stringsAsFactors=FALSE)
Delambre1=Delambre1[c(1:6318),] # Oligocene-Miocene
head(Delambre1)
plot(Delambre1, type="l", xlim = c(900, 1850), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Ermine1_interpolated <- linterp(Ermine1, dt = 0.2, genplot = F)
Delambre1_interpolated <- linterp(Delambre1, dt = 0.2, genplot = F)

# Scaling the data
Emean = Gmean(Ermine1_interpolated$GR)
Estd = Gsd(Ermine1_interpolated$GR)
Ermine1_scaled = (Ermine1_interpolated$GR - Emean)/Estd
Ermine1_rescaled = data.frame(Ermine1_interpolated$DEPT, Ermine1_scaled)

Dmean = Gmean(Delambre1_interpolated$GR)
Dstd = Gsd(Delambre1_interpolated$GR)
Delambre1_scaled = (Delambre1_interpolated$GR - Dmean)/Dstd
Delambre1_rescaled = data.frame(Delambre1_interpolated$DEPT, Delambre1_scaled)

# Resampling the data using moving window statistics
Ermine1_scaled = mwStats(Ermine1_rescaled, cols = 2, win=3, ends = T)
Ermine1_standardized = data.frame(Ermine1_scaled$Center_win, Ermine1_scaled$Average)

Delambre1_scaled = mwStats(Delambre1_rescaled, cols = 2, win=3, ends = T)
Delambre1_standardized = data.frame(Delambre1_scaled$Center_win, Delambre1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Ermine1_standardized, type="l", xlim = c(600, 2000), ylim = c(-20, 20), xlab = "Ermine1 Resampled Depth", ylab = "Normalized GR")
plot(Delambre1_standardized, type="l", xlim = c(900, 1900), ylim = c(-20, 20), xlab = "Delambre1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_d1_e1_ap1 <- dtw(Delambre1_standardized$Delambre1_scaled.Average, Ermine1_standardized$Ermine1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_d1_e1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Delambre1_on_Ermine1_depth = tune(Delambre1_standardized, cbind(Delambre1_standardized$Delambre1_scaled.Center_win[al_d1_e1_ap1$index1s], Ermine1_standardized$Ermine1_scaled.Center_win[al_d1_e1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Ermine1_standardized, type = "l", ylim = c(-20, 20), xlim = c(600, 2000), xlab = "Ermine1 Resampled Depth", ylab = "Normalized GR")
lines(Delambre1_on_Ermine1_depth, col = "red")

# DTW Distance
al_d1_e1_ap1$normalizedDistance
al_d1_e1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Delambre1_standardized),ncol=nrow(Ermine1_standardized))
image(x=Ermine1_standardized[,1],y=Delambre1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Ermine1_standardized[,1] - 1050))
base_1_y <- which.min(abs(Delambre1_standardized[,1] - 1240))

# Depth values for second datum
base_2_x <- which.min(abs(Ermine1_standardized[,1] - 1150))
base_2_y <- which.min(abs(Delambre1_standardized[,1] - 1470))

# Depth values for third datum
base_3_x <- which.min(abs(Ermine1_standardized[,1] - 1530))
base_3_y <- which.min(abs(Delambre1_standardized[,1] - 1720))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Delambre1_standardized), ncol = nrow(Ermine1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(Delambre1_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-300),(base_1_x+300):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(Delambre1_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-300),(base_2_x+300):ncol(compare.window)] <- 0

compare.window[(base_3_y+350):nrow(Delambre1_standardized),1:(base_3_x-350)] <- 0
compare.window[1:(base_3_y-350),(base_3_x+350):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Ermine1_standardized[,1],y=Delambre1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Ermine1_standardized[,1],y=Delambre1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_d1_e1_ap1 <- dtw(Delambre1_standardized$Delambre1_scaled.Average, Ermine1_standardized$Ermine1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_d1_e1_ap1, type = "threeway")

# DTW Distance measure
al_d1_e1_ap1$normalizedDistance
al_d1_e1_ap1$distance

image(y = Ermine1_standardized[,1], x = Delambre1_standardized[,1], z = compare.window, useRaster = T)
lines(Delambre1_standardized$Delambre1_scaled.Center_win[al_d1_e1_ap1$index1], Ermine1_standardized$Ermine1_scaled.Center_win[al_d1_e1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Delambre1_on_Ermine1_depth = tune(Delambre1_standardized, cbind(Delambre1_standardized$Delambre1_scaled.Center_win[al_d1_e1_ap1$index1s], Ermine1_standardized$Ermine1_scaled.Center_win[al_d1_e1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Ermine1_standardized, type = "l", ylim = c(-20, 30), xlim = c(600, 2000), xlab = "Ermine1 Resampled Depth", ylab = "Normalized GR (Delambre-1)")
lines(Delambre1_on_Ermine1_depth, col = "red")

# Tuning Delambre1 data on Hyde1 depth
Delambre1_on_Hyde1_depth = tune(Delambre1_on_Ermine1_depth, cbind(Ermine1_standardized$Ermine1_scaled.Center_win[al_e1_h1_ap1$index1s], Hyde1_standardized$Hyde1_scaled.Center_win[al_e1_h1_ap1$index2s]), extrapolate = F)

# Tuning Delambre1 data on Andromeda1 depth
Delambre1_on_Andromeda1_depth = tune(Delambre1_on_Hyde1_depth, cbind(Hyde1_standardized$Hyde1_scaled.Center_win[al_h1_am1_ap1$index1s], Andromeda1_standardized$Andromeda1_scaled.Center_win[al_h1_am1_ap1$index2s]), extrapolate = F)

# Tuning Delambre1 data on Fisher1 depth
Delambre1_on_Fisher1_depth = tune(Delambre1_on_Andromeda1_depth, cbind(Andromeda1_standardized$Andromeda1_scaled.Center_win[al_am1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_am1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Delambre1 data on Finucane1 depth
Delambre1_on_Finucane1_depth = tune(Delambre1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Delambre1 data on Picard1 depth
Delambre1_on_Picard1_depth = tune(Delambre1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Delambre-1)")
lines(Delambre1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Delambre1_originalGR_on_Picard1_depth = data.frame(Delambre1_on_Picard1_depth$X1, Delambre1_interpolated[1:4815,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 70), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Delambre-1)")
lines(Delambre1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Delambre1 

U1463Age_on_Delambre1_depth = tune(Delambre1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Delambre1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Delambre1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Delambre1_depth) <- new_column_names
write.csv(U1463Age_on_Delambre1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Delambre 1.csv", row.names = FALSE)
