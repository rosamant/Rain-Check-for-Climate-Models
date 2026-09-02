install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Plymouth1 and Lacerta1 datasets

Plymouth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Plymouth 1.csv", header=TRUE, stringsAsFactors=FALSE)
Plymouth1=Plymouth1[c(1:7444),] # Oligocene-Miocene
head(Plymouth1)
plot(Plymouth1, type="l", xlim = c(200, 1650), ylim = c(0, 70))

Lacerta1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Lacerta 1.csv", header=TRUE, stringsAsFactors=FALSE)
Lacerta1=Lacerta1[c(1:7429),] # Oligocene-Miocene
head(Lacerta1)
plot(Lacerta1, type="l", xlim = c(200, 1650), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Plymouth1_interpolated <- linterp(Plymouth1, dt = 0.2, genplot = F)
Lacerta1_interpolated <- linterp(Lacerta1, dt = 0.2, genplot = F)

# Scaling the data
Plmean = Gmean(Plymouth1_interpolated$GR)
Plstd = Gsd(Plymouth1_interpolated$GR)
Plymouth1_scaled = (Plymouth1_interpolated$GR - Plmean)/Plstd
Plymouth1_rescaled = data.frame(Plymouth1_interpolated$DEPT, Plymouth1_scaled)

Lamean = Gmean(Lacerta1_interpolated$GR)
Lastd = Gsd(Lacerta1_interpolated$GR)
Lacerta1_scaled = (Lacerta1_interpolated$GR - Lamean)/Lastd
Lacerta1_rescaled = data.frame(Lacerta1_interpolated$DEPT, Lacerta1_scaled)

# Resampling the data using moving window statistics
Plymouth1_scaled = mwStats(Plymouth1_rescaled, cols = 2, win=3, ends = T)
Plymouth1_standardized = data.frame(Plymouth1_scaled$Center_win, Plymouth1_scaled$Average)

Lacerta1_scaled = mwStats(Lacerta1_rescaled, cols = 2, win=3, ends = T)
Lacerta1_standardized = data.frame(Lacerta1_scaled$Center_win, Lacerta1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Plymouth1_standardized, type="l", xlim = c(200, 1650), ylim = c(-20, 20), xlab = "Plymouth1 Resampled Depth", ylab = "Normalized GR")
plot(Lacerta1_standardized, type="l", xlim = c(200, 1650), ylim = c(-20, 20), xlab = "Lacerta1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_la1_pl1_ap1 <- dtw(Lacerta1_standardized$Lacerta1_scaled.Average, Plymouth1_standardized$Plymouth1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_la1_pl1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Lacerta1_on_Plymouth1_depth = tune(Lacerta1_standardized, cbind(Lacerta1_standardized$Lacerta1_scaled.Center_win[al_la1_pl1_ap1$index1s], Plymouth1_standardized$Plymouth1_scaled.Center_win[al_la1_pl1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Plymouth1_standardized, type = "l", ylim = c(-20, 20), xlim = c(200, 1700), xlab = "Plymouth1 Resampled Depth", ylab = "Normalized GR")
lines(Lacerta1_on_Plymouth1_depth, col = "red")

# DTW Distance 
al_la1_pl1_ap1$normalizedDistance
al_la1_pl1_ap1$distance

# Tuning Lacerta1 data on Finucane1 depth
Lacerta1_on_Finucane1_depth = tune(Lacerta1_on_Plymouth1_depth, cbind(Plymouth1_standardized$Plymouth1_scaled.Center_win[al_pl1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_pl1_f1_ap1$index2s]), extrapolate = F)

# Tuning Lacerta1 data on Picard1 depth
Lacerta1_on_Picard1_depth = tune(Lacerta1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Lacerta1_originalGR_on_Picard1_depth = data.frame(Lacerta1_on_Picard1_depth$X1, Lacerta1_interpolated$GR)

dev.off()
plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Plymouth-1)")
lines(Lacerta1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Lacerta1 

U1463Age_on_Lacerta1_depth = tune(Lacerta1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Lacerta1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Lacerta1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Lacerta1_depth) <- new_column_names
write.csv(U1463Age_on_Lacerta1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Lacerta 1.csv", row.names = FALSE)
