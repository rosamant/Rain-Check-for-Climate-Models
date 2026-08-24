install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Andromeda1 and Hyde1 datasets

Andromeda1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Andromeda1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting signal
An1 = Gmean(Andromeda1[c(5997:9740),2])
An2 = Gmean(Andromeda1[c(3000:5000),2])
SD1 = Gsd(Andromeda1[c(5997:9740),2])
SD2 = Gsd(Andromeda1[c(3000:7000),2])
Andromeda1[c(5997:9740),2]=(Andromeda1[c(5997:9740),2]+(An2-An1))*(SD1/SD2)

# Recorrecting signal
An1 = Gmean(Andromeda1[c(7099:9740),2])
An2 = Gmean(Andromeda1[c(3000:7000),2])
SD1 = Gsd(Andromeda1[c(7099:9740),2])
SD2 = Gsd(Andromeda1[c(3000:7000),2])
Andromeda1[c(7099:9740),2]=(Andromeda1[c(7099:9740),2]+(An2-An1))*(SD1/SD2)

Andromeda1=Andromeda1[c(1:9573),] # Oligocene-Miocene
head(Andromeda1)
plot(Andromeda1, type="l", xlim = c(400, 2200), ylim = c(0, 60))

Hyde1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Hyde1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
H1 = Gmean(Hyde1[c(1:300),2])
H2 = Gmean(Hyde1[c(300:700),2])
SD1 = Gsd(Hyde1[c(1:300),2])
SD2 = Gsd(Hyde1[c(300:700),2])
Hyde1[c(1:300),2]=(Hyde1[c(1:300),2]+(H2-H1))*(SD1/SD2)

Hyde1=Hyde1[c(1:8335),] # Oligocene-Miocene
head(Hyde1)
plot(Hyde1, type="l", xlim = c(400, 2150), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Andromeda1_interpolated <- linterp(Andromeda1, dt = 0.2, genplot = F)
Hyde1_interpolated <- linterp(Hyde1, dt = 0.2, genplot = F)

# Scaling the data
Ammean = Gmean(Andromeda1_interpolated$GR)
Amstd = Gsd(Andromeda1_interpolated$GR)
Andromeda1_scaled = (Andromeda1_interpolated$GR - Ammean)/Amstd
Andromeda1_rescaled = data.frame(Andromeda1_interpolated$DEPT, Andromeda1_scaled)

Hmean = Gmean(Hyde1_interpolated$GR)
Hstd = Gsd(Hyde1_interpolated$GR)
Hyde1_scaled = (Hyde1_interpolated$GR - Hmean)/Hstd
Hyde1_rescaled = data.frame(Hyde1_interpolated$DEPT, Hyde1_scaled)

# Resampling the data using moving window statistics
Andromeda1_scaled = mwStats(Andromeda1_rescaled, cols = 2, win=3, ends = T)
Andromeda1_standardized = data.frame(Andromeda1_scaled$Center_win, Andromeda1_scaled$Average)

Hyde1_scaled = mwStats(Hyde1_rescaled, cols = 2, win=3, ends = T)
Hyde1_standardized = data.frame(Hyde1_scaled$Center_win, Hyde1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Andromeda1_standardized, type="l", xlim = c(600, 2000), ylim = c(-20, 20), xlab = "Andromeda1 Resampled Depth", ylab = "Normalized GR")
plot(Hyde1_standardized, type="l", xlim = c(400, 2100), ylim = c(-20, 20), xlab = "Hyde1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_h1_am1_ap1 <- dtw(Hyde1_standardized$Hyde1_scaled.Average, Andromeda1_standardized$Andromeda1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_h1_am1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Hyde1_on_Andromeda1_depth = tune(Hyde1_standardized, cbind(Hyde1_standardized$Hyde1_scaled.Center_win[al_h1_am1_ap1$index1s], Andromeda1_standardized$Andromeda1_scaled.Center_win[al_h1_am1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Andromeda1_standardized, type = "l", ylim = c(-20, 30), xlim = c(400, 2100), xlab = "Andromeda1 Resampled Depth", ylab = "Normalized GR")
lines(Hyde1_on_Andromeda1_depth, col = "red")

# DTW Distance
al_h1_am1_ap1$normalizedDistance
al_h1_am1_ap1$distance

# Tuning Hyde1 data on Fisher1 depth
Hyde1_on_Fisher1_depth = tune(Hyde1_on_Andromeda1_depth, cbind(Andromeda1_standardized$Andromeda1_scaled.Center_win[al_am1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_am1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Hyde1 data on Finucane1 depth
Hyde1_on_Finucane1_depth = tune(Hyde1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Hyde1 data on Picard1 depth
Hyde1_on_Picard1_depth = tune(Hyde1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Hyde1 Resampled Depth", ylab = "Normalized GR (Hyde-1)")
lines(Hyde1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Hyde1_originalGR_on_Picard1_depth = data.frame(Hyde1_on_Picard1_depth$X1, Hyde1_interpolated[1:8379,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 60), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Hyde-1)")
lines(Hyde1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Hyde1 

U1463Age_on_Hyde1_depth = tune(Hyde1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Hyde1_depth, type = "l", ylim = c(0, 60), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Hyde1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Hyde1_depth) <- new_column_names
write.csv(U1463Age_on_Hyde1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Hyde 1.csv", row.names = FALSE)
