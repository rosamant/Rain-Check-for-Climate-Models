install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Hadrian 1 and BankaBanka 1 datasets

Hadrian1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Hadrian 1.csv", header=TRUE, stringsAsFactors=FALSE)
Hadrian1=Hadrian1[c(28:6786),] # Oligocene-Miocene
head(Hadrian1)
plot(Hadrian1, type="l", xlim = c(300, 1350), ylim = c(0, 60))

BankaBanka1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Banka Banka 1.csv", header=TRUE, stringsAsFactors=FALSE)
BankaBanka1=BankaBanka1[c(1:4793),] # Oligocene-Miocene
head(BankaBanka1)
plot(BankaBanka1, type="l", xlim = c(450, 1400), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Hadrian1_interpolated <- linterp(Hadrian1, dt = 0.2, genplot = F)
BankaBanka1_interpolated <- linterp(BankaBanka1, dt = 0.2, genplot = F)

# Scaling the data
Hmean = Gmean(Hadrian1_interpolated$GR)
Hstd = Gsd(Hadrian1_interpolated$GR)
Hadrian1_scaled = (Hadrian1_interpolated$GR - Hmean)/Hstd
Hadrian1_rescaled = data.frame(Hadrian1_interpolated$DEPT, Hadrian1_scaled)

BBmean = Gmean(BankaBanka1_interpolated$GR)
BBstd = Gsd(BankaBanka1_interpolated$GR)
BankaBanka1_scaled = (BankaBanka1_interpolated$GR - BBmean)/BBstd
BankaBanka1_rescaled = data.frame(BankaBanka1_interpolated$DEPT, BankaBanka1_scaled)

# Resampling the data using moving window statistics
Hadrian1_scaled = mwStats(Hadrian1_rescaled, cols = 2, win=3, ends = T)
Hadrian1_standardized = data.frame(Hadrian1_scaled$Center_win, Hadrian1_scaled$Average)

BankaBanka1_scaled = mwStats(BankaBanka1_rescaled, cols = 2, win=3, ends = T)
BankaBanka1_standardized = data.frame(BankaBanka1_scaled$Center_win, BankaBanka1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Hadrian1_standardized, type="l", xlim = c(300, 1350), ylim = c(-20, 20), xlab = "Hadrian1 Resampled Depth", ylab = "Normalized GR")
plot(BankaBanka1_standardized, type="l", xlim = c(450, 1400), ylim = c(-10, 20), xlab = "BankaBanka 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_bb1_h1_ap1 <- dtw(BankaBanka1_standardized$BankaBanka1_scaled.Average, Hadrian1_standardized$Hadrian1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_bb1_h1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
BankaBanka1_on_Hadrian1_depth = tune(BankaBanka1_standardized, cbind(BankaBanka1_standardized$BankaBanka1_scaled.Center_win[al_bb1_h1_ap1$index1s], Hadrian1_standardized$Hadrian1_scaled.Center_win[al_bb1_h1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Hadrian1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1350), xlab = "Hadrian1 Resampled Depth", ylab = "Normalized GR")
lines(BankaBanka1_on_Hadrian1_depth, col = "red")

# DTW Distance
al_bb1_h1_ap1$normalizedDistance
al_bb1_h1_ap1$distance

# Tuning BankaBanka1 data on SG1 depth
BankaBanka1_on_Brontosaurus1_depth = tune(BankaBanka1_on_Hadrian1_depth, cbind(Hadrian1_standardized$Hadrian1_scaled.Center_win[al_h1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_h1_b1_ap1$index2s]), extrapolate = F)

# Tuning BankaBanka1 data on SG1 depth
BankaBanka1_on_SG1_depth = tune(BankaBanka1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (BankaBanka-1)")
lines(BankaBanka1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
BankaBanka1_originalGR_on_SG1_depth = data.frame(BankaBanka1_on_SG1_depth$X1, BankaBanka1_interpolated[1023:4776,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (BankaBanka-1)")
lines(BankaBanka1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to BankaBanka1 

SG1Age_on_BankaBanka1_depth = tune(BankaBanka1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_BankaBanka1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "BankaBanka-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_BankaBanka1_depth) <- new_column_names
SG1Age_on_BankaBanka1_depth[,1] = SG1Age_on_BankaBanka1_depth[,1] * 1000
write.csv(SG1Age_on_BankaBanka1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Banka Banka 1.csv", row.names = FALSE)
