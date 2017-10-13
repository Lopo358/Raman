#FHNW 2017 Sept 15
#Native NC191 and PC30 were also measured
#all tissues were fixed in formalin

#WiTec .spc file format - single spectra
library ("hyperSpec")
collapseHyperspec = function(path) {
  data <- Sys.glob (path) 
  data <- data [seq (1, length (data))] # import low wavenumber region only???
  spc <- lapply (data, read.spc)
  #read.txt.long imports long format ASCII files, i.e. one intensity value per row
  length (spc)  
  spc <- collapse (spc) #into a single hyperspec object
  return (spc)
}

#Native cartilage
NC = collapseHyperspec("./Raman-spectroscopy/2017_FHNW/engineered cartilage/NC191-native-Large area scan_40x*.spc")
PC = collapseHyperspec("./Raman-spectroscopy/2017_FHNW/engineered cartilage/PC30-native-Large area scan_40x*.spc")
#write csv for claudia
write.csv(df, file = "./Raman-spectroscopy/2017_FHNW/engineered cartilage/NC-PC/PC.txt")
wavenums = wl(PC)
df = PC$spc
colnames(df) = wavenums

plot(mean(PC), col = "red", add = TRUE)
plot(mean(NC[1:5]), col = "blue")

legend("top",c("NC","PC", "background"), lty = 1, col = c("blue", "red", "black"))
title(main = "Mean P0 NC, PC, and background spectra")

#median filter to remove cosmic rays
library(FBN)
#apply the median filter with specified windowSize to the rows (1) of the matrix
mfNC = apply(NC[], 1, medianFilter, windowSize = 3) 
mfPC = apply(PC[], 1, medianFilter, windowSize = 3)

plot(mean(mfPC), col = "red", add = TRUE)
plot(mean(mfNC[1:5]), col = "blue")

#smooth cell signal
smoothNC = spc.loess(mfNC, seq(0,3600,4))
smoothPC = spc.loess(mfPC, seq(0,3600,4))

plot(mean(smoothPC[,,500~2000]), col = "red", add = TRUE)
plot(mean(smoothNC[1:5,,]), col = "blue") 

plot(mean(smoothPC), col = "red", add = TRUE)
plot(mean(smoothNC), col = "blue", add = TRUE) 

legend("topleft",c("good d0", "good d7",  "bad d0","bad d7", "NC","PC"), lty = 1, col = c("green","dark green",  "violet", "purple", "blue", "red"))
title(main = "Engineered cartilage at 532 nm, median filtered, smoothed, baseline corrected")

polyBaseline = function(spectra) {
  baselines = spc.fit.poly.below(spectra[,,500~3050], poly.order = 6) #setting npts.min = 40 to a higher number doesn't seem to help remove the fluorescence arc
  subbed = spectra[,,500~3050] - baselines
  return(subbed)
}

#polynomaial baseline (least squares fit)
NC_bl = polyBaseline(smoothNC)
PC_bl = polyBaseline(smoothPC)

plot(mean(PC_bl), col = "red", add = TRUE)
plot(mean(NC_bl), col = "blue", add = TRUE) 

#investigate polybaseline subtraction
lowbaselines = spc.fit.poly.below(mfNC[1:5,,600~1800], poly.order = 3) #poly.order = 1 defualt
highbaselines = spc.fit.poly.below(mfNC[1:5,,2800~3040], poly.order = 6) #setting npts.min = 40 to a higher number doesn't seem to help remove the fluorescence arc
polybaselines = spc.fit.poly.below(smoothNC[,,500~3050], poly.order = 6)
plot(smoothNC[6000:6005,,], col = rainbow(6))
plot(mean(polybaselines[6000:6005,,]), col = rainbow(6), add = TRUE) #poly.order = 6 500~3050 looks good! Fitting with npts.min =  32
sub1 = smoothd0good[1:5,,600~1800] - lowbaselines
sub2 = smoothd0good[1:5,,2800~3600] - highbaselines
subbed = cbind(sub1, sub2)
plot(subbed[1:5])
plotspc (subbed, wl.range = c (600 ~ 1800, 2800 ~ 3600), xoffset = 750) #plot low and high regions with an x-offset
#investigation results = a separate baseline is subtracted from each spectra :D

#normalize to a specific wavenumber
normfnc = function(spectra) {
  factors = 1 / apply (spectra[, , 650~1800], 1, mean) 
  normed <- sweep (spectra, 1, factors, "*")
  return(normed)
}

NC_n <- normfnc(NC_bl)
PC_n <- normfnc(PC_bl)

plot(mean(PC_n), col = "red", add = TRUE)
plot(mean(NC_n[,,600~1800]), col = "blue", add = TRUE) 

plotspc(mean(PC_n), wl.range = c (600 ~ 1800, 2800 ~ 3600), xoffset = 750, col = "red", add = TRUE)
plotspc(mean(NC_n), wl.range = c (600 ~ 1800, 2800 ~ 3600), xoffset = 750, col = "blue", add = TRUE)

a = 1000
b = 1100
plot(mean(NC_n[,,a ~ b]), col = "blue") 
plot(mean(PC_n[,,a ~ b]), col = "red", add = TRUE)

plot(NC_n[,,a ~ b], "spcmeansd", col = "blue") 
plot(PC_n[,,a ~ b], "spcmeansd", col = "red", add = TRUE)

plot(PC_n[,,600~1800], "spcmeansd", col = "red")
plot(NC_n[,,600~1800], "spcmeansd", col = "blue", add = TRUE) 
plot(PC_n[,,2800~3100], "spcmeansd", col = "red")
plot(NC_n[,,2800~3100], "spcmeansd", col = "blue", add = TRUE)

#perhaps there is a lot of background, do cluster analysis to remove these spectra
NC_clusters = kmeans(NC_bl[,,600~1800], 6)
NC_clusters$cluster
NC_bl$cluster <- NC_clusters$cluster #transfer cluster parameter to the original hyperspec object
cluster.means = aggregate(NC_bl, NC_bl$cluster, mean_pm_sd)
plot(cluster.means[,,600~1800], stacked = ".aggregate", fill = ".aggregate")
#yes, clusters 1 and 3 do not have cell signal

legend("topleft",c("Day 0", "Day 7",  "Day 14", "NC"), lty = 1, col = c("light green","green",  "dark green","blue"))
legend("topleft",c("Day 0", "Day 7",  "Day 14", "NC"), lty = 1, col = c("sandy brown","brown1",  "brown4","blue"))
legend("topleft",c( "good 0wk","good 1wk", "good 2wk", "bad 0wk", "bad 1wk", "bad 2wk", "NC"), lty = 1, col = c("light green","green",  "dark green","sandy brown","brown1",  "brown4","blue"))
title(main = " ")

#choose which Raman shifts demonstrate cartilage quality - 1050-1150 cm-1 ?
legend("topleft",c("Native NC", "Native PC"), lty = 1, col = c("blue","red"))

#GAG at 1410 cm-1 (Bergholt 2017) 
#Common GAG peak at 1061 cm-1 (Gamsjaeger 2014 --> Bansil 1978)
peakintensity.fnc = function() {
  f1 = max(mean(NC_n[, ,peakstart~peakend]))
  f3 = max(mean(d0good_n[, ,peakstart~peakend]))
  f4 = max(mean(d7good_n[, ,peakstart~peakend]))
  f5 = max(mean(d14good_n[, ,peakstart~peakend]))
  f6 = max(mean(d0bad_n[, ,peakstart~peakend]))
  f7 = max(mean(d7bad_n[, ,peakstart~peakend]))
  f8 = max(mean(d14bad_n[, ,peakstart~peakend]))
  fmax = max(c(f1 ,f3 ,f4 ,f5 ,f6 , f7, f8))
  f1sd = sd(NC_n[, ,peakstart~peakend]$spc)/fmax
  f3sd = sd(d0good_n[, ,peakstart~peakend]$spc)/fmax
  f4sd = sd(d7good_n[, ,peakstart~peakend]$spc)/fmax
  f5sd = sd(d14good_n[, ,peakstart~peakend]$spc)/fmax
  f6sd = sd(d0bad_n[, ,peakstart~peakend]$spc)/fmax
  f7sd = sd(d7bad_n[, ,peakstart~peakend]$spc)/fmax
  f8sd = sd(d14bad_n[, ,peakstart~peakend]$spc)/fmax
  bardata = c(f1/fmax, f3/fmax, f4/fmax, f5/fmax, f6/fmax, f7/fmax, f8/fmax)
  peak.sd = c(f1sd/fmax, f3sd/fmax, f4sd/fmax, f5sd/fmax, f6sd/fmax, f7sd/fmax, f8sd/fmax)
  df = data.frame(bardata, peak.sd)
  row.names(df) =  c("NC", "good 0wk","good 1wk", "good 2wk", "bad 0wk", "bad 1wk", "bad 2wk")
  return(df)
}
GAGstart = 1060
GAGend = 1065
DNAstart = 780
DNAend = 790
GAGtoDNA.fnc = function() {
  gd1 = mean(rowMeans(NC_n[,,GAGstart~GAGend])$spc / rowMeans(NC_n[,,DNAstart~DNAend])$spc)
  gd1sd = sd(rowMeans(NC_n[,,GAGstart~GAGend])$spc / rowMeans(NC_n[,,DNAstart~DNAend])$spc)
  gd3 = mean(rowMeans(d0good_n[,,GAGstart~GAGend])$spc / rowMeans(d0good_n[,,DNAstart~DNAend])$spc)
  gd3sd = sd(rowMeans(d0good_n[,,GAGstart~GAGend])$spc / rowMeans(d0good_n[,,DNAstart~DNAend])$spc)
  gd4 = mean(rowMeans(d7good_n[,,GAGstart~GAGend])$spc / rowMeans(d7good_n[,,DNAstart~DNAend])$spc)
  gd4sd = sd(rowMeans(d7good_n[,,GAGstart~GAGend])$spc / rowMeans(d7good_n[,,DNAstart~DNAend])$spc)
  gd5 = mean(rowMeans(d14good_n[,,GAGstart~GAGend])$spc / rowMeans(d14good_n[,,DNAstart~DNAend])$spc)
  gd5sd = sd(rowMeans(d14good_n[,,GAGstart~GAGend])$spc / rowMeans(d14good_n[,,DNAstart~DNAend])$spc)
  gd6 = mean(rowMeans(d0bad_n[,,GAGstart~GAGend])$spc / rowMeans(d0bad_n[,,DNAstart~DNAend])$spc)
  gd6sd = sd(rowMeans(d0bad_n[,,GAGstart~GAGend])$spc / rowMeans(d0bad_n[,,DNAstart~DNAend])$spc)
  gd7 = mean(rowMeans(d7bad_n[,,GAGstart~GAGend])$spc / rowMeans(d7bad_n[,,DNAstart~DNAend])$spc)
  gd7sd = sd(rowMeans(d7bad_n[,,GAGstart~GAGend])$spc / rowMeans(d7bad_n[,,DNAstart~DNAend])$spc)
  gd8 = mean(rowMeans(d14bad_n[,,GAGstart~GAGend])$spc / rowMeans(d14bad_n[,,DNAstart~DNAend])$spc)
  gd8sd = sd(rowMeans(d14bad_n[,,GAGstart~GAGend])$spc / rowMeans(d14bad_n[,,DNAstart~DNAend])$spc)
  bardata = c(gd1, gd3, gd4, gd5, gd6, gd7, gd8)
  peak.sd = c(gd1sd, gd3sd, gd4sd, gd5sd, gd6sd, gd7sd, gd8sd)
  df = data.frame(bardata, peak.sd)
  row.names(df) =  c("NC", "good 0wk","good 1wk", "good 2wk", "bad 0wk", "bad 1wk", "bad 2wk")
  return(df)
}

#not ideal, but function input is entered as a global environment variable here:
peakstart = 1060
peakend = 1065
df2 = peakintensity.fnc() 
write.csv(df2, file = "./Raman-spectroscopy/2017_FHNW/engineered cartilage/good-bad GAG to DNA intensity.csv")
df2 = GAGtoDNA.fnc()

p <- ggplot(df2, aes(row.names(df1), bardata)) + 
  geom_col() +
  geom_errorbar(aes(ymin = bardata - peak.sd, ymax = bardata + peak.sd), width=0.2) +
  labs(y="Relative intensity (arb. u.)", x = " ") +
  ggtitle("GAG to DNA peak")
print(p)

#PCA nc and pc
pcafnc2 = function(nc, pc) {
  df = rbind(nc[,,600~1800]$spc, pc[,,600~1800]$spc)
  row.names(df) = c(1:nrow(df))
  pca = prcomp (df, center = TRUE, scale = TRUE)  
  pca$type = c(rep("Native NC",length(nc)), rep("Native PC",length(pc))) #ggbiplot groups
  return(pca)
}
pca = pcafnc2(NC_n, PC_n)

library(devtools)
library(ggbiplot)
#choices are the PCs
g <- ggbiplot(pca, choices = c(1,3), obs.scale = 1, var.scale = 1, 
              groups = pca$type, ellipse = TRUE, 
              circle = TRUE, var.axes=FALSE)
g = g + theme(aspect.ratio = 1)
g = g + scale_color_manual(values = c("blue", "red"))
g <- g + theme(legend.direction = 'vertical', legend.position = 'right')
g = g + ggtitle("Native NC&PC")
print(g)

#principal component weights
wavenums = seq(648, 1800, length.out = length(pca$rotation[,3]))
plot(x = wavenums, y = as.vector(pca$rotation[,3]), type = "l", col = "green", xlab = "Raman shift 1/cm", ylab = "Intensity (Arb. u.)") #first principal component
lines(x = wavenums, y = pca$rotation[,2], type = "l", col = "red")
lines(x = wavenums, y = pca$rotation[,1], type = "l", col = "blue")

legend("topright",c("PC 1", "PC 2", "PC 3"), lty = 1, col = c("blue", "red", "green"))
title("Loadings of first 3 principal components")
