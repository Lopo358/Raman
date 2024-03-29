# Pre-process Raman spectra using the hyperSpec package

# https://r-hyperspec.github.io/hyperSpec/

#WiTec .spc file format - single spectra
install.packages("hyperSpec")
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

#cells on CaF2
NC = collapseHyperspec("./Data/Single_cells/NC*.spc")
PC = collapseHyperspec("./Data/Single_cells/PC*.spc")

plot(mean(PC), col = "red")
plot(mean(NC), col = "blue", add = TRUE)
plot(PC, col = "red")
plot(NC, col = "blue", add = TRUE)
legend("top",c("NC","PC", "background"), lty = 1, col = c("blue", "red", "black"))
title(main = "Mean P0 NC, PC, and background spectra")

# smooth spectra with Loess
smoothNC = spc.loess(NC, seq(0,3600,4))
smoothPC = spc.loess(PC, seq(0,3600,4))
png("smooth_NC_PC_P0.png")
plot(smoothNC[,,500~3050], col = "blue") 
plot(smoothPC[,,500~3050], col = "red", add = TRUE) 
dev.off()

#polynomaial baseline subtraction (least squares fit)
polyBaseline = function(spectra) {
  plot(spectra[1:4,,], col = rainbow(4))
  baselines = spc.fit.poly.below(spectra[,,500~3050], poly.order = 6) #setting npts.min = 40 to a higher number doesn't seem to help remove the fluorescence arc
  print(spectra)
  plot(baselines[1:4], add = TRUE, col = rainbow(4))
  subbed = spectra[,,500~3050] - baselines
  return(subbed)
}
NC_bl = polyBaseline(smoothNC)
PC_bl = polyBaseline(smoothPC)

#area normalization using the mean spectra
normfnc = function(spectra) {
  factors = 1 / apply (spectra[, , 600~1800], 1, mean) 
  normed <- sweep (spectra, 1, factors, "*")
  plot(mean(normed))
  return(normed)
}
NC_n = normfnc(NC_bl)
PC_n = normfnc(PC_bl)

plot(NC_n, wl.range = c(600~1800,2800~3050), "spcmeansd", col = "blue", xoffset = 750)
plot(PC_n, "spcmeansd", col= "red", add = TRUE)
plot(mean(NC_n), col = "blue")
plot(mean(PC_n), col= "red", add = TRUE)

legend("topleft",c("NC","PC"), lty = 1, col = c("blue", "red"))
title(main = "Mean of processed spectra measured at 532 nm")

#save processed spectra
saveRDS(NC_n, file = "./Output/NC-processed-spectra.rds")
saveRDS(PC_n, file = "./Output/PC-processed-spectra.rds")
