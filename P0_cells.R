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
NC = collapseHyperspec("./Data/NC*.spc")
PC = collapseHyperspec("./Data/PC*.spc")

plot(mean(PC), col = "red")
plot(mean(NC), col = "blue", add = TRUE)
plot(PC, col = "red")
plot(NC, col = "blue", add = TRUE)
legend("top",c("NC","PC", "background"), lty = 1, col = c("blue", "red", "black"))
title(main = "Mean P0 NC, PC, and background spectra")

smoothNC = spc.loess(NC, seq(0,3600,4))
smoothPC = spc.loess(PC, seq(0,3600,4))
png("smooth_NC_PC_P0.png")
plot(smoothNC[,,500~3050], col = "blue") 
plot(smoothPC[,,500~3050], col = "red", add = TRUE) 
dev.off()
