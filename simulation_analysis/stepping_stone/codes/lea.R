start <- Sys.time()
library(LEA)

max_idx <- 15

output = geno2lfmm("genofile.geno", "genofile.lfmm")

project = snmf("genofile.geno",
               K = 1:max_idx, 
               entropy = TRUE, 
               repetitions = 10,
               project = "new")

knee_chord <- function(y, idx = seq_along(y), normalize = TRUE) {
  xx <- idx; yy <- y
  if (normalize) {
    xx <- (xx - min(xx)) / (max(xx) - min(xx))
    yy <- (yy - min(yy)) / (max(yy) - min(yy))
  }
  
  x1 <- xx[1]; y1 <- yy[1]
  x2 <- xx[length(xx)]; y2 <- yy[length(yy)]
  
  a <- y2 - y1
  b <- x1 - x2
  c <- x2*y1 - x1*y2
  
  d <- abs(a*xx + b*yy + c) / sqrt(a^2 + b^2)
  idx[which.max(d)]
}

write.csv(summary(project)$crossEntropy[2,], file="cross_ent.txt", quote=F)
Kmin <- which.min(summary(project)$crossEntropy[2,])

if(Kmin == max_idx){
  Kmin <- knee_chord(summary(project)$crossEntropy[2,])
}

start2 <- Sys.time()

mod <- lfmm2(input="genofile.lfmm",
               env="envfile.env",
               K = Kmin)

pv <- lfmm2.test(object = mod,
                 input="genofile.lfmm",
                 env="envfile.env",
                 linear = TRUE)

write.table(pv$pvalues, file="p-values.txt", row.names=F, col.names=F)
write.table(Kmin, file="K.txt", row.names=F, col.names=F)


end <- Sys.time()
elapsed <- as.numeric(difftime(end, start, units = "secs"))
elapsed2 <- as.numeric(difftime(end, start2, units = "secs"))
cat("Elapsed seconds:", elapsed, "\nElapsed seconds (only second part):", elapsed2, "\n", file = "runtime.log", append = TRUE)
