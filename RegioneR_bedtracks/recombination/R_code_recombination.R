

dat <- gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.500kb
#dat <- gilamonster_rptmsk_F4.PREDICT.BSCORRECTED

dat <- subset(gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.500kb, V1 == c("chr5", "chr6"))
chromosomes <- unique(dat$V1)

#pdf("gilamonster_18_chromosomes.pdf", width = 12, height = 7)

par(
  mfrow = c(3, 6),
  mar = c(2.5, 2.5, 2, 0.5),
  oma = c(3, 3, 1, 1)
)

for (chr in chromosomes) {
  
  chr_dat <- dat[dat$V1 == chr, ]
  
  plot(
    log(as.numeric(chr_dat$V4)),
    type = "l",
    ylim = c(-27, -23),
    main = chr,
    xlab = "",
    ylab = "",
    xaxt = "n"
  )
}

mtext("500-kb windows", side = 1, outer = TRUE, line = 1)
mtext("log(value)", side = 2, outer = TRUE, line = 1)

#dev.off()














dat <- subset(gilamonster_rptmsk_F4.PREDICT.BSCORRECTED.500kb, V1 == c("chr5", "chr6"))
chromosomes <- unique(dat$V1)

par(
  mfrow = c(1, 2),
)

for (chr in chromosomes) {
  
  chr_dat <- dat[dat$V1 == chr, ]
  
  plot(
    log(as.numeric(chr_dat$V4)),
    type = "l",
    ylim = c(-27, -23),
    main = chr,
    xlab = "",
    ylab = "",
  )
  
}







