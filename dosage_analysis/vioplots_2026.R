#males
#10, 16, KI01
#females
#30, 35, L007

suppressWarnings(library("vioplot"))
suppressWarnings(library("dplyr"))

'%!in%' <- function(x,y)!('%in%'(x,y))

############Heloderma suspectum
#import males
male_import <- read.delim("Gila_males.gff.tsv")
#import females
female_import <- read.delim("Gila_females.gff.tsv")

#add average expression per sex
male_import$Mmean <- rowMeans(male_import[,14:16], na.rm = T)
female_import$Fmean <- rowMeans(female_import[,14:16], na.rm = T)

################################F vs. M gene expression

#remove genes missing expression
male_import2 <- subset(male_import, Mmean %!in% c('Inf','NaN', "0"))
male_import2 <- male_import2 %>% mutate(value = log2(Mmean))
female_import2 <- subset(female_import, Fmean %!in% c('Inf','NaN', "0"))
female_import2 <- female_import2 %>% mutate(value = log2(Fmean))

#remove genes missing expression in either sex
female_sub <- subset(female_import2, GeneID %in% male_import2$GeneID)
male_sub <- subset(male_import2, GeneID %in% female_sub$GeneID)

#extract chrZ
male_chrZ <- subset(male_sub, Chromosome %in% 'chrZ')
female_chrZ <- subset(female_sub, Chromosome %in% 'chrZ')

#extract autosomes
male_auto <- subset(male_sub, Chromosome %!in% 'chrZ')
female_auto <- subset(female_sub, Chromosome %!in% 'chrZ')

#plot 2-halved violins
F_M_list <- list(
  F_auto = log2(female_auto$Fmean),
  F_Z = log2(female_chrZ$Fmean),
  M_Z = log2(male_chrZ$Mmean),
  M_auto = log2(male_auto$Mmean))

vioplot(F_M_list$F_auto, F_M_list$F_Z, side = "left", plotCentre = "line", col = "grey40", names = NULL,
        ylab = "log2(TPM)", main = "Heloderma F + M expression", ylim = range(-7:15))
vioplot(F_M_list$M_auto, F_M_list$M_Z, side = "right", plotCentre = "line", col = "grey80", add = T)
legend("topright", col = c("grey40", "grey80"), pch = 15,
       bg = "white", legend = c("Female", "Male"), cex = 1.5)
mtext("Heloderma auto                                                Heloderma Z", side=1, line=1)

#M vs. F, chrZ
x <- log2(female_chrZ$Fmean)
y <- log2(male_chrZ$Mmean)
wilcox.test(x, y)

#M vs. F, auto
x <- log2(male_auto$Mmean)
y <- log2(female_auto$Fmean)
wilcox.test(x, y)

#M (Z) vs. M (auto)
x <- log2(male_auto$Mmean)
y <- log2(male_chrZ$Mmean)
wilcox.test(x, y)

#F (Z) vs. F (auto)
x <- log2(female_auto$Fmean)
y <- log2(female_chrZ$Fmean)
wilcox.test(x, y)

write.csv(male_sub, "Heloderma_M_expression.csv", row.names = FALSE)
write.csv(female_sub, "Heloderma_F_expression.csv", row.names = FALSE)


#########################
#plot F/M gene expression
F_M_expression <- female_sub[c("Chromosome","Start","Stop","Name","BioType","Function","Length","Fmean")]
F_M_expression$Mmean <- male_sub$Mmean

#add average expression per sex
F_M_expression <- F_M_expression %>% mutate(quotient = Fmean / Mmean)
F_M_expression <- F_M_expression %>% mutate(value = log2(quotient))

#extract chrZ
F_M_chrZ <- subset(F_M_expression, Chromosome %in% 'chrZ')

#extract autosomes
F_M_auto <- subset(F_M_expression, Chromosome %!in% 'chrZ')

F_M_list2 <- list(
  F_M_auto = F_M_auto$value,
  F_M_chrZ = F_M_chrZ$value)

vioplot(F_M_list2, ylim = range(-12:16), 
        ylab = "log2(F/M expression)", 
        main = "Heloderma F/M expression", names = NULL)
abline(h = 0)
mtext("Heloderma auto                                                Heloderma Z", side=1, line=1)
wilcox.test(F_M_list2$F_M_auto, F_M_list2$F_M_chrZ)
t.test(F_M_list2$F_M_chrZ, F_M_list2$F_M_auto)

write.csv(F_M_expression, "Heloderma_F_M_expression.csv", row.names = FALSE)
write.csv(F_M_auto, "Heloderma_F_M_autosome-expression.csv", row.names = FALSE)
write.csv(F_M_chrZ, "Heloderma_F_M_chrZ-expression.csv", row.names = FALSE)

plot(F_M_chrZ$Start, log2(F_M_chrZ$quotient), 
     pch = 16,
     xlab = "chrZ",
     ylab = "log2(F/M expression)",
     main = "Expression by Position")
abline(h = 0)



#############Gallus gallus
#import males
male_import <- read.delim("Gallus_males.gff.tsv")
#import females
female_import <- read.delim("Gallus_females.gff.tsv")

male_import$Mmean <- rowMeans(male_import[,14:16], na.rm = T)
female_import$Fmean <- rowMeans(female_import[,14:16], na.rm = T)

################################F vs. M gene expression

#remove genes missing expression
male_import2 <- subset(male_import, Mmean %!in% c('Inf','NaN', "0"))
male_import2 <- male_import2 %>% mutate(value = log2(Mmean))
female_import2 <- subset(female_import, Fmean %!in% c('Inf','NaN', "0"))
female_import2 <- female_import2 %>% mutate(value = log2(Fmean))

#remove genes missing expression in either sex
female_sub <- subset(female_import2, GeneID %in% male_import2$GeneID)
male_sub <- subset(male_import2, GeneID %in% female_sub$GeneID)

#extract chrZ
male_chrZ <- subset(male_sub, Chromosome %in% 'NC_052572.1')
female_chrZ <- subset(female_sub, Chromosome %in% 'NC_052572.1')

#extract autosomes
male_auto <- subset(male_sub, Chromosome %!in% 'NC_052572.1')
female_auto <- subset(female_sub, Chromosome %!in% 'NC_052572.1')

F_M_list <- list(
  F_auto = log2(female_auto$Fmean),
  F_Z = log2(female_chrZ$Fmean),
  M_Z = log2(male_chrZ$Mmean),
  M_auto = log2(male_auto$Mmean))

vioplot(F_M_list$F_auto, F_M_list$F_Z, side = "left", plotCentre = "line", col = "grey40", xlab = "",
        ylab = "log2(TPM)", main = "X", ylim = range(-7:15))
vioplot(F_M_list$M_auto, F_M_list$M_Z, side = "right", plotCentre = "line", col = "grey80", add = T)
legend("topright", col = c("grey40", "grey80"), pch = 15,
       bg = "white", legend = c("Female", "Male"), cex = 1.5)
mtext("Gallus auto                                                Gallus Z", side=1, line=1)

#M vs. F, chrZ
x <- log2(female_chrZ$Fmean)
y <- log2(male_chrZ$Mmean)
wilcox.test(x, y)

#M vs. F, auto
x <- log2(male_auto$Mmean)
y <- log2(female_auto$Fmean)
wilcox.test(x, y)

#M (Z) vs. M (auto)
x <- log2(male_auto$Mmean)
y <- log2(male_chrZ$Mmean)
wilcox.test(x, y)

#F (Z) vs. F (auto)
x <- log2(female_auto$Fmean)
y <- log2(female_chrZ$Fmean)
wilcox.test(x, y)

write.csv(male_sub, "Gallus_M_expression.csv", row.names = FALSE)
write.csv(female_sub, "Gallus_F_expression.csv", row.names = FALSE)

#########################
#plot F/M gene expression
F_M_expression <- female_sub[c("Chromosome","Start","Stop","Name","BioType","Function","Length","Fmean")]
F_M_expression$Mmean <- male_sub$Mmean

#add average expression per sex
F_M_expression <- F_M_expression %>% mutate(quotient = Fmean / Mmean)
F_M_expression <- F_M_expression %>% mutate(value = log2(quotient))

#extract chrZ
F_M_chrZ <- subset(F_M_expression, Chromosome %in% 'NC_052572.1')

#extract autosomes
F_M_auto <- subset(F_M_expression, Chromosome %!in% 'NC_052572.1')

F_M_list2 <- list(
  F_M_auto = F_M_auto$value,
  F_M_chrZ = F_M_chrZ$value)

vioplot(F_M_list2, ylim = range(-6:8), 
        ylab = "log2(F/M expression)", 
        main = "Gallus F/M expression", names = NULL)
abline(h = 0)
#mtext("Gallus auto                                              Gallus Z", side=1, line=1)
wilcox.test(F_M_list2$F_M_auto, F_M_list2$F_M_chrZ)
t.test(F_M_list2$F_M_chrZ, F_M_list2$F_M_auto)

write.csv(F_M_expression, "Gallus_F_M_expression.csv", row.names = FALSE)
write.csv(F_M_auto, "Gallus_F_M_autosome-expression.csv", row.names = FALSE)
write.csv(F_M_chrZ, "Gallus_F_M_chrZ-expression.csv", row.names = FALSE)

plot(F_M_chrZ$Start, log2(F_M_chrZ$quotient), 
     pch = 16,
     xlab = "chrZ",
     ylab = "log2(F/M expression)",
     main = "Expression by Position")
abline(h = 0)
