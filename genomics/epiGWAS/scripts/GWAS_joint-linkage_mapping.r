# Title: More GWAS from PhD work

# Clean workspace by removing any previously saved session data
unlink(".RData")

# Record description of the script purpose
info <- paste("Script to perform ExtFlexGWAS()")

# Load required libraries
library(NAM)
library(data.table)

# Load genotype and phenotype dataset
load("CHH_genoPhen.RData")

# Extract phenotype vector and family information
y = test_2$blup
fam = test_2$famno

# Extract genotype matrix (exclude first 3 columns)
gen = test_2[,4:ncol(test_2)]

# Load list of overlapping genes with DMRs (differentially methylated regions)
ovlps <- fread("CHG_DMRovGenes.csv", header = TRUE, sep = ",")
ovlps <- ovlps[,-1]  # Remove first column (assumed to be row numbers or IDs)

# Subset genotype matrix to retain only columns matching overlapping DMRs
t <- which(colnames(gen) %in% ovlps$DMR)
gen <- data.matrix(gen)

# Save intermediate genotype matrix
save(gen, file = "CHG_genOVLPS.Rdata")

# Quality control: filter markers with >50% missing data
MissPerc = colMeans(is.na(gen))
gen = gen[, MissPerc < 0.5]

# Save filtered genotype matrix
save(gen, file = "CHG_genOVLPSminusMiss.Rdata")

# Remove raw data to free up memory
rm(test_2)

# ----------------------------------------------------
# Define Extended Flexible GWAS function
# ----------------------------------------------------
ExtFlexGwas = function(y, gen, fam, rescale = TRUE) {
  tag = 0
  require(lme4)
  cat('Running marker: ')
  sma = apply(gen, 2, function(x, y, fam, rescale) {
    # Print current marker being processed
    tag <<- tag + 1
    cat(paste(tag, '. ', sep=''))
    
    # Build data frame for modeling
    TmpDta = data.frame(y = y, f = factor(fam), x = x)
    lvl = levels(TmpDta$f)
    
    # Drop levels with missing data
    TmpDta = droplevels.data.frame(TmpDta[rowMeans(!is.na(TmpDta)) == 1,])
    
    # Optionally rescale phenotype and genotype
    if(rescale) {
      TmpDta$y = scale(TmpDta$y)
      TmpDta$x = scale(TmpDta$x)
    }

    Vy = c(var(TmpDta$y))  # Total variance

    # Model 0: Null model with family as random effect
    fit0 = lmer(y ~ (1|f), TmpDta)
    ll0 = logLik(fit0)

    # Model 1: Within-family marker effect (random slope)
    fit1 = suppressMessages(lmer(y ~ (1|f) + (1|f):x, TmpDta))
    eff1 = ranef(fit1)$f[,2]
    names(eff1) = rownames(ranef(fit1)$f)
    eff1 = eff1[lvl]
    ll1 = logLik(fit1)
    LRT1 = ll1 - ll0
    PVAL1 = -log10(1 - pchisq(LRT1, 1))

    # Model 2: Across-family fixed effect
    fit2 = suppressMessages(lmer(y ~ x + (1|f), TmpDta))
    eff2 = fit2@beta[2]
    ll2 = logLik(fit2)
    LRT2 = ll2 - ll0
    PVAL2 = -log10(1 - pchisq(LRT2, 1))

    # R-squared values for each model
    R2 = 1 - c(fit0@devcomp$cmp['sigmaREML'],
               fit1@devcomp$cmp['sigmaREML'],
               fit2@devcomp$cmp['sigmaREML']) / Vy
    names(R2) = paste0('R2.model', 0:2)

    # Output vector
    NumObs = nrow(TmpDta)
    out = c(NumObs, P1 = PVAL1, P2 = PVAL2, R2, Fxd = eff2, eff1)
    return(out)
  }, fam = fam, y = y, rescale = rescale)

  cat('\n')
  
  # Set output row names
  rownames(sma) = c('NumObs', 'PVAL1', 'PVAL2', 'NHR2', 'WFR2', 'AFR2', 'EffOverall',
                    paste('Eff', sort(unique(fam)), sep=''))
  sma[is.na(sma)] = 0
  return(data.frame(t(sma)))
}

# Debug: test function with a few markers
if(F){
  test = ExtFlexGwas(y, gen[,1:5], fam)
}

# Create a reduced subset of markers for quicker tests
x = dim(gen)[2]
sml <- seq(1, x, 2)
sml_gen <- subset(gen, select = sml)
sml_y <- y[sml]
sml_fam <- fam[sml]
save(sml_gen, sml_y, sml_fam, file = "sml_CHH_gen.Rdata")

# Remove columns with all missing data before analysis
ex <- -1 * which(sapply(as.data.frame(gen), function(x) all(is.na(x))))
gen <- data.matrix(gen[, ex])

# Perform GWAS using the full genotype matrix
Fit = ExtFlexGwas(y, gen, fam, rescale = TRUE)

# Save session image after GWAS
save.image(file = "GWAS_CHG_geneovlps-ExtFlexGWAS_rsclT.RData")

# Filter for markers with significant associations (p ≥ 10^-3)
signif <- Fit[which(Fit$PVAL1 >= 3 | Fit$PVAL2 >= 3),]

# Save results to file
write.csv(signif, "CHG-INgeneONLY-ExtFlexGWAS_signif_scores.csv")
write.csv(Fit, "CHG-INgeneONLY-ExtFlexGWAS_scores.csv")

