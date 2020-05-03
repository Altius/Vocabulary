source("../../../code/general.R")
library(gplots)

### Load sample & NMF metadata
load("../../WM20180608_733_sample_metadata/metadata.RData") # metadata
comp_ord <- rev(comp_ord) # reverse for plotting prettiness
ncomp <- length(comp_ord)

labs <- c("Tissue invariant", "Stromal A", "Primitive / embryonic", "Stromal B",
          "Lymphoid", "Renal / cancer", "Placental / trophoblast", "Neural",
          "Cardiac", "Organ devel. / renal", "Pulmonary devel.", "Musculoskeletal",
          "Digestive", "Vascular / endothelial", "Myeloid / erythroid", "Cancer / epithelial")
#comp_ord <- c(7,3,10,11,5,15,9,12,14,8,13,6,16,2,4,1)

nam_dat <- read.delim("../../WM20180608_construct_masterlist_733samples/sampnams_733.txt", header=FALSE, as.is=TRUE)
nam_dat <- nam_dat[metadata$system != "Hematopoietic" | (metadata$system == "Hematopoietic" & metadata$Unique.cellular.condition != 1),] ### Select subset
nams <- paste(nam_dat[,3], nam_dat[,2], sep=".")
nams_short <- nam_dat[,3]

outdir <- "."

###################################################################################################################

#dat_raw <- as.matrix(read.delim("data/2020-04-15NC16_NNDSVD_hemoUniqOnly_Basis.tsv.gz", sep=",", as.is=TRUE, row.names=1))
#dat_raw <- as.matrix(read.delim("data/2020-04-15NC16_NNDSVD_uniqOnly_Basis.tsv.gz", sep=",", as.is=TRUE, row.names=1))
dat_raw <- as.matrix(read.delim("data/2020-04-15NC16_NNDSVD_hemoNonUniqOnly_Basis.tsv.gz", sep=",", as.is=TRUE, row.names=1))
dat_raw <- dat_raw[,c(0,  1,  5,  3,  4,  15,  6,  7,  8, 10,  9, 11, 12, 13,  2, 14)+1]

## Plot (maximally) a single winning component per sample
dat_single <- t(apply(dat_raw, 1, function(x) { x[x != max(x)] <- 0; x }))

plotfile(paste(outdir, "/barplot_comp_ALLincl1_comp_ordered_pruned_raw_single_top15", sep=""), type="pdf", height=4*4, width=4*7)
par(mfrow=c(4,4), oma=c(3,6,0,0))
for (i in rev(comp_ord)) {
  ord <- order(dat_single[,i], decreasing=FALSE)
  ord <- tail(ord, 15)
  #if (i == 1) next; # skip component 1 (sink)
  par(mar=c(2,7,2,2))
  barplot(dat_raw[ord,i], col=cols[i], border=NA, las=2, horiz=T, cex.names=min(1, 10/length(ord)),
          names.arg=gsub(".DS.*$", "", nams[ord]), cex.axis=0.75, xaxs="i", yaxs="i", xlab="",
          xlim=c(0, max(dat_raw)), xpd=FALSE, main=labs[i])
  #abline(v=p95_val, lwd=2, col="black", lty=2)
}      
mtext("NMF loadings", side=1, outer=TRUE, line=1)
dev.off()

