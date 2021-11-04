library(greedyBAPs)
library(ngBap)
plot.dat <- read.csv("../data/grace_plot_level.csv")
site.dat <- read.csv("../data/grace_site_level.csv")



plot.sem.dat <- with(plot.dat, data.frame(psitecode))
names(plot.sem.dat)[names(plot.sem.dat)=="psitecode"] <- "PlotSiteCode"
plot.sem.dat$PlotRichRaw    <- plot.dat$prich
plot.sem.dat$PlotRich       <- plot.dat$ln.prich
plot.sem.dat$SiteRich       <- plot.dat$ln.site.rich
plot.sem.dat$PlotProdRaw    <- plot.dat$pprod
plot.sem.dat$PlotProd       <- plot.dat$ln.pprod
plot.sem.dat$SiteProd       <- plot.dat$ln.site.prod
plot.sem.dat$PlotBiomass    <- plot.dat$ln.ptotmass
plot.sem.dat$SiteBiomass    <- plot.dat$ln.site.totmass
plot.sem.dat$PlotShade      <- plot.dat$ln.pshade
plot.sem.dat$PlotSoilSuit   <- plot.dat$SoilSuitability
plot.sem.dat$SoilWithShade  <- plot.dat$SoilWithShade #control variable for soil influences on effectiveness of shading
plot.sem.dat$PlotN          <- plot.dat$pn
plot.sem.dat$PlotN.sqr      <- plot.dat$sqr.pn
plot.sem.dat$PlotN.ln       <- plot.dat$ln.pn
plot.sem.dat$PlotP.ln       <- plot.dat$ln.pp
plot.sem.dat$PlotC.ln       <- plot.dat$ln.pc

Y <- cbind(plot.sem.dat[,c("PlotRich", "PlotSoilSuit", "SiteRich", "PlotShade", "PlotBiomass",
                      "PlotProd","SiteBiomass", "SiteProd")])#, SiteSoilSuit = site.dat$SoilSuitability[match(plot.sem.dat$PlotSiteCode, site.dat$site.code)])


Y <- scale(Y)
out_01 <- ngBap::bang(Y, 4, level = .01, verbose = F)


grace_modB <- matrix(0, nrow = dim(Y)[2], ncol = dim(Y)[2])
grace_modOm <- matrix(0, nrow = dim(Y)[2], ncol = dim(Y)[2])
colnames(grace_modOm) <- rownames(grace_modOm) <-
  colnames(grace_modB) <- rownames(grace_modB) <- colnames(Y)


grace_modB[matrix(c("SiteProd", "PlotProd",
                    "SiteProd", "PlotBiomass",
                    "SiteProd", "PlotRich",
                    "SiteBiomass", "PlotRich",
                    "SiteBiomass", "PlotBiomass",
                    "SiteBiomass", "PlotProd",
                    "SiteRich", "PlotRich",
                    "SiteRich", "PlotBiomass",
                    "SiteRich", "PlotProd",
                    "PlotBiomass", "PlotShade",
                    "PlotShade", "PlotRich",
                    "PlotSoilSuit", "PlotRich",
                    "PlotRich", "PlotProd"), byrow = T, ncol = 2)] <- 1
grace_modB <- t(grace_modB)


grace_modOm[matrix(c("SiteProd", "SiteRich",
                     "SiteRich", "SiteBiomass",
                     "SiteProd", "SiteBiomass",
                     "SiteRich", "PlotSoilSuit"), byrow = T, ncol = 2)] <- 1
grace_modOm <- grace_modOm + t(grace_modOm) + diag(rep(1, dim(Y)[2]))

colnames(out_01$dEdge) <- rownames(out_01$dEdge) <- colnames(grace_modB)
colnames(out_01$bEdge) <- rownames(out_01$bEdge) <- colnames(grace_modB)

sum(grace_modB == 1 & out_01$dEdge == 1)
(sum(grace_modOm == 1 & out_01$bEdge == 1) - dim(Y)[2]) / 2
sum((grace_modB + t(grace_modB) + grace_modOm) == 0 & (out_01$dEdge + t(out_01$dEdge) +  out_01$bEdge) == 0) / 2

sum(dbinom(16:28, 28, p = .25))






#### Greedy BAP Comparison ###
set.seed(123)
sim.size <- 500
rec <- matrix(0,nrow =  sim.size, ncol = 2)
for(i in 1:sim.size){
  greedy_output <- greedyBAPs::greedySearch(cov(Y),n = dim(Y)[1], n.restarts = 1)
  greedyB <- t(ifelse(greedy_output$final.bap == 1, 1, 0))
  greedyOm <- ifelse(greedy_output$final.bap == 100, 1, 0) + diag(rep(1, dim(Y)[2]))
  rec[i] <-   sum(grace_modB == 1 & greedyB == 1) + (sum(grace_modOm == 1 & greedyOm == 1) - dim(Y)[2]) / 2 +sum((grace_modB + t(grace_modB) + grace_modOm) == 0 & (greedyB + t(greedyB) +  greedyOm) == 0) / 2
  cat(i)
  cat(": ")
  cat(rec[i])
  cat("\n")
}





# dat <- read.csv("graceOutput.csv")
#
# dat[which(abs(dat[,3] - max(dat[,3])) < 1e-6), 2]
# max(dat)
# postscript("graceGBS.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 4)
# par(mfrow = c(1,2), oma = c(0, 0, 0,0))
# plot(dat[, -1], xlab = "Correct Edges", ylab = "Penalized Log-Likelihood")
# hist(read.csv("graceOutput.csv")[,2], main = "", xlab = "Correct Edges")
# abline(v = 16, col = "red", lwd = 2)
# abline(v = 12, col = "blue", lwd = 2)
# mtext("Greedy BAP Search", outer = T, line = -2, cex = 2)
# dev.off()
#
