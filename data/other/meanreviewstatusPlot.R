
library(ggplot2)
lightgray <- "#cccccc"; gray <- "#999999"; orange <- "#E69F00"; skyblue <- "#56B4E9"; blueishgreen <- "#009E73"
yellow <- "#F0E442"; blue <-"#0072B2"; vermillion <- "#D55E00"; reddishpurple <- "#CC79A7"

#data based on: gavin r0.3 calibrations, and ClinVar of 1/12/2016 "variant_summary.txt"
#Gene	Category	ReviewStatusMean	NumberSubmittersMean	SubmitterCategoriesMean
## C1 = CADD scores highly significantly predictive for pathogenicity (pval < 0.01).
## C2 = CADD scores significantly predictive for pathogenicity (pval < 0.05).
## C3 = CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).
## C4 = CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).
pathToGavinGitRepo <- "/Users/joeri/github/gavin"
mrd <- read.table(paste(pathToGavinGitRepo,"/data/other/meanreviewstatus.tsv",sep=""), sep="\t", header=T)
mrdSub <- subset(mrd, Category == 'C1' | Category == 'C2' | Category == 'C3' | Category == 'C4')

plot(mrd$Pvalue ~ mrd$NrOfVariants)
lmfit <- lm(ReviewStatusMean ~ NrOfVariants, data=mrd)

ggplot() +
  geom_point(data = mrdSub, aes(x = NrOfVariants, y = ReviewStatusMean, color = Pvalue), alpha=0.5) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_abline(intercept = lmfit$coefficients[1], slope = lmfit$coefficients[2], color="blue") +
  #scale_colour_manual(values=rainbow(16)) +
  ylab("Gene means of ClinVar review status") + xlab("Number of ClinVar variants for gene") +
  scale_x_continuous(lim=c(0,250)) +
  scale_y_continuous(lim=c(0,1.5)) +
  scale_colour_gradient(low="red", high="green", limits=c(0, 0.1))
#ggsave("NrVariantsVsReviewStatus.pdf", width = 12, height = 12, units = "cm")

#showing: more variants means better p-value for calibration
lmfit <- lm(Pvalue ~ NrOfVariants, data=mrdSub)
ggplot() +
  geom_point(data = mrdSub, aes(x = NrOfVariants, y = Pvalue), alpha=0.5) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_abline(intercept = lmfit$coefficients[1], slope = lmfit$coefficients[2], color="red", size=1) +
  ylab("GAVIN gene calibration p-value") +
  xlab("Number of ClinVar variants for gene") +
  scale_x_continuous(lim=c(0,2200)) +
  scale_y_continuous(lim=c(0,1))
ggsave("NrVariantsVsPvalue.pdf", width = 12, height = 12, units = "cm")

#showing: better review quality means better p-value for calibration
lmfit <- lm(Pvalue ~ ReviewStatusMean, data=mrdSub)
ggplot() +
  geom_point(data = mrdSub, aes(x = ReviewStatusMean, y = Pvalue), alpha=0.2) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  geom_abline(intercept = lmfit$coefficients[1], slope = lmfit$coefficients[2], color="red", size=1) +
  scale_x_continuous(lim=c(0,3)) +
  ylab("GAVIN gene calibration p-value") +
  xlab("Mean of variant ClinVar review status")
ggsave("ReviewStatusMeanVsPvalue.pdf", width = 12, height = 12, units = "cm")


ggplot() +
  geom_boxplot(data = mrdSub, aes(x = Category, y = ReviewStatusMean), fill=orange, size=1.1) +
  theme_bw() + theme(panel.grid.major = element_line(colour = "black"), axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) +
  ylab("Gene means of ClinVar review status") + xlab("GAVIN gene calibration category")
ggsave("CalibVsReviewStatus.pdf", width = 12, height = 12, units = "cm")



### 3D plots


#install.packages("scatterplot3d")
library(scatterplot3d)
s3d <- scatterplot3d(mrdSub$ReviewStatusMean, mrdSub$Pvalue, mrdSub$NrOfVariants,
              highlight.3d = TRUE, col.axis = "blue",
              col.grid = "lightblue", main = "Title", pch = 20, zlim = c(0,250))
s3d$plane3d(my.lm <- lm(mrdSub$Pvalue ~ mrdSub$ReviewStatusMean + mrdSub$NrOfVariants), lty = "dotted")

scatterplot3d(mrdSub$Pvalue, mrdSub$NrOfVariants, mrdSub$ReviewStatusMean,
              highlight.3d = TRUE, col.axis = "blue",
              col.grid = "lightblue", main = "Title", pch = 20, ylim = c(0,250))
scatterplot3d(mrdSub$NrOfVariants, mrdSub$ReviewStatusMean, mrdSub$Pvalue,
              highlight.3d = TRUE, col.axis = "blue",
              col.grid = "lightblue", main = "Title", pch = 20, xlim = c(0,250))




#install.packages("rgl")
library(rgl)

plot3d(mrdSub$ReviewStatusMean, mrdSub$Pvalue, mrdSub$NrOfVariants, zlim=c(0,250), xlim=c(0,1))
plot3d(mrdSub$NrOfVariants, mrdSub$ReviewStatusMean, mrdSub$Pvalue, xlim=c(0,250), ylim=c(0,1))


fit1 <- lm(mrdSub$Pvalue ~ mrdSub$NrOfVariants)
fit2 <- lm(mrdSub$Pvalue ~ mrdSub$ReviewStatusMean)
plot(mrdSub$Pvalue ~ mrdSub$NrOfVariants)
abline(fit1)
plot(mrdSub$Pvalue ~ mrdSub$ReviewStatusMean)
abline(fit2)


