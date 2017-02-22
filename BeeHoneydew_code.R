#####################################################################
######## Bees Without Flowers: Honeydew Use Manuscript Code #########
#####################################################################
############# By Joan Meiners, Last updated Sept 2016 ###############
#####################################################################

###########################################################################################################
## This is code for the Honeydew model version used in the Meiners 2015 thesis and the pre-print versions of the manuscript, as well as for the branch temperature model and some of the Figures used for publication with The American Naturalist. For code on the final model and model selection procedure used in the AmNat publication, as well as additional figures, please see either the Supplementary Materials associated with the publication or the mixed-models.R file in Beecycles github repository (https://github.com/beecycles/Bees-without-flowers-project2016). Treatment results did not vary between models. Supplementary files just provide additional detail on model seletion and refine the model slightly. ##
############################################################################################################

## Import dataset and Rename Chamise_summ
Chamise_summ <- read.csv("Meiners_BeeHoneydew_data.csv", header = TRUE)
#View(Chamise_summ)
dim(Chamise_summ)

## Make sure Site and Treatment values are categorical
Chamise_summ$Plant_Code = as.factor(Chamise_summ$Plant_Code)
Chamise_summ$Site = as.factor(Chamise_summ$Site)
Chamise_summ$Treatment_Code = as.factor(Chamise_summ$Treatment_Code)


########### Main Model: GLMM ############ 
###### NOTE: (This is the model used in M.S. Thesis, Meiners 2016, Utah State University)
###### (For the model used in the publication please see the mixed-models.R file in beecycles github repository.(Treatment results do not vary between models.))
library(lme4)
library(car)

hist(Chamise_summ$Bee_Count)   ### Determine negative binomial distribution of response variable

Honeydew <- glmer.nb(Bee_Count ~ Mold * Insecticide + Sugar * Paint + (1|Plant_Code:Site) + (1|min_day:julDate), data=Chamise_summ) 

summary(Honeydew)
anova(Honeydew)

Anova(Honeydew, type=3, test.statistic="F")

########### Branch Temperature Model: GLMM ############
library(lme4)
library(car)

BranchTemps <- read.csv("Chamise_BranchTemps.csv", header = TRUE)
dim(BranchTemps)
BranchTemps = na.omit(BranchTemps)

## Make sure Site and Treatment values are categorical
BranchTemps$Plant_Code = as.factor(BranchTemps$Plant_Code)
BranchTemps$Site = as.factor(BranchTemps$Site)
BranchTemps$Treatment_Code = as.factor(BranchTemps$Treatment_Code)
BranchTemps$Black = as.factor(BranchTemps$Black)

hist(BranchTemps$BranchTempC)  ### Determine normal distribution of response variable
qqnorm(BranchTemps$BranchTempC)  ### Distribution looks pretty normal

#### Test of difference in branch temps between blackened and not blackened branches
# Old Model: BranchTemp <- lmer(BranchTempC ~ Black * julDate + (1|Plant_Code:Site), data=BranchTemps)
BranchTemp <- lmer(BranchTempC ~ Black * julDate + Site + (1|Plant_Code), data=BranchTemps)

summary(BranchTemp)
anova(BranchTemp)

Anova(BranchTemp, type=3, test.statistic="F")


#######################################################
################ Code for Four Figures ################
#######################################################


####### (Figure 1 is a colletion of photographs) #######


############### Figure 2: Pilot Study ###############
quartz(width = 7, height = 7)
plot.new()
par(mar=c(18,4,2,2)+0.1)
bee_count = c(124, 30, 0, 6)
barplot(bee_count, names.arg = c("Moldy", "Mold-free", "Moldy", "Mold-free"), col = c("grey15", "grey75", "grey15", "grey75"),  las=1)
axis(side=1, tick= FALSE, at=c(1.3, 2.5, 3.7), labels= c("Pre-bloom", "|", "Flowering"), line=2, font=1)
mtext("Total pilot study bee count", side=2, line=2.5, font=1)

tiff(filename = "ChamisePilot3.tiff", units = "in", compression = "lzw", res = 300, width = 7, height = 7)
dev.off()

#### Former Figure 2 plot (with pictures as background) ####
#library(ggplot2)
#library(jpeg)
#library(grid)

#plot.new()
#quartz(width = 8, height = 5)
#Adeno_nonfl = readJPEG("/Users/joanmeiners/Dropbox/Chamise_manuscript/Pilot_picture_Shebs.jpg")
#pilot.df = data.frame(plant = factor(c("Moldy ", "Mold-free ", "Moldy", "Mold-free"), levels = c("Moldy ", "Mold-free ", "Moldy", "Mold-free")), bee_count = c(124, 30, 0, 6))
#pilot_1 = ggplot(data = pilot.df, aes(x=plant, y=bee_count)) + geom_bar(stat="identity") +
#  annotation_custom(rasterGrob(Adeno_nonfl, width=unit(1,"npc"), height=unit(1,"npc")), -Inf, Inf, -Inf, Inf) +
#  scale_y_continuous(expand=c(0,0), limits = c(0,max(pilot.df$bee_count)*1.05)) +
#  geom_bar(stat="identity", colour=c("white", "black", "white", "black"), fill=c("black", "white", "black", "white"), width=0.5, size = 1) + 
#  labs(x="Pre-bloom                 |                  Flowering", y="Total pilot study bee count", font = 2, las = 1, cex = 1.8) +
#  theme(axis.title.y = element_text(size = rel(1.5), colour = "black", face = "bold", angle = 90))+
#  theme(axis.title.x = element_text(size = rel(1.5), colour = "black", face = "bold", angle = 00 ,vjust = 2))+
#  theme(text=element_text(family=NULL),
#        panel.border = element_rect(colour = "black", fill = NA, size = 2),
#        axis.text.y=element_text(size = rel(1.5), colour = "black", face = "bold", vjust = 1.5),
#        axis.text.x=element_text(size = rel(1.5), colour = "black", face = "bold", vjust = 10))
#pilot_1

#tiff(filename = "ChamisePilot.tiff", units = "in", compression = "lzw", res = 300, width = 8, height = 5)
#dev.off()


####### Previous Figure 3: Effect size vs. Total Bee Count (Removed from AmNat version for publication)#######
#### Cohens_D calculated for use with Master's thesis version of project, Dec. 2015 ########
#Cohens_D = c(0.508291, 0.040563, 0.222611, 0.175647, 1.032771, 0.662503)
#sorted_D = sort(Cohens_D, decreasing = FALSE)

plot.new()
quartz(width = 10, height = 6)
par(mar=c(5,5,2,1))
par(oma=c(0,0,0,0) )
plot(sorted_D, c(12, 7, 17, 41, 101, 119), type = "p", axes = FALSE, col = c("gray", "gray", "gray", "orange", "orange", "red"), pch = 19, cex = 3, xlab = "Cohen's d effect size compared to control", ylab = "Total bee count", cex.lab = 2, xlim = c(-0.07, 1.05), ylim = c(-15, 135))
axis(2, at=seq(0 , 120, by=40), cex.axis = 1.5)
axis(1, c(0, 0.2, 0.5, 0.8, 1.0), cex.axis = 1.5)
text(c(0.2, 0.5, 0.8), c(-14, -14, -15), c("Small", "Medium", "Large"), 
     col = c("gray", "orange", "red"), cex = 1.6, font = 2)
text(sorted_D, c(21, -1, 29, 50, 110, 128), c("Mold + Insecticide", "Paint", "Insecticide", "Mold", "Sugar + Paint", "Sugar"), cex = 1.4, font = 2, col = "blue")

tiff(filename = "ChamiseEffectSize.tiff", units = "in", compression = "lzw", res = 200, width = 10, height = 6)
dev.off()


################# Figure 3 (previously Figure 4 in pre-print): Interaction plot with CIs #################

## Color version
library(sciplot)

plot.new()
quartz(width = 10, height = 6)
par(mfrow=c(1,2))
par(mfrow=c(1,2), oma = c(0,0,0,0), mar=c(6,6,2,2)+0.1)
lineplot.CI(as.factor(Mold), Bee_Count, group = as.factor(Insecticide), data = Mold_trt, type = c("b"), x.cont=FALSE, legend = FALSE, xlab = "Plant natural condition", ylab = "Mean bee abundance", cex.lab = 1.7, ylim=c(0,3), font=2, lty = c(2,1), pch = c(16, 19), col = c("slateblue", "red"), cex = 1.7, lwd = 3, xaxt = 'n', axes = FALSE) 

axis(side = 1, at = c(1,2), labels = c("Mold-free", "Moldy"), cex.axis = 1.5)
axis(2, at=seq(0 , 3, by=1), font = 2, las=1)

legend("top", legend = c("No insecticide", "Insecticide"), col = c("slateblue", "red"), lwd = 3, pch = 19, lty = c(2,1), title = "", cex = 1.7, bty = "n")
#~~
lineplot.CI(as.factor(Sugar), Bee_Count, group = as.factor(Paint), data = Sugar_trt, type = c("b"),  x.cont=FALSE, legend = FALSE, xlab = "On mold-free plants", cex.lab = 1.7, ylab = "", ylim=c(0,3), lty = c(2,1), pch = c(16, 19), col = c("chartreuse4", "black"), cex = 1.7, lwd = 3, xaxt = 'n', yaxt = 'n', axes = FALSE) 

legend("top", legend = c("No paint", "Paint"), col = c("chartreuse4", "black"), lwd = 3, pch = 19, lty = c(2,1), title = "", cex = 1.7, bty = "n")

axis(side = 1, at = c(1,2), labels = c("No sugar spray", "Sugar spray"), cex.axis = 1.5)
axis(2, at=seq(0 , 3, by=1), font = 2, las=1, labels = NULL)


## Grey scale version of Figure 3
library(sciplot)

quartz(width = 10, height = 6)
plot.new()
par(mfrow=c(1,2))
par(mfrow=c(1,2), oma = c(0,0,0,0), mar=c(6,6,2,2)+0.1)
lineplot.CI(as.factor(Mold), Bee_Count, group = as.factor(Insecticide), data = Mold_trt, type = c("b"), x.cont=FALSE, legend = FALSE, xlab = "Plant natural condition", ylab = "Mean sample bee abundance", cex.lab = 1.6, ylim=c(0,1.5), font=2, lty = c(2,1), pch = c(16, 19), col = c("darkgrey", "black"), cex = 1.6, lwd = 3, xaxt = 'n', axes = FALSE) 

axis(side = 1, at = c(1,2), labels = c("Mold-free", "Moldy"), cex.axis = 1.5)
axis(2, at=seq(0 , 1.5, by=0.5), font = 2, las=1)

legend("top", legend = c("No insecticide", "Insecticide"), col = c("darkgrey", "black"), lwd = 3, pch = 19, lty = c(2,1), title = "", cex = 1.6, bty = "n")
#~~
lineplot.CI(as.factor(Sugar), Bee_Count, group = as.factor(Paint), data = Sugar_trt, type = c("b"),  x.cont=FALSE, legend = FALSE, xlab = "On mold-free plants", cex.lab = 1.6, ylab = "", ylim=c(0,3), lty = c(2,1), pch = c(16, 19), col = c("darkgrey", "black"), cex = 1.6, lwd = 3, xaxt = 'n', yaxt = 'n', axes = FALSE) 

legend("top", legend = c("No paint", "Paint"), col = c("darkgrey", "black"), lwd = 3, pch = 19, lty = c(2,1), title = "", cex = 1.6, bty = "n")

axis(side = 1, at = c(1,2), labels = c("No sugar spray", "Sugar spray"), cex.axis = 1.5)
axis(2, at=seq(0 , 3, by=1), font = 2, las=1, labels = NULL)


tiff(filename = "ChamiseInteractionGrey1.tiff" , units = "in", compression = "lzw", res = 300, width = 10, height = 6)
dev.off()


####### Appendix Figure 1 Branch Temp Diffs ########
plot.new()
quartz(width = 7, height = 5)
boxplot(BranchTemps$BranchTempC ~ BranchTemps$Black, xlab = "Branches", ylab = "Average branch temperature (C)", cex.lab = 1.3, axes = FALSE)
axis(side = 1, at = c(1,2), labels = c("Not blackened", "Blackened"), cex.axis = 1.2)
axis(2, at=seq(15, 35, by=5))

tiff(filename = "ChamiseBranchTemps.tiff", units = "in", compression = "lzw", res = 200, width = 7, height = 5)
dev.off()


##### Appendix Figure 2: Response Variable Distribution ######
plot.new()
quartz(width = 8, height = 6)
par(mfrow = c(1,1))
hist(Chamise_summ$Bee_Count, breaks = 22, main = "", 
     xlab = "Bee abundance per plant sample", ylab = "Frequency of samples")

tiff(filename = "Chamise_hist_BeeAbund.tiff" , units = "in", compression = "lzw", res = 150, width = 8, height = 6)
dev.off()


###########################################################################################################
## This is code for the Honeydew model version used in the Meiners 2015 thesis and the pre-print versions of the manuscript, as well as for the branch temperature model and some of the Figures used for publication with The American Naturalist. For code on the final model and model selection procedure used in the AmNat publication, as well as additional figures, please see either the Supplementary Materials associated with the publication or the mixed-models.R file in Beecycles github repository (https://github.com/beecycles/Bees-without-flowers-project2016). Treatment results did not vary between models. Supplementary files just provide additional detail on model seletion and refine the model slightly. ##
############################################################################################################