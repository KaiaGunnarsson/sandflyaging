setwd("/Users/jamescotton/Work_newlaptop/Glasgow_Teaching/L4_project_KaiaGunnarson/MIRS_analysis")
data <- read.delim("kaia_notrim_matrix.txt.dat")

library(tidyverse)

#replace intensity values with 'normalised' ones - subtract rowmean and divide by sd.
data.rowmeans <- rowMeans(data[,8:1807])
data.rowSD <-  apply(data[,8:1807],1,sd)
new_data_bits <- (data[,8:1807] - data.rowmeans) / data.rowSD
new_data <- cbind(data[,1:7],new_data_bits)
data <- new_data
#getting the data into some sensible formats, and adding columns with species, dpe etc for each group
#
data.MIRS <- data %>% pivot_longer(cols=starts_with("X"),names_to = "wavelength",values_to="intensity")
data.MIRS <- data.MIRS %>% mutate(wavelength=as.numeric(str_remove(wavelength,"X")))


data.MIRS$species <- "control"
data.MIRS$species[str_detect(data.MIRS$Cat1,"RAIN")] <- "sergenti"
data.MIRS$species[str_detect(data.MIRS$Cat1,"PURP")] <- "papatasi"
data.MIRS$group <- str_extract(data.MIRS$Cat2,"G([0-9])",group=TRUE)
#make sure group is treated as a categorical factor
data.MIRS$group <- as.factor(data.MIRS$group)

#make 'days post eclosion' variable
data.MIRS$dpe <- NA
data.MIRS$dpe[data.MIRS$group == "1"] <- NA
data.MIRS$dpe[data.MIRS$group == "2"] <- 1
data.MIRS$dpe[data.MIRS$group %in% c("3","4","8")] <- 3
data.MIRS$dpe[data.MIRS$group == "5"] <- 10
data.MIRS$dpe[data.MIRS$group %in% c("6","7")] <- 13
data.MIRS$dpe <- as.numeric(data.MIRS$dpe)

#sex variable
data.MIRS$sex <- NA
data.MIRS$sex[data.MIRS$group %in% c("2","3","4","5","6","7")] <- "female"
data.MIRS$sex[data.MIRS$group == "8"] <- "male"
data.MIRS$sex[data.MIRS$group == "1"] <- "larval"

#blood or no blood
data.MIRS$blood <- NA
data.MIRS$blood[data.MIRS$group %in% c("1","2","4","7","8")] <- FALSE
data.MIRS$blood[data.MIRS$group %in% c("3","5","6")] <- TRUE
table(data.MIRS$blood,data.MIRS$group)
data.MIRS$blood <- as.factor(data.MIRS$blood)

#other factors
data.MIRS$species <- as.factor(data.MIRS$species)
data.MIRS$sex <- as.factor(data.MIRS$sex)

#remove columns I don't need now:
data.MIRS <- data.MIRS %>% select(!(Cat1:Cat6))

#-------
#example plots of some spectra
#get average intensity across replicates for each group
data.MIRS.replicate_average <- data.MIRS %>% group_by(species,group,wavelength,dpe,sex,blood) %>% summarize(mean_intensity=mean(intensity))
#plot these
blood.labs <- c("blood-fed", "no blood")
names(blood.labs) <- c("TRUE","FALSE")
ggplot(data.MIRS.replicate_average,aes(y=mean_intensity,x=wavelength,colour=as.factor(dpe))) + geom_line() + facet_grid(species ~ blood+sex, labeller=labeller(blood=blood.labs)) + 
theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(color="dpe")

#---------
#principal components analysis
#need to start from original matrix
data.matrix <- pivot_wider(data.MIRS,names_from = "wavelength",values_from="intensity")
data.matrix.2 <- data.matrix %>% select(!(Cat7:blood))
pca <- prcomp(data.matrix.2,center=TRUE,scale=TRUE)
pca_points <- data.matrix %>% select(species:blood) %>% bind_cols(as_tibble(pca$x))

library(ggnewscale)
pca_var_frac <- pca$sdev^2 / sum(pca$sdev^2)
#amount of variation explained by PC1 and PC2
ggplot(pca_points,aes(x=PC1,y=PC2,color=dpe,shape=blood))+ geom_point() +
  scale_colour_distiller(palette="Spectral") +
  new_scale("shape") + geom_point(aes(x=PC1,y=PC2,shape=sex,size=3)) + 
  scale_shape(solid=FALSE) + facet_grid(.~ species) + xlab(paste0("PC1 (",round(pca_var_frac[1],4)*100,"%)")) + ylab(paste0("PC2 (",round(pca_var_frac[2],4)*100,"%)")) + guides(size = "none") + theme_bw()
#remove '3' from legend??

#fit a model to PCs
m1.PC1 <- lm(PC1 ~ species*dpe+sex+blood, data=pca_points)
m1.PC2 <- lm(PC2 ~ species*dpe+sex+blood, data=pca_points)
summary(m1.PC1)
summary(m1.PC2)


#------
#fit a model for each wavelength
library(broom)
#this code does a statistical test for the effect of dpe on intensity for 
#every wavelength in the data, then returns a p-value for each wavelength

model.fitter <- function(d) {
  m1 <- lm(intensity ~ species*dpe+sex+blood,data=d)
  t <- tidy(m1)
  as.numeric(t[t$term=="dpe","p.value"])
}

model.fitter2 <- function(d) {
  m1 <- lm(intensity ~ species*dpe+sex+blood,data=d)
  t <- tidy(m1)
  list(p=as.numeric(t[t$term=="dpe","p.value"]),pap_slope=as.numeric(t[t$term=="dpe","estimate"]),
       serg_slope=as.numeric(t[t$term=="dpe","estimate"])+as.numeric(t[t$term=="speciessergenti:dpe","estimate"]))
}


data.MIRS_pvals <- data.MIRS %>% group_by(wavelength) %>% nest() %>% mutate(dpe_pval = map_dbl(data, model.fitter)) %>% mutate(neglogp=-1*log10(dpe_pval))

signif.level <- 0.05 / nrow(data.MIRS)
ggplot(data.MIRS_pvals,aes(x=wavelength,y=neglogp)) + geom_point(size=1) + theme_bw() + ylab("-ve log(p) for days post eclosion") + geom_hline(yintercept =-1*log10(signif.level),lty=2)
ggplot(data.MIRS_pvals,aes(x=wavelength,y=neglogp,color=neglogp < -1*log10(signif.level))) + scale_color_manual(values=c("black","lightgrey")) + geom_point(size=1) + theme_bw() + ylab("-ve log(p) for days post eclosion") + geom_hline(yintercept =-1*log10(signif.level),lty=2)

data.MIRS_slopes <- data.MIRS %>% group_by(wavelength) %>% nest() %>% mutate(dpe_model = map(data, model.fitter2)) %>% hoist(dpe_model,"p","pap_slope","serg_slope")
ggplot(data.MIRS_slopes,aes(x=wavelength,y=pap_slope))  + geom_point(size=1) + theme_bw() + ylab("slope with dpe") 
data.MIRS_slopes2 <- data.MIRS_slopes %>% pivot_longer(cols=c(pap_slope,serg_slope),names_to="spp",values_to="slope")
ggplot(data.MIRS_slopes2,aes(x=wavelength,y=slope,color=spp))  + geom_point(size=1) + theme_bw() + ylab("slope with dpe") 

g1 <- ggplot(data.MIRS_pvals,aes(x=wavelength,y=neglogp,color=neglogp < -1*log10(signif.level))) + scale_color_manual(name="significant",values=c("black","lightgrey")) + geom_point(size=1) + theme_bw() + ylab("-ve log(p) for days post eclosion") + geom_hline(yintercept =-1*log10(signif.level),lty=2)
g2 <- ggplot(data.MIRS_slopes2,aes(x=wavelength,y=slope,color=spp))  + geom_point(size=1) + theme_bw() + ylab("slope with dpe") + scale_color_discrete(labels=c("papatasi","segenti"))
g1 / g2 + plot_annotation(tag_levels = 'A')
