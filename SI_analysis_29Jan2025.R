setwd("/Users/jamescotton/Work_newlaptop/Glasgow_Teaching/L4_project_KaiaGunnarson")
data <- read.delim("JC241216NCS_StableIsotopeData.txt")
colnames(data) <- c("identifier","amount","blank1","d15N","d13C","d34S","blank2","Nug","Cug","percentN","percentC","CNratio")
#library for string matching and manipulation
library(stringr)
#create species column to distinguish sandfly samples from controls
data$control <- "control"
data$control[str_detect(data$identifier,"purple|rainbow")] <- "sample"
#create species column based on identifiers
data$species <- "control"
data$species[str_detect(data$identifier,"rainbow")] <- "sergenti"
data$species[str_detect(data$identifier,"purple")] <- "papatasi"
data$group <- str_extract(data$identifier,"_([0-9])_",group=TRUE)
#make sure group is treated as a categorical factor
data$group <- as.factor(data$group)
data$replicate <- str_extract(data$identifier,"_N([0-9])",group=TRUE)
#make sure replicate is treated as a categorical factor
data$replicate <- as.factor(data$replicate)
#remove some blank lines from data frame
data <- subset(data,data$identifier != "tests" & data$identifier != "blank")
#set numerical columns to be numeric
data$d15N <- as.numeric(data$d15N)
data$d13C <- as.numeric(data$d13C)
data$d34S <- as.numeric(data$d34S)
data$Nug <- as.numeric(data$Nug)
data$Cug <- as.numeric(data$Cug)
data$percentN <- as.numeric(data$percentN)
data$percentC <- as.numeric(data$percentC)
data$CNratio <- as.numeric(data$CNratio)
#still have blank columns in there but that shouldn't break anything 

#exploratory data analysis
library(ggplot2)
ggplot(data,aes(x=d15N,y=d13C,color=species)) + geom_point() + theme_bw()
dotplot_C <- ggplot(subset(data,control=="sample"), aes(x=group, y=d13C,color=species)) + geom_jitter(width=0.2) + theme_bw() + theme(legend.position = "none")
dotplot_N <- ggplot(subset(data,control=="sample"), aes(x=group, y=d15N,color=species)) + geom_jitter(width=0.2) + theme_bw() + theme(legend.position = "none")
dotplot_S <- ggplot(subset(data,control=="sample"), aes(x=group, y=d34S,color=species)) + geom_jitter(width=0.2) + theme_bw()
#example of the easiest way to export ggplot figures "nicely"
ggsave("C_per_group.pdf",dotplot_C)


library(patchwork)
dotplot_C + dotplot_N + dotplot_S + plot_annotation(tag_levels = 'A')

#make 'days post eclosion' variable
data$dpe <- NA
data$dpe[data$group == "1"] <- NA
data$dpe[data$group == "2"] <- 1
data$dpe[data$group %in% c("3","4","8")] <- 3
data$dpe[data$group == "5"] <- 10
data$dpe[data$group %in% c("6","7")] <- 13
data$dpe <- as.numeric(data$dpe)


#sex variable
data$sex <- NA
data$sex[data$group %in% c("2","3","4","5","6","7")] <- "female"
data$sex[data$group == "8"] <- "male"
data$sex[data$group == "1"] <- "larval"

#blood or no blood
data$blood <- NA
data$blood[data$group %in% c("1","2","4","7","8")] <- FALSE
data$blood[data$group %in% c("3","5","6")] <- TRUE
table(data$blood,data$group)
data$blood <- as.factor(data$blood)

#amount -> 
data$amount <- as.numeric(data$amount)

#other factors
data$species <- as.factor(data$species)
data$sex <- as.factor(data$sex)

sample.data <- subset(data,control == "sample")

#simple model
m1 <- lm(d13C ~ species,data=sample.data)
m2 <- lm(d13C ~ species+dpe+sex+blood,data=sample.data)
m3 <- lm(d13C ~ species+dpe+sex+blood+amount,data=sample.data)
m4 <- lm(d13C ~ species*dpe+sex+blood+amount,data=sample.data)

#in fact, I think we can drop 'amount' from the model..
m2.interaction <- lm(d13C ~ species*dpe+sex+blood,data=sample.data)

#to test the effect of 'blood'
#build the simpler model without blood included
m4.noblood <- lm(d13C ~ species*dpe+sex+amount,data=sample.data)
#this is the test based on F-statistics that "summary" is using:
anova(m4,m4.noblood)
#pvalue will match pv-alue for blood in summary(m4)
library(lmtest)
lrtest(m4,m4.noblood)

#in fact we'd only test this if we've already tested the interaction term and found it to be not significant
#so you probably don't want to do this.
#larval samples don't have a dpe (they are 'NA' or not applicable), so this model includes additional data..
m3.nodpe <- lm(d13C ~ species+sex+blood+amount,data=subset(sample.data,!is.na(sample.data$dpe)))
lrtest(m3,m3.nodpe)
# I don't think there is any reason to prefer the LRT for these models!
#


#is a standard 'normal' linear model appropriate or do we need something fancy like a GLM
#are the model assumptions violated?
#do you need to include interaction terms?

#model 'criticism' - checking assumptions
#plot residuals against fitted values --> look for 'a nice rectangle'
#plot QQplot or histogram of residuals --> looks roughtly normal.
# for this model we are happy with our normal GLM approach --
plot(m4)
#but not very pretty
#pretty plot of residuals.. like ch16 of Dan's book
#make a dataframe that includes the residuals and fitted values..
m4.plots <- data.frame(subset(sample.data,!is.na(sample.data$dpe)),residuals=residuals(m4),fitted=fitted.values(m4))
#pretty residuals vs fitted values
#just an example
rvf1 <- ggplot(m4.plots,aes(x=fitted,y=residuals,color=group)) + geom_point() + xlab("fitted values") + ylab("residuals") + theme_bw()
#histogram of residuals
rhist1 <- ggplot(m4.plots,aes(x=residuals)) + geom_histogram(fill="grey") + xlab("residuals") + theme_bw()
#residuals QQplot
qq1 <- ggplot(m4.plots,aes(sample=residuals)) + geom_qq(size=1,color="darkgrey") + geom_qq_line(lty=2,lwd=0.5) + xlab("thoretical quantiles") + ylab("standardised residuals") + theme_bw()
ggsave("residuals_plots_1_m4_C.pdf",rvf1 + rhist1 + qq1)
#residuals against explanatory variables
g1 <- ggplot(m4.plots,aes(y=residuals,x=dpe)) + geom_point() + theme_bw()
g2 <- ggplot(m4.plots,aes(y=residuals,x=amount)) + geom_point() + theme_bw()
g3 <- ggplot(m4.plots,aes(y=residuals,x=species)) + geom_boxplot(fill="lightgrey") + theme_bw()
g4 <- ggplot(m4.plots,aes(y=residuals,x=sex)) + geom_boxplot(fill="lightgrey") + theme_bw()
g5 <- ggplot(m4.plots,aes(y=residuals,x=blood)) + geom_boxplot(fill="lightgrey") + theme_bw()
ggsave("residuals_plots_2_m4_C.pdf",(g1 + g2) / (g3 + g4 + g5))

#there is some hint here (in plot g1) that our straight line model isn't fitting the data perfectly.
#its underestaimting the 10-day samples and overestimating for the 1-day samples. 
#We could try with a quadratic term in dpe (so include dpe^2 in the model..)
#I think that would look like this..
m6 <- lm(formula = d13C ~ species*dpe + I(dpe^2)*species + sex + blood + amount, data = sample.data)
m6.plots <- data.frame(subset(sample.data,!is.na(sample.data$dpe)),residuals=residuals(m6),fitted=fitted.values(m6))
g1 <- ggplot(m6.plots,aes(y=residuals,x=dpe)) + geom_point() + theme_bw()
#this looks a bit better, but itsn't perfect in that for high dpe the fit curves downwards.
#you might want to stick with the simpler m4 model.
#You should re-check all of the other plots too for this new model if you go with this....

#plot your model results

#the hard way..
#make a dataframe with all possible combinations of input variables..
#I think the same code will work if you go with the linear fit in model m4.
#obviously you'll need to update hte 'm6' below.
#points along dpe axis, 1 value of amount
newdata1 <- expand.grid(dpe=pretty(m6.plots$dpe,n=10),amount=mean(m6.plots$amount),species=unique(m6.plots$species),sex=c("male","female"),blood=as.factor(c(TRUE,FALSE)))
#predicted values and CIs for all combinations
m6.predict <- predict(m6,newdata1,type="response",interval="predict")
m6.predict.data <- cbind(newdata1,m6.predict)
#now can plot stuff
g1 <- ggplot(m6.predict.data,aes(x=dpe,y=fit,color=interaction(species,blood),lty=sex)) + geom_line() + geom_point(data=sample.data,aes(x=dpe,y=d13C)) + scale_color_manual(values=c("lightblue","darkblue","pink","red"),labels=c("P. sergenti, no blood meal","P. papatasi, no bloodmeal","P. sergenti, blood meal","P. papatasi, blood meal"),name="") + theme_bw()
#points for different dpe and amount of material, fixed sex and bloodmeal
newdata2 <- expand.grid(dpe=pretty(m6.plots$dpe,n=10),amount=pretty(m6.plots$amount,n=5),species=unique(m6.plots$species),sex=c("female"),blood=as.factor(c(FALSE)))
m6.predict.2 <- predict(m6,newdata2,type="response",interval="predict")
m6.predict.2.data <- cbind(newdata2,m6.predict.2)
g2 <- ggplot(m6.predict.2.data,aes(x=dpe,y=fit,color=species,lty=as.factor(amount))) + geom_line() + theme_bw() + scale_linetype(name="amount of material (mg)")
ggsave("fitted_models_d15C_m6.pdf",g1+g2)

#or use ggpredict to make plotting fitted values very easy!
#this example is for m4
install.packages("see", dependencies =
                   TRUE)
library(ggeffects)
gp1 <- ggpredict(m4, terms = c("dpe","amount","species","sex"))
gp1.g <- plot(gp1)

#look at covariation between your explanatory variables
#
#look at covariation between C, N and S
#check missing value for N for sample purple_3_N4 is in original spreadsheet

#check blood is coded properly in the data frame.


#try this for the other isotopes


#is there more variation in the biolgoy samples than the controls?
#how does control variation vary with sample size, doese this expalin how variable the samples are?
#is there any effect of 'order' in the experiment?
sample.data$order <- as.numeric(rownames(sample.data))
m5 <- lm(d13C ~ species*dpe+sex+blood+amount+order,data=sample.data)
summary(m5)
#order is not significant so ignored in all other analyses
#effect size is also tiny (1/100th of the effect of dpe)


#control data---
control.data <- subset(data,control == "control")
control.data$identifier[control.data$identifier == "cns1.1"] <- "CNS1.1"
CNS.control.data <- control.data[str_detect(control.data$identifier,"CNS"),]

#does "amount" have an effect?
CNS.m1 <- lm(d13C ~ amount+identifier, data=CNS.control.data)
#no
#does order matter?
CNS.control.data$order <- as.numeric(rownames(CNS.control.data))
CNS.m2 <- lm(d13C ~ amount+identifier+order, data=CNS.control.data)
#no


#OTHER THINGS we could try and do:
#treat dpe as categorical, see if still significant
#use post-hoc tests to see if difference is due to some particular step in the ageing

#how precisely can we estimate age from this model?
#does this need all 3 isotopes?
#might need a multivariate model.

#does adding MIRS data help?

#make pretty plots.
#see if you can get MIRS data
#presentation
#intro
  #Leishmaniasis
  #sandflies
  #stable isotopes
  

