##################### Example from Best Practice for Linear Mixed Effect Models in Psychological Science ############################
# original script run in R Studio with R version 3.3.3
# script created by LM, checked and annotated by RD (-- RD --), further comments by LM as -- LM --
# NB: script threw up errors for RD, LM diagnosed this as caused by RD using tidyverse, whereas LM had used
# ggplot2 / plyr alone. See lines 102-105.

##################### Experiment structure ############################

# Complete reference
# Meteyard, L., & Bose, A. (2018). What does a cue do? comparing phonological and semantic cues for picture naming in aphasia. 
# Journal of Speech, Language, and Hearing Research, 61(3), 658-674.

# Picture naming experiment with 10 individuals with aphasia

# Four conditions / two sessions
# Semantic: related (associated) and unrelated (not associated) word presented with picture
# Phonological: phonological cue (word onset) and tone presented with picture
# Counterbalanced so that in each session, half items cued and half uncued
# 2 phonological sessions, 2 semantic sessions

# Pictures taken from the Philadelphia Naming Test, 175 pictures
# Vary along different lexical and picture variables.
# For this tutorial, we will look at effect of Length and Frequency

# Dependent variable: Accuracy of naming (coded as 0/1)
# Generalised LMM used with logistic link function / logistic regression

# Guide to data file: NamingData
# Word - word to be named from the picture
# Subject - ID for participant
# Trial - trial number during naming session
# ItemNumber - Number code for picture 
# Session - testing session (1-4)
# Condition - Phonological or Semantic condition of session
# Acc - naming accuracy (0 or 1)
# Cue - whether item was cued (associated word or word onset) or not cued (unassociated word or tone)
# CueCondition - full label for condition - Phonological Cue (PhCue), Associated word (SemCue), Unassociated word (Unassoc) or Tone (tone)
# Length.Phon - Length of word to be named in phonemes
# FreqLogLemma - Frequency of the word, log lemma value from Celex

####### Notes on models ####

# For each model, ONLY the random effect under consideration is fit (see manuscript)
# This allowed the figures to be generated more easily, and each one to be considered on its own
# In reality, a model was fit with multiple random effects (e.g. for the effects of frequency and length and Cue type)


####### Data import and checking ####

# -- RD -- request use of here library, and run here::here() function call to locate the data file
# https://cran.r-project.org/web/packages/here/
library(here)
here::here()

# -- RD -- read-in the data file for use in the following analyses
NamingData <- read.delim("NamingData.txt")

# check responese - expect variation across groups / participants
# check number responses across participants
table(NamingData$Acc,NamingData$Subject,NamingData$Condition)
# check performance across patients and four testing sessions
table(NamingData$Acc,NamingData$Subject,NamingData$Session)
# data looks OK re: no. responses across conditions and sessions


#### Mean centering and standardize continuous IVs, check coding of predictors #######

NamingData$zLengthPh <- scale(NamingData$Length.Phon, scale = TRUE, center = TRUE)
NamingData$zFreq <- scale(NamingData$FreqLogLemma, scale = TRUE, center = TRUE)
#check
head(NamingData$zFreq)
head(NamingData$zLengthPh)

# -- RD -- ensure that cue condition is treated by R as a categorical variable, factor
NamingData$CueCondition <- as.factor(NamingData$CueCondition)


#### Response variable coding #######

#coding of response variable needs to be factor with 0=error and 1=correct, currently integer
NamingData$Acc <- as.factor(NamingData$Acc)
x <- as.integer(NamingData$Acc) 
x #coded internally as 2 and 1
#code as factor so correct = true, incorrect = false
NamingData$ACC=factor(x)
# From here on, use ACC as dependent variable
rm(x)



##################### Packages ############################

require(lmerTest)
require(lme4)
require(effects)
require(multcomp)
require(merTools)
require(tidyverse) # -- RD -- subsumes ggplot
require(lattice)
# require(reshape2) # -- RD -- subsumed under tidyverse
# require(plyr) # -- RD -- may need to uncomment to get revalue() but then get warning about dplyr-plyr clash

# -- LM -- Rob completed check with tidyverse and I think this caused errors to flag up in code
# -- LM -- Lotte has not yet transitioned to the tidyverse! For simplicity, have reverted to plyr, gglplot2 and reshape2 etc
# -- LM -- this set is below. To try running with tidyverse, use above

require(lmerTest)
require(lme4)
require(effects)
require(multcomp)
require(merTools)
require(ggplot2) 
require(lattice)
require(reshape2)
require(plyr)

###### FIGURE 1 Subjects only model to visualise random effects for subjects #########

# -- RD -- as preparation for later plot construction, we relabel subjects according to accuracy (based on analysis below)
# -- using forcats library function fct_recode()
# -- note that forcats is part of the tidyverse suite of data processing libraries:
# https://forcats.tidyverse.org/reference/fct_recode.html
# first check the level names for the original participant ID coding variable, Subject
levels(NamingData$Subject)

# then recode the participant ID (levels), creating a new ID coding variable, Participant, in the same dataset
NamingData$Participant <- fct_recode(NamingData$Subject, A = "CB", B = "AW", C = "EW", D = "JV",
                                                        E = "WR", F = "JK", G = "FF", H = "DH",
                                                        I = "MH", J = "AM")

# -- LM -- if above fails then use plyr to get revalue
require(plyr)
NamingData$Participant <- revalue(NamingData$Subject, c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                                        "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                                        ,"MH" = "I", "AM" = "J"))

# to illustrate the estimation of random deviations, we look at the random intercepts model
# -- RD -- the model is fitted using the glmer() function from the lme4 library
# -- RD -- we specify a random effect of subjects on intercepts (coded 1) with (1|Subject)

# -- RD -- random intercepts only model, with accuracy (ACC) predicted by fixed effects of cue type condition (CueCondition), plus 
# standardised target name length in phonemes and target name lexical frequency, plus the random effect due to between-subject
# variation in intercepts
Cue.Subj.lmer = glmer(ACC ~ CueCondition + zLengthPh + zFreq + (1|Subject), 
                      
                      data = NamingData, family = "binomial",
                      control = glmerControl(optimizer="bobyqa"))

summary(Cue.Subj.lmer)

#####  OUTPUT  #######

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
                                                                                  ]
#Family: binomial  ( logit )
#Formula: ACC ~ CueCondition + zLengthPh + zFreq + (1 | Subject)
#Data: NamingData
#Control: glmerControl(optimizer = "bobyqa")

#AIC      BIC     logLik    deviance  df.resid 
#6600.7   6648.0  -3293.4   6586.7     6348 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-5.0036 -0.6876  0.3306  0.5796  4.0567 

#Random effects:
#  Groups  Name        Variance Std.Dev.
#  Subject (Intercept) 1.58     1.257   
#Number of obs: 6355, groups:  Subject, 10

#Fixed effects:
#                       Estimate  Std. Error z value Pr(>|z|)    
#  (Intercept)          0.95583    0.40251   2.375   0.0176 *  
#  CueConditionSemCue  -0.62757    0.08667  -7.241 4.46e-13 ***
#  CueConditionTone    -0.59457    0.08720  -6.818 9.22e-12 ***
#  CueConditionUnAssoc -0.62037    0.08768  -7.075 1.49e-12 ***
#  zLengthPh           -0.34702    0.03439 -10.092  < 2e-16 ***
#  zFreq                0.21299    0.03432   6.205 5.46e-10 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#             (Intr) CCndSC CCndtT CCndUA zLngtP
# CueCndtnSmC -0.114                            
# CueCondtnTn -0.114  0.525                     
# CCndtnUnAss -0.113  0.523  0.520              
# zLengthPh   -0.007  0.031  0.029  0.028       
# zFreq        0.005 -0.019 -0.016 -0.021  0.407



# Plot random intercepts/variation for participants

#Extract fitted values and align with original data for plotting
a<-fitted.values(Cue.Subj.lmer)
length(a)
Dat<-na.omit(NamingData) #remove NAs from original data so aligns with model values
Dat$fitted<-a
names(Dat)


# -- RD -- run into problems with revalue in the following
# -- tried to use fct_recode() instead but that didn't work
# -- LM -- to fix didn't use tidyverse, reverted back to plyr and revalue

Dat$Accuracy <- as.integer(as.character(revalue(Dat$Acc, c("0" = 0,"1"= 1))))
temp <- ddply(Dat, "Participant", summarise,
      MeanAcc = mean(Accuracy))
temp #-- LM --check output
str(temp)#-- LM --check structure, note 'Participant' is factor
temp$MeanAcc <- round(temp$MeanAcc, digits = 2) #-- LM -- round numbers for plotting
temp$Participant<-as.character(temp$Participant)
temp[11,] <- c("Mean", round(mean(Dat$Accuracy),digits=2)) #-- LM -- add labels and mean text for graph
temp$MeanAcc <- as.numeric(temp$MeanAcc)  #-- LM -- ensure numeric data so plotting works

##### Figure 1a Plot accuracy by subject and with grand mean ######
#-- LM -- I used ggplot2 alone here not as part of tidyverse - see lines 102-105

p <- ggplot(data=temp, aes(x=as.factor(Participant), y=MeanAcc)) + geom_bar(stat="identity") 
p <- p + xlab("Participant and Grand Mean Accuracy") + ylab("Accuracy")
p <- p + ggtitle("1a) Raw participant accuracy and group/grand mean")
p <- p + expand_limits(y = c(0,1))
p <- p + theme_bw() + theme(text=element_text(size=14))
p

# -- RD -- following runs OK

# produce a quick plot of random effect of subject deviations in accuracy - subject deviations from grand mean
lattice::dotplot(ranef(Cue.Subj.lmer))

# extract random effect variances
randomeffects <- ranef(Cue.Subj.lmer)
str(randomeffects)
randomeffects$Subject
rownames(randomeffects[[1]]) <- revalue(rownames(randomeffects[[1]]), c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                                                        "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                                                        ,"MH" = "I", "AM" = "J"))
names(randomeffects$Subject) <- "Intercept"
randomeffects$Subject$Participant <- rownames(randomeffects[[1]])

##### Figure 1b Plot random effects for subjects ######
p <- ggplot(data=randomeffects$Subject, aes(x=as.factor(Participant), y=Intercept)) + geom_point(stat="identity") 
p <- p + ylab("Deviation from Grand Mean Intercept (0)") + xlab("Participant")
p <- p + geom_hline(yintercept=0)
p <- p + ggtitle("1b) Participant accuracy as deviations from the grand mean")
p <- p + theme_bw() + theme(text=element_text(size=14))
p



###### FIGURE 2 Cue Conditions over Subjects  ######

# -- RD -- fit a mixed-effects model with random effect of subjects on slopes for cue condition effects
# -- RD -- note that here we specify the random effect of subjects on the slope of the CueCondition effect
# -- RD -- by specifying (0 + CueCondition|Subject) we are asking the glmer() function to ignore potential
# deviations between subjects in the intercepts, i.e. we are excluding the random effect of subjects on intercepts
Cue.Subj.lmer = glmer(ACC ~ CueCondition + zLengthPh + zFreq + (0 + CueCondition|Subject), 
                      data = NamingData, family = "binomial",
                      control = glmerControl(optimizer="bobyqa"))

summary(Cue.Subj.lmer)

# -- RD -- note boundary singular fit warning for this model

# -- RD -- do we need to explain dummy or reference coding, and how the fixed effects and random effects estimate differ in
# results presentation from glmer()
# -- some would prefer the use of effects coding, can we flag that?
# -- LM -- I think keep it simple to illustrate random effects for the example
# -- LM -- we can't address everything here.

#####  OUTPUT  #######

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
                                                                                  ]
#Family: binomial  ( logit )
#Formula: ACC ~ CueCondition + zLengthPh + zFreq + (0 + CueCondition |      Subject)
#Data: NamingData
#Control: glmerControl(optimizer = "bobyqa")

#AIC      BIC     logLik    deviance df.resid 
#6601.2   6709.3  -3284.6   6569.2     6339 

#Scaled residuals: 
#  Min      1Q    Median      3Q     Max 
#-5.3847 -0.6704  0.3238  0.5808  3.9320 

#Random effects:
#Groups  Name                Variance Std.Dev. Corr          
#Subject CueConditionPhCue   1.608    1.268                  
#CueConditionSemCue          1.585    1.259    0.94          
#CueConditionTone            1.715    1.310    0.99 0.97     
#CueConditionUnAssoc         1.572    1.254    0.97 0.99 0.99
#Number of obs: 6355, groups:  Subject, 10

#Fixed effects:
#                       Estimate  Std. Error z value Pr(>|z|)    
#  (Intercept)          0.95128    0.40659   2.340 0.019303 *  
#  CueConditionSemCue  -0.62145    0.16964  -3.663 0.000249 ***
#  CueConditionTone    -0.58364    0.10478  -5.570 2.55e-08 ***
#  CueConditionUnAssoc -0.61682    0.13412  -4.599 4.24e-06 ***
#  zLengthPh           -0.34961    0.03448 -10.141  < 2e-16 ***
#  zFreq                0.21370    0.03442   6.209 5.34e-10 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#             (Intr) CCndSC CCndtT CCndUA zLngtP
#CueCndtnSmC -0.231                            
#CueCondtnTn -0.011  0.647                     
#CCndtnUnAss -0.204  0.821  0.665              
#zLengthPh   -0.008  0.018  0.026  0.020       
#zFreq        0.004 -0.008 -0.012 -0.012  0.407


# extract fixed effect fitted values
eff<-allEffects(Cue.Subj.lmer)
data.summary<-eff$CueCondition$fit
data.summary<-as.data.frame(data.summary)
data.summary[,2]<-eff$CueCondition$lower
data.summary[,3]<-eff$CueCondition$upper
data.summary[,4]<-eff$CueCondition$se
eff$CueCondition$variables
rownames(data.summary)<-c("PhCue","SemCue","Tone","UnAssoc")
data.summary[,5]<-c("Shared onset","Associated word","Tone","Non-associated word")
names(data.summary)[1:5]<-c("fit","lower","upper","se","Cue")
data.summary


###### Figure 2a mean accuracy across four cue conditions ######

# set up data and plot
plot <- ggplot(data.summary,aes(x=as.factor(Cue),y=fit))
plot <- plot + geom_bar(position=position_dodge(),stat="identity")
# add error bars with standard error
plot <- plot + geom_errorbar(aes(ymin=lower,ymax=upper),width=.2)
plot <- plot + ggtitle("2a) Average accuracy across four Cue Type conditions")
# add axis labels
plot<-plot + xlab("Cue Type")
plot<-plot + ylab("Accuracy (model fit)")
# use theme
plot <- plot + theme_bw() + theme(text=element_text(size=14))
plot


###### Figure 2b Plot individual subject data across cue conditions ######

# extract fitted values and align with original data for plotting
a <- fitted.values(Cue.Subj.lmer)
length(a)
Dat <- na.omit(NamingData) #remove NAs from original data so aligns with model values
Dat$fitted <- a
names(Dat)

# rename cue conditions for plot
Dat$CueType <- revalue(Dat$CueCondition, c("PhCue" = "Sh.Ons", "SemCue" = "Assoc",
                                           "Tone" = "Tone", "UnAssoc" = "NonAssoc"))

# change to black and white/grey scale 
trellis.par.set(canonical.theme(color = FALSE))

unique(Dat$Participant) #--LM-- check order of participants for plotting

# box plot of patient by cue type data
lattice::bwplot(fitted(Dat) ~ CueType|Participant, data=Dat,
                index.cond=list(c(9,1,6,4,10,7,5,8,3,2)), #panel order A-J
                ylab="Accuracy (model fit)",xlab="Cue Type",
                pch ="|", notch = TRUE, layout=c(2, 5, 1),
                main="2b) Effect of Cue Type by Participant")




###### Figure 2c Plot cue condition by subject random effects ######

# -- RD -- get error for following, in trying to produce 2c:
#   Error in `$<-.data.frame`(`*tmp*`, "se", value = c(0.156386150052519,  : 
#                                                        replacement has 40 rows, data has 50
# -- maybe something to do with what I was unable to run earlier
# -- do get a plot but panels are wrong, i.e. shows cue type within subjects rather than subjects within cue type
# --LM-- I don't get above error - think it is to do with tidyverse vs ggplot2 and plyr, see lines 102-105

                                                     
#Extract random effect variances
randomeffects <- ranef(Cue.Subj.lmer)
str(randomeffects) #look at what has been extracted
rownames(randomeffects[[1]]) <- revalue(rownames(randomeffects[[1]]), c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                                                        "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                                                        ,"MH" = "I", "AM" = "J"))
names(randomeffects[[1]]) <- c("Shared Onset", "Associated Word", "Tone", "Non associated word")
randomeffects$Subject$Participant <- rownames(randomeffects[[1]])

str(randomeffects) # --LM-- check structure of data before plot

#Quick plot of random effects (will order participants by deviation, and panel by Cue Type)
lattice::dotplot(ranef(Cue.Subj.lmer))

#Prep data for ggplot, put into long format
x<-melt(randomeffects$Subject,id=c("Participant"))
names(x) <- c("Participant","CueType","Deviation")

#Plot to match - panel by cue type as for Figure 2c
p <- ggplot(data=x, aes(x=as.factor(Participant), y=Deviation)) + geom_point(stat="identity") 
p <- p + ylab("Deviation from Condition Mean (0)") + xlab("Participant")
p <- p + geom_hline(yintercept=0)
p <- p + facet_grid(. ~ CueType)
p <- p + ggtitle("2c) Effect of Cue Type by Participant, as deviations from Condition Mean")
p <- p + theme_bw() + theme(text=element_text(size=14))
p

#Plot to panel by participant rather than cue type (just to show alternative)
p <- ggplot(data=x, aes(x=as.factor(CueType), y=Deviation)) + geom_point(stat="identity") 
p <- p + ylab("Deviation from Condition Mean (0)") + xlab("Cue Type")
p <- p + geom_hline(yintercept=0)
p <- p + facet_grid(. ~ Participant)
p <- p + theme_bw() + theme(text=element_text(size=12))
p


#### FIGURE 3 Length in phonemes #######

# -- RD -- fit a model with fixed effects and with the random effect of subjects on intercepts and the slopes for length 
Cue.Subj.lmer = glmer(ACC ~ CueCondition + zLengthPh + zFreq + (1+zLengthPh|Subject), 
                      data = NamingData, family = "binomial",
                      control = glmerControl(optimizer="bobyqa"))

summary(Cue.Subj.lmer)

#####  OUTPUT  #######

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
                                                                                  ]
#Family: binomial  ( logit )
#Formula: ACC ~ CueCondition + zLengthPh + zFreq + (1 + zLengthPh | Subject)
#Data: NamingData
#Control: glmerControl(optimizer = "bobyqa")

# AIC      BIC      logLik   deviance df.resid 
# 6595.9   6656.7  -3288.9   6577.9     6346 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-4.7441 -0.6793  0.3361  0.5689  4.2515 

#Random effects:
#  Groups  Name        Variance Std.Dev. Corr
#  Subject (Intercept) 1.58402  1.2586       
#  zLengthPh          0.01816  0.1348   0.34
#Number of obs: 6355, groups:  Subject, 10

#Fixed effects:
#                      Estimate  Std. Error z value Pr(>|z|)    
# (Intercept)          0.94763    0.40310   2.351   0.0187 *  
#  CueConditionSemCue  -0.62891    0.08686  -7.241 4.47e-13 ***
#  CueConditionTone    -0.59569    0.08742  -6.814 9.47e-12 ***
#  CueConditionUnAssoc -0.62221    0.08789  -7.080 1.44e-12 ***
#  zLengthPh           -0.34817    0.05531  -6.295 3.07e-10 ***
#  zFreq                0.21124    0.03432   6.154 7.55e-10 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) CCndSC CCndtT CCndUA zLngtP
#CueCndtnSmC -0.114                            
#CueCondtnTn -0.114  0.526                     
#CCndtnUnAss -0.113  0.524  0.521              
#zLengthPh    0.256  0.019  0.018  0.017       
#zFreq        0.005 -0.018 -0.016 -0.020  0.254


# quick plot of random effects from model
lattice::dotplot(ranef(Cue.Subj.lmer))

# extract fitted values and align with original data for plotting
a<-fitted.values(Cue.Subj.lmer)
length(a)
Dat<-na.omit(NamingData) #remove NAs from original data so aligns with model values
Dat$fitted<-a
names(Dat)
#relabel participants according to initial severity from random intercepts (see above)
Dat$Participant <- revalue(Dat$Subject, c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                          "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                          ,"MH" = "I", "AM" = "J"))

###### Figure 3a Plot average slope for effect of length ######
p <- ggplot(data=Dat,aes(x=zLengthPh, y=fitted))
p <- p + geom_smooth(method = "glm")
p <- p + xlab("Length in Phonemes (z score)") + ylab("Accuracy (model fit)")
p <- p + ggtitle("3a) Average (group) effect of Length in Phonemes")
p <- p + ylim(0,1)
p <- p + theme_bw() + theme(text=element_text(size=14))
p

#--LM-- appears to be ignoring call to black and white, gave up trying to fix
#will make greyscale when saving final image.

#### Figure 3b Plot individual slopes for effect of length #######
#--LM-- need to fiddle with participant labelling so order is correct for plot
str(Dat) #check format of participant coding
names(Dat) 
#going to rename column with participant labelled as factor, and create new column with text/character participant labels
names(Dat)[15]<-"PartCodeFact" #participant labels as factor in this column
Dat$Participant<-as.character(Dat$PartCodeFact) #convert participant labels to text

p <- ggplot(data=Dat, aes(x=zLengthPh, y=fitted, linetype=Participant))
p <- p + geom_smooth(se = FALSE, method = "glm")
p <- p + xlab("Length in Phonemes (z score)") + ylab("Accuracy (model fit)")
p <- p + ggtitle("3b) Effect of Length in Phonemes by Participant")
p <- p + ylim(0,1)
p <- p + theme(text=element_text(size=14))
p <- p + scale_colour_grey() + theme_bw()
p  

# -- RD -- ggplot appears to be refusing to recognise the call for grey scale as slopes shown in blue when I run this
# --LM-- yes, gave up trying to fix within R

#### Figure 3c random effects for Length in phonemes #######

# extract random effect variances
randomeffects <- ranef(Cue.Subj.lmer)
str(randomeffects) #look at what has been extracted
randomeffects[[1]]
names(randomeffects[[1]]) <- c("Intercept", "Length in Phonemes")
rownames(randomeffects[[1]]) <- revalue(rownames(randomeffects[[1]]), c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                                                        "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                                                        ,"MH" = "I", "AM" = "J"))
randomeffects$Subject$Participant <- rownames(randomeffects[[1]])

# quick plot of random effects (will order participants by deviation, and using original participants codes from model)
lattice::dotplot(ranef(Cue.Subj.lmer))

# prepare data for ggplot, put into long format
x<-melt(randomeffects$Subject,id=c("Participant"))
names(x) <- c("Participant","variable","Deviation")

# produce plot to match previous using ggplot2 
p <- ggplot(data=x, aes(x=as.factor(Participant), y=Deviation)) + geom_point(stat="identity") 
p <- p + ylab("Deviation from Mean (0)") + xlab("Participant")
p <- p + ggtitle("3c) Participant intercepts and slopes for Length")
p <- p + geom_hline(yintercept=0)
p <- p + facet_grid(. ~ variable)
p <- p + theme_bw() + theme(text=element_text(size=14))
p


#### FIGURE 4 Frequency #######

# -- RD -- slopes and intercepts for frequency
Cue.Subj.lmer = glmer(ACC ~ CueCondition + zLengthPh + zFreq + (1+zFreq|Subject), 
                      data = NamingData, family = "binomial",
                      control = glmerControl(optimizer="bobyqa"))

summary(Cue.Subj.lmer)

#####  OUTPUT  #######

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
                                                                                  ]
#Family: binomial  ( logit )
#Formula: ACC ~ CueCondition + zLengthPh + zFreq + (1 + zFreq | Subject)
#Data: NamingData
#Control: glmerControl(optimizer = "bobyqa")

#AIC      BIC   logLik deviance df.resid 
#6577.3   6638.1  -3279.6   6559.3     6346 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-3.9959 -0.6727  0.3477  0.5622  5.7565 

#Random effects:
#  Groups  Name        Variance Std.Dev. Corr 
#Subject (Intercept)   1.60700  1.2677        
#zFreq                 0.03981  0.1995   -0.84
#Number of obs: 6355, groups:  Subject, 10

#Fixed effects:
#                      Estimate Std. Error z value Pr(>|z|)    
#(Intercept)            0.93783    0.40600   2.310   0.0209 *  
#  CueConditionSemCue  -0.63507    0.08714  -7.288 3.15e-13 ***
#  CueConditionTone    -0.59970    0.08767  -6.841 7.88e-12 ***
#  CueConditionUnAssoc -0.62466    0.08813  -7.088 1.36e-12 ***
#  zLengthPh           -0.36034    0.03513 -10.256  < 2e-16 ***
#  zFreq                0.17644    0.07239   2.437   0.0148 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#             (Intr) CCndSC CCndtT CCndUA zLngtP
#CueCndtnSmC -0.114                            
#CueCondtnTn -0.113  0.528                     
#CCndtnUnAss -0.113  0.525  0.522              
#zLengthPh   -0.006  0.031  0.029  0.029       
#zFreq       -0.726 -0.008 -0.007 -0.009  0.205



# quick plot of random effects from model
lattice::dotplot(ranef(Cue.Subj.lmer))

# extract fitted values and align with original data for plotting
a<-fitted.values(Cue.Subj.lmer)
length(a)
Dat<-na.omit(NamingData) #remove NAs from original data so aligns with model values
Dat$fitted<-a
names(Dat)
#relabel participants according to initial severity from random intercepts (see above)
Dat$Participant <- revalue(Dat$Subject, c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                          "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                          ,"MH" = "I", "AM" = "J"))

#### FIGURE 4a Plot average slope for effect of Frequency  #######

p <- ggplot(data=Dat,aes(x=zFreq, y=fitted))
p <- p + geom_smooth(method = "glm")
p <- p + xlab("Frequency (z score)") + ylab("Accuracy (model fit)")
p <- p + ggtitle("4a) Average (group) effect of Frequency")
p <- p + ylim(0,1)
p <- p + theme_bw() + theme(text=element_text(size=12))
p

#### FIGURE 4b Plot individual slopes for effect of Frequency  #######
#--LM-- need to fiddle with participant labelling so order is correct for plot
str(Dat) #check format of participant coding
names(Dat) 
#going to rename column with participant labelled as factor, and create new column with text/character participant labels
names(Dat)[15]<-"PartCodeFact" #participant labels as factor in this column
Dat$Participant<-as.character(Dat$PartCodeFact) #convert participant labels to text

p <- ggplot(data=Dat, aes(x=zFreq, y=fitted, linetype=Participant))
p <- p + geom_smooth(se = FALSE, method = "glm")
p <- p + xlab("Frequency (z score)") + ylab("Accuracy (model fit)")
p <- p + ylim(0,1)
p <- p + ggtitle("4b) Effect of Frequency by Participant")
p <- p + theme(text=element_text(size=14))
p <- p + scale_colour_grey() + theme_bw()
p  

# extract random effect variances
randomeffects <- ranef(Cue.Subj.lmer)
str(randomeffects) #look at what has been extracted
randomeffects[[1]]
names(randomeffects[[1]]) <- c("Intercept", "Frequency")
rownames(randomeffects[[1]]) <- revalue(rownames(randomeffects[[1]]), c("CB" = "A", "AW" = "B","EW" = "C","JV" = "D",
                                                                        "WR" = "E", "FF" = "G", "JK" = "F", "DH" = "H"
                                                                        ,"MH" = "I", "AM" = "J"))
randomeffects$Subject$Participant <- rownames(randomeffects[[1]])

#Quick plot of random effects (will order participants by deviation)
lattice::dotplot(ranef(Cue.Subj.lmer))

# prepare data for ggplot, put into long format
x<-melt(randomeffects$Subject,id=c("Participant"))
names(x) <- c("Participant","variable","Deviation")

#### FIGURE 4c Random effects for frequency #######

p <- ggplot(data=x, aes(x=as.factor(Participant), y=Deviation)) + geom_point(stat="identity") 
p <- p + ylab("Deviation from Mean (0)") + xlab("Participant")
p <- p + ggtitle("4c) Participant intercepts and slopes for Frequency")
p <- p + geom_hline(yintercept=0)
p <- p + facet_grid(. ~ variable)
p <- p + theme_bw() + theme(text=element_text(size=14))
p

#### FIGURE 5 Pubmed citations per year #######

number <- c(11,15,20,22,36,51,66,88,105,173,219,272,346,423,590,685,849,979,1067)
year <- c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018)
dat <- as.data.frame(cbind(number,year))

p <- ggplot(data=dat, aes(x=as.factor(year), y=number)) + geom_bar(stat="identity") 
p <- p + ylab("Number of Pubmed citations") + xlab("Year")
p <- p + theme_bw() + theme(text=element_text(size=14))
p
