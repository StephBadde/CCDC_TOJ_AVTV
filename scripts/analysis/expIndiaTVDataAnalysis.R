# clears workspace:  
rm(list=ls(all=TRUE)) 

#### functions ####

#function to detach all packages (to avoid interferences)
detachAllPackages <- function() { 
  #store basic packages names in a list
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")  
  #make list of all loaded packages
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  #remove basic packages from the list
  package.list <- setdiff(package.list,basic.packages)
  #remove all packages from the list
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
#detach all packages
detachAllPackages()

#negative log-likelihood of cumulative Gaussian fit, takes parameter vector and data as input
nll <- function(p, d) { 
  #read out expected cumulative gaussian probabilities for each distance
  phi <- p[1] + (1-2*p[1])*pnorm(d$SOA, p[2], p[3]) 
  #ensure not to hit boundaries, otherwise log will cause problems
  phi <- pmax(pmin(phi, .99999), .00001)
  #negative log-likelihood of probability to get the number of right responses in the number of trials given probability phi
  -sum(d$noFurtherRight * log(phi) + (d$noTrials - d$noFurtherRight) * log(1 - phi) )
}

#### load requiered libraries ####
#load plyr
library(plyr)
#need car package for recode function
library(car)
#load R.utils
library(R.utils)
#need gdata package for reorder function
library(gdata)
#package for plotting
library(ggplot2)
#needed for unit
library(grid)
#load ez package
library(ez)
#for string manipulations
library(stringr)
#needed for glmm
library(lme4)
#needed for interaction testing
library(phia)
#better data.frames
library(data.table)
#partial correlation
library(ppcor)


###### general settings #######
#save path to main directory
path="~/Documents/Experimente/expXV"

#store colours
myColours.TV <- c('#B68F8B','#e0e1e0', '#8c3d35')

#indicate whether single subjects should be plotted
plotSingleSubjects <- T

###### data collection and transformation #######
print("####################### data collection and transformation #############################")
#change working directory to data directory
setwd(paste(path,"/data/expIndia_TV/prepared",sep=""))
#list all txt files in data directory
listOfFiles = list.files(pattern= "bimodalTOJ")
#read every .txt file in and append it to data frame
dataSet.old  <-  do.call(rbind.fill, lapply(listOfFiles, function(fname) {
  #print to console the name of the file
  print(fname)
  #read in the file
  dum  <-  read.table(fname, header = TRUE, sep="\t")
  #return data.frame to calling function
  return(dum)
}))

#list all txt files in data directory
listOfFiles = list.files(pattern= "expXV")
#read every .txt file in and append it to data frame
dataSet.newer  <-  do.call(rbind, lapply(listOfFiles, function(fname) {
  #print to console the name of the file
  print(fname)
  #read in the file
  dum  <-  read.table(fname, header = TRUE, sep="\t")
  #return data.frame to calling function
  return(dum)
}))

#read in subject list
subjectList.India <- read.table("TOJ_VT.txt", header = TRUE, sep="\t", as.is = T)

#change working directory to data directory
setwd(paste(path,"/data/expIndia_TV/prepared_HH",sep=""))
#list all txt files in data directory
listOfFiles.HH = list.files(pattern= "_VT.txt")
#read every .txt file in and append it to data frame
dataSet.newer.HH  <-  do.call(rbind, lapply(listOfFiles.HH, function(fname) {
  #print to console the name of the file
  print(fname)
  #read in the file
  dum  <-  read.table(fname, header = TRUE, sep="\t")
  #return data.frame to calling function
  return(dum)
}))

#add suffix to subject IDs to distinguish from India
dataSet.newer.HH$subNo <- paste(dataSet.newer.HH$subNo, '_HH', sep='')
#read in subject list
subjectList.HH <- read.table("TOJ_VT_HH.txt", header = TRUE, sep="\t", as.is = T)

#change working directory to data directory
setwd(paste(path,"/data/expIndia_TV/prepared_HH_control",sep=""))
#list all txt files in data directory
listOfFiles.HH.control = list.files(pattern= "_VT.txt")
#read every .txt file in and append it to data frame
dataSet.newer.HH.control  <-  do.call(rbind, lapply(listOfFiles.HH.control, function(fname) {
  #print to console the name of the file
  print(fname)
  #read in the file
  dum  <-  read.table(fname, header = TRUE, sep="\t")
  #return data.frame to calling function
  return(dum)
}))

#add suffix to subject IDs to distinguish from India
dataSet.newer.HH.control$subNo <- paste(dataSet.newer.HH.control$subNo, '_HH_control', sep='')
#read in subject list
subjectList.HH.control <- read.table("TOJ_VT_HH_control.txt", header = TRUE, sep="\t", as.is = T)

#change working directory to data directory
setwd(paste(path,"/data/expIndia_TV/prepared_control",sep=""))
#list all txt files in data directory
listOfFiles.control = list.files(pattern= "_VT.txt")
#read every .txt file in and append it to data frame
dataSet.newer.control  <-  do.call(rbind, lapply(listOfFiles.control, function(fname) {
  #print to console the name of the file
  print(fname)
  #read in the file
  dum  <-  read.table(fname, header = TRUE, sep="\t")
  #return data.frame to calling function
  return(dum)
}))

#add suffix to subject IDs to distinguish from India
dataSet.newer.control$subNo <- paste(dataSet.newer.control$subNo, '_control', sep='')
#read in subject list
subjectList.control <- read.table("TOJ_VT_control.txt", header = TRUE, sep="\t", as.is = T)


###### rename variables #######
print("####################### rename variables #############################")
#rename some variables
dataSet.old <- rename(dataSet.old, c("responsetimeTOJ" = "RT", "accuracyTOJ" = "accuracy",
                                     "Hand" = "handPosition", "firstHand" = "firstSide", "secondHand" = "secondSide",
                                     "trialCondition" = "modalityCombination"))
dataSet.newer <- rename(dataSet.newer, c("responsetimeTOJ" = "RT", "accuracyTOJ" = "accuracy"))
dataSet.newer.HH <- rename(dataSet.newer.HH, c("responsetimeTOJ" = "RT", "accuracyTOJ" = "accuracy"))
dataSet.newer.control <- rename(dataSet.newer.control, c("responsetimeTOJ" = "RT", "accuracyTOJ" = "accuracy"))
dataSet.newer.HH.control <- rename(dataSet.newer.HH.control, c("responsetimeTOJ" = "RT", "accuracyTOJ" = "accuracy"))

###### recode variables #######
print("####################### recode variables #############################")
#somebody forgot to code the posture for no 31
dataSet.old$handPosition[is.na(dataSet.old$handPosition)] <- 1
#recode posture
dataSet.old$handPosition <- factor(dataSet.old$handPosition, levels = c(1, 2), labels = c("handsUncrossed", "handsCrossed"))
#recode hand
dataSet.old$firstSide <- factor(dataSet.old$firstSide, levels = c(1, 2), labels = c("left","right"))
#recode hand
dataSet.old$secondSide <- factor(dataSet.old$secondSide, levels = c(1, 2), labels = c("left","right"))
#determine modality
dataSet.old$firstModality <- substr(dataSet.old$subversion, dataSet.old$firstModality, dataSet.old$firstModality)
#determine modality
dataSet.old$secondModality <- substr(dataSet.old$subversion, dataSet.old$secondModality, dataSet.old$secondModality)
#re-do modality combination (coding the same for AT and AV :( )
dataSet.old$modalityCombination <- paste(dataSet.old$firstModality, dataSet.old$secondModality, sep="")
#recode modality
dataSet.old$firstModality <- factor(dataSet.old$firstModality, levels=c('A', 'T', 'V'), labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.old$secondModality <- factor(dataSet.old$secondModality, levels=c('A', 'T', 'V'), labels = c("audio", "tactile", "visual"))

#recode modality
dataSet.newer$firstModality <- factor(dataSet.newer$firstModality, 
                                      levels = c("audioFirst", "tactileFirst", "visualFirst"), 
                                      labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer$secondModality <- factor(dataSet.newer$secondModality, 
                                       levels = c("audioSecond", "tactileSecond", "visualSecond"), 
                                       labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.HH$firstModality <- factor(dataSet.newer.HH$firstModality, 
                                         levels = c("audioFirst", "tactileFirst", "visualFirst"), 
                                         labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.HH$secondModality <- factor(dataSet.newer.HH$secondModality, 
                                          levels = c("audioSecond", "tactileSecond", "visualSecond"), 
                                          labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.control$firstModality <- factor(dataSet.newer.control$firstModality, 
                                              levels = c("audioFirst", "tactileFirst", "visualFirst"), 
                                              labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.control$secondModality <- factor(dataSet.newer.control$secondModality, 
                                               levels = c("audioSecond", "tactileSecond", "visualSecond"), 
                                               labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.HH.control$firstModality <- factor(dataSet.newer.HH.control$firstModality, 
                                                 levels = c("audioFirst", "tactileFirst", "visualFirst"), 
                                                 labels = c("audio", "tactile", "visual"))
#recode modality
dataSet.newer.HH.control$secondModality <- factor(dataSet.newer.HH.control$secondModality, 
                                                  levels = c("audioSecond", "tactileSecond", "visualSecond"), 
                                                  labels = c("audio", "tactile", "visual"))

#recode side
dataSet.newer$firstSide <- factor(dataSet.newer$firstSide, 
                                  levels = c("leftSideFirst", "rightSideFirst"), 
                                  labels = c("left", "right"))
#recode side
dataSet.newer.HH$firstSide <- factor(dataSet.newer.HH$firstSide, 
                                     levels = c("leftSideFirst", "rightSideFirst"), 
                                     labels = c("left", "right"))
#recode side
dataSet.newer.control$firstSide <- factor(dataSet.newer.control$firstSide, 
                                          levels = c("leftSideFirst", "rightSideFirst"), 
                                          labels = c("left", "right"))
#recode side
dataSet.newer.HH.control$firstSide <- factor(dataSet.newer.HH.control$firstSide, 
                                             levels = c("leftSideFirst", "rightSideFirst"), 
                                             labels = c("left", "right"))

###### combine dataSets #######
print("####################### combine dataSets #############################")
#missing columns will be filled with NA
dataSet <- rbind.fill(dataSet.old, dataSet.newer)
#missing columns will be filled with NA
dataSet <- rbind.fill(dataSet, dataSet.newer.HH)
#missing columns will be filled with NA
dataSet <- rbind.fill(dataSet, dataSet.newer.control)
#missing columns will be filled with NA
dataSet <- rbind.fill(dataSet, dataSet.newer.HH.control)
#add subject list
dataSet <- merge(dataSet, rbind.fill(subjectList.India, subjectList.HH, subjectList.control, subjectList.HH.control), by="subNo")

#remove single datasets
rm(dataSet.old, dataSet.newer, dataSet.newer.control, dataSet.newer.HH, dataSet.newer.HH.control)
rm(subjectList.India, subjectList.HH, subjectList.HH.control, subjectList.control)
#rm directory inventories
rm(listOfFiles, listOfFiles.control, listOfFiles.HH, listOfFiles.HH.control)

#turn dataSet into a data.table
dataSet <- data.table(dataSet)

###### add variables #######
print("####################### add variables #############################")

#add SOA as factor for analysis
dataSet[, SOAFactor := as.factor(SOA)]
#add catch trial accuracy
dataSet[, catchTrialAccuracy := ifelse(accuracy == 6, 0, ifelse(accuracy == 7, 1, NA))]
#add right first responses
dataSet[, rightSideResp := dataSet$buttonTOJ - 1]
#recode subNo
dataSet[, subNo := paste('ID', subNo)]
#change group acronyms
dataSet[,group := factor(group, 
                         levels = c("CC", "CC_C", "DC", "DC_C", "unclear"), 
                         labels = c("CC", "MCC", "DC", "MDC", "unclear"))]

#add modalityMatch where it is missing
dataSet[, modalityMatch := ifelse(modalityCombination == 'VV' | modalityCombination == 'TT', 'unimodal', 'bimodal')]
#add right modality as it is missing in old dataset
dataSet[, rightModality := ifelse(firstSide == 'right', as.character(firstModality), as.character(secondModality))]
#code SOA with respect to modality
dataSet[, SOA.Modality := ifelse(firstModality == 'visual', abs(SOA), - abs(SOA))]
#store modality condition
dataSet[, modalityCondition := ifelse(modalityMatch == 'unimodal', 
                                      ifelse(firstModality == 'tactile', 
                                             'tactile', 
                                             'visual'),
                                      'bimodal')]
#set order of levels
dataSet$modalityCondition <- factor(dataSet$modalityCondition, levels = c('bimodal', 'visual', 'tactile'))
#reorder posture
dataSet$handPosition <- factor(dataSet$handPosition, levels=c("handsUncrossed", "handsCrossed"), labels=c("uncrossed", "crossed"))
#code different response
dataSet[, visualFirstResp := ifelse(modalityMatch == 'bimodal',
                                    ifelse(rightModality == 'visual' & rightSideResp == 1, 1,
                                           ifelse(rightModality == 'visual' & dataSet$rightSideResp == 0, 0,
                                                  ifelse(rightModality != 'visual' & rightSideResp == 1, 0,
                                                         ifelse(rightModality != 'visual' & rightSideResp == 0, 1, NA)))), NA)]
#calculate duration of sight restored
dataSet[,durationSight := age * 12 - surgery]

#translate Snellen-20 into logMar
dataSet[, logMar.OD := -round(log10(20/Snellen20OD),2)]
dataSet[, logMar.OS := -round(log10(20/Snellen20OS),2)]
dataSet[, logMar.betterEye := ifelse(!is.na(logMar.OD) & !is.na(logMar.OS), 
                                     pmin(logMar.OD, logMar.OS),
                                     ifelse(!is.na(logMar.OD), logMar.OD,
                                            ifelse(!is.na(logMar.OS), logMar.OS,
                                                   NA)))]
dataSet[, logMar.worseEye := ifelse(!is.na(logMar.OD) & !is.na(logMar.OS), 
                                    pmax(logMar.OD, logMar.OS),
                                    ifelse(!is.na(logMar.OD), 5,
                                           ifelse(!is.na(logMar.OS), 5,
                                                  NA)))]

###### data cleaning #######
print("####################### data cleaning #############################")

#filter for invalid trials
dataSet <- dataSet[dataSet$accuracy <= 1 | is.na(dataSet$accuracy),]

#filter for practice trials
dataSet <- dataSet[dataSet$practiceLevel == 0,]

#filter SOAs not used in all subjects
dataSet <- dataSet[is.element(abs(dataSet$SOA), c(400, 135, 90, 30)),]

#function to clean by standard deviations
clean <- function(x) {replace(x, x > mean(x, na.rm=T) + 2.5 * sd(x, na.rm=T)  | x < 100  | x > 40000, NA)}
#clean the data for each subject and condition (reorders the data by conditions as a side product)
dataSet[, RTclean := clean(RT), by = c('subNo', 'handPosition')]
#calculate percentage of excluded trials
print(paste('trials excluded: ', round(100 - 100*sum(!is.na(dataSet$RTclean))/sum(!is.na(dataSet$RT)), 2), '%', sep=''))

#remove these trials 
dataSet <- dataSet[!is.na(RTclean),]

###### data preparation #######
print("####################### data preparation #############################")

#change type of each variable
dataSet <- data.table(data.frame(mapply(function(dataSet, class) match.fun(paste("as", class, sep="."))(dataSet),
                                        dataSet, colClasses("fffffffdfffffdddddddfffffffdffdfdddfdddfddddddd"), SIMPLIFY=FALSE),
                                 stringsAsFactors=FALSE))

#determine list of subjects
subjList <- unique(dataSet$subNo)
#determine number of subjects
n <- length(subjList)

#divide into CC and DC
n.CC <- length(unique(dataSet$subNo[dataSet$action == 'keep' & dataSet$group == 'CC']))
n.DC <- length(unique(dataSet$subNo[dataSet$action == 'keep' & dataSet$group == 'DC']))
n.MCC <- length(unique(dataSet$subNo[dataSet$action == 'keep' & dataSet$group == 'MCC']))
n.MDC <- length(unique(dataSet$subNo[dataSet$action == 'keep' & dataSet$group == 'MDC']))

# store mean values in a txt file
write.table(subset(dataSet, action == 'keep' & handPosition == 'uncrossed'),
            file = paste(path, '/data/expIndia_TV/aggregated/expXV_India_TV_dataSet.txt', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = T)


#### data aggregation for graphix ####
print("####################### data aggregation for graphix #############################")

# sample stats
my.summary = function(x) list(mean = mean(x, na.rm = T), median = median(x, na.rm = T),
                              min = min(x, na.rm = T), max = max(x, na.rm = T))
# extract only sample data
sampleData <- unique(subset(dataSet, action == 'keep', select = c('subNo', 'group', 'age', 'surgery', 'durationSight', 'logMar.betterEye', 'Snellen20OD','Snellen20OD')))
# do sample stats
sampleData[, as.list(unlist(lapply(.SD, my.summary))), .SDcols = c('age', 'surgery', 'durationSight', 'logMar.betterEye'), by = c('group')]


#calculate mean values for each subject and factor level (excludes factor SOA)
meanValuesSOASubj <- dataSet[, list(propFurtherRight = mean(rightSideResp, na.rm = T),
                                    noFurtherRight = sum(rightSideResp),
                                    noTrials = sum(!is.na(rightSideResp))), 
                             by = c('SOA', 'modalityMatch', 'rightModality', "subversion", "group",
                                    "age", "action", 'handPosition', "subNo", 'surgery', 'durationSight', 
                                    'logMar.OD', 'logMar.OS', 'logMar.betterEye')]

#aggregate across subjects
meanValuesSOAGroup <- meanValuesSOASubj[action == 'keep', list(propFurtherRight.M = mean(propFurtherRight),
                                                               propFurtherRight.SE = sd(propFurtherRight)/sqrt(.N),
                                                               noFurtherRight = mean(noFurtherRight), 
                                                               noTrials = mean(noTrials)),
                                        by = c('SOA', 'modalityMatch', 'rightModality', "subversion", "group", 'action', 'handPosition')]

#function to wrap call to fitting function
fitCumulativeGaussian <- function(x){
  #fit sigmoid
  parameters <- optim(c(0.1, 0, 20), nll, d = x, lower = c(0.0, -400, 5), upper = c(0.2, 400, 2000), method='L-BFGS-B')$par
  #return parameters as list
  list(sigma=parameters[3], lambda=parameters[1], mu=parameters[2])
}

#add parameter estimates per condition
meanValuesSOAGroup[,  c('sigma', 'lambda', 'mu') := fitCumulativeGaussian(list(SOA=SOA,noTrials=noTrials, noFurtherRight=noFurtherRight)),
                   by = c('modalityMatch', 'rightModality', "subversion", "group", 'action', 'handPosition')]
#recode group
meanValuesSOAGroup$group2 <- factor(meanValuesSOAGroup$group, levels = c('CC', 'MCC', 'DC', 'MDC'), 
                                   labels = c(paste('CC (N = ',n.CC,')', sep=''), paste('MCC (N = ',n.MCC,')', sep=''), 
                                              paste('DC (N = ',n.DC,')', sep=''), paste('MDC (N = ',n.MDC,')', sep=''))) 
#generate datapoints for estimated distribution
curvesSOAGroup <- merge(unique(meanValuesSOAGroup[,c("group", "group2", 'modalityMatch', "rightModality", 'action', 'handPosition', 'sigma', 'lambda', 'mu')]),
                        CJ(SOA = seq(-420, 420, 1), group = unique(meanValuesSOAGroup$group)), by = 'group', all = T, allow.cartesian=TRUE)
#calculate function values
curvesSOAGroup[, propFurtherRight.M := lambda + (1-2*lambda)*pnorm(SOA, mu, sigma)]

#calculate mean values for each subject and factor level (excludes factor SOA)
meanValuesSubj <- dataSet[, list(visualFirstResp = mean(visualFirstResp, na.rm = T),
                                 accuracy = mean(accuracy, na.rm = T)), 
                          by = c('modalityCondition', "subversion", "group",
                                 "age", "action", "subNo", 'surgery', 'durationSight', 
                                 'logMar.OD', 'logMar.OS', 'logMar.betterEye', 'handPosition')]

#aggregate across subjects
meanValuesGroup <- meanValuesSubj[action == 'keep', list(visualFirstResp.M = mean(visualFirstResp),
                                                         visualFirstResp.SE = sd(visualFirstResp)/sqrt(.N),
                                                         accuracy.M = mean(accuracy), 
                                                         accuracy.SE = sd(accuracy)/sqrt(.N)),
                                  by = c('modalityCondition', "subversion", "group", 'handPosition')]


# store mean values in a txt file
write.table(subset(meanValuesSubj, action == 'keep' & handPosition == 'uncrossed'),
            file = paste(path, '/data/expIndia_TV/aggregated/expXV_India_TV_meanValues.txt', sep = ''),
            quote = F, sep = '\t', row.names = F, col.names = T)


###### graphix - mean values SOA ######
print('###### graphix - mean values SOA ######')

#set basic layout of graph
gp <- ggplot(data = subset(meanValuesSOAGroup, handPosition == 'uncrossed'), 
             aes(x=SOA, y=propFurtherRight.M, colour=rightModality))
#set theme
gp <- gp + theme_bw(8)
#plot mean values per condition as point
gp <- gp + geom_errorbar(aes(ymin=propFurtherRight.M-propFurtherRight.SE, 
                             ymax=propFurtherRight.M+propFurtherRight.SE), 
                         width=50, size = 0.2)
#plot mean values per condition as line
gp <- gp + geom_line(data = subset(curvesSOAGroup, handPosition == 'uncrossed'), size = 0.3)
#plot mean values per condition as point
gp <- gp + geom_point(size=0.75)
#divide plot into condition panels
gp <- gp + facet_grid(modalityMatch ~ group2)
#set x-axis label
gp <- gp + xlab("SOA (negative numbers indicate left side first stimulation)")
#set y-axis label
gp <- gp + ylab("p('right side first' - response)")
#axis limits
gp <- gp + coord_cartesian(xlim=c(-450,450), ylim=c(-0.1,1.1))
#axis ticks
gp <- gp + scale_y_continuous(breaks=seq(0.0,1.0,0.5))
gp <- gp + scale_x_continuous(lim=c(-450,450), breaks=seq(-300,300,300))
#set colours
gp <- gp + scale_colour_manual(values=c('#848354', '#C3A3A4'), name = "Modality at right side")
gp <- gp + scale_colour_manual(values=myColours.TV[c(3,2)], name = "Modality at right side")
#legend
gp <- gp + guides(colour=guide_legend(title='Modality presented at right side'))
#remove lines and shades
gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank(), panel.grid.minor = element_blank())
gp <- gp + theme(axis.ticks = element_blank())
gp <- gp + theme(strip.text = element_text(size = 7))
gp <- gp + theme(legend.position = "bottom",
                 legend.box.background = element_blank(),
                 legend.box.margin = margin(0, 0, 0, 0),
                 legend.box.spacing = unit(0.2, 'lines'),
                 legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = 'lines'),
                 legend.key = element_blank(),
                 legend.key.height = unit(.95, 'lines'))
#show plot
gp
#set filename
filename="/results/expIndia_TV/expIndia_TV_TOJ_RightSide_TV_meanValues_Fig_1.pdf"
#save picture
ggsave(paste(path,filename,sep=""), width = 9, height = 6, unit="cm")


#set basic layout of graph
gp <- ggplot(data = subset(meanValuesGroup, handPosition == 'uncrossed' & modalityCondition == 'bimodal'),
             aes(x=group, y=visualFirstResp.M - 0.5))
#set theme
gp <- gp + theme_bw(9)
#plot mean values per condition as point
gp <- gp + geom_errorbar(aes(ymax=visualFirstResp.M + visualFirstResp.SE - 0.5,
                             ymin=visualFirstResp.M - visualFirstResp.SE - 0.5),
                         width=0.5, size = 0.15,  alpha = 1, colour = 'grey20')
#plot mean values per condition as bar
gp <- gp + geom_bar(stat="identity", position="dodge", colour = '#B68F8B', fill = '#B68F8B')
#set y-axis label
gp <- gp + ylab("P('visual 1st' - response)")
#set x-axis label
gp <- gp + theme(axis.title.x = element_blank())
gp <- gp + theme(axis.title.y = element_text(size = 9, colour = 'black'))
gp <- gp + theme(axis.text.x = element_text(size = 8))
gp <- gp + theme(strip.text = element_text(size = 8))
#axis limits
gp <- gp + coord_cartesian(ylim=c(0.37,0.63)-0.5)
#axis ticks
gp <- gp + scale_y_continuous(breaks=seq(-0.1,0.1,0.1), labels = seq(-0.1,0.1,0.1) + 0.5)
#remove lines and shades
gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank())
gp <- gp + theme(axis.ticks = element_blank(), panel.grid.major.x = element_blank())
gp <- gp + theme(panel.spacing = unit(0.75, "lines"))
#show plot
gp
#set filename
filename <- "/results/expIndia_TV/expIndia_TV_TOJ_visualFirst_Fig_2.pdf"
#save picture
ggsave(paste(path,filename,sep=""), width = 5.25, height = 5, unit="cm")

#set basic layout of graph
gp <- ggplot(data = subset(meanValuesGroup, handPosition == 'uncrossed'),
             aes(x=group, y=accuracy.M, group = modalityCondition, fill = modalityCondition, colour = modalityCondition))
#set theme
gp <- gp + theme_bw(9)
#plot mean values per condition as point
gp <- gp + geom_errorbar(aes(ymax=accuracy.M+accuracy.SE,
                             ymin=accuracy.M-accuracy.SE),
                         width=0.5, size = 0.15,
                         alpha = 1, colour = 'grey20',
                         position = position_dodge(width = 0.9))
#split by modality condition
gp <- gp + facet_grid(modalityCondition ~.)
#plot mean values per condition as point
#gp <- gp + geom_bar(stat="identity", position="dodge", colour = '#A5A4AE', fill = '#A5A4AE')
gp <- gp + geom_bar(stat="identity", position="dodge", size = 0.25)
gp <- gp + scale_fill_manual(values=myColours.TV, name = "Modality Condition")
gp <- gp + scale_colour_manual(values=myColours.TV, name = "Modality Condition")
#set y-axis label
gp <- gp + ylab("P(correct response)")
#set x-axis label
gp <- gp + theme(axis.title.x = element_blank())
gp <- gp + theme(axis.title.y = element_text(size = 9, colour = 'black'))
gp <- gp + theme(axis.text.x = element_text(size = 8))
gp <- gp + theme(strip.text = element_text(size = 8))
#axis limits
gp <- gp + coord_cartesian(ylim=c(0.5,0.975))
#axis ticks
gp <- gp + scale_y_continuous(breaks=seq(0.5,0.9,0.25))
#remove lines and shades
gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank())
gp <- gp + theme(axis.ticks = element_blank(), panel.grid.major.x = element_blank())
gp <- gp + theme(panel.spacing = unit(0.5, "lines"))
gp <- gp + theme(legend.position = 'none')
#show plot
gp
#set filename
filename <- "/results/expIndia_TV/expIndia_TV_TOJ_accuracy_Fig_3.pdf"
#save picture
ggsave(paste(path,filename,sep=""), width = 5.25, height = 5, unit="cm")


for (ii in 1:length(subjList)){ 
  #set basic layout of graph
  gp <- ggplot(data = subset(meanValuesSOASubj, subNo == subjList[ii] & handPosition == 'uncrossed'), 
               aes(x=SOA, y=propFurtherRight, colour=rightModality))
  
  #plot mean values per condition as line
  gp <- gp + geom_line(data = subset(curvesSOAGroup, group == meanValuesSOASubj[subNo == subjList[ii], group]  & handPosition == 'uncrossed'), 
                       aes(x=SOA, y=propFurtherRight.M, colour=rightModality), size = 0.5, alpha = 0.3)
  #set theme
  gp <- gp + theme_bw(7)
  #plot mean values per condition as line
  gp <- gp + geom_line(size = 0.5)
  #plot mean values per condition as point
  gp <- gp + geom_point(size=0.75)
  #divide plot into condition panels
  gp <- gp + facet_grid(modalityMatch ~ group)
  #set x-axis label
  gp <- gp + xlab("SOA")
  #set y-axis label
  gp <- gp + ylab("p('right side first' - response)")
  #axis limits
  gp <- gp + coord_cartesian(xlim=c(-450,450), ylim=c(-0.1,1.1))
  #axis ticks
  gp <- gp + scale_y_continuous(breaks=seq(0.0,1.0,0.5))
  gp <- gp + scale_x_continuous(lim=c(-450,450), breaks=seq(-300,300,300))
  #set colours
  gp <- gp + scale_colour_manual(values=myColours.TV[c(3,2)], name = "Modality at right side")
  #legend
  gp <- gp + guides(colour=guide_legend(title='Modality presented at right side'))
  #remove lines and shades
  gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank(), panel.grid.minor = element_blank())
  gp <- gp + theme(axis.ticks = element_blank())
  gp <- gp + theme(strip.text = element_text(size = rel(1.0)))
  gp <- gp + theme(legend.position = "none",
                   legend.box.background = element_blank(),
                   legend.box.margin = margin(0, 0, 0, 0),
                   legend.box.spacing = unit(0.2, 'lines'),
                   legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = 'lines'),
                   legend.key = element_blank(),
                   legend.key.height = unit(.95, 'lines'))
  #show plot
  gp
  #set filename
  filename=paste("/results/expIndia_TV/singleSubjects/expIndia_TV_TOJ_PsychPhysFunct_", subjList[ii], ".pdf", sep='')
  #save picture
  ggsave(paste(path,filename,sep=""), width = 6, height = 5, unit="cm")
  
  #set basic layout of graph
  gp <- ggplot(data = subset(meanValuesSubj, modalityCondition == 'bimodal' & subNo == subjList[ii]  & handPosition == 'uncrossed'),
               aes(x=group, y=visualFirstResp - 0.5))
  #set theme
  gp <- gp + theme_bw(8)
  #plot mean values per condition as point
  gp <- gp + geom_bar(stat="identity", position="dodge", colour = '#B68F8B', fill = '#B68F8B')
  #set y-axis label
  gp <- gp + ylab("P('visual 1st' - response)")
  #set x-axis label
  gp <- gp + theme(axis.title.x = element_blank())
  gp <- gp + theme(axis.title.y = element_text(size = 8, colour = 'black'))
  gp <- gp + theme(axis.text.x = element_text(size = 7))
  #axis limits
  gp <- gp + coord_cartesian(ylim=c(0.37,0.63)-0.5)
  #axis ticks
  gp <- gp + scale_y_continuous(breaks=seq(-0.1,0.1,0.1), labels = seq(-0.1,0.1,0.1) + 0.5)
  #remove lines and shades
  gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank())
  gp <- gp + theme(axis.ticks = element_blank(), panel.grid.major.x = element_blank())
  gp <- gp + theme(panel.spacing = unit(0.75, "lines"))
  #show plot
  gp
  #set filename
  filename <- paste("/results/expIndia_TV/singleSubjects/expIndia_TV_TOJ_visualFirst_", subjList[ii], ".pdf", sep='')
  #save picture
  ggsave(paste(path,filename,sep=""), width = 4.2, height = 3.8, unit="cm")
  
  #set basic layout of graph
  gp <- ggplot(data = subset(meanValuesSubj, subNo == subjList[ii] & handPosition == 'uncrossed'),
               aes(x=group, y=accuracy, group = modalityCondition, fill = modalityCondition, colour = modalityCondition))
  #set theme
  gp <- gp + theme_bw(8)
  #split by modality condition
  gp <- gp + facet_grid(modalityCondition ~.)
  #plot mean values per condition as point
  #gp <- gp + geom_bar(stat="identity", position="dodge", colour = '#A5A4AE', fill = '#A5A4AE')
  gp <- gp + geom_bar(stat="identity", position="dodge", size = 0.25)
  gp <- gp + scale_fill_manual(values=myColours.TV, name = "Modality Condition")
  gp <- gp + scale_colour_manual(values=myColours.TV, name = "Modality Condition")
  #set y-axis label
  gp <- gp + ylab("P(correct response)")
  #set x-axis label
  gp <- gp + theme(axis.title.x = element_blank())
  gp <- gp + theme(axis.title.y = element_text(size = 8, colour = 'black'))
  gp <- gp + theme(axis.text.x = element_text(size = rel(1.0)))
  #axis limits
  gp <- gp + coord_cartesian(ylim=c(0.5,0.95))
  #axis ticks
  gp <- gp + scale_y_continuous(breaks=seq(0.5,1.,0.25))
  #remove lines and shades
  gp <- gp + theme(panel.border = element_blank(), strip.background = element_blank())
  gp <- gp + theme(axis.ticks = element_blank(), panel.grid.major.x = element_blank())
  gp <- gp + theme(panel.spacing = unit(0.05, "lines"))
  gp <- gp + theme(legend.position = 'none')
  #show plot
  gp
  #set filename
  filename <- paste("/results/expIndia_TV/singleSubjects/expIndia_TV_TOJ_accuracy", subjList[ii], ".pdf", sep='')
  #save picture
  ggsave(paste(path,filename,sep=""), width = 4.2, height = 3.8, unit="cm")

  }


##### data analysis #####
print("##### data analysis #####")

#calculate glmm on accuracy data
visualFirst.GLMM <- glmer(visualFirstResp ~ group + (1|subNo),
                          dataSet[modalityMatch == 'bimodal' & action == 'keep' &
                                    !is.na(visualFirstResp) & abs(SOA) < 500 & handPosition == 'uncrossed',],
                          family="binomial", control=glmerControl(optimizer="bobyqa"))
print(Anova(visualFirst.GLMM))
#do post-hoc tests to resolve group differences
visualFirst.GLMM.interactions.I <- as.data.table(
  testInteractions(visualFirst.GLMM, pairwise=c("group"), adjust="none"),
  keep.rownames = 'comparison')
#split comparison column
visualFirst.GLMM.interactions.I[,c("groupDifference") := tstrsplit(comparison, ":", fixed=TRUE)]
#trim the group variable, round others
visualFirst.GLMM.interactions.I[, ':=' (groupDifference = trim(groupDifference), Value = round(Value, 2),
                                        Chisq = round(Chisq, 2), pValue = round(`Pr(>Chisq)`, 3)), key = groupDifference]
#add significance stars
visualFirst.GLMM.interactions.I[, ':=' (signif = ifelse(pValue < 0.001, '***',
                                                        ifelse(pValue < 0.01, '**',
                                                               ifelse(pValue < 0.05, '*',
                                                                      ifelse(pValue < 0.1, '.',
                                                                             ' ')))))]
#loose unintended comparisons and select variables
print(visualFirst.GLMM.interactions.I[is.element(groupDifference, c("CC-MCC", "DC-MDC")),
                                      c('groupDifference', 'Value',
                                        'Chisq', 'pValue', 'signif')])
#do post-hoc tests to resolve group differences
visualFirst.GLMM.interactions.II <- as.data.table(
  testInteractions(visualFirst.GLMM, fixed = c("group"), adjust="none"),
  keep.rownames = 'comparison')
#split comparison column
visualFirst.GLMM.interactions.II[,c("group") := tstrsplit(comparison, ":", fixed=TRUE)]
#trim the group variable, round others
visualFirst.GLMM.interactions.II[, ':=' (group = trim(group), Value = round(Value, 2),
                                         Chisq = round(Chisq, 2), pValue = round(`Pr(>Chisq)`, 3)), key = group]
#add significance stars
visualFirst.GLMM.interactions.II[, ':=' (signif = ifelse(pValue < 0.001, '***',
                                                         ifelse(pValue < 0.01, '**',
                                                                ifelse(pValue < 0.05, '*',
                                                                       ifelse(pValue < 0.1, '.',
                                                                              ' ')))))]
#loose unintended comparisons and select variables
print(visualFirst.GLMM.interactions.II[, c('group', 'Value', 'Chisq', 'pValue', 'signif')])

#calculate glmm on accuracy data
accuracy.GLMM <- glmer(accuracy ~ group * modalityCondition + (1|subNo),
                       droplevels(dataSet[action == 'keep' & !is.na(accuracy) & abs(SOA) < 500 & handPosition == 'uncrossed',]),
                       family="binomial", control=glmerControl(optimizer="bobyqa"))
#print results per factor
print(Anova(accuracy.GLMM))
#do post-hoc tests to resolve interaction
accuracy.GLMM.interactions.I <- as.data.table(
  testInteractions(accuracy.GLMM, pairwise=c("group","modalityCondition"), adjust="none"),
  keep.rownames = 'comparison')
#split comparison column
accuracy.GLMM.interactions.I[,c("groupDifference", "modalityDifference") := tstrsplit(comparison, ":", fixed=TRUE)]
#trim the group variable, round others
accuracy.GLMM.interactions.I[, ':=' (groupDifference = trim(groupDifference), Value = round(Value, 2),
                                     Chisq = round(Chisq, 2), pValue = round(`Pr(>Chisq)`, 3)), key = groupDifference]
#add significance stars
accuracy.GLMM.interactions.I[, ':=' (signif = ifelse(pValue < 0.001, '***',
                                                     ifelse(pValue < 0.01, '**',
                                                            ifelse(pValue < 0.05, '*',
                                                                   ifelse(pValue < 0.1, '.',
                                                                          ' ')))))]
#loose unintended comparisons and select variables
print(accuracy.GLMM.interactions.I[is.element(groupDifference, c("CC-MCC", "DC-MDC")),
                                   c('groupDifference', 'modalityDifference', 'Value',
                                     'Chisq', 'pValue', 'signif')])
#do post-hoc tests to resolve single group differences
accuracy.GLMM.interactions.II <- as.data.table(
  testInteractions(accuracy.GLMM, pairwise=c("group"), fixed = c("modalityCondition"), adjust="none"),
  keep.rownames = 'comparison')
#split comparison column
accuracy.GLMM.interactions.II[,c("groupDifference", "modality") := tstrsplit(comparison, ":", fixed=TRUE)]
#trim the group variable, round others
accuracy.GLMM.interactions.II[, ':=' (groupDifference = trim(groupDifference), Value = round(Value, 2),
                                      Chisq = round(Chisq, 2), pValue = round(`Pr(>Chisq)`, 3)), key = groupDifference]
#add significance stars
accuracy.GLMM.interactions.II[, ':=' (signif = ifelse(pValue < 0.001, '***',
                                                      ifelse(pValue < 0.01, '**',
                                                             ifelse(pValue < 0.05, '*',
                                                                    ifelse(pValue < 0.1, '.',
                                                                           ' ')))))]
#loose unintended comparisons and select variables
print(accuracy.GLMM.interactions.II[is.element(groupDifference, c("CC-MCC")),
                                    c('groupDifference', 'modality', 'Value',
                                      'Chisq', 'pValue', 'signif')])

#do post-hoc tests to resolve single group differences
accuracy.GLMM.interactions.III <- as.data.table(
  testInteractions(accuracy.GLMM, pairwise=c("group"), adjust="none"),
  keep.rownames = 'comparison')
#rename comparison column
accuracy.GLMM.interactions.III[,"groupDifference":= comparison]
#trim the group variable, round others
accuracy.GLMM.interactions.III[, ':=' (groupDifference = trim(groupDifference), Value = round(Value, 2),
                                      Chisq = round(Chisq, 2), pValue = round(`Pr(>Chisq)`, 3)), key = groupDifference]
#add significance stars
accuracy.GLMM.interactions.III[, ':=' (signif = ifelse(pValue < 0.001, '***',
                                                      ifelse(pValue < 0.01, '**',
                                                             ifelse(pValue < 0.05, '*',
                                                                    ifelse(pValue < 0.1, '.',
                                                                           ' ')))))]
#loose unintended comparisons and select variables
print(accuracy.GLMM.interactions.III[is.element(groupDifference, c("DC-MDC")),
                                    c('groupDifference', 'Value',
                                      'Chisq', 'pValue', 'signif')])


