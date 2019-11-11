# clears workspace:  
rm(list=ls(all=TRUE)) 

#### functions ####

#function to detach all packages (to avoid interferences)
detachAllPackages <- function() { 
  #store basic packages names in a list
  basic.packages <- c("package:stats","package:graphics","package:grDevices",
                      "package:utils","package:datasets","package:methods","package:base")  
  #make list of all loaded packages
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  #remove basic packages from the list
  package.list <- setdiff(package.list,basic.packages)
  #remove all packages from the list
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
#detach all packages
detachAllPackages()

#### load requiered libraries ####
#better data.frames
library(data.table)
#partial correlation
library(ppcor)
#pretty plots
library(ggplot2)
#for arrangeGrob
library(gridExtra)

###### general settings #######
#save path to main directory
path="~/Documents/Experimente/expXV"

#store colours
myColours <- c('#e0e1e0', '#7897AA','#104c74','#e0e1e0', '#B68F8B','#8c3d35')

# read mean values from a txt file
dataSet.AV <- fread(file = paste(path, '/data/expIndia_AV/aggregated/expXV_India_AV_dataSet.txt', sep = ''),
                           sep = '\t', header = T)

# read mean values from a txt file
dataSet.TV <- fread(file = paste(path, '/data/expIndia_TV/aggregated/expXV_India_TV_dataSet.txt', sep = ''),
                           sep = '\t', header = T)

# combine datasets
dataSet <- rbind(dataSet.AV, dataSet.TV, fill = T)


# read mean values from a txt file
meanValuesSubj.AV <- fread(file = paste(path, '/data/expIndia_AV/aggregated/expXV_India_AV_meanValues.txt', sep = ''),
            sep = '\t', header = T)

# read mean values from a txt file
meanValuesSubj.TV <- fread(file = paste(path, '/data/expIndia_TV/aggregated/expXV_India_TV_meanValues.txt', sep = ''),
                           sep = '\t', header = T)

# combine datasets
meanValuesSubj <- rbind(meanValuesSubj.AV, meanValuesSubj.TV, fill = T)

# set new factors
meanValuesSubj[, modalityCondition2 := paste(modalityCondition, subversion, sep='-')]

# set order of factor levels
meanValuesSubj[, modalityCondition2 := factor(modalityCondition2, 
                                              levels = c("visual-VA", "bimodal-VA", "auditory-VA", 
                                                         "visual-VT", "bimodal-VT", "tactile-VT"))]

# set new factors
meanValuesSubj[is.element(modalityCondition, c('visual', 'bimodal')), modalityCondition3 := modalityCondition]
meanValuesSubj[!is.element(modalityCondition, c('visual', 'bimodal')), modalityCondition3 := 'otherModality']
# set order of factor levels
meanValuesSubj[, modalityCondition3 := factor(modalityCondition3, 
                                              levels = c("visual", "bimodal", "otherModality"))]
# sort data table
setkey(meanValuesSubj, modalityCondition3, subversion)

# set order of factor levels
meanValuesSubj[, subversion := factor(subversion, 
                                              levels = c("VA", "VT"), labels = c('Expt. 1 - VA', 'Expt. 2 - VT'))]

#calculate accuracy and predictor residuals (for partial regression plots)
meanValuesSubj[(group == 'CC') & action == 'keep' & !is.na(surgery), 
               ':=' (accuracy.age = lm(accuracy ~ age)$residuals,
                     surgery.age = lm(surgery ~ age)$residuals,
                     durationSight.age = lm(durationSight ~ age)$residuals),
               by = .(modalityCondition, subversion)]
#calculate accuracy and predictor residuals (for partial regression plots)
meanValuesSubj[(group == 'CC') & action == 'keep' & !is.na(surgery) & modalityCondition == 'bimodal', 
               ':=' (visualFirstResp.age = lm(visualFirstResp ~ age)$residuals),
               by = .(modalityCondition, subversion)]

# function for partial correlation changed to work well with data tables
correlation <- function(x, y){
  tmp <- cor.test(x, y, method = 'spearman')
  return(list(rValue = round(tmp$estimate,2), pValue = round(tmp$p.value,3)))
}


# calculate partial correlation and store in big data.table
correlations.accuracy <- rbind(
  cbind(medVar = 'logMar.betterEye', depVar = 'accuracy',
        meanValuesSubj[(group == 'CC' | group == 'CC') & action == 'keep' & !is.na(surgery), 
                       correlation(logMar.betterEye, accuracy), by = .(modalityCondition, subversion)]),
  cbind(medVar = 'surgery.age', depVar = 'accuracy.age',
        meanValuesSubj[group == 'CC' & action == 'keep' & !is.na(surgery), 
                       correlation(surgery.age, accuracy.age), by = .(modalityCondition, subversion)]),
  cbind(medVar = 'durationSight.age', depVar = 'accuracy.age',
        meanValuesSubj[group == 'CC' & action == 'keep' & !is.na(surgery), 
                       correlation(durationSight.age, accuracy.age), by = .(modalityCondition, subversion)])
)
# adjust p-values
correlations.accuracy[, pAdjust := round(p.adjust(pValue, method = 'BH'), 3), by = .(modalityCondition, subversion)]
print(correlations.accuracy)



# calculate partial correlation and store in big data.table
correlations.visualFirst <- rbind(
  cbind(medVar = 'logMar.betterEye', depVar = 'visualFirstResp',
        meanValuesSubj[(group == 'CC' | group == 'CC') & action == 'keep' & modalityCondition == 'bimodal' & !is.na(surgery), 
                       correlation(logMar.betterEye, visualFirstResp), 
                       by = .(modalityCondition, subversion, group)]),
  cbind(medVar = 'surgery.age', depVar = 'visualFirstResp.age',
        meanValuesSubj[group == 'CC' & action == 'keep' & modalityCondition == 'bimodal', 
                       correlation(surgery.age, visualFirstResp.age), 
                       by = .(modalityCondition, subversion, group)]),
  cbind(medVar = 'durationSight.age', depVar = 'visualFirstResp.age', 
        meanValuesSubj[group == 'CC' & action == 'keep' & modalityCondition == 'bimodal', 
                       correlation(durationSight.age, visualFirstResp.age), 
                       by = .(modalityCondition, subversion, group)])
)
# adjust p-values
correlations.visualFirst[, pAdjust := round(p.adjust(pValue, method = 'BH'), 3), by = .(modalityCondition, subversion, group)]
print(correlations.visualFirst)


#set basic layout of graph
gp.surgery.visualFirstResp.cor <- ggplot(droplevels(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal')), 
                                  aes(x = surgery, 
                                      y = visualFirstResp,
                                      colour = modalityCondition2,
                                      group = modalityCondition2))
#set theme and font size
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + geom_point(size = 2)
#set y-axis label
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + ylab("P('visual first')")
#set x-axis label
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + xlab("Duration visual deprivation (months)")
# #adjust scaling and breaks
# gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(axis.ticks = element_blank())
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(panel.border = element_blank(), 
                                                           panel.grid.minor = element_blank())
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(strip.background = element_blank())
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(panel.spacing = unit(2, 'mm'))
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(plot.margin = margin(1.5,2,2,0.5, 'mm'))
#change legend
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + theme(legend.position = 'none')
#split into panels
gp.surgery.visualFirstResp.cor <- gp.surgery.visualFirstResp.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.surgery.visualFirstResp.cor

#set basic layout of graph
gp.durationSight.visualFirstResp.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal'), 
                                        aes(x = durationSight, 
                                            y = visualFirstResp,
                                            colour = modalityCondition2,
                                            group = modalityCondition2))
#set theme and font size
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + geom_point(size = 2)
#set y-axis label
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + ylab("P('visual first')")
#set x-axis label
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + xlab("Time passed since surgery (months)")
# #adjust scaling and breaks
# gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(axis.ticks = element_blank())
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(panel.border = element_blank(), 
                                                                       panel.grid.minor = element_blank())
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(strip.background = element_blank())
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(panel.spacing = unit(2, 'mm'))
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + theme(legend.position = 'none')
#split into panels
gp.durationSight.visualFirstResp.cor <- gp.durationSight.visualFirstResp.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.durationSight.visualFirstResp.cor


#set basic layout of graph
gp.age.visualFirstResp.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal'), 
                              aes(x = age, 
                                  y = visualFirstResp,
                                  colour = modalityCondition2,
                                  group = modalityCondition2))
#set theme and font size
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + geom_point(size = 2)
#set y-axis label
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + ylab("P('visual first')")
#set x-axis label
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + xlab("Age (years)")
# #adjust scaling and breaks
# gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(axis.ticks = element_blank())
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(panel.border = element_blank(), 
                                                   panel.grid.minor = element_blank())
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(strip.background = element_blank())
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(panel.spacing = unit(2, 'mm'))
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + theme(legend.position = 'none')
#split into panels
gp.age.visualFirstResp.cor <- gp.age.visualFirstResp.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.age.visualFirstResp.cor

#set basic layout of graph
gp.surgery.age.visualFirstResp.age.cor <- ggplot(droplevels(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal')), 
                                         aes(x = surgery.age, 
                                             y = visualFirstResp.age,
                                             colour = modalityCondition2,
                                             group = modalityCondition2))
#set theme and font size
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + geom_point(size = 2)
#set y-axis label
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + ylab("P('visual first') | age")
#set x-axis label
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + xlab("Duration visual deprivation | age")
# #adjust scaling and breaks
# gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
# gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(axis.ticks = element_blank())
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(panel.border = element_blank(), 
                                                                         panel.grid.minor = element_blank())
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(strip.background = element_blank())
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(panel.spacing = unit(2, 'mm'))
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(plot.margin = margin(1.5,2,2,0.5, 'mm'))
#change legend
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + theme(legend.position = 'none')
#split into panels
gp.surgery.age.visualFirstResp.age.cor <- gp.surgery.age.visualFirstResp.age.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.surgery.age.visualFirstResp.age.cor

#set basic layout of graph
gp.durationSight.age.visualFirstResp.age.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal'), 
                                               aes(x = durationSight.age, 
                                                   y = visualFirstResp.age,
                                                   colour = modalityCondition2,
                                                   group = modalityCondition2))
#set theme and font size
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + geom_point(size = 2)
#set y-axis label
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + ylab("P('visual first') | age")
#set x-axis label
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + xlab("Time passed since surgery | age")
# #adjust scaling and breaks
# # gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
# gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(axis.ticks = element_blank())
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(panel.border = element_blank(), 
                                                                                     panel.grid.minor = element_blank())
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(strip.background = element_blank())
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(panel.spacing = unit(2, 'mm'))
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + theme(legend.position = 'none')
#split into panels
gp.durationSight.age.visualFirstResp.age.cor <- gp.durationSight.age.visualFirstResp.age.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.durationSight.age.visualFirstResp.age.cor

#set basic layout of graph
gp.logMar.visualFirstResp.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC' & modalityCondition == 'bimodal'), 
                                 aes(x = logMar.betterEye, 
                                     y = visualFirstResp,
                                     colour = modalityCondition2,
                                     group = modalityCondition2))
#set theme and font size
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + geom_point(size = 2)
#set y-axis label
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + ylab("P('visual first')")
#set x-axis label
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + xlab("Visual acuity of the better eye (logMar)")
# #adjust scaling and breaks
# gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(axis.ticks = element_blank())
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(panel.border = element_blank(), 
                                                         panel.grid.minor = element_blank())
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(strip.background = element_blank())
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(panel.spacing = unit(2, 'mm'))
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + theme(legend.position = 'none')
#split into panels
gp.logMar.visualFirstResp.cor <- gp.logMar.visualFirstResp.cor + facet_grid(.~subversion)
#show plot
gp.logMar.visualFirstResp.cor

#set basic layout of graph
gp.DC.logMar.visualFirstResp.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'DC' & modalityCondition == 'bimodal'), 
                                    aes(x = logMar.betterEye, 
                                        y = visualFirstResp,
                                        colour = modalityCondition2,
                                        group = modalityCondition2))
#set theme and font size
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + geom_point(size = 2)
#set y-axis label
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + ylab("P('visual first')")
#set x-axis label
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + xlab("Visual acuity of the better eye (logMar)")
# #adjust scaling and breaks
# gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + scale_colour_manual(values=myColours[c(2,5)], name = "Modality Condition")
#remove ugly lines and shades
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(axis.ticks = element_blank())
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(panel.border = element_blank(), 
                                                               panel.grid.minor = element_blank())
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(strip.background = element_blank())
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(strip.text.y = element_blank(), strip.text.x= element_text(size = rel(1.1)))
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(panel.spacing = unit(2, 'mm'))
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(plot.margin = margin(2,0.5,0.5,2, 'mm'))
#change legend
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + theme(legend.position = 'none')
#split into panels
gp.DC.logMar.visualFirstResp.cor <- gp.DC.logMar.visualFirstResp.cor + facet_grid(.~subversion)
#show plot
gp.DC.logMar.visualFirstResp.cor



#set basic layout of graph
gp.surgery.accuracy.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                   aes(x = surgery, 
                                       y = accuracy,
                                       colour = modalityCondition2,
                                       group = modalityCondition2))
#set theme and font size
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + geom_point(size = 2)
#set y-axis label
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + ylab("P('correct')")
#set x-axis label
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + xlab("Duration visual deprivation (months)")
# #adjust scaling and breaks
# gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(axis.ticks = element_blank())
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(panel.border = element_blank(), 
                                                           panel.grid.minor = element_blank())
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(strip.background = element_blank())
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(panel.spacing = unit(2, 'mm'))
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(plot.margin = margin(1.5,2,2,0.5, 'mm'))
#change legend
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + theme(legend.position = 'none')
#split into panels
gp.surgery.accuracy.cor <- gp.surgery.accuracy.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.surgery.accuracy.cor


#set basic layout of graph
gp.durationSight.accuracy.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                  aes(x = durationSight, 
                                      y = accuracy,
                                      colour = modalityCondition2,
                                      group = modalityCondition2))
#set theme and font size
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + geom_point(size = 2)
#set y-axis label
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + ylab("P('correct')")
#set x-axis label
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + xlab("Time passed since surgery (months)")
# #adjust scaling and breaks
# gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(axis.ticks = element_blank())
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(panel.border = element_blank(), 
                                                           panel.grid.minor = element_blank())
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(strip.background = element_blank())
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(panel.spacing = unit(2, 'mm'))
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + theme(legend.position = 'none')
#split into panels
gp.durationSight.accuracy.cor <- gp.durationSight.accuracy.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.durationSight.accuracy.cor



#set basic layout of graph
gp.age.accuracy.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                        aes(x = age, 
                                            y = accuracy,
                                            colour = modalityCondition2,
                                            group = modalityCondition2))
#set theme and font size
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.age.accuracy.cor <- gp.age.accuracy.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.age.accuracy.cor <- gp.age.accuracy.cor + geom_point(size = 2)
#set y-axis label
gp.age.accuracy.cor <- gp.age.accuracy.cor + ylab("P('correct')")
#set x-axis label
gp.age.accuracy.cor <- gp.age.accuracy.cor + xlab("Age (years)")
# #adjust scaling and breaks
# gp.age.accuracy.cor <- gp.age.accuracy.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.age.accuracy.cor <- gp.age.accuracy.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.age.accuracy.cor <- gp.age.accuracy.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(axis.ticks = element_blank())
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(panel.border = element_blank(), 
                                                                       panel.grid.minor = element_blank())
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(strip.background = element_blank())
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(panel.spacing = unit(2, 'mm'))
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.age.accuracy.cor <- gp.age.accuracy.cor + theme(legend.position = 'none')
#split into panels
gp.age.accuracy.cor <- gp.age.accuracy.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.age.accuracy.cor

#set basic layout of graph
gp.surgery.age.accuracy.age.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                  aes(x = surgery.age, 
                                      y = accuracy.age,
                                      colour = modalityCondition2,
                                      group = modalityCondition2))
#set theme and font size
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + geom_point(size = 2)
#set y-axis label
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + ylab("P('correct') | age")
#set x-axis label
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + xlab("Duration visual deprivation | age")
# #adjust scaling and breaks
# gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
# gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(axis.ticks = element_blank())
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(panel.border = element_blank(), 
                                                           panel.grid.minor = element_blank())
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(strip.background = element_blank())
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(panel.spacing = unit(2, 'mm'))
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(plot.margin = margin(1.5,2,2,0.5, 'mm'))
#change legend
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + theme(legend.position = 'none')
#split into panels
gp.surgery.age.accuracy.age.cor <- gp.surgery.age.accuracy.age.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.surgery.age.accuracy.age.cor


#set basic layout of graph
gp.durationSight.age.accuracy.age.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                        aes(x = durationSight.age, 
                                            y = accuracy.age,
                                            colour = modalityCondition2,
                                            group = modalityCondition2))
#set theme and font size
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + geom_point(size = 2)
#set y-axis label
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + ylab("P('correct') | age")
#set x-axis label
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + xlab("Time passed since surgery | age")
# #adjust scaling and breaks
# gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
# gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(axis.ticks = element_blank())
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(panel.border = element_blank(), 
                                                                       panel.grid.minor = element_blank())
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(strip.background = element_blank())
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(panel.spacing = unit(2, 'mm'))
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + theme(legend.position = 'none')
#split into panels
gp.durationSight.age.accuracy.age.cor <- gp.durationSight.age.accuracy.age.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.durationSight.age.accuracy.age.cor


#set basic layout of graph
gp.logMar.accuracy.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'CC'), 
                                  aes(x = logMar.betterEye, 
                                      y = accuracy,
                                      colour = modalityCondition2,
                                      group = modalityCondition2))
#set theme and font size
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + geom_point(size = 2)
#set y-axis label
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + ylab("P('correct')")
#set x-axis label
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + xlab("Visual acuity of the better eye (logMar)")
# #adjust scaling and breaks
# gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(axis.ticks = element_blank())
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(panel.border = element_blank(), 
                                                           panel.grid.minor = element_blank())
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(strip.background = element_blank())
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(strip.text.y = element_blank(), strip.text.x = element_text(size = rel(1.1)))
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(panel.spacing = unit(2, 'mm'))
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(plot.margin = margin(1.5,2,2,2, 'mm'))
#change legend
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + theme(legend.position = 'none')
#split into panels
gp.logMar.accuracy.cor <- gp.logMar.accuracy.cor + facet_grid(modalityCondition3~subversion)
#show plot
gp.logMar.accuracy.cor

#set basic layout of graph
gp.DC.logMar.accuracy.cor <- ggplot(subset(meanValuesSubj, action == 'keep' & group == 'DC'), 
                                 aes(x = logMar.betterEye, 
                                     y = accuracy,
                                     colour = modalityCondition2,
                                     group = modalityCondition2))
#set theme and font size
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme_bw(7)
#plot regression line indicating the correlation
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + geom_smooth(method = 'lm', se = T, alpha = 0.05)
#plot mean values per condition as point
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + geom_point(size = 2)
#set y-axis label
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + ylab("P('correct')")
#set x-axis label
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + xlab("Visual acuity of the better eye (logMar)")
# #adjust scaling and breaks
# gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + scale_x_continuous(limits = c(-38,38), breaks = seq(-30,30,30))
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + scale_y_continuous(limits = c(0.375,1.025), breaks = seq(0.4,1.0,0.2))
#set colours
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + scale_colour_manual(values=myColours, name = "Modality Condition")
#remove ugly lines and shades
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(axis.ticks = element_blank())
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(panel.border = element_blank(), 
                                                         panel.grid.minor = element_blank())
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(strip.background = element_blank())
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.1)))
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(panel.spacing = unit(2, 'mm'))
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(plot.margin = margin(2,0.5,0.5,2, 'mm'))
#change legend
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + theme(legend.position = 'none')
#split into panels
gp.DC.logMar.accuracy.cor <- gp.DC.logMar.accuracy.cor + facet_grid(subversion~modalityCondition3)
#show plot
gp.DC.logMar.accuracy.cor

#set filename
filename <- paste(path, "/results/expIndia_AV/expIndia_AVTV_correlations.pdf", sep = '')


# save the plot, combine both plots (without legends) using arrangeGrob and then add legend using grid.arrange
ggsave(filename, width = 18.00, height = 18, unit='cm',
       plot = grid.arrange(
         arrangeGrob(
           gp.surgery.visualFirstResp.cor + labs(tag = "A") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           gp.durationSight.visualFirstResp.cor + labs(tag = "B") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           gp.age.visualFirstResp.cor + labs(tag = "C") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           nrow=1),
         arrangeGrob(
           gp.surgery.age.visualFirstResp.age.cor + labs(tag = "D") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           gp.durationSight.age.visualFirstResp.age.cor + labs(tag = "E") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           gp.logMar.visualFirstResp.cor + labs(tag = "F") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.95)),
           nrow=1),
         arrangeGrob(
         gp.surgery.accuracy.cor + labs(tag = "G") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         gp.durationSight.accuracy.cor + labs(tag = "H") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         gp.age.accuracy.cor + labs(tag = "I") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         nrow=1),
       arrangeGrob(
         gp.surgery.age.accuracy.age.cor + labs(tag = "J") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         gp.durationSight.age.accuracy.age.cor + labs(tag = "K") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         gp.logMar.accuracy.cor + labs(tag = "L") + theme(plot.tag = element_text(face = 'bold', size = 8), plot.tag.position = c(0.99,0.975)),
         nrow=1),
       heights=c(0.45,0.45,1,1)
       )
)



