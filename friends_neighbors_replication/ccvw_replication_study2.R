# Replication code for Study 2 in 
# Campbell, R., Cowley, P., Vivyan, N & Wagner, M. 
# 'Why friends and neighbors? Explaining the electoral appeal of local roots'.
# Journal of Politics


### set working directory and options

rm(list=ls()); gc()
setwd("~/Dropbox/MP background & focus/replication_materials")

### load packages

library(ggplot2)
library(foreign)
library(rms)
library(reshape2)
library(plyr)
library(arm)



### create useful functions

add.top <- function(df, new.level){
  to.add <- data.frame(mean = c(NA,NA, 0), ci.lo = c(NA,NA, 0), ci.hi = c(NA,NA, 0),
                       category = rep("", 3), attribute = rep(df$attribute[1],3),
                       level = c("", " ", new.level))
  return(rbind(to.add, df))
}
add.justify <- function(df){
  df$left.justify <- rep(0, nrow(df))
  df$left.justify[2] <- 1
  return(df)
}

## Function to get regression-based AMCE estimates for each attribute level using 
## OLS estimator with clustered SEs (Hainmueller, Hopkins and Yammamoto 2014)

# function called within amce.tab to get the actual amce's 

get.amcetab <- function(data, variables, J = 2){    
  
  # important to have the list of variables length to know how many AMCE's to
  # calculate
  
  Nvar <- length(variables)
  amce.list <- vector("list", length = Nvar)
  
  for(i in 1:Nvar){ # get AMCE for each variable attribute
    
    # regression on every variables with mp.preferred individually. fmla builds
    # each regression through each iteration of the loop
    
    fmla <- as.formula(paste("mp.preferred ~ ",variables[i], sep = ""))
    
    # fmla is placed into ols. Runs a linear regression based on the formula.
    # Gets same results for the regression as just calling lm but there are
    # other things in the output as well
    
    model <- ols(fmla, data = data, x = T, y = T) 
    # NOTE: The data for the model has to have no NAs on any model variables 
    # for the robcov(cluster()) function to work 
    
    #adjusts the variance-covariance matrix of a fit to correct for
    #heteroscedasticity and for correlated responses from cluster samples.
    #Clusters by ID, which is each individual person tested as they were tested
    #10 times.
    
    model.clus <- robcov(model, cluster = data$ID, method = "efron")
    
    # list of the coefficients from the model. taking out the baseline because
    # they do comparisons to the baseline
    
    coef <- model.clus$coef[names(model.clus$coef)!="Intercept"]
    
    # calculates as well as gets the standard error for each of the coefficients 
    
    se.clus <- sqrt(diag(model.clus$var))
    se.clus <- se.clus[names(se.clus)!="Intercept"]   
    
    # creating a table from the info drawn from the model 
    
    sub.tab <- data.frame("AMCE" = coef, 
                          "ci.lo" = coef - (1.96*se.clus),
                          "ci.hi" = coef + (1.96*se.clus),
                          "cluster.se" = se.clus)
    
    # making the name of each of the coefficients from the model a column in the
    # created dataframe
    
    sub.tab$category <- names(coef)
    
    # splits into two columns. Since these are categorical variables, the
    # rgression output is the name of the variable being regressed on, and equal
    # sign, and then the category. This separates the variable and the category
    # of the coefficient into two columns
    
    sub.tab <- cbind(sub.tab, colsplit(sub.tab$category, "=", c("attribute","level")))
    sub.tab$level <- as.character(sub.tab$level)    
    
    # no longer needed since made the two columns out of the variable names
    
    row.names(sub.tab) <- NULL
    
    ## add in gaps and baselines
    
    # these had been taken out when we did not want coefficient for the
    # intercept. Now it is added back into table for reference as 0's
    
    to.add <- data.frame(AMCE = c(NA,NA, 0), ci.lo = c(NA,NA, 0), ci.hi = c(NA,NA, 0),
                         cluster.se = c(NA,NA, 0),
                         category = rep("", 3), attribute = rep(sub.tab$attribute[1],3),
                         level = c("", " ", "baseline"))
    
    # since this is a for loop will do the same for each regression
    # individually. All are saved into amce.list
    
    amce.list [[i]] <- rbind(to.add, sub.tab)
  } 
  
  # create one large table from the list of individual subtables for each regression 
  
  amce.tab <- do.call("rbind", amce.list)
  
  # re-make initial labels column
  
  amce.tab$category <- paste(amce.tab$attribute, amce.tab$level, sep = ": ")
  
  # make this into ordered factor
  amce.tab$category <- factor(amce.tab$category, levels = rev(amce.tab$category), order =T)    
  
  return(amce.tab)
}


## author: Function that calls get.amcetab for multiple predictors and combines results

# inputs for the function are a dataframe, variables from dataframe wanting to
# be used to get their amce and whether there are multiple datasets, and if they
# test subjects are in same political party.

# amce.tab only function is to prepare data so as to call get_amcetab. That is
# why there are many if statements for different number of datasets and parties.

amce.tab <- function(data, variables, multi = F, same.party = F){
  # data must be a single data frame or a list of data frames (if multi = T)
  # with named elements
  # Also relies on specific ordering of explanatory variables
  
  # first try when multi is true. What is going on here is iterating over the
  # list of dataframes all calling get_amcetab for each of the dataframes in the
  # list. Necessary that all dataframes have the variables inputed into the
  # orignal function
  
  if(multi == T & same.party == F){
    amce.tab.list <- list(NA, length = length(data))
    for(i in 1:length(data)){
      tmp <- get.amcetab(data[[i]], variables = variables)
      tmp$set <- rep(names(data)[i], nrow(tmp))
      amce.tab.list[[i]] <- tmp
    }
    
    # combining all the dataframes together to create one large dataset with set
    # column to be able to tell which amce belongs to which dataset from the list.
    
    amce.tab <- do.call("rbind",amce.tab.list)
    amce.tab$set <- factor(amce.tab$set)
    return(amce.tab)
  }
  
  
  # same as before in needing to iterate over every dataset in the list of
  # datasets and call get_amcetab for all of them
  
  if(multi == T & same.party == T){
    amce.tab.list <- list(NA, length = length(data))
    for(i in 1:length(data)){
      
      # check for same party. If they are the same party. If true, the first
      # variable is taken out of the variable list. This is why author stated
      # that the order of variables is important.
      
      vars <- if(grepl("same party", names(data)[i])==T|grepl("Same Party", names(data)[i])==T) variables[2:length(variables)] else variables
      
      # regardless now calls get.amcetab for every element of the list
      
      tmp <- get.amcetab(data[[i]], variables = vars)
      tmp$set <- rep(names(data)[i], nrow(tmp))
      amce.tab.list[[i]] <- tmp
    }
    names(amce.tab.list) <- names(data)
    
    # separtes the data where they are same party and datasets where they
    # are not in same party
    diff.party <- amce.tab.list[grepl("same party", names(amce.tab.list))==F&
                                  grepl("Same Party", names(amce.tab.list))==F]
    same.party <- amce.tab.list[grepl("same party", names(amce.tab.list))==T |
                                  grepl("Same Party", names(amce.tab.list))==T]
    
    # dataframe to fill in when they are same party 
    
    to.add <- data.frame(AMCE = rep(NA, 4), ci.lo = rep(NA, 4), ci.hi =  rep(NA, 4),
                         cluster.se =  rep(NA, 4),
                         category =  diff.party[[1]]$category[1:4], 
                         attribute =  diff.party[[1]]$attribute[1:4], 
                         level = diff.party[[1]]$level[1:4],
                         set = rep(NA, 4))
    
    for(i in 1:length(data)){ 
      
      #replacing results from amcetab when same party is true.  
      
      if(grepl("same party", names(data)[i])==T|grepl("Same Party", names(data)[i])==T){
        amce.tab.list[[i]] <- rbind(to.add,amce.tab.list[[i]])
        amce.tab.list[[i]]$set[1:4] <-  amce.tab.list[[i]]$set[5:8]
      }
      else amce.tab.list[[i]] <- amce.tab.list[[i]]
    }    
    
    # finally combines all the datasets from the list of datasets together and
    # outputs the results together.
    
    amce.tab <- do.call("rbind",amce.tab.list)
    amce.tab$set <- factor(amce.tab$set)
    return(amce.tab)
  } 
  
  # last if statement and it is the easiest just when multi is false. All
  # amce.tab needs does is just call get.amcetab
  
  if(multi == F)   {
    
    get.amcetab(data, variables)
  }
  
}



### Load in Study 2 data 

long.dat <- readRDS("./friends_neighbors_replication/study2data_long.rds")
wide.dat <- readRDS("./friends_neighbors_replication/study2data_wide.rds")


### Create labels for plotting

# for full results

# more informative labels 

labels.full <- rev(expression(
  "", italic("Party & position (baseline = Labour left-wing)"), "Labour centre", "Conservative centre", "Conservative right-wing",
  "", italic("Local roots (baseline = lives elsewhere)"),"5 years in area", "20 years in area", "Grew up and lives in area", 
  "", italic("Constituency work (baseline = 1 day)"), "2 days", "3 days", "4 days",
  "", italic("Main policy influence (baseline = party)"), "constituents' views","own personal views",
  "", italic("Policy interests (baseline = economy and tax)"), "education and health",
  "", italic("MP sex (baseline = female)"), "male"))

# for analysis of local roots only
labels.sub <- expression("Grew up and lives in area", "20 years in area", "5 years in area", 
                         italic("Local roots (baseline = lives elsewhere)"))


# label for x axis

effect.label <- "Change in probability of MP being preferred,\n relative to baseline"


### Figure 3: AMCEs for all attributes

# calling amce.tab to get the amce values for all these variables when used as
# explanatory variable for mp.preferred

res <- amce.tab(data = long.dat, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = F)

# reversing the order of the factors for the various categories of amce, unsure
# why

res$category <- factor(as.character(res$category), levels = rev(as.character(res$category)), order =T)
#write.csv(res, "amce-all.csv")# write results to csv file

# Full plot for all attributes

# since each had two spaces in between, the top one is not needed.

res <- res[2:nrow(res),] # chop off top empty layer
res <- subset(res, level != "baseline")# remove artificial 'baseline' rows

# labels created above

labels <- labels.full

# initiates the plot. The attribute column is each variable which is why every
# category within a variable has same color.

ggplot(res, aes(x = category, y = AMCE, color = attribute)) + 
  
  # baseline of 0 to show how significant each coefficient is 
  
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  
  # gives the standard error range and a point at the amce
  
  geom_pointrange(aes(ymin = ci.lo, ymax = ci.hi), size = 0.75) +
  
  # inputs better label for y axis (Will become x axis)
  labs(y = effect.label, x = "") + 
  
  # With the categories on x axis, everything is jumbled, thus flipping the
  # coordinates so the categories are legible
  
  coord_flip() + 
  
  # adding nice formating 
  
  theme_bw() + 
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 15)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1), # remove ticks and justify
        axis.title.x = element_text(size = 13, vjust = 0)) + 
  
  # change labels from those from the dataset that were abbreivated to something
  # more infomrative
  
  scale_x_discrete(labels=labels)
ggsave("figure3.eps", height = 7, width = 8)
#ggsave("figure3.png", dpi = 600, height = 7, width = 8)



### Appendix F: conjoint analysis diagnostics


### Test for Carryover effects (i.e. results by choice task)

long.dat$choicetask <- factor(long.dat$comparison)

# F-test

# all of the variables 

the.xvars <- c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")

# F-test for each variable listed above 

for(i in 1:length(the.xvars)){
  x.var <- the.xvars[i]
  cat("##", x.var, "\n\n")
  fmla <- as.formula(paste("mp.preferred ~ ", x.var, "*choicetask", sep = ""))
  model <- ols(fmla, data = long.dat, x = T, y = T)
  model.clus <- robcov(model, cluster = long.dat$ID, method = "efron")
  print(anova(model.clus, main.effect = T, ss = FALSE ), which = "subscripts"); cat("\n","\n")
}

# combining into one table 

cat("\n","\n", "## Investigate any F-tests where p<.1:", "\n", "\n")
model <- ols(mp.preferred ~ mp.policy*choicetask, data = long.dat, x = T, y = T)
model.clus <- robcov(model, cluster = long.dat$ID, method = "efron")
print(model.clus)



### Profile order effects

long.dat$mpnum.factor <- factor(long.dat$mp.num)

# F-test
the.xvars <- c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
cat("\n","\n", "## F-tests for interaction between profile order and MP attributes:", "\n", "\n")
for(i in 1:length(the.xvars)){
  x.var <- the.xvars[i]
  cat("##", x.var, "\n\n")
  fmla <- as.formula(paste("mp.preferred ~ ", x.var, "*mpnum.factor", sep = ""))
  model <- ols(fmla, data = long.dat, x = T, y = T)
  model.clus <- robcov(model, cluster = long.dat$ID, method = "efron")
  print(anova(model.clus, main.effect = T, ss = FALSE ), which = "subscripts"); cat("\n","\n")
}

## Randomization checks

# F-test
the.xvars <- c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
the.yvars <- c("gender", "age_grouped", "socialgrade4", "qualdegree")
rchecks <- expand.grid(xvar = the.xvars, yvar = the.yvars)
rchecks$df <- NA
rchecks$Fstat <- NA
rchecks$pvalue <- NA
for(i in 1:nrow(rchecks)){
  xvar <- rchecks$xvar[i]
  yvar <- rchecks$yvar[i]
  fmla <- as.formula(paste(" ~ ", xvar, " + ", yvar, sep = ""))
  mytable <- xtabs(fmla, data = long.dat)
  #ftable(mytable)
  rchecks$df[i] <- summary(mytable)$parameter
  rchecks$Fstat[i] <- round(summary(mytable)$statistic,2)
  rchecks$pvalue[i] <- round(summary(mytable)$p.value,2)
  if(summary(mytable)$p.value >= 0.05) rchecks$sig[i] <- " "
  if(summary(mytable)$p.value < 0.1 & summary(mytable)$p.value >= 0.05) rchecks$sig[i] <- "*"
  if(summary(mytable)$p.value < 0.05 & summary(mytable)$p.value >= 0.01) rchecks$sig[i] <- "**"
  if(summary(mytable)$p.value < 0.01) rchecks$sig[i] <- "***"
}

cat("\n","\n", "## F-tests for interaction between profile order and MP attributes:", "\n", "\n")
print(rchecks)

cat("\n","\n", "## Investigate any cases where chi-square p-value < 0.1:", "\n", "\n")
for(i in 1:nrow(rchecks)){
  if(rchecks$pvalue[i] < 0.1){
    xvar <- rchecks$xvar[i]
    yvar <- rchecks$yvar[i]
    fmla <- as.formula(paste(" ~ ", xvar, " + ", yvar, sep = ""))
    mytable <- xtabs(fmla, data = long.dat)
    print(round(prop.table(ftable(mytable), 1)*100 , 2))
    cat("\n\n")
  }
}






### Appendix G: local identificaton and friends & neighbors voting


### Treatment effects on measures of local identity (see p.20 of appendix)

## Create variables summarizing average local-ness of MPs viewed in experiment

# Assign 1 if MP lives locally and zero otherwise, then average over MPs a respondent observes
wide.dat$mplocalpercent <-  1 - rowMeans(wide.dat[,grep("^local_", names(wide.dat))]=="lives in <b>another part of the country</b>")

# Assign 1 if MP born and bread, zero otherwise, then average over MPs a respondent observes
wide.dat$mpbornbredpercent <-  rowMeans(wide.dat[,grep("^local_", names(wide.dat))]=="<b>grew up and lives in</b> your local area (the area within a 15-20 minute walk from your home)")

# Assign 1-4 depending on MP localism level, then average over the MPs a respondent observes
localscores.df <- wide.dat[,grep("^local_", names(wide.dat))]
for(j in 1:ncol(localscores.df)){
  localscores.df[,j] <- as.numeric(factor(localscores.df[,j], levels = c(
    "lives in <b>another part of the country</b>",
    "originally lived in another part of the country but <b>five years ago</b> moved to live in your local area (the area within a 15-20 minute walk from your home)",
    "originally lived in another part of the country but <b>twenty years ago</b> moved to live in your local area (the area within a 15-20 minute walk from your home)",
    "<b>grew up and lives in</b> your local area (the area within a 15-20 minute walk from your home)")))
}
wide.dat$mplocalscore <- rowMeans(localscores.df)


summary(mideff1a <- lm(localfeel.num ~ mplocalpercent, data = wide.dat))
summary(mideff1b <- lm(localfeel.num ~ mpbornbredpercent, data = wide.dat))
summary(mideff1c <- lm(localfeel.num ~ mplocalscore, data = wide.dat))

summary(mideff2a <- lm(localcare.num ~ mplocalpercent, data = wide.dat))
summary(mideff2b <- lm(localcare.num ~ mpbornbredpercent, data = wide.dat))
summary(mideff2c <- lm(localcare.num ~ mplocalscore, data = wide.dat))

summary(mideff4a <- lm(localyears ~ mplocalpercent, data = wide.dat))
summary(mideff4b <- lm(localyears ~ mpbornbredpercent, data = wide.dat))
summary(mideff4c <- lm(localyears ~ mplocalscore, data = wide.dat))




### Figure G1 Panel (a)

## Conditional effects by strength of local feelings (note: use 3-category rather than 4-category coding due to small number of cases in lowest local feeling group)
data.list <- list(NULL, length = length(levels(long.dat$localfeel)))
for(i in 1:length(levels(long.dat$localfeel))){
  data.list[[i]] <- long.dat[!is.na(long.dat$localfeel) & long.dat$localfeel == levels(long.dat$localfeel)[i],]
}
names(data.list) <- levels(long.dat$localfeel)
res <- amce.tab(data = data.list, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = T)
res$set <- factor(res$set, levels = rev(levels(long.dat$localfeel)))

# Plot
subres <- subset(res, !is.na(AMCE) & attribute == "mp.localroots")# remove artificial NA rows and subset to local roots attribute
labels <- labels.sub
ggplot(subres, aes(x = category, y = AMCE, fill = set)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x = category, xend = category, ymin = ci.lo, ymax = ci.hi),
                 position = position_dodge(width = - 0.6), size = 0.7) + 
  geom_point(position = position_dodge(width = - 0.6), 
             show_guide = TRUE, size = 3.5, shape = 21) +
  labs(y = effect.label, x = "") + 
  theme(axis.title = element_text(size = 12, face = "bold")) +
  coord_flip() + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1), # remove ticks and justify
        axis.title.x = element_text(vjust = 0)) +
  scale_x_discrete(labels=labels) +
  scale_fill_manual(name = "Strength of feelings about local area",  
                    values = c("White","Gray50","Black"),
                    breaks = rev(levels(res$set))) +
  theme(legend.position="bottom")
ggsave("figureG1a.png", dpi = 600, height = 5, width = 10)


### Figure G1 Panel (b)

## Conditional effects by strength of local focus
data.list <- list(NULL, length = length(levels(long.dat$localcare)))
for(i in 1:length(levels(long.dat$localcare))){
  data.list[[i]] <- long.dat[!is.na(long.dat$localcare) & long.dat$localcare == levels(long.dat$localcare)[i],]
}
names(data.list) <- levels(long.dat$localcare)
res <- amce.tab(data = data.list, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = T)
res$set <- factor(res$set, levels = rev(levels(long.dat$localcare)))

# Sub plot for mp local roots only
subres <- subset(res, !is.na(AMCE) & attribute == "mp.localroots")# remove artificial NA rows and subset to local roots attribute
labels <- labels.sub
ggplot(subres, aes(x = category, y = AMCE, fill = set)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x = category, xend = category, ymin = ci.lo, ymax = ci.hi),
                 position = position_dodge(width = - 0.6), size = 0.7) + 
  geom_point(position = position_dodge(width = - 0.6), 
             show_guide = TRUE, size = 3.5, shape = 21) +
  labs(y = effect.label, x = "") + 
  theme(axis.title = element_text(size = 12, face = "bold")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1), # remove ticks and justify
        axis.title.x = element_text(vjust = 0)) +
  scale_x_discrete(labels=labels) +
  scale_fill_manual(name = "Voter local vs national focus",  
                    values = c("White","Gray50","Black"),
                    breaks = rev(levels(res$set))) +
  theme(legend.position="bottom")
ggsave("figureG1b.png", dpi = 600, height = 5, width = 10)




### Figure G1 Panel (c)


## Conditional effects by years lived in area
data.list <- list(NULL, length = length(levels(long.dat$localyears4cat)))
for(i in 1:length(levels(long.dat$localyears4cat))){
  data.list[[i]] <- long.dat[!is.na(long.dat$localyears4cat) & long.dat$localyears4cat == levels(long.dat$localyears4cat)[i],]
}
names(data.list) <- levels(long.dat$localyears4cat)
res <- amce.tab(data = data.list, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = T)
res$set <- factor(res$set, levels = rev(levels(long.dat$localyears4cat)))

# Sub plot for mp local roots only
subres <- subset(res, !is.na(AMCE) & attribute == "mp.localroots")# remove artificial NA rows and subset to local roots attribute
labels <- labels.sub
ggplot(subres, aes(x = category, y = AMCE, fill = set)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x = category, xend = category, ymin = ci.lo, ymax = ci.hi),
                 position = position_dodge(width = - 0.6), size = 0.7) + 
  geom_point(position = position_dodge(width = - 0.6), 
             show_guide = TRUE, size = 3.5, shape = 21) +
  labs(y = effect.label, x = "") + 
  theme(axis.title = element_text(size = 12, face = "bold")) +
  coord_flip() + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1), # remove ticks and justify
        axis.title.x = element_text(vjust = 0)) +
  scale_x_discrete(labels=labels) +
  scale_fill_manual(name = "Years lived in local area",  
                    values = c("White","Gray30","Gray70","Black"),
                    breaks = rev(levels(res$set))) +
  theme(legend.position="bottom")
ggsave("figureG1c.png", dpi = 600, height = 5, width = 10)


### F-tests quoted in main text
### These test interactions between MP local roots attribute and voter characteristics

# list all variables to test
the.zvars <- c("localfeel",
               "localcare",
               "localyears4cat", "localyears10yrblocks", "localyears")

## extract F-test results and put in table
Ftable <- data.frame(zvar = the.zvars, F = NA, d.f. = NA, P = NA)
for(i in 1:length(the.zvars)){
  z <- Ftable$zvar[i]
  fmla <- as.formula(paste("mp.preferred ~ mp.localroots*", z, sep = ""))
  model <- ols(fmla, data = long.dat, x = T, y = T)  
  model.clus <- robcov(model, cluster = long.dat$ID, method = "efron")
  Ftable[i, 2:4] <- anova(model.clus, main.effect = T, ss = FALSE )[7,]
}
print(Ftable, digits= 2)


### Inspect conditioning effect of years living in local area controlling for age, social grade and qualifications

z <- "localyears4cat"
fmla <- as.formula(paste("mp.preferred ~ mp.localroots*age_grouped + 
                         mp.localroots*gender +
                         mp.localroots*qualdegree +
                         mp.localroots*socialgrade2 +
                         mp.localroots*", z, sep = ""))
model <- ols(fmla, data = long.dat, x = T, y = T)  
model.clus <- robcov(model, cluster = long.dat$ID, method = "efron")

print(model.clus)
anova(model.clus, main.effect = T, ss = FALSE )





### Appendix H1: Same- and different-party comparisons

### Figure H1

data.list <- list(NULL, length = length(unique(long.dat$same.party)))
for(i in 1:length(unique(long.dat$same.party))){
  data.list[[i]] <- long.dat[long.dat$same.party == unique(long.dat$same.party)[i],]
}
names(data.list) <- c("MPs from same party","MPs from different parties")

res <- amce.tab(data = data.list, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = T, same.party = T)
res$Comparison <- res$set

# Sub plot for mp local roots only
subres <- subset(res, !is.na(AMCE) & attribute == "mp.localroots")# remove artificial NA rows and subset to local roots attribute
labels <- labels.sub
ggplot(subres, aes(x = category, y = AMCE)) + #, color = set, shape = set)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x=category, xend=category, ymin=ci.lo, ymax=ci.hi), size = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = "white") +
  labs(y = effect.label, x = "") + 
  facet_wrap( ~ set, ncol = 1) +
  coord_flip() + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels=labels)
ggsave("figureH1.png", dpi = 600, height = 5, width = 8)
#ggsave("figureH1.eps",height = 5, width = 8)



### Appendix I: Further analysis of ideological proximity in Study 2

### Figure I1

## Conditional effects by lrgroup

long.dat$lrgroup <- factor(long.dat$lrgroup, levels = c("Labour: left of", "Labour: other", "Indifferent", "Conservatives: other", "Conservatives: right of"))

data.list <- list(NULL, length = length(levels(long.dat$lrgroup)))
for(i in 1:length(levels(long.dat$lrgroup))){
  data.list[[i]] <- long.dat[!is.na(long.dat$lrgroup) & long.dat$lrgroup == levels(long.dat$lrgroup)[i],]
}
names(data.list) <- levels(long.dat$lrgroup)
res <- amce.tab(data = data.list, 
                variables = c("mp.partypos", "mp.localroots", "mp.const", "mp.influence", "mp.policy", "mp.gender")
                , multi = T)
#res$category <- factor(as.character(res$category), levels = rev(as.character(res$category)), order =T)
res$set <- factor(res$set, levels = levels(long.dat$lrgroup))

# Sub plot for MP partypos only
subres <- subset(res, !is.na(AMCE) & attribute == "mp.partypos")# remove artificial NA rows and subset to local roots attribute
labels <- expression("Conservative right-wing", "Conservative centre",
                     "Labour centre", italic("Party & position (baseline = Labour left-wing)"))

subres$set.neat <- car:::recode(sub(":", "", as.character(subres$set)), 
                                '"Labour left of" = "Lab closest, left of Lab";
                                "Labour other" = "Lab closest, not left of Lab";
                                "Indifferent" = "Equally close to Lab and Con";
                                "Conservatives other" = "Con closest, not to right of Con";
                                "Conservatives right of" = "Con closest, right of Con"',
                                as.factor.result = TRUE, 
                                levels = c("Lab closest, left of Lab",
                                           "Lab closest, not left of Lab",
                                           "Equally close to Lab and Con",
                                           "Con closest, not to right of Con",
                                           "Con closest, right of Con"
                                ))  

ggplot(subres, aes(x = category, y = AMCE)) + #, color = set, shape = set)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x = category, xend = category, ymin = ci.lo, ymax = ci.hi),
                 position = position_dodge(width = - 0.6), size = 0.7) + 
  geom_point(position = position_dodge(width = - 0.6), 
             show_guide = TRUE, size = 3.5, shape = 21, fill = "white") +
  labs(y = effect.label, x = "") + 
  facet_wrap( ~ set.neat, ncol = 1) +
  coord_flip() + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels=labels)
ggsave("figureI1.png", dpi = 600, height = 9, width = 8)


### Figure I2

## Now get AMCE of local roots controlling for lrgroup*MP party/position interaction

the.dat <- subset(long.dat, !is.na(lrgroup))
m <- ols(mp.preferred ~ mp.localroots + mp.partypos*lrgroup, data = the.dat, x = T, y = T)
mclus <- robcov(m, cluster = the.dat$ID, method = "efron")

coef <- mclus$coef[names(mclus$coef)!="Intercept"]
se.clus <- sqrt(diag(mclus$var))
se.clus <- se.clus[names(se.clus)!="Intercept"]       
outdf <- data.frame("AMCE" = coef, 
                    "ci.lo" = coef - (1.96*se.clus),
                    "ci.hi" = coef + (1.96*se.clus),
                    "cluster.se" = se.clus)
outdf$label <- rownames(outdf)


plotdf <- outdf[grep("mp.localroots",outdf$label),] # subset to local roots effects
plotdf <- rbind(plotdf, data.frame(AMCE = 0, ci.lo= 0,ci.hi = 0, cluster.se = NA, label = "baseline"))# add empty row
plotdf$localroots <- factor(as.character(plotdf$label), 
                            levels = rev(c("baseline", "mp.localroots=5 years in area"
                                           ,"mp.localroots=20 years in area", 
                                           "mp.localroots=Grew up and lives in area")))
labels <- labels.sub

ggplot(plotdf, aes(x = localroots, y = AMCE)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x = localroots, xend = localroots, ymin = ci.lo, ymax = ci.hi),
                 position = position_dodge(width = - 0.6), size = 0.7) + 
  geom_point(position = position_dodge(width = - 0.6), 
             show_guide = TRUE, size = 3.5, shape = 21, fill = "white") +
  
  labs(y = effect.label, x = "") + 
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels=labels)
ggsave("figureI2.png", dpi = 600, height = 4, width = 8)





