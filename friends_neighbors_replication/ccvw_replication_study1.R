# Replication code for Study 1 in 
# Campbell, R., Cowley, P., Vivyan, N & Wagner, M. 
# 'Why friends and neighbors? Explaining the electoral appeal of local roots'.
# Journal of Politics



### set working directory and options

rm(list=ls()); gc()
setwd("~/Dropbox/MP background & focus/replication_materials")

### load packages

library(ggplot2)
library(sandwich)
library(lmtest)
library(stargazer)
library(cowplot)
library(cobalt)
library(xtable)
library(devtools)
library(margins)
library(car)



### Create useful functions

## Function to get predicted levels for different treatment combinations
predict.rob <- function(object, vcov,newdata){
  tt <- terms(object)
  if(missing(newdata)){ newdata <- x$model }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  #m.mat <- model.matrix(x$terms,data=newdata)
  m.coef <- object$coef
  fit <- as.vector(X %*% object$coef)
  se.fit <- sqrt(diag(X%*%vcov%*%t(X)))
  return(list(fit=fit,se.fit=se.fit))
}




### Load in Study 1 data

d <- readRDS("./friends_neighbors_replication/study1data.rds")



### Balance and randomization checks (Appendix C)

### Check balance across treatment groups

covariates <- subset(d, select = c(gender, agegrp, qualdegree, socgrade))

bp1 <- bal.plot(treat = d$treatgroup,  obj = covariates, var.name = "gender") + ggtitle("Distribution balance for gender") + theme_gray() + xlab("")
bp2 <- bal.plot(treat = d$treatgroup,  obj = covariates, var.name = "agegrp") + ggtitle("Distribution balance for age") + theme_gray() + xlab("")
bp3 <- bal.plot(treat = d$treatgroup,  obj = covariates, var.name = "qualdegree") + ggtitle("Distribution balance for education") + theme_gray() + xlab("")
bp4 <- bal.plot(treat = d$treatgroup,  obj = covariates, var.name = "socgrade") + ggtitle("Distribution balance for social grade") + theme_gray() + xlab("")


pcomb <- plot_grid(bp1, bp2, bp3, bp4, align = "l", nrow = 4)

pcomb
save_plot("figureC1.png", pcomb, dpi = 600, base_width = 8, base_height = 10)


### Randomization checks



checkvars <- c("gender", "agegrp", "qualdegree", "socgrade")
neatvars <- c("Gender", "Age group", "Education", "Social grade")
chisqs <- vector("list", length = length(checkvars))
for(i in 1:length(checkvars)){
  chisqs[[i]] <- data.frame("Background variable" = neatvars[i], 
                            "Chi sq." = unlist(chisq.test(table(d[,checkvars[i]], d$treatgroup))[1]),
                            "df" = unlist(chisq.test(table(d[,checkvars[i]], d$treatgroup))[2]),
                            "P-value" = unlist(chisq.test(table(d[,checkvars[i]], d$treatgroup))[3]))
}
resdf <- do.call("rbind", chisqs)
names(resdf) <- c("Respondent attribute", "Chi.sq", "df", "P-value")
rownames(resdf) <- NULL

### Now output as table
print(xtable(resdf, align = c("l", "l", "c", "c", "c"),
             digits = c(0,0, 2, 0,2),
             label = "tab:rchecks", caption = "Randomization checks"), 
      type = "html", file = "tableC1.html", include.rownames = FALSE)





### Table 2

m1 <- lm(nickminusphil ~ localtreat*behtreatsimple, data = d)

se1 <- sqrt(diag(vcovHC(m1)))

m2 <- lm(nickminusphil ~ localtreat*behtreatsimple +
           gender + agegrp + socgrade + qual
         , data = d)
se2 <- sqrt(diag(vcovHC(m2)))

m3 <- lm(nickminusphil ~ localtreat*behtreat, data = d)
se3 <- sqrt(diag(vcovHC(m3)))

m4 <- lm(nickminusphil ~ localtreat*behtreat +
           gender + agegrp + socgrade + qual
         , data = d)
se4 <- sqrt(diag(vcovHC(m4)))

## Output to table
filename <- paste0("table2.html")
stargazer(mget(paste0("m",1:4)),
          se = mget(paste0("se",1:4)),
          type = "html",
          float = FALSE,
          out = filename,
          dep.var.caption = "",
          dep.var.labels.include = FALSE,
          title = "", 
          #intercept.bottom = FALSE, intercept.top = TRUE,
          keep.stat = c("n", "rsq", "adj.rsq"),
          omit = "socgrade|agegrp|gender|qual",
          order = c("Constant", "^localtreatLocal roots$", 
                    "^behtreatsimpleBehavioural info$", 
                    "^behtreatConst. focus$", 
                    "^behtreatWestmin. focus$",
                    "^localtreatLocal roots:behtreatsimpleBehavioural info$",
                    "^localtreatLocal roots:behtreatConst. focus$",
                    "^localtreatLocal roots:behtreatWestmin. focus$"
          ),
          add.lines = list(c("Controls for voter characteristics?", 
                             rep(c("No","Yes"), 2))),
          covariate.labels = c("Intercept", "Local roots", "Behavioral localism information",
                               "Behavioral localism: High (vs. no info)",
                               "Behavioral localism: Low (vs. no info)",
                               "Local roots X Behavioral info.", 
                               "Local roots X High behavioral localism",
                               "Local roots X Low behavioral localism")
)



### Figure 1


## Avg treatment effect of local roots with const and westmihn. info
out <- summary(margins(m4, vcov = vcovHC(m4), 
                       at = list(behtreat = c("No behavioural info", "Const. focus", "Westmin. focus"))))
out <- subset(out, factor == "localtreatLocal roots")
margins.m4 <- out 


## make margins plot
margins.comb <- margins.m4
margins.comb$behtreat.neat <- car:::recode(margins.comb$behtreat, 
                                           '"No behavioural info" = "Behavioral information treatment --\\nNo information (Vignettes 1-2)";
                                           "Westmin. focus" = "Behavioral information treatment --\\nLow behavioral localism (Vignettes 5-6)";
                                           "Const. focus" = "Behavioral information treatment --\\nHigh behavioral localism (Vignettes 3-4)"')
                                           #as.factor.result = TRUE)
margins.comb$behtreat.neat <- factor(sub(" --", ":", margins.comb$behtreat.neat),
                                     levels = c("Behavioral information treatment:\nNo information (Vignettes 1-2)", 
                                                "Behavioral information treatment:\nLow behavioral localism (Vignettes 5-6)",
                                                "Behavioral information treatment:\nHigh behavioral localism (Vignettes 3-4)"
                                     ))
p1 <- ggplot(margins.comb, aes(x = factor, y = AME)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x=factor, xend=factor, ymin=lower, ymax=upper), size = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = "white") +
  labs(x = "", y = "") + 
  coord_flip() +
  facet_wrap( ~ behtreat.neat, ncol = 1) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ggtitle("(b) Effect of MP local roots treatment") 


p1


## predlevels for m4
newdf <- data.frame(expand.grid(localtreat = factor(c("No local roots", "Local roots"), 
                                                    levels = levels(d$localtreat)),
                                behtreat = factor(levels(d$behtreat), levels = levels(d$behtreat))
),
gender = factor(levels(d$gender)[which.max(table(d$gender))], 
                levels = levels(d$gender)),
agegrp = factor(levels(d$agegrp)[which.max(table(d$agegrp))], 
                levels = levels(d$agegrp)),
socgrade = factor(levels(d$socgrade)[which.max(table(d$socgrade))], 
                  levels = levels(d$socgrade)),
qual = factor(levels(d$qual)[which.max(table(d$qual))], 
              levels = levels(d$qual))
)
preds <- predict.rob(m4, vcov = vcovHC(m4), newdata = newdf)
newdf$yhat <- preds$fit
newdf$se.yhat <- preds$se.fit
newdf$lo <- newdf$yhat - (1.96*newdf$se.yhat)
newdf$hi <- newdf$yhat + (1.96*newdf$se.yhat)
predlevels.m4 <- newdf


## make predlevels plot
predlevels.comb <- predlevels.m4
predlevels.comb$behtreat.neat <- car:::recode(predlevels.comb$behtreat, 
                                              '"No behavioural info" = "Behavioral information treatment --\\nNo information (Vignettes 1-2)";
                                              "Westmin. focus" = "Behavioral information treatment --\\nLow behavioral localism (Vignettes 5-6)";
                                              "Const. focus" = "Behavioral information treatment --\\nHigh behavioral localism (Vignettes 3-4)"')
                                              #as.factor.result = TRUE)
predlevels.comb$behtreat.neat <- factor(sub(" --", ":", predlevels.comb$behtreat.neat),
                                        levels = c("Behavioral information treatment:\nNo information (Vignettes 1-2)", 
                                                   "Behavioral information treatment:\nLow behavioral localism (Vignettes 5-6)",
                                                   "Behavioral information treatment:\nHigh behavioral localism (Vignettes 3-4)"
                                        ))

p2 <- ggplot(predlevels.comb, aes(x = localtreat, y = yhat)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x=localtreat, xend=localtreat, ymin=lo, ymax=hi), size = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = "white") +
  labs(x = "", y = "") + 
  coord_flip() +
  facet_wrap( ~ behtreat.neat, ncol = 1) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  ggtitle("(a) Predicted relative rating")


## Now combine aspects of above plots into 3 x 2 plot
pcomb <- plot_grid(p2, p1, align = "h", rel_widths = c(1,0.8))
#save_plot("figure1.png", pcomb, dpi = 600, base_width = 8, base_height = 4)
save_plot("figure1.eps", pcomb, base_width = 8, base_height = 4)




### Table B1

## Now do models where DV is MP Nick rating

m5 <- lm(nickscore ~ localtreat*behtreatsimple, data = d)
se5 <- sqrt(diag(vcovHC(m5)))


m6 <- lm(nickscore ~ localtreat*behtreatsimple +
           gender + agegrp + socgrade + qual
         , data = d)
se6 <- sqrt(diag(vcovHC(m6)))


m7 <- lm(nickscore ~ localtreat*behtreat, data = d)
se7 <- sqrt(diag(vcovHC(m7)))

m8 <- lm(nickscore ~ localtreat*behtreat +
           gender + agegrp + socgrade + qual
         , data = d)
se8 <- sqrt(diag(vcovHC(m8)))

filename <- paste0("tableB1.html")
stargazer(mget(paste0("m",5:8)),
          se = mget(paste0("se",5:8)),
          type = "html",
          float = FALSE,
          out = filename,
          dep.var.caption = "",
          dep.var.labels.include = FALSE,
          title = "", 
          #intercept.bottom = FALSE, intercept.top = TRUE,
          keep.stat = c("n", "rsq", "adj.rsq"),
          omit = "socgrade|agegrp|gender|qual",
          order = c("Constant", "^localtreatLocal roots$", 
                    "^behtreatsimpleBehavioural info$", 
                    "^behtreatConst. focus$", 
                    "^behtreatWestmin. focus$",
                    "^localtreatLocal roots:behtreatsimpleBehavioural info$",
                    "^localtreatLocal roots:behtreatConst. focus$",
                    "^localtreatLocal roots:behtreatWestmin. focus$"
          ),
          add.lines = list(c("Controls for voter characteristics?", 
                             rep(c("No","Yes"), 2))),
          covariate.labels = c("Intercept", "Local roots", "Behavioral localism information",
                               "Behavioral localism: High (vs. no info)",
                               "Behavioral localism: Low (vs. no info)",
                               "Local roots X Behavioral info.", 
                               "Local roots X High behavioral localism",
                               "Local roots X Low behavioral localism")
)




### Figure B1: summarise models where DV is MP Nick rating


## Avg treatment effect of local roots with const and westmin. info
out <- summary(margins(m8, vcov = vcovHC(m8), 
                       at = list(behtreat = c("No behavioural info", "Const. focus", "Westmin. focus"))))
out <- subset(out, factor == "localtreatLocal roots")
margins.m8 <- out 


## predlevels for m8
newdf <- data.frame(expand.grid(localtreat = factor(c("No local roots", "Local roots"), 
                                                    levels = levels(d$localtreat)),
                                behtreat = factor(levels(d$behtreat), levels = levels(d$behtreat))
),
gender = factor(levels(d$gender)[which.max(table(d$gender))], 
                levels = levels(d$gender)),
agegrp = factor(levels(d$agegrp)[which.max(table(d$agegrp))], 
                levels = levels(d$agegrp)),
socgrade = factor(levels(d$socgrade)[which.max(table(d$socgrade))], 
                  levels = levels(d$socgrade)),
qual = factor(levels(d$qual)[which.max(table(d$qual))], 
              levels = levels(d$qual))
)
preds <- predict.rob(m8, vcov = vcovHC(m8), newdata = newdf)
newdf$yhat <- preds$fit
newdf$se.yhat <- preds$se.fit
newdf$lo <- newdf$yhat - (1.96*newdf$se.yhat)
newdf$hi <- newdf$yhat + (1.96*newdf$se.yhat)
predlevels.m8 <- newdf


## Now combine marginal effects and predicted levels into 3 x 2 plot (column 1 = margins, column 2 = levels)

# make margins plot

margins.comb <- margins.m8
margins.comb$behtreat.neat <- car:::recode(margins.comb$behtreat, 
                                           '"No behavioural info" = "Behavioral information treatment --\\nNo information (Vignettes 1-2)";
                                           "Westmin. focus" = "Behavioral information treatment --\\nLow behavioral localism (Vignettes 5-6)";
                                           "Const. focus" = "Behavioral information treatment --\\nHigh behavioral localism (Vignettes 3-4)"')
                                           #as.factor.result = TRUE)
margins.comb$behtreat.neat <- factor(sub(" --", ":", margins.comb$behtreat.neat),
                                     levels = c("Behavioral information treatment:\nNo information (Vignettes 1-2)", 
                                                "Behavioral information treatment:\nLow behavioral localism (Vignettes 5-6)",
                                                "Behavioral information treatment:\nHigh behavioral localism (Vignettes 3-4)"
                                     ))

p1 <- ggplot(margins.comb, aes(x = factor, y = AME)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x=factor, xend=factor, ymin=lower, ymax=upper), size = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = "white") +
  labs(x = "", y = "") + 
  coord_flip() +
  facet_wrap( ~ behtreat.neat, ncol = 1) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ggtitle("(b) Effect of MP local roots treatment") 


# make predlevels plot

predlevels.comb <- predlevels.m8
predlevels.comb$behtreat.neat <- car:::recode(predlevels.comb$behtreat, 
                                              '"No behavioural info" = "Behavioral information treatment --\\nNo information (Vignettes 1-2)";
                                              "Westmin. focus" = "Behavioral information treatment --\\nLow behavioral localism (Vignettes 5-6)";
                                              "Const. focus" = "Behavioral information treatment --\\nHigh behavioral localism (Vignettes 3-4)"')
                                              #as.factor.result = TRUE)
predlevels.comb$behtreat.neat <- factor(sub(" --", ":", predlevels.comb$behtreat.neat),
                                        levels = c("Behavioral information treatment:\nNo information (Vignettes 1-2)", 
                                                   "Behavioral information treatment:\nLow behavioral localism (Vignettes 5-6)",
                                                   "Behavioral information treatment:\nHigh behavioral localism (Vignettes 3-4)"
                                        ))

p2 <- ggplot(predlevels.comb, aes(x = localtreat, y = yhat)) + 
  #geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(x=localtreat, xend=localtreat, ymin=lo, ymax=hi), size = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = "white") +
  labs(x = "", y = "") + 
  coord_flip() +
  facet_wrap( ~ behtreat.neat, ncol = 1) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  ggtitle("(a) Predicted rating")


pcomb <- plot_grid(p2, p1, align = "h", rel_widths = c(1,0.8))
save_plot("figureB1.png", pcomb, dpi = 600, base_width = 8, base_height = 4)
#save_plot("figureB1.eps", pcomb, base_width = 8, base_height = 4)




### Table B2

## Now do models on binary DV

m9 <- lm(nickgt ~ localtreat*behtreatsimple, data = d)
se9 <- sqrt(diag(vcovHC(m9)))
coeftest(m9, vcovHC(m9))


m10 <- lm(nickgt ~ localtreat*behtreatsimple +
            gender + agegrp + socgrade + qual
          , data = d)
se10 <- sqrt(diag(vcovHC(m10)))
coeftest(m10, vcovHC(m10))



m11 <- lm(nickgt ~ localtreat*behtreat, data = d)
se11 <- sqrt(diag(vcovHC(m11)))
coeftest(m11, vcovHC(m11))



m12 <- lm(nickgt ~ localtreat*behtreat +
            gender + agegrp + socgrade + qual
          , data = d)
se12 <- sqrt(diag(vcovHC(m12)))
coeftest(m12, vcovHC(m12))



filename <- paste0("tableB2.html")
stargazer(mget(paste0("m",9:12)),
          se = mget(paste0("se",9:12)),
          type = "html",
          float = FALSE,
          out = filename,
          dep.var.caption = "",
          dep.var.labels.include = FALSE,
          title = "", 
          #intercept.bottom = FALSE, intercept.top = TRUE,
          keep.stat = c("n", "rsq", "adj.rsq"),
          omit = "socgrade|agegrp|gender|qual",
          order = c("Constant", "^localtreatLocal roots$", 
                    "^behtreatsimpleBehavioural info$", 
                    "^behtreatConst. focus$", 
                    "^behtreatWestmin. focus$",
                    "^localtreatLocal roots:behtreatsimpleBehavioural info$",
                    "^localtreatLocal roots:behtreatConst. focus$",
                    "^localtreatLocal roots:behtreatWestmin. focus$"
          ),
          add.lines = list(c("Controls for voter characteristics?", 
                             rep(c("No","Yes"), 2))),
          covariate.labels = c("Intercept", "Local roots", "Behavioral localism information",
                               "Behavioral localism: High (vs. no info)",
                               "Behavioral localism: Low (vs. no info)",
                               "Local roots X Behavioral info.", 
                               "Local roots X High behavioral localism",
                               "Local roots X Low behavioral localism")
)





### Table D1

## Now do models assessing whether effect of local roots is conditional on local political context


m13 <- lm(nickminusphil ~ localtreat*behtreat + localtreat*vote15share + behtreat*vote15share, data = d)
se13 <- sqrt(diag(vcovHC(m13)))
coeftest(m13, vcovHC(m13))


m14 <- lm(nickminusphil ~ localtreat*behtreat + localtreat*vote15share + behtreat*vote15share +
            gender + agegrp + socgrade + qual
          , data = d)
se14 <- sqrt(diag(vcovHC(m14)))
coeftest(m14, vcovHC(m14))



m15 <- lm(nickminusphil ~ localtreat*behtreat + localtreat*votewinner15 + behtreat*votewinner15, data = d)
se15 <- sqrt(diag(vcovHC(m15)))
coeftest(m15, vcovHC(m15))


m16 <- lm(nickminusphil ~ localtreat*behtreat + localtreat*votewinner15 + behtreat*votewinner15 +
            gender + agegrp + socgrade + qual
          , data = d)
se16 <- sqrt(diag(vcovHC(m16)))
coeftest(m16, vcovHC(m16))


filename <- paste0("tableD1.html")
stargazer(mget(paste0("m",c(15:16, 13:14))),
          se = mget(paste0("se",c(15:16, 13:14))),
          type = "html",
          float = FALSE,
          out = filename,
          dep.var.caption = "",
          dep.var.labels.include = FALSE,
          title = "", 
          #intercept.bottom = FALSE, intercept.top = TRUE,
          keep.stat = c("n", "rsq", "adj.rsq"),
          omit = "socgrade|agegrp|gender|qual",
          order = c("Constant", "^localtreatLocal roots$", 
                    "^behtreatConst. focus$", 
                    "^behtreatWestmin. focus$",
                    "^votewinner15",
                    "^vote15share",
                    "^localtreatLocal roots:behtreatConst. focus$",
                    "^localtreatLocal roots:behtreatWestmin. focus$",
                    "^localtreatLocal roots:votewinner15$",
                    "^behtreatConst. focus:votewinner15$",
                    "^behtreatWestmin. focus:votewinner15$"
          ),
          add.lines = list(c("Controls for voter characteristics?", 
                             rep(c("No","Yes"), 2))),
          covariate.labels = c("Intercept", "Local roots",
                               "Behavioral localism: High (vs. no info)",
                               "Behavioral localism: Low (vs. no info)",
                               "Supported party winning local seat 2015",
                               "Local vote share for 2015 vote choice",
                               "Local roots X High behavioral localism",
                               "Local roots X Low behavioral localism",
                               "Local roots X Supported winning party 2015",
                               "High behavioral localism X Supported winning party 2015",
                               "High behavioral localism X Supported winning party 2015",
                               "Local roots X Local vote for 2015 vote choice",
                               "Low behavioral localism X Local share for 2015 vote choice",
                               "Low behavioral localism X Local share for 2015 vote choice")
)

## Now print F-tests

print(anova(m15))

print(anova(m16))

print(anova(m13))

print(anova(m14))


