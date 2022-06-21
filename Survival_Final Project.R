library(survival)
library(readxl)
library(survminer)
library(survMisc)
library(dplyr)
library(survminer)
library(ggplot2)

data.co <- colon[seq(1,1858,2),]
data.co <- na.omit(data.co)

data.co$sex <-factor(data.co$sex)
data.co$obstruct <-factor(data.co$obstruct)
data.co$perfor   <- factor(data.co$perfor)
data.co$adhere   <- factor(data.co$adhere)
data.co$node4   <- factor(data.co$node4)
data.co$differ   <- factor(data.co$differ)
data.co$extent   <- factor(data.co$extent)
data.co$surg     <- factor(data.co$surg)
data.co$rx     <- factor(data.co$rx)

#------------------#
###Raw Analysis###
#------------------#
data.age<- 
  data.co %>% 
  mutate(age_group <- case_when (age <40 ~ "0-40",
                               age >= 40 & age < 60 ~ "40-60",
                               age >= 60 & age < 80 ~ "60-80",
                               age >= 80  ~ ">80" ))
  
#Colon cancer death
table(data.co$status)
prop.table(table(data.co$status))

#distribution of categorical variables
sex.table<-cbind(table(data.co$sex),aggregate(x=data.co$time, by=list(data.co$sex),FUN=mean)[2])
rx.table<-cbind(table(data.co$rx),aggregate(x=data.co$time, by=list(data.co$rx),FUN=mean)[2])
differ.table<-cbind(table(data.co$differ),aggregate(x=data.co$time, by=list(data.co$differ),FUN=mean)[2])
node4.table<-cbind(table(data.co$node4),aggregate(x=data.co$time, by=list(data.co$node4),FUN=mean)[2])
age.table<-cbind(table(data.age$`age_group <- ...`),aggregate(x=data.age$time, by=list(data.age$`age_group <- ...`),FUN=mean)[2])

b=boxplot(data.co$age~data.co$status,main="Statistics 1:Cancer death")

#subset original data with non-death(status=0) and death(status=1)
d_death<-subset(data.co,data.co$status==1)
d_ndeath<-subset(data.co,data.co$status==0)
table(d_death$rx)
table(d_ndeath$rx)
table(d_death$node4)
table(d_ndeath$node4)

##Descriptive Statistics sex/surg/diff
ggplot(data.co, aes(x = rx, fill = sex)) +
  geom_bar(stat="count",width=0.5,position = "dodge") +
  geom_text(stat='count',aes(label=..count..), color="black")
ggplot(data.co, aes(x = rx, fill = node4)) +
  geom_bar(stat="count",width=0.5,position = "dodge") +
  geom_text(stat='count',aes(label=..count..), color="black")
ggplot(data.co, aes(x = rx, fill = differ)) +
  geom_bar(stat="count",width=0.5,position = "dodge") +
  geom_text(stat='count',aes(label=..count..), color="black")
ggplot(data.age, aes(x = rx, fill = data.age$`age_group <- ...`)) +
  geom_bar(stat="count",width=0.5,position = "dodge") +
  geom_text(stat='count',aes(label=..count..), color="black")

#-------------------------------#
#Survival Analysis--KM Survival
#-------------------------------#

##Overall Time to recurrence KM survival and Nelson-aalen cummulative hazard 
all_surv<-survfit(Surv(time,status)~1, data = data.co)
ggsurv <- ggsurvplot(all_surv, risk.table = TRUE,ylim=c(0.4,1),xlim=c(0,3500),
                     tables.theme = clean_theme())
ggsurv

ggsurvplot(all_surv,fun="cumhaz")

##plot KM for different treatment Variable
fit.rx <- survfit(Surv(time,status)~rx, data = data.co, type = "kaplan-meier")
fit.rx
surv.rx<- summary(fit.rx)
surv.rx
par("mgp")
ggsurvplot(fit.rx,
           title="K-M Survival Estimates",
           ylab =c(expression(paste(hat(S),"(t)"))),
           xlab="Days",
           surv.median.line = "hv",  # Add median survival times (horizontal and vertical)
           legend.title="Treatment Group",
           legend.labs = c("Obs", "Lev", "Lev+5FU"),
           pval = TRUE,   # Add p-value and 95% confidence intervals
           conf.int = TRUE,
           # Add risk table
           risk.table = TRUE,
           tables.theme = theme_cleantable(),
           palette = c("green4", "orange", "red"),
           ggtheme = theme_bw()
           )
##Facet 
ggsurvplot_facet(fit.rx, 
           data = data.co,
           ylab =c(expression(paste(hat(S),"(t)"))),
           xlab="Days",
           facet.by = "sex",
           surv.median.line = "hv",  
           legend.title="Treatment Group",
           legend.labs = c("Obs", "Lev", "Lev+5FU"),
           pval = TRUE,   
           conf.int = TRUE,
           
)

##-----Transform KM for different treatment Variable

#“log”: log transformation of the survivor function
ggsurvplot(fit.rx,
           conf.int = TRUE,
           ggtheme = theme_bw(), 
           fun = "log")

#“cumhaz” plots the cumulative hazard function (f(y) = -log(y))
ggsurvplot(fit.rx,
           conf.int = TRUE,
           ggtheme = theme_bw(), 
           fun = "cumhaz")

##plot KM for sex
fit.sex <- survfit(Surv(time,status)~sex, data = data.co, type = "kaplan-meier")
fit.sex
surv.sex<- summary(fit.rx)
surv.sex
ggsurvplot(fit.sex,ylim=c(0.1,1),
           ylab =c(expression(paste(hat(S),"(t)"))),
           xlab="Days",
           surv.median.line = "hv",  
           pval = TRUE,   
           conf.int = TRUE,   
           tables.theme = theme_cleantable(),
           palette = c("green4", "orange"),
           ggtheme = theme_bw())


##plot KM for differ
fit.differ <- survfit(Surv(time,status)~differ, data = data.co, type = "kaplan-meier")
fit.differ
surv.differ<- summary(fit.differ)
surv.differ
ggsurvplot(fit.differ,
           ylab =c(expression(paste(hat(S),"(t)"))),
           xlab="Days",
           surv.median.line = "hv",  
           pval = TRUE,   
           conf.int = TRUE,   
           tables.theme = theme_cleantable(),
           palette = c("green4", "orange","red"),
           ggtheme = theme_bw())

##plot KM for node4
fit.node4 <- survfit(Surv(time,status)~node4, data = data.co, type = "kaplan-meier")
fit.node4
surv.node4<- summary(fit.node4)
surv.node4
ggsurvplot(fit.node4,
           ylab =c(expression(paste(hat(S),"(t)"))),
           xlab="Days",
           surv.median.line = "hv",  
           pval = TRUE,   
           conf.int = TRUE,   
           tables.theme = theme_cleantable(),
           palette = c("green4", "orange"),
           ggtheme = theme_bw())


#-------------------------------#
#Survival Analysis--Logrank Test
#-------------------------------#

#logrankTest on TreatmentGroup
rx.diff<-survdiff(Surv(time,status)~rx, data = data.co, rho = 0)
rx.diff
rx.diff1<-survdiff(Surv(time,status)~rx, data = data.co[-which(data.co$rx=="Lev+5FU"),],rho = 0)
rx.diff1
rx.diff2<-survdiff(Surv(time,status)~rx, data = data.co[-which(data.co$rx=="Obs"),],rho = 0)
rx.diff2
rx.diff3<-survdiff(Surv(time,status)~rx, data = data.co[-which(data.co$rx=="Lev"),],rho = 0)
rx.diff3

#logrankTest on sex /node4
sex.diff<-survdiff(Surv(time,status)~sex, data = data.co, rho = 0)
sex.diff
node4.diff<-survdiff(Surv(time,status)~node4, data = data.co, rho = 0)
node4.diff

#logrankTest on Differ Group
differ.diff<-survdiff(Surv(time,status)~differ, data = data.co, rho = 0)
differ.diff
differ.diff1<-survdiff(Surv(time,status)~differ, data = data.co[-which(data.co$differ=="1"),],rho = 0)
differ.diff1
differ.diff2<-survdiff(Surv(time,status)~differ, data = data.co[-which(data.co$differ=="2"),],rho = 0)
differ.diff2
differ.diff3<-survdiff(Surv(time,status)~differ, data = data.co[-which(data.co$differ=="3"),],rho = 0)
differ.diff3

#-----------Cox Regression------------#
##############################

##model1: Cox single cat variable:node4
model.PHnode4 <-coxph(Surv(time,status)~node4,
                      data = data.co)
summary(model.PHnode4)

##model2: Cox 2 cat variable:node4$rx
model.PH2cv <-coxph(Surv(time,status)~node4 + rx,
                   data = data.co)
summary(model.PH2cv)

##compare model 1 and model 2
anova(model.PHnode4,model.PH2cv, test = "LRT")

##model3 Cox single num variable:nodes 
model.PHnodes <-coxph(Surv(time,status)~nodes,
                      data = data.co)
summary(model.PHnodes)

##model4: Cox 2 num variable:nodes & age
model.PH2nv <-coxph(Surv(time,status)~age + nodes,
                    data = data.co)
summary(model.PH2nv)

##compare model 3 and model 4
anova(model.PHnodes,model.PH2nv, test = "LRT")

#Full cox-PH model Model 5
model.PHfull <-coxph(Surv(time,status)~rx + sex + age + obstruct + perfor + adhere 
                     +nodes + differ + extent + surg + node4,
                  data = data.co)
summary(model.PHfull)

#reduced cox model model 6
model.PHreduced <-coxph(Surv(time,status)~rx + obstruct + +nodes + extent + surg + node4,
                        data = data.co)

summary(model.PHreduced)

##Compare model 5 and model 6
anova(model.PHfull,model.PHreduced, test = "LRT")

plot(survfit(Surv(time,status)~1,data=data.co)$time, 
     survfit(Surv(time,status)~1,data=data.co)$surv,
     type="s",
     xlab="Time", 
     ylab=c(expression(paste(hat(S),"(t)"))), 
     main="KM Estimate vs Cox PH Estimate of Survival Function")

lines(survfit(model.PHreduced,type="aalen"), col="red")
legend("topright", 
       c("KM Estimate", "Cox PH Estimate"), 
       lty=1,
       col=c("black","red"), 
       bty = "n")

########Model Check#######
#Reduced model PH checking
temp <- cox.zph(model.PHreduced) 
temp
plot(temp[7])     

plot(data.co$time[data.co$status==1], 
     residuals(model.PHreduced, "schoenfeld")[,3], 
     main = "Schoenfeld residuals over time of node4"
     )
lines(lowess(data.co$time[data.co$status==1], 
             residuals(model.PHreduced, type="schoenfeld")[,3]))


cox.zph(model.PH4)

########Stratified#######

## Solution : Stratify

data.co$SurvObj <- with(data.co, Surv(time, status == 1))

res.stra <- coxph(SurvObj ~ rx + strata(obstruct) + +nodes + extent + surg + node4,
                  data =  data.co)
cox.zph(res.stra)
summary(res.stra)


