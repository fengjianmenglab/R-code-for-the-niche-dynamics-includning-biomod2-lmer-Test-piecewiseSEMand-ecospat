setwd('C:/Users/Administrator/Desktop/feng-paper2020/20210506')
data=read.csv('2021plan0506.csv',header = T)
str(data)
colnames(data)
data=data[,1:15]
data=data[1:91,]

###
data$wildtimeln=scale(data$wildtimeln)
data$cultimeln=scale(data$cultimeln)
data$culsrod=scale(data$culsrod)
data$wildsrod=scale(data$wildsrod)
data$LATcul=scale(data$LATcul)
data$LATwild=scale(data$LATwild)
data$E=scale(data$E)
data$S=scale(data$S)
data$U=scale(data$U)
data$BR=scale(data$BR)
data$SI=scale(data$SI)

library(lmerTest)
#####E#####
m2=lmer(data=data, formula = E~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m2)
AIC(m2)
m3=lmer(data=data, formula = E~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m3)
AIC(m3,m2)

m5=lmer(data=data, formula = E~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm))
summary(m5)
AIC(m5,m3)

m6=lmer(data=data, formula = E~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m6)
AIC(m5,m6)

m7=lmer(data=data, formula = E~cultimeln+culsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m7)
AIC(m7,m6)

m8=lmer(data=data, formula = E~cultimeln+LATcul+LATwild+gf+(1|biogeorealm))
summary(m8)
AIC(m7,m8)

m9=lmer(data=data, formula = E~cultimeln+LATcul+LATwild+(1|biogeorealm))
summary(m9)
AIC(m9,m8)

m10=lmer(data=data, formula = E~LATcul+LATwild+(1|biogeorealm))
summary(m10)
AIC(m9,m10)

m12=lmer(data=data, formula = E~LATwild+(1|biogeorealm))  ###best model
summary(m12)
AIC(m10,m12)
###E best models
library(piecewiseSEM)
library(lme4)
library(nlme)
data(shipley)
data$E=as.numeric(data$E)
data$LATwild=as.numeric(data$LATwild)
data$biogeorealm=as.factor(data$biogeorealm)
shipley.list=list(lmer(data=data, formula = E~LATwild+(1|biogeorealm)))
shipley.psem=as.psem(shipley.list)
summary(shipley.psem, .progresssBar=F)

#####S#####
m2=lmer(data=data, formula = S~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m2)
AIC(m2)
m3=lmer(data=data, formula = S~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m3)
AIC(m3,m2)

m5=lmer(data=data, formula = S~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm))
summary(m5)
AIC(m5,m3)

m6=lmer(data=data, formula = S~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m6)
AIC(m5,m6)

m7=lmer(data=data, formula = S~cultimeln+culsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m7)
AIC(m7,m6)

m8=lmer(data=data, formula = S~cultimeln+LATcul+LATwild+gf+(1|biogeorealm))
summary(m8)
AIC(m7,m8)

m9=lmer(data=data, formula = S~cultimeln+LATcul+LATwild+(1|biogeorealm))
summary(m9)
AIC(m9,m8)

m10=lmer(data=data, formula = S~LATcul+LATwild+(1|biogeorealm))
summary(m10)
AIC(m9,m10)

m12=lmer(data=data, formula = S~LATwild+(1|biogeorealm))  ###best model
summary(m12)
AIC(m10,m12)
###S best models
shipley.list=list(lmer(data=data, formula = S~LATwild+(1|biogeorealm)))
shipley.psem=as.psem(shipley.list)
summary(shipley.psem, .progresssBar=F)

#####U#####
m1=lmer(data=data, formula = U~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m1)
AIC(m1)
m2=lmer(data=data, formula = U~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m2)
AIC(m2,m1)

m3=lmer(data=data, formula = U~cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m3)
AIC(m3,m2)

m4=lmer(data=data, formula = U~cultimeln+culsrod+LATcul+gf+use+(1|biogeorealm)+LFAO)
summary(m4)
AIC(m3,m4)

m5=lmer(data=data, formula = U~culsrod+LATcul+gf+use+(1|biogeorealm)+LFAO)
summary(m5)
AIC(m5,m4)

m6=lmer(data=data, formula = U~culsrod+LATcul+gf+(1|biogeorealm)+LFAO)
summary(m6)
AIC(m5,m6)

m7=lmer(data=data, formula = U~LATcul+gf+(1|biogeorealm)+LFAO)
summary(m7)
AIC(m7,m6)

m8=lmer(data=data, formula = U~LATcul+gf+(1|biogeorealm))
summary(m8)
AIC(m7,m8)

m9=lmer(data=data, formula = U~LATcul+(1|biogeorealm)) ##best models
summary(m9)
AIC(m9,m8)


###U best models
data$U=as.numeric(data$U)
shipley.list=list(lmer(data=data, formula = U~LATcul+(1|biogeorealm)))
shipley.psem=as.psem(shipley.list)
summary(shipley.psem, .progresssBar=F)

#####BR#####
m1=lmer(data=data, formula = BR~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m1)
AIC(m1)
m2=lmer(data=data, formula = BR~wildtimeln+cultimeln+culsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m2)
AIC(m1,m2)

m3=lmer(data=data, formula = BR~wildtimeln+cultimeln+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m3)
AIC(m3,m2)

m4=lmer(data=data, formula = BR~wildtimeln+cultimeln+LATcul+LATwild+use+(1|biogeorealm)+LFAO)
summary(m4)
AIC(m3,m4)

m5=lmer(data=data, formula = BR~cultimeln+LATcul+LATwild+use+(1|biogeorealm)+LFAO)
summary(m5)
AIC(m5,m4)

m6=lmer(data=data, formula = BR~LATcul+LATwild+use+(1|biogeorealm)+LFAO)
summary(m6)
AIC(m5,m6)

m7=lmer(data=data, formula = BR~LATcul+LATwild+(1|biogeorealm)+LFAO)
summary(m7)
AIC(m7,m6)

m8=lmer(data=data, formula = BR~LATcul+LATwild+(1|biogeorealm)) ##best
summary(m8)
AIC(m7,m8)


###BR best models
data$BR=as.numeric(data$BR)
shipley.list=list(lmer(data=data, formula = BR~LATcul+LATwild+(1|biogeorealm)))
shipley.psem=as.psem(shipley.list)
summary(shipley.psem, .progresssBar=F)

#####SI#####
m2=lmer(data=data, formula = SI~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+use+(1|biogeorealm)+LFAO)
summary(m2)
AIC(m2)
m3=lmer(data=data, formula = SI~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+(1|biogeorealm)+LFAO)
summary(m3)
AIC(m3,m2)

m4=lmer(data=data, formula = SI~wildtimeln+cultimeln+culsrod+wildsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m4)
AIC(m3,m4)

m5=lmer(data=data, formula = SI~wildtimeln+culsrod+wildsrod+LATcul+LATwild+gf+(1|biogeorealm))
summary(m5)
AIC(m5,m4)

m6=lmer(data=data, formula = SI~wildtimeln+culsrod+wildsrod+LATwild+gf+(1|biogeorealm))
summary(m6)
AIC(m5,m6)

m7=lmer(data=data, formula = SI~wildtimeln+culsrod+LATwild+gf+(1|biogeorealm))
summary(m7)
AIC(m7,m6)

m8=lmer(data=data, formula = SI~culsrod+LATwild+gf+(1|biogeorealm))
summary(m8)
AIC(m7,m8)

m9=lmer(data=data, formula = SI~LATwild+gf+(1|biogeorealm))
summary(m9)
AIC(m9,m8)

m10=lmer(data=data, formula = SI~LATwild+(1|biogeorealm)) ##best
summary(m10)
AIC(m9,m10)
###SI best models
data$SI=as.numeric(data$SI)
shipley.list=list(lmer(data=data, formula = SI~LATwild+(1|biogeorealm)))
shipley.psem=as.psem(shipley.list)
summary(shipley.psem, .progresssBar=F)
