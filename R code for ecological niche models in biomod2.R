
library(sp)
library(raster)
library(parallel)
library(reshape)
library(ggplot2)
library(biomod2)
library(rgdal)


setwd("J:/biomod2_zaipei")
DataSpecies <- read.csv("specise_data/60_rarefied_points.csv")
head(DataSpecies)
myRespName <-'Dioscorea.alata.L.' 
myResp <- as.numeric(DataSpecies[,myRespName])
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]




myExpl = stack(
                raster( "current/bio2.tif"),
              
                raster( "current/bio4.tif"),
                
                raster( "current/bio8.tif"),
               
                raster( "current/bio10.tif"),
               
                raster( "current/bio12.tif"),
               
                raster( "current/bio15.tif"),
                
                raster( "current/bio19.tif"),
              
                raster( "current/aspect.tif"),
                raster( "current/awcclass.tif"),
                raster( "current/bulkdensity3.tif"),
                raster( "current/clay.tif"),
                raster( "current/cropland.tif"),
                raster( "current/drainage.tif"),
               
                raster( "current/gravel.tif"),
                raster( "current/hii_v2ge.tif"),
                raster( "current/organicc2.tif"),
               
                raster( "current/refdepth.tif"),
                raster( "current/sand.tif"),
                raster( "current/silt.tif"),
                raster( "current/slope.tif"),
               
                raster( "current/t_bs3.tif"),
                raster( "current/tesp.tif"),
               
                raster( "current/topece.tif"))

setwd("J:/biomod2_zaipei/biomod2/60")
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 5,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'random')

write.csv(myBiomodData@coord,file = "./PA_60.csv")



myBiomodData

plot(myBiomodData)
myBiomodOption <- BIOMOD_ModelingOptions()


getwd()
myBiomodModelOut <- BIOMOD_Modeling( 
                      myBiomodData, 
                      models = c("GLM", "GBM",  "CTA", "ANN", "FDA",  "RF","MAXENT.Phillips"),
                      models.options = myBiomodOption, 
                      NbRunEval=5,
                      DataSplit=70,
                      Prevalence=0.5,
                      VarImport=3, 
                      models.eval.meth = c('TSS','ROC','KAPPA'), 
                      SaveObj = TRUE, 
                      rescal.all.models = TRUE,
                      do.full.models = FALSE,
                      modeling.id = paste(myRespName,"FirstModeling",sep=""))



myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
write.csv(myBiomodModelEval,file = "./myBiomodModelEval_60.csv")
KAPPA=myBiomodModelEval["KAPPA","Testing.data",,,]
write.csv(KAPPA,"./KAPPA_60.CSV")
ROC=myBiomodModelEval["ROC","Testing.data",,,]
write.csv(ROC,"./ROC_60.CSV")
TSS=myBiomodModelEval["TSS","Testing.data",,,]
write.csv(TSS,"./TSS_60.CSV")
varimportan=get_variables_importance(myBiomodModelOut)
write.csv(varimportan,"./varimportan_60.csv")


myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut, 
                                       chosen.models = 'all', 
                                       em.by='all', 
                                       eval.metric = c('TSS','ROC'), 
                                       eval.metric.quality.threshold = c(0.6,0.8),
                                       prob.mean = T, 
                                       prob.cv = T, 
                                       prob.ci = T,
                                       prob.ci.alpha = 0.05, 
                                       prob.median = T, 
                                       committee.averaging = T, 
                                       prob.mean.weight = T,
                                       prob.mean.weight.decay = 'proportional' )

myBiomodEM




EM_evaluations=get_evaluations(myBiomodEM)
write.csv(EM_evaluations,"/.EM_evaluations_60.csv")
myBiomodProj <- BIOMOD_Projection( modeling.output = myBiomodModelOut, 
                                   new.env = myExpl, 
                                   proj.name = 'current', 
                                   selected.models = 'all', 
                                   binary.meth = 'TSS', 
                                   compress = 'xz', 
                                   build.clamping.mask = FALSE,
                                   output.format='.img')



myBiomodProj

plot(myBiomodProj)
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj


myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                         EM.output = myBiomodEM)


myBiomodEF
plot(myBiomodEF)








setwd("J:/biomod2_zaipei/biomod2/60/Dioscorea.alata.L./proj_current")
raster_<- stack('proj_current_Dioscorea.alata.L._ensemble.grd')
writeRaster(raster_, file="proj_current.asc", format="ascii", overwrite=TRUE, bylayer=TRUE, suffix=names(raster_))


setwd("J:/biomod2_zaipei/biomod2/60")

