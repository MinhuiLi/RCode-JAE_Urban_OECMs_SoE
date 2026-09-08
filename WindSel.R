###### Select Observation covariates ###########################################
################################################################################
WindSel <- function(var_names) {
  # create unmarked data frame
  pandatestD <- unmarkedFrameOccu(y=y, siteCovs=xraw, obsCovs=obsCovs)
  
  # fit models
  vars <- paste(names(pandatestD@siteCovs), collapse = " + ")
  occuw1 <- occu(as.formula(paste("~wind1~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw2 <- occu(as.formula(paste("~wind2~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw3 <- occu(as.formula(paste("~wind3~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw4 <- occu(as.formula(paste("~wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw12 <- occu(as.formula(paste("~wind1 + wind2~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw13 <- occu(as.formula(paste("~wind1 + wind3~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw14 <- occu(as.formula(paste("~wind1 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw23 <- occu(as.formula(paste("~wind2 + wind3~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw24 <- occu(as.formula(paste("~wind2 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw34 <- occu(as.formula(paste("~wind3 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw123 <- occu(as.formula(paste("~wind1 + wind2 + wind3~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw124 <- occu(as.formula(paste("~wind1 + wind2 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw234 <- occu(as.formula(paste("~wind2 + wind3 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  occuw1234 <- occu(as.formula(paste("~wind1 + wind2 + wind3 + wind4~", vars)), data = pandatestD, control = list(maxit = 1000))
  null <- occu(as.formula(paste("~1~", vars)), data = pandatestD, control = list(maxit = 1000))
  
  # model selection
  modList <- fitList(DWind1 = occuw1, DWind2 = occuw2, DWind3 = occuw3, DWind4 = occuw4,
                     DWind12 = occuw12, DWind13 = occuw13, DWind14 = occuw14, 
                     DWind23 = occuw23, DWind24 = occuw24, DWind34 = occuw34,
                     DWind123 = occuw123, DWind124 = occuw124, DWind234 = occuw234, 
                     DWind1234 = occuw1234, Null=null)
  modSel(modList, nullmod = 'Null')
}