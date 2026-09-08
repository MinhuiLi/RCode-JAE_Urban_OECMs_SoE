###### Select Site covariates ##################################################
################################################################################
SiteCovSel <- function(x, WindBest) {
  
  nx <- dim(x)[2]
  
  for (i in 1:nx) { 
    index <- combn(1:nx,i) 
    
    list.out <- list() # 存储每个组合生成的模型
    model.names <- character(dim(index)[2]) # 存储模型的名称
    
    for (j in 1:dim(index)[2]) {
      predSCS        <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x[,index[,j]]), obsCovs=obsCovs)
      formula_str    <- paste(names(predSCS@siteCovs), collapse = " + ") 
      list.out[j]    <- occu(as.formula(paste(WindBest, formula_str)), predSCS, se=F, control=list(maxit=1000)) 
      model.names[j] <- paste(paste('s', index[,j], sep=''), collapse=",") 
    }
    
    if (i==1) {aic.table <- aictab(list.out, modnames=model.names)}
    else {aic.table <- rbind(aic.table, aictab(list.out, modnames=model.names))}
  }
  
  aic.table2 <- aic.table[order(aic.table$AICc),]
  aic.table2 <- data.frame(aic.table2)
  names(aic.table2)[1] <- 'model'
  aic.table3 <- aic.table2[,c('model','K','LL','AICc')]
  aic.table.out <- transform(aic.table3, delta.AICc=AICc-min(AICc))
  
  top_models <- aic.table.out[aic.table.out$delta.AICc <= 2, ]$model
  predUMF <- unmarkedFrameOccu(y = y, siteCovs = x, obsCovs = obsCovs)
  preds <- list()
  for (top_model in top_models) {
    formula_str <- paste(unlist(strsplit(top_model, ",")), collapse = " + ")
    print(formula_str)
    full_formula_str = paste(WindBest, formula_str)
    pred <- occu(as.formula(full_formula_str), data = predUMF)
    preds <- append(preds, list(pred))
  }
  # weighted average model
  aic_avg <- model.avg(preds)
  avg_model <- summary(aic_avg)
  best_model <- aic.table.out[aic.table.out$delta.AICc <= 2, ]$model[1]
  
  return(list(top_models = top_models, avg_model = avg_model, best_model = best_model))
}