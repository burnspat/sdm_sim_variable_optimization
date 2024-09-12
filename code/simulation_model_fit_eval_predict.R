##################
# About
# fit, evaluate, and predict simulated SDMs
# arguments are passed from shell script


##################
# Source functions
# may need to change the path below
source("~/simulation_funcs.R")

# Libraries specific to this script
library(terra)



##### Import data

# get arguments passed from shell script
args <- commandArgs(TRUE)

# the list of arguments from the shell script
# $save_dir $dt $runid $pa_i $seed $method $vif $mrmr $boruta $dredge $parsimony $rast_path

# where to save the output
save_dir <- args[1]

# pre-processed pa-predictor table
DT_orig_f <- args[2]

# unique runid
runid <- args[3]

# the simulation presence-absence iteration number
pa_i <- args[4]

# the random seed
seed <- as.numeric(args[5])

# which method to use
# options are: "GLM-uni", "RF-uni", "RF-MIR", "RF-BAF", "GLM-net", "GLM", "RF"
method <- args[6]

# whether or not to use VIF for screening (TRUE/FALSE)
use_vif <- as.logical(args[7])

# whether or not to use MRMR for screening (TRUE/FALSE)
use_mrmr <- as.logical(args[8])

# whether or not to use boruta feature selection in RF models (TRUE/FALSE)
use_boruta <- as.logical(args[9])

# whether or not to use dredge in GLM
use_dredge <- as.logical(args[10])

# whether or not to use parsimony argument in RF models (TRUE/FALSE)
use_parsimony <- as.logical(args[11])

# the raster(s) to use for prediction
# should be 'NULL' if predictor rasters are not available
rast_path <- args[12]



##### Processing
cat("\n")

# need to specify model prediction type (different for RFs and GLMs)
base_meth <- unlist(strsplit(x = method, split = "-"))[1]
if (base_meth == "RF"){
  pred_type <- "prob"
} else if (base_meth == "GLM"){
  pred_type <- "response"
} else {
  cat("Check prediction type for this method...\n")
  pred_type <- 'NULL'
}

# check if fitted model exists already
save_base <- paste0(pa_i, "_run_", runid, "_seed_", seed)
save_mod <- paste0(save_dir, 'models/', save_base, ".rds")
save_eval <- paste0(save_dir, 'eval/', save_base, ".csv")


# if the model doesn't exist, go through the fitting process, else load the model
if (!file.exists(save_mod)){
  
  # set the seed
  set.seed(seed)
  
  ### prep the PA and predictors data table
  DT_orig <- read.csv(DT_orig_f, check.names = FALSE)
  pred_var_names <- colnames(DT_orig)[grepl(glob2rx("*_mean_*m|*_stdev_*m"),colnames(DT_orig))] #[85:328]
  DT <- DT_orig[,c(pa_i, pred_var_names)]
  DT <- na.omit(DT)
  DT$id<-1:nrow(DT)
  
  # reset the name of the response variable
  response_var <- "pa"
  names(DT)[1]<-response_var
  
  cat("Loaded PA + predictors table for", pa_i, "\n")
  cat("\n")
  
  
  
  ### Sampling for training and validation 25:75% -->
  cat("Training / testing split \n")
  #TRAINING dataset 75% pres and 75% abs
  pres.t<-DT[DT$pa==1,]#sampling presence
  n_pres.t<-round((0.75*nrow(pres.t)),digits=0)#number of rows
  smpl_p.t<-pres.t[sample(nrow(pres.t),size=n_pres.t,replace=F),]
  abs.t<-DT[DT$pa==0,]#sampling presence
  n_abs.t<-round((0.75*nrow(abs.t)),digits=0)#number of rows
  smpl_a.t<-abs.t[sample(nrow(abs.t),size=n_abs.t,replace=F),]
  smpl.t<-rbind(smpl_p.t,smpl_a.t) #final sampled dataset
  cat(nrow(smpl.t), "rows for training (75%) \n")
  
  #VALIDATION dataset 25% pres and 25% abs
  smpl.v<-DT[!(DT$id %in% smpl.t$id),] #delete rows used for training
  cat(nrow(smpl.v), "rows for testing (25%) \n")
  cat("\n")
  
  # drop id column
  smpl.t <- smpl.t[,!(names(smpl.t) %in% c("id"))]
  smpl.v <- smpl.v[,!(names(smpl.v) %in% c("id"))]
  
  
  ### Preliminary filtering
  # use VIF (only RF and GLM)
  if(use_vif){
    vif_label <- 'vif'
    
    # select only the predictors and remove those with little variance
    preds_df <- smpl.t[, pred_var_names]
    small_var_cols <- colnames(preds_df)[round(apply(preds_df, 2, var), 4) == 0]
    preds_df_f <- preds_df %>% select(-one_of(small_var_cols))
    
    # Apply VIF reduction
    cat("Starting feature selection with VIF...\n")
    preds_vif <- spatialRF::auto_vif(x = preds_df_f, vif.threshold = 5)
    
    vars_sel <- preds_vif$selected.variables
    
    cat("Filtered predictors using VIF. Reduced from", length(pred_var_names),
        "to", length(vars_sel), "predictors \n")
    cat("\n")
  } else {
    vif_label <- 'novif'
  }
  
  
  # use MRMR
  # possible to limit number of output variables (we choose 20)
  if (use_mrmr) {
    mrmr_label <- 'mrmr'
    
    nsv <- 20
    dd <- smpl.t
    dd$pa <- as.numeric(dd$pa)
    dd$id <- NULL
    
    cat("Starting feature selection with MRMR...\n")
    dd <- mRMR.data(data = dd)
    mi1 <- mRMR.classic(data = dd, target_indices = c(1), feature_count = nsv)
    
    vars_sel <- mi1@feature_names[as.numeric(mi1@filters[[1]])]
    
    cat("Filtered predictors using MRMR. Reduced from", length(pred_var_names),
        "to", length(vars_sel), "predictors \n")
    cat("\n")
  } else {
    mrmr_label <- 'nomrmr'
  }
  
  
  if (use_boruta) {
    bor_label <- 'bor'
    
    cat("Starting feature selection with Boruta...\n")
    boruta_output <- Boruta(pa ~ ., data=na.omit(smpl.t))
    vars_sel <- getSelectedAttributes(boruta_output)
    
    cat("Filtered predictors using Boruta. Reduced from", length(pred_var_names),
        "to", length(vars_sel), "predictors \n")
    cat("\n")
  } else {
    bor_label <- 'nobor'
  }
  
  # Specify the selected variables (should be NULL if neither screening method is used)
  if (!any(use_vif, use_mrmr, use_boruta)){
    vars_sel <- NULL
    cat("No screening applied to predictors. \n")
    cat("\n")
  }
  
  
  
  ### model fitting and evaluation
  # apply the methods - fit models and get the selected predictor variables
  
  if(method == "GLM-uni"){
    cat("Running univariate GLM with dredge... \n")
    mod_fit <- fit_glm_uni(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names)
    
  } else if (method == "RF-uni"){
    cat("Running univariate RF... \n")
    mod_fit <- fit_rf_uni(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names,
                          use_parsimony = use_parsimony)
    
  } else if (method == "RF-MIR"){
    cat("Running MIR RF... \n")
    mod_fit <- fit_rf_mir(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names,
                          use_parsimony = use_parsimony)
    
  } else if (method == "RF-BAF"){
    cat("Running MIR BAF... \n")
    mod_fit <- fit_rf_baf(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names)
    
  } else if(method == "GLM-net"){
    cat("Running GLM net... \n")
    mod_fit <- fit_glmnet_lasso(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names,
                                pred_vars_sel = vars_sel)
    
  } else if (method == "RF"){
    cat("Running RF... \n")
    mod_fit <- fit_rf(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names,
                      pred_vars_sel = vars_sel)
    
  } else if (method == "GLM"){
    cat("Running GLM... \n")
    mod_fit <- fit_glm(smpl.t = smpl.t, response_var = "pa", pred_vars = pred_var_names,
                       pred_vars_sel = vars_sel, use_dredge = use_dredge)
    
  } else{
    cat("No functionality for this method. Check the input argument. \n")
  }
  
  
  # the basename to using for saving different things
  if (use_dredge){
    dredge_label <- "dredge"
  } else {
    dredge_label <- "nodredge"
  }
  
  if (use_parsimony){
    pars_label <- "pars"
  } else {
    pars_label <- "nopars"
  }
  
  saveRDS(object = mod_fit[[1]], file = save_mod)
  cat("Saved model to", save_mod, "\n")
  cat("\n")
  
  
  ### Model evaluation
  # independent (held-out) first
  # initialize empty dfs
  df.acc<-data.frame(pa_i = pa_i,
                     method = method,
                     parsimony = use_parsimony,
                     vif = use_vif,
                     mrmr = use_mrmr,
                     dredge = use_dredge,
                     seed = seed)
  
  #select final variables from validation datasets
  if (method != "GLM-net"){
    df <- mod_fit[[2]]
    sel.var <- colnames(df[colSums(df > 0) > 0])
    df.v <- smpl.v[,c(response_var, sel.var)]
    df.t <- smpl.t[,c(response_var, sel.var)]
  } else{
    sel.var <- mod_fit[[3]]
    df.v <- as.matrix(smpl.v[,sel.var])
    df.t <- as.matrix(smpl.t[,sel.var])
  }
  
  # predict the model on the held-out test set
  m_pred<-terra::predict(mod_fit[[1]], df.v, type = pred_type)
  if (base_meth == "RF"){
    m_pred<-m_pred[,2]
  } else if (base_meth == "GLM"){
    if (method == "GLM-net"){
      # predictions are in a matrix and named
      m_pred<-unname(m_pred[,1])
    } else {
      # predictions are named
      m_pred<-unname(m_pred)
    }
  } else {
    cat("Check prediction index for this method...\n")
    m_pred <- 'NULL'
  }
  m_pred<-as.data.frame(m_pred)
  names(m_pred)[1]<-"Predicted"
  # get eval metrics
  m_val<-cbind(smpl.v[,response_var],m_pred)
  names(m_val)[1]<-"Observed"
  m_val$ID<-seq_len(nrow(m_val))# add column ID otherwise pres absence won't run
  m_val<-m_val[,c(3,1,2)]
  acc<-presence.absence.accuracy(m_val,threshold=500,st.dev=FALSE)
  maxKappa<-acc[acc$Kappa == max(acc$Kappa),] #extract max Kappa
  maxKappa<-maxKappa[1,2:7]
  maxKappa$TSS<-maxKappa$sensitivity+maxKappa$specificity-1
  names(maxKappa) <- c("threshold.v", "PCC.v", "sensitivity.v", "specificity.v", "Kappa.v", "AUC.v", "TSS.v")
  df.acc.v<-cbind(df.acc,maxKappa)
  
  
  # non-independent (training set) next
  # predict the model on the training set
  m_pred.t<-terra::predict(mod_fit[[1]], df.t, type = pred_type)
  if (base_meth == "RF"){
    m_pred.t<-m_pred.t[,2]
  } else if (base_meth == "GLM"){
    if (method == "GLM-net"){
      # predictions are in a matrix and named
      m_pred.t<-unname(m_pred.t[,1])
    } else {
      # predictions are named
      m_pred.t<-unname(m_pred.t)
    }
  } else {
    cat("Check prediction index for this method...\n")
    m_pred <- 'NULL'
  }
  m_pred.t<-as.data.frame(m_pred.t)
  names(m_pred.t)[1]<-"Predicted"
  # get eval metrics
  m_val.t<-cbind(smpl.t[,response_var],m_pred.t)
  names(m_val.t)[1]<-"Observed"
  m_val.t$ID<-seq_len(nrow(m_val.t))# need to add column ID otherwise pres absence won't run
  m_val.t<-m_val.t[,c(3,1,2)]
  acc.t<-presence.absence.accuracy(m_val.t,threshold=500,st.dev=FALSE)
  maxKappa.t<-acc.t[acc.t$Kappa == max(acc.t$Kappa),] #extract max Kappa
  maxKappa.t<-maxKappa.t[1,2:7]
  maxKappa.t$TSS.t<-maxKappa.t$sensitivity+maxKappa.t$specificity-1
  names(maxKappa.t) <- c("threshold.t", "PCC.t", "sensitivity.t", "specificity.t", "Kappa.t", "AUC.t", "TSS.t")
  
  # final results table with evaluation metrics and selected predictor vars
  # combine the independent eval, non-independent eval, and predictors together
  DT_eval <- cbind(df.acc.v, maxKappa.t, mod_fit[[2]])
  write.csv(x = DT_eval, file = save_eval, row.names = FALSE)
  cat("Saved evaluation and predictors to", save_eval, "\n")
  cat("\n")
} else {
  # load the saved model
  mod_fit <- list(readRDS(save_mod))
  
  # load the saved eval file so we can access the selected variable names for subsetting the raster stack for prediction
  save_eval_re <- read.csv(save_eval)
  save_eval_re <- save_eval_re[,grep("mean_scl", names(save_eval_re), value=TRUE)]
  sel.var <- names(save_eval_re[colSums(save_eval_re)==1])
}


### Prediction step
if (rast_path != 'NULL'){
  # load the raster stack for prediction
  rast_list <- list.files(path = rast_path, pattern=".tif$", full.names=T)
  cat("Predicting on", length(rast_list), "rasters... \n")
  
  
  # custom prediction function (mostly needed to deal with glmnet)
  predfun <- function(model, data, ...){
    
    if (base_meth == "RF"){
      m_pred <- predict(model, data, type = pred_type)
      m_pred <- m_pred[,2]
    } else if (base_meth == "GLM"){
      if (method == "GLM-net"){
        m_pred <- predict(model, newx=as.matrix(data), s="lambda.min", type = pred_type, ...)
        # predictions are in a matrix and named
        m_pred <- unname(m_pred[,1])
      } else {
        # predictions are named
        m_pred <- predict(model, data, type = pred_type)
        m_pred <- unname(m_pred)
      }
    }
    return(p=as.vector(m_pred))
  }
  
  if (length(rast_list) == 1){
    rast <- terra::rast(x = rast_list[[1]])
    new_band_names <- gsub(pattern = "mean-", replacement = "mean_scl", x = names(rast), fixed = TRUE)
    names(rast) <- new_band_names
    names(rast) <- gsub("tree-plant","treeplant",names(rast))
    names(rast) <- gsub("elevation_mean", "Topo_elev",names(rast))
    rast <- terra::subset(rast, sel.var)
    rast_pred <- terra::predict(rast, model = mod_fit[[1]], fun = predfun, na.rm = TRUE)
    cat("...predicted to raster 1 \n")
    
  } else if (length(rast_list) > 1){
    # multiple rasters to predict over
    rast_pred_list <- list()
    
    for (i in 1:length(rast_list)){
      # load the raster tile
      rast_tile <- terra::rast(x = rast_list[[i]])
      new_band_names <- gsub(pattern = "mean-", replacement = "mean_scl", x = names(rast_tile), fixed = TRUE)
      names(rast_tile) <- new_band_names
      names(rast_tile) <- gsub("tree-plant","treeplant",names(rast_tile))
      names(rast_tile) <- gsub("elevation_mean", "Topo_elev_mean",names(rast_tile))
      rast_tile <- terra::subset(rast_tile, sel.var)
      cat(length(names(rast_tile)), "bands selected \n")
      rast_pred_tile <- terra::predict(rast_tile, model = mod_fit[[1]], fun = predfun, na.rm = TRUE)
      rast_pred_list[[i]] <- rast_pred_tile
      rm(rast_pred_tile)
      cat("...predicted to raster", i, "\n")
    }
    
    cat("\n")
    rast_pred <- terra::mosaic(terra::sprc(rast_pred_list))
    cat("Mosaiced raster tiles. \n")
    
  } else {
    cat("check input rasters \n")
  }
  
  # save the prediction rasters
  save_rast_pred <- paste0(save_dir, 'pred_maps/', save_base, ".tif")
  terra::writeRaster(x = rast_pred, filename = save_rast_pred,
                     overwrite = TRUE, NAflag = -9999, gdal=c("COMPRESS=LZW", "TILED=YES"))
  cat("Saved prediction raster to", save_rast_pred, "\n")
  cat("\n")
} else {
  cat("Not predicting to rasters.\n")
  cat("\n")
}
