##################
# Packages
.libPaths(c("/projects/above_gedi/users/pburns/envs/biodiv_mod2/lib/R/library"))
library(randomForest, warn.conflicts = F, quietly = T)
library(rfUtilities, warn.conflicts = F, quietly = T)
library(spatialRF, warn.conflicts = F, quietly = T)
library(glmnet, warn.conflicts = F, quietly = T)
library(mRMRe, warn.conflicts = F, quietly = T)
library(plyr, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(PresenceAbsence, warn.conflicts = F, quietly = T)
library(MuMIn, warn.conflicts = F, quietly = T)
library(caret, warn.conflicts = F, quietly = T)
library(data.table, warn.conflicts = F, quietly = T)
library(stringr, warn.conflicts = F, quietly = T)
library(Boruta, warn.conflicts = F, quietly = T)



##################
# Functions

# fit univariate GLM models and choose the scale with the lowest AIC for each variable
# also applies correlation screening retaining variables with lowest AIC
# also dredge the global model (all combinations) and select the top model (lowest AIC)
# and extract the final list of variables (0 = not used, 1 = used)
fit_glm_uni <- function(smpl.t, response_var, pred_vars){
  v<-strsplit(pred_vars,'_scl',2)#spliting name into two columns
  v1<-data.frame(do.call(rbind,v)) #combine into one dataframe and change format
  var_names<-unique(v1[,1])
  scales<-unique(v1[,2])

  # 2. Univariate glms
  df.aic<-NULL
  for (j in pred_vars){
    m.glm <- glm(smpl.t[,1] ~ smpl.t[,j], family="binomial")
    m.aic <- data.frame(cbind(j, AIC(m.glm))) #change row number
    df.aic<-rbind(df.aic,m.aic)

  }

  # needed for saving intermediate results
  names(df.aic)[1]<-"variable_sc"
  names(df.aic)[2]<-"AIC"

  # 3. SO: Select best scales for each variable
  best.sc_fin<-as.character(NULL)

  for (k in var_names){
    m.var<-df.aic %>% dplyr::filter(grepl(k,variable_sc))
    best.sc <- m.var[,1] [m.var[,2] %in% min(m.var[,2])]
    #best.sc<-best.sc[1]##in case there are 2 scales of same OOB it will take lower scale
    best.sc_fin<-append(best.sc_fin,best.sc)
  }

  var.SO_aic<-df.aic[df.aic$variable_sc %in% best.sc_fin,] #selecting best scales for each variable

  #order variables by AIC from high to low
  var.SO_aic1<-var.SO_aic[order(var.SO_aic$AIC, decreasing=T),]
  

  ###############################################
  # 4. CORRELATION SCREENING
  ###############################################

  sel.var<-var.SO_aic1$variable_sc

  cor_full<-cor(smpl.t[,sel.var],method="pearson")
  cor<-cor(smpl.t[,sel.var],method="pearson")

  #replace 1 with 0
  for(l in 1:ncol(cor)){
    cor[l,l]<-0
  }

  cor <- data.frame(cor)
  threshold<-0.7 #correlation threshold


  while (max(abs(cor)) > threshold) {
    breakloop <- FALSE
    for (l in 1:ncol(cor)) {
      for (m in 1:ncol(cor)) {
        if (abs(cor)[l,m] >= threshold) {
          imp.x <- var.SO_aic1[m,2] #AIC
          imp.y <-var.SO_aic1[l,2]
          if (imp.x > imp.y) { #Adapt this depending if AIC or impr
            drop <- colnames(cor)[m]
          }
          else {
            drop <- rownames(cor)[l]
          }
          breakloop <- TRUE
          break
        }
      }
      if (breakloop) break
    }
    cor <- cor[row.names(cor) != drop, colnames(cor) != drop]
  }

  var.red<-names(cor) #final set of variables after correlation reduction


  #####################################################
  # 5. Dredge the global GLM
  ####################################################

  var_fml<-paste(var.red,collapse="+")
  m_fml<-as.formula(paste0(response_var,"~ ",var_fml))

  options(na.action="na.fail") #important for dredge

  glm.fit <- glm(m_fml,data=smpl.t,family="binomial") #global model

  glm.d<-dredge(glm.fit)#all subsets from global model

  glm.best<-subset(glm.d, delta==0) #best model (top)

  m.sum<-summary(get.models(glm.d,1)[[1]])$coefficient #this is a full table for best model!
  coef.p<-as.data.frame(m.sum[,c(1,4)])
  names(coef.p)[1]<-"Coefficient"
  names(coef.p)[2]<-"p"

  var.d<-row.names(coef.p)#names of final variables with Intercept
  var.d.fin<-var.d[-1]
  
  #refit the model with var.d.fin
  var_fml.fin<-paste(var.d.fin,collapse="+")
  m_fml.fin<-as.formula(paste0(response_var, "~ ",var_fml.fin))
  glm.fin <- glm(m_fml.fin,data=smpl.t,family="binomial") #final model

  df.var<-setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)
  v.bin<-(ifelse(colnames(df.var) %in% var.d.fin,"1","0")) #filling up table with model variables
  df.var[1,] <- v.bin

  return(list(glm.fin, df.var))
} # end fit_glm_uni



# fit univariate RF
fit_rf_uni <- function(smpl.t, response_var, pred_vars, ntree = 1000, use_parsimony = FALSE){

  if (use_parsimony){
    parsimony <- 0.04
  } else {
    parsimony <- NULL
  }
  
  v<-strsplit(pred_vars,'_scl',2)#spliting name into two columns
  v1<-data.frame(do.call(rbind,v)) #combine into one dataframe and change format
  var_names<-unique(v1[,1])
  scales<-unique(v1[,2])

  # 2. Univariate rfs
  df.OOB<-NULL
  for (j in pred_vars){
    rf <- randomForest(as.factor(smpl.t[, response_var]) ~ ., data=smpl.t[, c(response_var, j)], importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
    rf.oob <- data.frame(cbind(j, rf$err.rate[ntree,1]))
    df.OOB<-rbind(df.OOB,rf.oob)
  }
  names(df.OOB)[1]<-"variable_sc"
  names(df.OOB)[2]<-"OOB"


  # 3. SO: Select best scales for each variable
  best.sc_fin<-as.character(NULL)

  for (k in var_names){
    OOB.var<-df.OOB %>% dplyr::filter(grepl(k,variable_sc))
    best.sc <- OOB.var[,1] [OOB.var[,2] %in% min(OOB.var[,2])]
    #best.sc<-best.sc[1]##in case there are 2 scales of same OOB it will take lower scale
    best.sc_fin<-append(best.sc_fin,best.sc)
  }

  var.SO_oob<-df.OOB[df.OOB$variable_sc %in% best.sc_fin,] #selecting

  rf.modelSO <- rf.modelSel(smpl.t[,best.sc_fin], as.factor(smpl.t[,response_var]), imp.scale="mir",
                            ntree=ntree, parsimony=parsimony, proximity=FALSE) ###!!! CHANGED - REDEFINED RESPONSE VAR

  sel.var <- rf.modelSO$selvars

  # now fit the final model (consisting of best combination of variables defined bby modelSel)
  rf.fit <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,sel.var], importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)

  #Add variables to the binary table with final variables
  df.var<-setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)
  v.bin<-(ifelse(colnames(df.var) %in% sel.var,"1","0")) #filling up table with model variables
  df.var[1,] <- v.bin

  return(list(rf.fit, df.var))
} # end fit_rf_uni



# fit MIR RF
fit_rf_mir <- function(smpl.t, response_var, pred_vars, ntree = 1000, use_parsimony = FALSE){

  if (use_parsimony){
    parsimony <- 0.04
  } else {
    parsimony <- NULL
  }

  # Model selection and final RF
  rf.model <- rf.modelSel(smpl.t[,pred_vars], as.factor(smpl.t[,response_var]),
                          imp.scale="mir", ntree=ntree,
                          parsimony=parsimony, proximity=FALSE)

  MIR<-rf.model$importance
  MIR$pred_vars_selc<-rownames(MIR) #add another column with variables names

  nc<-strsplit(as.character(MIR$pred_vars_selc),'_scl',2)#spliting name into two columns
  do.call(rbind,nc) #change format
  MIR_new<-data.frame(MIR,do.call(rbind,nc)) #combine into one dataframe

  MIR_best<-MIR_new %>% #MIR only for best scale !!!!
    group_by(X1) %>%
    slice(which.max(imp))

  pred_vars_selO<-MIR_best %>%pull(pred_vars_selc)
  df_SO <- smpl.t[,c(response_var,pred_vars_selO)]

  # Model selection for SO variables
  rf.modelSO <- rf.modelSel(df_SO[,4:ncol(df_SO)], as.factor(df_SO[,response_var]),
                            imp.scale="mir", ntree=ntree,
                            parsimony=parsimony, proximity=FALSE)

  sel.var <- rf.modelSO$selvars

  # Final RF
  rf.fit <- randomForest(y=as.factor(df_SO[,response_var]), x=df_SO[,sel.var], ntree=ntree, importance=TRUE, proximity=F, na.action=na.omit)

  #Add variables to the binary table with final variables
  df.var<-setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)
  v.bin<-(ifelse(colnames(df.var) %in% sel.var,"1","0")) #filling up table with model variables
  df.var[1,] <- v.bin

  return(list(rf.fit, df.var))
} # end fit_rif_mir



# fit BAF RF
fit_rf_baf <- function(smpl.t, response_var, pred_vars, ntree = 1000, use_parsimony = FALSE){

  if (use_parsimony){
    parsimony <- 0.04
  } else {
    parsimony <- NULL
  }
      
    list_i <- 1  
    RF_cv_k_fold <- function(input,response_var,pred_vars,seedLength=20){
      tmp_df <- as.data.frame(input[,c(response_var,pred_vars)])
      out_list<-list()
  
      for (baf_seed in sort(as.integer(runif(seedLength,1,10000)))){
        set.seed(baf_seed)
        seed_sample <- createDataPartition(tmp_df[, response_var], p = 0.8, list = FALSE)
        trainSet  <- tmp_df[seed_sample, ]
        testSet <- tmp_df[-seed_sample, ]
        if (length(unique(testSet$response))==1){next}
        RFmodel <- randomForest(y=as.factor(trainSet[,response_var]), x=trainSet[,pred_vars],
                   importance=TRUE, ntree=ntree, seed=baf_seed,proximity=FALSE, na.action=na.omit)
        RFpred <- predict(RFmodel, testSet,type='prob')
        
          test<-cbind(RFpred, testSet[,response_var])
          test[,1]<-1:nrow(test)
          acc.t<-presence.absence.accuracy(test,threshold=.500,st.dev=FALSE)
          
        out_list[[list_i]] <- data.frame(PCC=acc.t$PCC,
                                         R2 = R2(RFpred[,2], testSet[,response_var]),
                                         RMSE = RMSE(RFpred[,2], testSet[,response_var]),
                                         MAE = MAE(RFpred[,2], testSet[,response_var]))
        list_i <- list_i + 1
      }
      out_mat<-do.call(rbind,out_list)
      output<-apply(out_mat,2,mean, na.rm = TRUE)
      return(output)
    }

  # 2 CH_back_and_forth_alg ---------

    pred_scales <- pred_vars[!pred_vars %in% c("id")]
    preds<-unique(word(pred_scales,1,sep = "_scl"))

    # 2.1 preds best scale
      preds_best_scale<-vector();index<-1;output_list<-list()
      for (pred in preds){
        scales_tmp <- pred_scales[pred_scales %like% pred] 
        tmp <- smpl.t[,c(response_var, scales_tmp)]
        rRF_scale <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,scales_tmp],
                                  importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
        best_scale <- rRF_scale$importance[,4]
        preds_best_scale <-c(preds_best_scale,names(which.max(best_scale)))
      }

    # 2.2 preds forward selection
      forward1 <- rf.modelSel(smpl.t[,c(preds_best_scale)],as.factor(smpl.t[,response_var]),imp.scale="mir",
                            ntree=ntree, parsimony=parsimony, proximity=FALSE)
      forward1_preds <- forward1$selvars

    # 2.3 pred vars backward selection
      back_sel_results<-vector()
      back_tmp <- smpl.t[,c(response_var, forward1_preds)]
      rRF_preds <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,forward1_preds],
                   importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
      back_sel_vars <- names(rRF_preds$importance[,4])[!names(rRF_preds$importance[,4])%in%names(which.min(rRF_preds$importance[,4]))]

      for (steps in length(back_sel_vars):2){
        tmp <- smpl.t[,c(response_var, back_sel_vars)]
        rRF_back_sel <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,back_sel_vars],
                        importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
        tmp_output <- c(RF_cv_k_fold(tmp,response_var,back_sel_vars,seedLength=10)[1],back_sel_vars)  
        output_list[[index]]<- as.data.frame(matrix(tmp_output,1,length(tmp_output)));index<-index+1
        tmp_VIs <- rRF_back_sel$importance[,4]
        back_sel_vars <- names(tmp_VIs)[!names(tmp_VIs) %in% names(which.min(tmp_VIs))]
      }
      preds_output <- as.data.frame(rbindlist(output_list,fill=TRUE))
      names(preds_output)[1] <- "PCC"
      preds_tmp<-preds_output[which.max(preds_output$PCC),-1]
      preds_tmp<- as.character(preds_tmp)[complete.cases(as.character(preds_tmp))]

    # 2.4 scale backward selection
      for (i in 1:length(preds_tmp)){
        tmp_pred_scales<- pred_scales [pred_scales %ilike% word(preds_tmp[i],1,sep = "_scl")]
        for (j in 1:length(tmp_pred_scales)){
          tmp <- smpl.t[,c(response_var,preds_tmp[-i],tmp_pred_scales[j])]
          rRF_scale2 <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,c(preds_tmp[-i],tmp_pred_scales[j])],
                                       importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
          tmp_output <- c(RF_cv_k_fold(tmp,response_var,names(tmp)[-1],seedLength=10)[1],preds_tmp[-i],tmp_pred_scales[j]) 
          output_list[[index]]<- as.data.frame(matrix(tmp_output,1,length(tmp_output)));index<-index+1
        }
      }
      preds_sim_final <- as.data.frame(rbindlist(output_list,fill=TRUE))
      names(preds_sim_final)[1] <- "PCC"
      preds_tmp_scaled<-preds_sim_final[which.max(preds_sim_final$PCC),-1]
      preds_tmp_scaled<- as.character(preds_tmp_scaled)[complete.cases(as.character(preds_tmp_scaled))]

    # 2.5 final backward selection
      for (i in length(preds_tmp_scaled):2){
        tmp <- smpl.t[,c(response_var,preds_tmp_scaled)]
        rRF_back_fin <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,preds_tmp_scaled],
                        importance=TRUE, ntree=ntree, proximity=FALSE, na.action=na.omit)
        tmp_output <- c(RF_cv_k_fold(tmp,response_var,names(tmp)[-1],seedLength=10)[1],preds_tmp_scaled)
        output_list[[index]]<- as.data.frame(matrix(tmp_output,1,length(tmp_output)));index<-index+1
        tmp_VIs <- rRF_back_fin$importance[,4]
        preds_tmp_scaled <- names(tmp_VIs)[!names(tmp_VIs) %in% names(which.min(tmp_VIs))]
      }

      preds_sim_final_mat <- as.data.frame(rbindlist(output_list,fill=TRUE))
      names(preds_sim_final_mat)[1] <- "PCC"
      preds_final_mat<-preds_sim_final_mat[which.max(preds_sim_final_mat$PCC),]
      preds_final<- as.character(preds_final_mat[-1])[complete.cases(as.character(preds_final_mat[-1]))]
      sel_var <- as.character(as.data.frame(matrix(preds_final,1,length(preds_final)))) 
      rf_fin <- randomForest(y=as.factor(smpl.t[,response_var]), x=smpl.t[,sel_var],
                             importance=TRUE, ntree=ntree, seed=seed,proximity=FALSE, na.action=na.omit)
      df.var <- setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars) 
      v.bin <- (ifelse(pred_vars %in% sel_var,"1","0"))
      df.var[1,] <- v.bin

    return(list(rf_fin, df.var))
      
} # end fit_rf_BAF



# fit GLMnet
fit_glmnet_lasso <- function(smpl.t, response_var, pred_vars, pred_vars_sel = NULL){
  
  # Set model decision threshold (either lambda.min or lambda.1se)
  mdec <- 'lambda.min'
  
  # Initialize a table that will indicate which variables were selected
  df.var <- setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)

  if (!is.null(pred_vars_sel)){
    pred_vars <- pred_vars_sel
  }

  # Set up crossval
  nfold <- 10
  foldid <- sample(1:nfold, size = length(as.factor(smpl.t[,response_var])), replace = TRUE)

  # If prevalence is >= 10% use lasso
  # Otherwise, use regular glm
  if (sum(smpl.t[,response_var])/dim(smpl.t[,pred_vars])[1] >= 0.1) { 
    # Set alpha to one for lasso
    a <- 1
    cv <- cv.glmnet(as.matrix(smpl.t[,pred_vars]), as.factor(smpl.t[,response_var]), relax=T, family="binomial", type.measure="default", foldid=foldid, alpha = a, keep=T)
    # Get coefs
    cs <- as.matrix(coef(cv, s = mdec, exact=T))
  } else {
    # Set lambda penalty to zero
    cv <- glmnet(as.matrix(smpl.t[,pred_vars]), as.factor(smpl.t[,response_var]), relax=T, family="binomial", type.measure="default", foldid=foldid, alpha = 1, lambda = 0, keep=T)
    # Get coefs
    cs <- as.matrix(coef(cv, exact=T))
  }

  # Remove intercept and convert to data.frame
  cs <- data.frame(estimate=cs[-1,1])
  # Get variable importance (scaled coefficients) if needed at some point
  cs <- cs[cs$estimate != 0, , drop=F]
  # Get names of selected (nonzero) variables
  sel.var <- rownames(cs)

  #Add variables to the binary table with final variables
  v.bin<-(ifelse(colnames(df.var) %in% sel.var,"1","0")) #filling up table with model variables
  df.var[1,] <- v.bin

  return(list(cv, df.var, pred_vars))
} # end fit_glmnet_lasso



# fit GLM
fit_glm <- function(smpl.t, response_var, pred_vars, pred_vars_sel = NULL, use_dredge = FALSE){

  # Initialize a table that will indicate which variables were selected
  df.var <- setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)

  if (!is.null(pred_vars_sel)){
    pred_vars <- pred_vars_sel
  }
  var_fml <- paste(pred_vars, collapse="+")
  m_fml <- as.formula(paste("pa~",var_fml))

  # dredge
  if (use_dredge){
    options(na.action="na.fail") #important for dredge
    glm.fit <- glm(m_fml, data=smpl.t, family="binomial") #global model
    glm.d<-dredge(glm.fit)#all subsets from global model
    glm.best<-subset(glm.d, delta==0) #best model (top)

    m.sum<-summary(get.models(glm.d,1)[[1]])$coefficient #this is a full table for best model!
    coef.p<-as.data.frame(m.sum[,c(1,4)])
    names(coef.p)[1]<-"Coefficient"
    names(coef.p)[2]<-"p"
    var.d<-row.names(coef.p)#names of final variables with Intercept
    var.d.fin<-var.d[-1]
    
    #refit the model with var.d.fin
    var_fml.fin<-paste(var.d.fin,collapse="+")
    m_fml.fin<-as.formula(paste0(response_var, "~ ",var_fml.fin))
    glm.fin <- glm(m_fml.fin,data=smpl.t,family="binomial") #final model
    
    # Add variables to the binary table with final variables
    v.bin<-(ifelse(colnames(df.var) %in% var.d.fin,"1","0"))

  } else {
    glm.fit <- glm(m_fml, data=smpl.t, family="binomial") #global model
    glm.fin <- glm.fit
    summ <- summary(glm.fin)
    summ_coef <- summ$coefficients
    dt_coef <- as.data.table(summ_coef)
    dt_coef <- as.data.table(dt_coef[,vars:=rownames(summ_coef)])
    
    # Add variables to the binary table with final variables
    if (is.null(pred_vars_sel)){
      # Note: no variable selection happens in this instance
      v.bin <- rep(1, length(pred_vars))
    } else {
      v.bin <- (ifelse(colnames(df.var) %in% pred_vars,"1","0"))
    }
  }
  df.var[1,] <- v.bin

  return(list(glm.fin, df.var, dt_coef))
}



# fit RF
fit_rf <- function(smpl.t, response_var, pred_vars, pred_vars_sel = NULL, ntree = 1000){

  # Initialize a table that will indicate which variables were selected
  df.var <- setNames(data.frame(matrix(ncol=length(pred_vars),nrow=0)), pred_vars)

  if (!is.null(pred_vars_sel)){
    pred_vars <- pred_vars_sel
  }
  pred_df <- smpl.t[,c(response_var, pred_vars)]

  #rf.fit <- ranger(as.factor(response_var)~., data=pred_df, importance="impurity", num.trees = ntree)
  rf.fit <- randomForest(y=as.factor(pred_df[,response_var]), x=pred_df[,pred_vars], ntree=ntree, importance=TRUE, proximity=F, na.action=na.omit)

  # Add variables to the binary table with final variables
  # Note: no variable selection happens in this instance
  if (is.null(pred_vars_sel)){
    v.bin <- rep(1, length(pred_vars))
  } else {
    v.bin <- (ifelse(colnames(df.var) %in% pred_vars,"1","0"))
  }
  df.var[1,] <- v.bin

  return(list(rf.fit, df.var))
}
