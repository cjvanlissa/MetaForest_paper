library(parallel)
library(data.table)


#
#Define functions
#


#Metaforest stripped version for speed
metaforest <-
  function(formula,
           data,
           vi,
           tau2 = NULL,
           whichweights = "random",
           ...) {
    if (whichweights == "unif") {
      metaweights <- rep(1, nrow(data))
    }
    if (whichweights == "fixed") {
      metaweights <- (1 / vi)
    }
    if (whichweights == "random") {
      metaweights <- 1 / (vi + tau2)
    }
    metaweights <-
      (metaweights / sum(metaweights)) * length(vi)
    outputmodel <- do.call(ranger, list(
      formula = formula,
      data = data,
      importance = "permutation",
      case.weights = metaweights,
      ...))
    return(outputmodel)
  }

#metaCART
metaCART <- function(formula,
                     mods,
                     vi,
                     data,
                     c = 0,
                     tau2 = 0,
                     weights = "none",
                     method = "reg") {
  
  wi <- (1 / (vi + tau2)) / sum(1 / (vi + tau2)) * length(vi)
  tree <-
    do.call(rpart, list(
      formula = formula,
      data = data,
      weights = wi
    ))
  
  prunedtree <- treepruner(tree, .5)
  
  variable.importance <- prunedtree$variable.importance
  
  partytree <- as.party(prunedtree)
  
  if (length(partytree) == 0) {
    nodes <- factor(rep(1, nrow(data)))
  } else{
    nodes <- factor(predict(partytree, newdata = data, type = "node"))
  }
  
  if (length(unique(nodes)) > 1) {
    res <- try(rma(
      yi = data$yi,
      vi = vi,
      mods = ~ nodes - 1,
      intercept = FALSE), silent = TRUE)
    if(is.atomic(res))
    {
      res <- try(rma(
        yi = data$yi,
        vi = vi,
        mods = ~ nodes - 1,
        intercept = FALSE, 
        control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
      if(is.atomic(res))
      {
        res <- try(rma(
          yi = data$yi,
          vi = vi,
          mods = ~ nodes - 1,
          intercept = FALSE, 
          method="EB",
          control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
        if(is.atomic(res))
        {
          res <- try(rma(
            yi = data$yi,
            vi = vi,
            mods = ~ nodes - 1,
            intercept = FALSE, 
            method="ML",
            control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
          if(is.atomic(res))
          {
            largevalues <- which(abs(scale(data$yi)) > 3) 
            data[largevalues, ] <- data[largevalues-1, ]
            res <- try(rma(
              yi = data$yi,
              vi = vi,
              mods = ~ nodes - 1,
              intercept = FALSE), silent = TRUE) 
          }
        }
      }
    }
    
    if(is.atomic(res)){
      stop("rma did not converge. Please add another method for rma.")
    }
    
  } else {
    res <- try(rma(
      yi = data$yi,
      vi = vi), silent = TRUE)
    if(is.atomic(res))
    {
      res <- try(rma(
        yi = data$yi,
        vi = vi,
        control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
      if(is.atomic(res))
      {
        res <- try(rma(
          yi = data$yi,
          vi = vi,
          method="EB",
          control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
        if(is.atomic(res))
        {
          res <- try(rma(
            yi = data$yi,
            vi = vi,
            method="ML",
            control=list(stepadj=0.01, maxiter=1000)), silent = TRUE)  
          if(is.atomic(res))
          {
            largevalues <- which(abs(scale(data$yi)) > 3) 
            data[largevalues, ] <- data[largevalues-1, ]
            res <- try(rma(
              yi = data$yi,
              vi = vi), silent = TRUE) 
          }
        }
      }
    }
    
    if(is.atomic(res)){
      stop("rma did not converge. Please add another method for rma.")
    }
  }
  
  result <-
    list(
      tree = partytree,
      SubgroupTest = res,
      variable.importance = variable.importance
    )
  class(result) <- "metaCART"
  result
}


treepruner <- function(tree, c) {
  tree <- tree
  c <- c
  mindex <-
    which.min(tree$cptable[, 4]) #find the row of the minimum x-error
  cp.minse <-
    tree$cptable[mindex, 4] + c * tree$cptable[mindex, 5] #the minimum x-error + c*SE
  cp.row <-
    min(which(tree$cptable[, 4] <= cp.minse)) # find the smallest tree within the minimum x-error + c*SE
  cp.take <-
    tree$cptable[cp.row, 1] #get the cp value for the smallest tree
  prune(tree, cp = cp.take) #prune the tree
}




ModelAccuracy <-
  function(type,
           fit,
           newdata = NULL,
           observed,
           ymean = NULL) {
    if (type == "rma") {
      if (is.null(newdata)) {
        predicted	  <- predict(fit)$pred
        if (length(predicted) == 1) {
          predicted <- rep(predicted, length(observed))
        }
      } else {
        predicted <- predict.rma(fit, newmods = as.matrix(newdata))$pred
      }
    }
    
    if (type == "metaforest") {
      if (is.null(newdata)) {
        predicted <- fit$predictions
      } else {
        predicted <- predict(fit, data = newdata)$predictions
      }
      if (anyNA(predicted))
        predicted[is.na(predicted)] <- 0
    }
    
    if (type == "metaCART") {
      if (is.null(newdata) || length(fit$tree) == 0) {
        predicted	  <- predict(fit$SubgroupTest)$pred
        if (length(predicted) == 1) {
          predicted <- rep(predicted, length(observed))
        }
      } else {
        nodes <-
          factor(predict(fit$tree, newdata = newdata, type = "node")) #newdata=data,
        if (length(levels(nodes)) != nrow(fit$SubgroupTest$b)) {
          r_names <- row.names(fit$SubgroupTest$b)
          len <- nchar(r_names)
          nodes <- factor(nodes, levels = substr(r_names, 6, len))
        }
        
        predicted <-
          predict(fit$SubgroupTest,
                  newmods = model.matrix( ~ -1 + nodes, nodes, contrasts.arg = diag(nlevels(nodes))))$pred
      }
    }
    if (is.null(ymean))
      ymean <- mean(observed)
    SS.total    <- sum((observed - ymean) ^ 2)
    SS.residual <- sum((observed - predicted) ^ 2)
    SS.regression <- sum((predicted - ymean) ^ 2)
    
    r.2	<- 1 - SS.residual / SS.total
    mse	<- SS.residual / length(observed)
    sd(predicted)
    if (sd(predicted) == 0) {
      r.actual.pred <- 0
    } else{
      r.actual.pred <- cor(observed, predicted)
    }
    
    return(c(r.2, mse, r.actual.pred))
  }

CreateDataset <-
  function(ndataset,
           k_studies,
           mean_study_n,
           es,
           residual_heterogeneity,
           moderators,
           model) {

    #Set number of moderators
    #Use if construction to speed up
    moderators <- moderators + model
    if (model == 5)
      moderators <- moderators - 4
    
    #Make label for training and test datasets
    training <- c(rep(1, k_studies), rep(0, 100))
    
    #Randomly determine n from a normal distribution with mean mean_study_n and sd mean_study_n/3
    n <- rnorm(length(training), mean = mean_study_n, sd = mean_study_n /
                 3)
    #Round to even number
    n <- ceiling(n) - ceiling(n) %% 2
    #Truncate n so the lowest value can be 8; cases where n=2 will result in errors
    n[n < 8] <- 8
    
    #Generate moderator matrix x:
    x <- matrix(rnorm(length(n) * moderators), ncol = moderators)
    
    #Sample true effect sizes theta.i from a normal distribution with mean mu, and variance residual_heterogeneity,
    #where mu is the average population effect size. The value of mu depends on the values of the moderators and
    #the true model
    #mu <- eval(model)
    if (model == 1)
      mu <- es * x[, 1]
    if (model == 2)
      mu <- es * x[, 1] + es * x[, 2] + es * (x[, 1] * x[, 2])
    if (model == 3)
      mu <- es * x[, 1] + es * x[, 2] + es * x[, 3] + es * (x[, 1] * x[, 2]) +
      es * (x[, 1] * x[, 3]) + es * (x[, 2] * x[, 3]) + es * (x[, 1] * x[, 2] *
                                                                x[, 3])
    if (model == 4)
      mu <- es * x[, 1] + es * x[, 2] + es * (x[, 1] * x[, 2]) +
      es * x[, 3] + es * x[, 4] + es * (x[, 3] * x[, 4])
    if (model == 5)
      mu <- es * x[, 1] + es * (x[, 1] ^ 2) + es * (x[, 1] ^ 3)
    
    #theta.i: true effect size of study i
    theta.i <- mu + rnorm(length(n), 0, sd = sqrt(residual_heterogeneity))
    
    #Then the observed effect size yi is sampled from a non-central t-distribution under the assumption that
    #the treatment group and control group are both the same size
    p_ntk <- .5 #Percentage of cases in the treatment group
    ntk <- p_ntk * n #n in the treatment group for study i
    nck <- (1 - p_ntk) * n #n in the control group for study i
    df <- n - 2 #degrees of freedom
    j <- 1 - 3 / (4 * df - 1) #correction for bias
    nk = (ntk * nck) / (ntk + nck)
    ncp = theta.i * sqrt(nk) #Non-centrality parameter, c(mk) in Li et al.
    
    #Standardized mean difference drawn from a non-central t-distribution
    SMD <- mapply(FUN = rt,
                  n = 1,
                  df = df,
                  ncp = ncp)
    
    #yi is Hedges' g for study i
    yi <- SMD / ((j ^ -1) * (nk ^ .5))
    
    #Calculate the variance of the effect size
    vi <- j ^ 2 * (((ntk + nck) / (ntk * nck)) + ((yi / j) ^ 2 / (2 * (ntk +
                                                                         nck))))
    
    #Dersimonian and Laird estimate of tau2
    Wi <- 1 / vi[1:k_studies]
    tau2 <-
      max(0, (sum(Wi * (yi[1:k_studies] - (
        sum(Wi * yi[1:k_studies]) / sum(Wi)
      )) ^ 2) - (k_studies - 1)) / (sum(Wi) - (sum(Wi ^ 2) / sum(Wi))))
    
    data <- data.table(training, vi, yi, x)
    
    list(
      training = subset(data, training == 1,-1),
      testing = subset(data, training == 0,-c(1, 2)),
      housekeeping = cbind(n, mu, theta.i),
      tau2 = tau2
    )
  }


#
#Specify parameters of the simulation
#

hyper_parameters <- list(
  #Number of datasets per condition
  ndataset = 1:100,
  #Number of studies per dataset, normally distributed with mean n and sd n/3
  k_studies = c(20, 40, 80, 120),
  #Average n per study (k)
  mean_study_n = c(40, 80, 160),
  #Effect size
  es = c(.2, .5, .8),
  #Residual heterogeneity
  residual_heterogeneity = c(0, .04, .28),
  #Study-level moderators
  moderators = c(1, 2, 5),
  model = c(1, 2, 3, 4, 5)
)

summarydata <- expand.grid(hyper_parameters)

file_name <-
  paste0(paste(c(
    "summarydata",
    substr(Sys.time(), 1, 10),
    substr(Sys.time(), 12, 13),
    substr(Sys.time(), 15, 16)
  ), collapse = "_"), ".RData")
saveRDS(summarydata, file = file_name)
#summarydata<-readRDS("summarydata_2017-09-23_21_57.RData")


#
#Start clusters and load relevant packages
#
#stopCluster(cl=cl)
no_cores <- detectCores()
cl <- makeCluster(no_cores)

clusterEvalQ(cl, library(ranger))
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(metafor))
clusterEvalQ(cl, library(rpart))
clusterEvalQ(cl, library(party))
clusterEvalQ(cl, library(partykit))

clusterExport(cl = cl, "metaforest")
clusterExport(cl = cl, "metaCART")
clusterExport(cl = cl, "treepruner")
clusterExport(cl = cl, "CreateDataset")
clusterExport(cl = cl, "ModelAccuracy")



set.seed(42)

#Specify number of chunks
n_chunks <- 400
chunk_length <- (nrow(summarydata) / n_chunks)
if ((nrow(summarydata) %% n_chunks) > 0)
  stop("Make sure the number of chunks is an even divisor of the total job length")

#Save seeds
seeds <- round(runif(n_chunks, 0, 10000000))
while (length(unique(seeds)) < n_chunks) {
  addseeds <- n_chunks - length(unique(seeds))
  seeds <- c(unique(seeds), round(runif(addseeds, 0, 10000000)))
}
saveRDS(seeds, file = "seeds.RData", compress = FALSE)
#seeds<-readRDS("seeds.RData")

for (chunk in 1:n_chunks) {
  #chunk=1
  set.seed(seeds[chunk])
  chunk_conditions <-
    as.list(summarydata[(1 + ((chunk - 1) * chunk_length)):(chunk_length * chunk), ])
  
  simdata_test <-
    do.call(clusterMap, append(list(cl = cl, fun = CreateDataset), chunk_conditions))
  
  file_name <-
    paste0(paste(c("simdata", chunk), collapse = "_"), ".RData")
  saveRDS(simdata_test, file = file_name, compress = FALSE)
  
  #
  #Conduct metaForest models
  #
  
  metaforest_models <- parLapply(cl, simdata_test, function(data) {
    metaforest(
      yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      tau2 = data[[4]]
    )
  })
  
  metaforest.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        c(
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[1]], select = -c(1, 2)),
            observed = data[[1]]$yi
          ),
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[2]], select = -1),
            observed = data[[2]]$yi,
            ymean = mean(data[[1]]$yi)
          )
        )
      },
      models = metaforest_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  
  metaforest.importance <-
    parLapply(
      cl = cl,
      metaforest_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  file_name <-
    paste0(paste(c("metaforest_models", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest_models, file = file_name, compress = TRUE)
  
  file_name <-
    paste0(paste(c("metaforest_fits", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("metaforest_importance", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.importance,
          file = file_name,
          compress = FALSE)
  
  rm(metaforest_models, metaforest.fits, metaforest.importance)
  
  #
  #Conduct mf_fixed models
  #
  mf_fixed_models <- parLapply(cl, simdata_test, function(data) {
    metaforest(
      yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      whichweights = "fixed"
    )
  })
  
  
  mf_fixed.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        c(
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[1]], select = -c(1, 2)),
            observed = data[[1]]$yi
          ),
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[2]], select = -1),
            observed = data[[2]]$yi,
            ymean = mean(data[[1]]$yi)
          )
        )
      },
      models = mf_fixed_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  
  mf_fixed.importance <-
    parLapply(
      cl = cl,
      mf_fixed_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  file_name <-
    paste0(paste(c("mf_fixed_models", chunk), collapse = "_"), ".RData")
  saveRDS(mf_fixed_models, file = file_name, compress = TRUE)
  
  file_name <-
    paste0(paste(c("mf_fixed_fits", chunk), collapse = "_"), ".RData")
  saveRDS(mf_fixed.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("mf_fixed_importance", chunk), collapse = "_"), ".RData")
  saveRDS(mf_fixed.importance, file = file_name, compress = FALSE)
  
  rm(mf_fixed_models, mf_fixed.fits, mf_fixed.importance)
  
  
  #
  #Conduct mf_unif models
  #
  
  mf_unif_models <- parLapply(cl, simdata_test, function(data) {
    metaforest(
      yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      whichweights = "unif"
    )
  })
  
  
  mf_unif.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        c(
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[1]], select = -c(1, 2)),
            observed = data[[1]]$yi
          ),
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[2]], select = -1),
            observed = data[[2]]$yi,
            ymean = mean(data[[1]]$yi)
          )
        )
      },
      models = mf_unif_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  
  mf_unif.importance <-
    parLapply(
      cl = cl,
      mf_unif_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  res_tau2s <- clusterMap(cl, fun = function(fit, data){
    residuals <- data$training$yi - fit$predictions
    Wi <- 1 / data$training$vi
    max(0, (sum(Wi * (residuals - (
      sum(Wi * residuals) / sum(Wi)
    )) ^ 2) - (length(residuals) - 1)) / (sum(Wi) - (sum(Wi ^ 2) / sum(Wi))))
  }, fit = mf_unif_models, data = simdata_test, SIMPLIFY = TRUE)
  
  
  file_name <-
    paste0(paste(c("res_tau2s", chunk), collapse = "_"), ".RData")
  saveRDS(res_tau2s, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("mf_unif_models", chunk), collapse = "_"), ".RData")
  saveRDS(mf_unif_models, file = file_name, compress = TRUE)
  
  file_name <-
    paste0(paste(c("mf_unif_fits", chunk), collapse = "_"), ".RData")
  saveRDS(mf_unif.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("mf_unif_importance", chunk), collapse = "_"), ".RData")
  saveRDS(mf_unif.importance, file = file_name, compress = FALSE)
  
  rm(mf_unif_models, mf_unif.fits, mf_unif.importance)
  
  #Use unweighted MF tau
  metaforest_models <- clusterMap(cl = cl, fun = function(data, tau2) {
    metaforest(
      yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      tau2 = tau2
    )
  }, data = simdata_test, tau2 = res_tau2s)
  
  metaforest.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        c(
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[1]], select = -c(1, 2)),
            observed = data[[1]]$yi
          ),
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[2]], select = -1),
            observed = data[[2]]$yi,
            ymean = mean(data[[1]]$yi)
          )
        )
      },
      models = metaforest_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  
  metaforest.importance <-
    parLapply(
      cl = cl,
      metaforest_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  file_name <-
    paste0(paste(c("MF_restau_fits", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("MF_restau_importance", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.importance,
          file = file_name,
          compress = FALSE)
  
  rm(metaforest_models, metaforest.fits, metaforest.importance, res_tau2s)
  
  #Use known tau
  metaforest_models <- clusterMap(cl = cl, fun = function(data, tau2) {
    metaforest(
      yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      tau2 = tau2
    )
  }, data = simdata_test, tau2 = chunk_conditions$residual_heterogeneity)
  
  metaforest.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        c(
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[1]], select = -c(1, 2)),
            observed = data[[1]]$yi
          ),
          ModelAccuracy(
            type = "metaforest",
            fit = models,
            newdata = subset(data[[2]], select = -1),
            observed = data[[2]]$yi,
            ymean = mean(data[[1]]$yi)
          )
        )
      },
      models = metaforest_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  
  metaforest.importance <-
    parLapply(
      cl = cl,
      metaforest_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  file_name <-
    paste0(paste(c("MF_knowntau_fits", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("MF_knowntau_importance", chunk), collapse = "_"), ".RData")
  saveRDS(metaforest.importance,
          file = file_name,
          compress = FALSE)
  
  rm(metaforest_models, metaforest.fits, metaforest.importance)
  
  #Meta-cart
  cart_models <- clusterMap(cl = cl, fun = function(data, tau2) {
    metaCART(
      formula = yi ~ .,
      data = subset(data[[1]], select = -1),
      vi = data[[1]]$vi,
      c = .5,
      tau = data[[4]],
      weights = "RE",
      method = "reg"
    )
  }, data = simdata_test, tau2 = chunk_conditions$residual_heterogeneity)
  
  
  cart.importance <-
    parLapply(
      cl = cl,
      cart_models,
      fun = function(models) {
        models$variable.importance
      }
    )
  
  cart.fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        if (is.atomic(models$SubgroupTest)) {
          rep(NA, 6)
        } else {
          return(c(
            ModelAccuracy(
              type = "metaCART",
              fit = models,
              observed = data[[1]]$yi
            ),
            ModelAccuracy(
              type = "metaCART",
              fit = models,
              newdata = subset(data[[2]], select = -1),
              observed = data[[2]]$yi,
              ymean = mean(data[[1]]$yi)
            )
          ))
        }
      },
      models = cart_models,
      data = simdata_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
  file_name <-
    paste0(paste(c("cart_models", chunk), collapse = "_"), ".RData")
  saveRDS(cart_models, file = file_name, compress = TRUE)
  
  file_name <-
    paste0(paste(c("cart_fits", chunk), collapse = "_"), ".RData")
  saveRDS(cart.fits, file = file_name, compress = FALSE)
  
  file_name <-
    paste0(paste(c("cart_importance", chunk), collapse = "_"), ".RData")
  saveRDS(cart.importance, file = file_name, compress = FALSE)
  
  rm(cart_models, cart.fits, cart.importance)
  
}

#Close cluster
stopCluster(cl)


###############
# Merge files #
###############

#summarydata<-readRDS("summarydata_2017-08-30_16_57.RData")
data<-as.data.table(summarydata)
head(data)

nfiles <- 400
chunk_length <- 405

model_names <- unique(gsub("_fits_\\d+.RData", "", list.files(pattern="fits", all.files=FALSE, full.names=FALSE)))

model_names <- model_names[c(2,4,5,3,6,1)]

data[, paste0(unlist(lapply(model_names, rep, 6)), c("_train_r2", "_train_mse", "_train_r", "_test_r2", "_test_mse", "_test_r")) := double() ]

for(modelname in model_names){
  for(i in 1:nfiles){
    fileName<-paste0(c(modelname, "_fits_", i, ".RData"), collapse="")
    if (file.exists(fileName)){
      set(x = data, 
          i = (1+((i-1)*chunk_length)):(i*chunk_length), 
          j = paste0(rep(modelname, 6), c("_train_r2", "_train_mse", "_train_r", "_test_r2", "_test_mse", "_test_r")), 
          value = as.data.table(readRDS(fileName)))
    }
  }
}
data
#Replace names with short two-letter codes
model_names
newnames <- c("mf", "mk", "mr", "fx", "un", "ca")
mapply(FUN = function(old, new){
  names(data) <<- gsub(old, new, names(data))
}, old = model_names, new = newnames)

tail(data)

saveRDS(data, "data.RData")
#data<-readRDS("data.RData")

library(ggplot2)
library(data.table)


vars <- c("_train_mse", "_train_mse", "_train_r", "_test_mse", "_test_mse", "_test_r")

###################################
# Predictive performance analyses #
###################################

analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(newnames, "_train_r2"))]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

sink("output.txt")
tmp<-rbind(round(analyzedat[, lapply(.SD, mean), .SDcols=c(7:ncol(analyzedat))], 2),
           round(analyzedat[, lapply(.SD, sd), .SDcols=c(7:ncol(analyzedat))], 2))
tmp<-data.frame(names(tmp), t(tmp))
apply(tmp, 1, function(x){paste(x, collapse=", ")})
sink()

analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(newnames, "_test_r2"))]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]


analyzedat[ , paste(newnames[1], newnames[-1], sep = "_") := 
              lapply(paste0(newnames[-1], "_test_r2"), function(x){
                analyzedat[[paste0(newnames[1], "_test_r2")]] - analyzedat[[x]]
              })]

tmp<-rbind(round(analyzedat[, lapply(.SD, mean), .SDcols=c(7:ncol(analyzedat))], 3),
           round(analyzedat[, lapply(.SD, sd), .SDcols=c(7:ncol(analyzedat))], 3))
tmp<-data.frame(names(tmp), t(tmp))
sink("output.txt", append = TRUE)
apply(tmp, 1, function(x){paste(x, collapse=", ")})
sink()


#Check for significant moderators
EtaSq<-function (x) 
{
  anovaResults <- summary.aov(x)[[1]]
  anovaResultsNames <- rownames(anovaResults)
  SS <- anovaResults[,2] #SS effects and residuals
  k <- length(SS) - 1  # Number of factors 
  ssResid <- SS[k + 1]  # Sum of Squares Residual
  ssTot <- sum(SS)  # Sum of Squares Total
  SS <- SS[1:k] # takes only the effect SS
  anovaResultsNames <- anovaResultsNames[1:k]
  etaSquared <- SS/ssTot # Should be the same as R^2 values
  partialEtaSquared <- SS/(SS + ssResid)
  res <- cbind(etaSquared, partialEtaSquared)
  colnames(res) <- c("Eta^2", "Partial Eta^2")
  rownames(res) <- anovaResultsNames
  return(res)
}

yvars<-c(paste0(newnames, "_test_r2"), paste(newnames[1], newnames[-1], sep = "_"))

anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, "(es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2", sep="~")
  thisaov<-aov(as.formula(form), data=analyzedat)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs
names(coef.table)<-yvars

table.names<-row.names(coef.table)
table.names<-gsub("_studies", "", table.names)
table.names<-gsub("moderators", "M", table.names)
table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
table.names<-gsub("es", "$\\beta$", table.names)
table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names
names(coef.table)<-c("MetaForest", "MF_Knowntau", "MF_Residualtau" ,"Fixed-effects",  "Unweighted", "metaCART", "$\\Delta_{MF-MK}$", "$\\Delta_{MF-MR}$", "$\\Delta_{MF-FX}$", "$\\Delta_{MF-UN}$",   "$\\Delta_{MF-CA}$")

study2.coef.table<-coef.table[, -c(2,3, 7, 8)]

library(xtable)

print(xtable(study2.coef.table), file="table_study2_coefs.tex",sanitize.text.function=function(x){x})

sink("output.txt", append = TRUE)
study2.coef.table
sink()

##Figure out which predictors are most important per variable
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
sink("output.txt", append = TRUE)
lapply(study2.coef.table, function(x){
  tmp <- as.numeric.factor(x)
  names(tmp)<- c("k", "n", "beta", "tau2", "m", "model", "k:n", "k:beta", "k:tau2", "k:m", "k:model", "n:beta", "n:tau2", "n:m", "n:model", "beta:tau2", "beta:m", "beta:model", "tau2:m", "tau2:model", "m:model")
  sort(tmp, decreasing = TRUE)
})
sink()

##################
# Power analysis #
##################
analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", "mf_test_r2", "ca_test_r2")]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

measure.vars <- names(analyzedat)[-c(1:6)]
grouping.vars <- quote(list(k_studies,
                            mean_study_n,
                            es,
                            residual_heterogeneity,
                            moderators, model))

power2<-analyzedat[,lapply(.SD, function(x){ifelse(quantile(x, probs = .2, na.rm = TRUE)>0, 1, 0)}),by=eval(grouping.vars), .SDcols=measure.vars]

plotdat<-power2

plotdat[, power := "Neither"]
plotdat[mf_test_r2 == 1 & ca_test_r2 == 0, power := "Only MetaForest"]
plotdat[mf_test_r2 == 0 & ca_test_r2 == 1, power := "Only metaCART"]
plotdat[mf_test_r2 == 1 & ca_test_r2 == 1, power := "Both"]
plotdat[, power := factor(power)]
plotdat$power <- ordered(plotdat$power, levels = c("Neither", "Only MetaForest", "Only metaCART", "Both"))
table(plotdat$power)
names(plotdat)[c(1:6)]<-c("k", "n", "es", "tau2", "M", "Model")

categorical<-c("k", "n", "es", "tau2", "M", "Model")

plotdat[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
plotdat$tau2 <- factor(plotdat$tau2, labels = c("tau^2:~.00", "tau^2:~.04", "tau^2:~.28"))
plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
levels(plotdat$Model)<-c("a", "b", "c", "d", "e")
table(plotdat$power)
thisplot <- ggplot(plotdat, aes(x = n, y = k, fill = power)) +
  geom_raster(hjust = 0, vjust = 0) +
  facet_grid(
    Model + M ~ es + tau2,
    labeller = labeller(
      es = label_parsed,
      tau2 = label_parsed,
      Model = label_both,
      M = label_both)) +
  scale_fill_manual(values = c("white", "grey50", "black")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = expression(bar(n)), fill = "Power > .80")
thisplot
ggsave("power_plot.pdf", plot=thisplot, device="pdf", width=210, height=297, units="mm")  

#Percentage table
power <- analyzedat[!(model==1),lapply(.SD, function(x){sum(x > 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
powertmp <- analyzedat[model==1,lapply(.SD, function(x){sum(x < 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
power <- rbindlist(list(powertmp, power))

power <- power
power <- melt(power, measure.vars = c("mf_test_r2", "ca_test_r2"))

tmp <- dcast(power, model+moderators+residual_heterogeneity~variable+ es + k_studies + mean_study_n)
write.csv(tmp, "study2 power.csv")

###################################################
#               R-squared plots                   #
###################################################
includevars <- c("mf", "fx", "un", "ca")
plotdat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(includevars, "_test_r2"))]
plotdat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

id.vars=names(plotdat)[1:6]
plotdat2<-melt.data.table(plotdat, id.vars=id.vars)
names(plotdat2)<-c("k", "n", "es", "tau2", "Moderators", "Model", "Algorithm", "R2")
categorical<-c("es", "tau2", "Moderators", "Model", "Algorithm")
plotdat2[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
levels(plotdat2$Algorithm)
plotdat2$Algorithm <- ordered(plotdat2$Algorithm, levels = c("mf_test_r2", "fx_test_r2", "un_test_r2", "ca_test_r2"))#, "mk_test_r2", "mr_test_r2"
levels(plotdat2$Algorithm)<-c("MetaForest", "MetaForest FX", "MetaForest UN", "metaCART")#, "MetaForest KT", "MetaForest RT"
plotdat2$tau2 <- factor(plotdat2$tau2, labels = c("tau^2:~.00", "tau^2:~.04", "tau^2:~.28"))
plotdat2$es <- factor(plotdat2$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat2$n <- factor(plotdat2$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))
levels(plotdat2$Model)<-c("a", "b", "c", "d", "e")
#plotdat2 <- plotdat2[!(Algorithm %in% c("MetaForest FX", "MetaForest UN")),]

#Explore interaction between beta:Model
plotdat3<-plotdat2
plotdat3$es <- factor(plotdat3$es, labels = c("0.2", "0.5", "0.8"))
plots <- ggplot(plotdat3, aes(x = es, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(0.2, .5, 0.8), limits=c(0.2,.8))+
  facet_grid(. ~ Model, labeller = labeller(Model=label_both))+ 
  theme_bw()+ labs(y=expression(R^2), x=expression(beta))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("int_beta_model.pdf", plot=plots, device="pdf")

#Explore interaction between tau2:Model
plotdat3<-plotdat2
plotdat3$tau2 <- factor(plotdat3$tau2, labels = c("0", "0.04", "0.28"))
plots <- ggplot(plotdat3, aes(x = tau2, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(0, .5, 1), limits=c(0,1))+
  facet_grid(. ~ Model, labeller = labeller(Model=label_both))+ 
  theme_bw()+ labs(y=expression(R^2), x=expression(tau^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("int_tau2_model.pdf", plot=plots, device="pdf")


#Explore interaction between k:Model
plotdat3<-plotdat2
#plotdat3$k <- as.numeric.factor(plotdat3$k)
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20,120))+
  facet_grid(. ~ Model, labeller = labeller(Model=label_both))+ 
  theme_bw()+ labs(y=expression(R^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("int_k_model.pdf", plot=plots, device="pdf")

#Explore interaction between k:tau2
plotdat3<-plotdat2
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20,120))+
  facet_grid(. ~ tau2, labeller = labeller(tau2 = label_parsed))+ 
  theme_bw()+ labs(y=expression(R^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("int_k_tau2.pdf", plot=plots, device="pdf")

#Marginal effect of k
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20,120))+
  theme_bw()+ labs(y=expression(R^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_k.pdf", plot=plots, device="pdf")

#Marginal effect of n
plotdat3<-plotdat2
plotdat3$n <- factor(plotdat3$n, labels = c("40", "80", "160"))
plots <- ggplot(plotdat3, aes(x = n, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(40, 80, 160), limits=c(40,160))+
  theme_bw()+ labs(y=expression(R^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_n.pdf", plot=plots, device="pdf")


#Marginal effect of M
plotdat3<-plotdat2
#plotdat3$Moderators <- as.numeric.factor(factor(plotdat3$Moderators))
plots <- ggplot(plotdat3, aes(x = Moderators, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(1, 2, 5), limits=c(1, 5))+
  theme_bw()+ labs(y=expression(R^2))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_M.pdf", plot=plots, device="pdf")

#Marginal effect of beta
plotdat3<-plotdat2
plotdat3$es <- factor(plotdat3$es, labels=c(".2", ".5", ".8"))
plots <- ggplot(plotdat3, aes(x = es, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(.2, .5, .8), limits=c(.2,.8))+
  theme_bw()+ labs(y=expression(beta))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_beta.pdf", plot=plots, device="pdf")

#Marginal effect of tau2
plotdat3<-plotdat2
plotdat3$tau2 <- factor(plotdat3$tau2, labels=c("0", "0.04", "0.28"))
plots <- ggplot(plotdat3, aes(x = tau2, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(0, .5, 1), limits=c(0,1))+
  theme_bw()+ labs(y=expression(beta))+ 
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_beta.pdf", plot=plots, device="pdf")



###################################################
#               Big R-squared plots               #
###################################################
includevars <- c("mf", "fx", "un", "ca")
plotdat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(includevars, "_test_r2"))]
plotdat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

#Summarize several columns:
measure.vars <- names(plotdat)[-c(1:6)]
grouping.vars <- quote(list(k_studies,
                            mean_study_n,
                            es,
                            residual_heterogeneity,
                            moderators, model))

plotdat2<-plotdat[,lapply(.SD,mean, na.rm=TRUE),by=eval(grouping.vars), .SDcols=measure.vars]

id.vars=names(plotdat2)[1:6]
plotdat3<-melt.data.table(plotdat2, id.vars=id.vars)#measure.vars = measure.vars)
plotdat3[, "Algorithm" := substr(variable, 1,2)]
plotdat3[, variable := gsub("\\w+r2", "r2", variable)]

plotdat4<-dcast.data.table(plotdat3, k_studies + mean_study_n + es + residual_heterogeneity + moderators + model+Algorithm~variable, value.var="value")

names(plotdat4)[1:8]<-c("k", "n", "es", "tau2", "Moderators", "Model", "Algorithm", "R2")

categorical<-c("es", "tau2", "Moderators", "Model", "Algorithm")
plotdat4[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
levels(plotdat4$Algorithm)
plotdat4$Algorithm <- ordered(plotdat4$Algorithm, levels = c("mf", "fx", "un", "ca"))
levels(plotdat4$Algorithm)<-c("MetaForest RE", "MetaForest FX", "MetaForest UN", "metaCART")

plotdat4$tau2 <- factor(plotdat4$tau2, labels = c("tau^2:~0", "tau^2:~0.04", "tau^2:~0.28"))
plotdat4$es <- factor(plotdat4$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat4$n <- factor(plotdat4$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))
plotdat4$Model <- factor(plotdat4$Model, labels = c("a", "b", "c", "d", "e"))

plots<-lapply(levels(plotdat4$Model), function(x){
  tmp<-plotdat4[Model==x, -6]
  #tmp<-plotdat4[Model==x&!(Algorithm %in% c("MetaForest FX", "MetaForest UN")), -6]
  ggplot(tmp, aes(x = k, y = R2, linetype = Algorithm, shape = Algorithm, group=Algorithm)) +
    geom_point() +
    geom_path() + 
    facet_grid(es + tau2 ~ Moderators+n, labeller = labeller(es = label_parsed, tau2 = label_parsed, n=label_parsed, Moderators=label_both))+
    theme_bw()+ labs(y=expression(R^2))+ 
    geom_hline(yintercept = 0)
})

for(i in 1:length(plots)){
  ggsave(paste0("Rsquare_plot_", i,".pdf"), plot=plots[[i]], device="pdf", width=297, height=210, units="mm")  
}


##########################################
# Variable importance analyses and plots #
##########################################

#
#Merge files
#
data<-data.table(readRDS(list.files(pattern = "summarydata")))

nfiles<-400
chunk_length<-405

###Read files 
#Metaforest
data[,c("mf_V1", "mf_V2", "mf_V3", "mf_V4", "mf_V5", "mf_V6", "mf_V7", "mf_V8", "mf_V9"):=double() ]

for(i in 1:nfiles){
  #i=279
  fileName<-paste0(c("metaforest_importance_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-readRDS(fileName)
    imp.values<-lapply(imp.values, function(x){c(x, rep(NA, 9-length(x)))})
    imp.values<-as.data.table(t(do.call(cbind, imp.values)))
    set(x=data, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c("mf_V1", "mf_V2", "mf_V3", "mf_V4", "mf_V5", "mf_V6", "mf_V7", "mf_V8", "mf_V9"), value=imp.values)
  }
}

#Metacart
data[,c("ca_V1", "ca_V2", "ca_V3", "ca_V4", "ca_V5", "ca_V6", "ca_V7", "ca_V8", "ca_V9"):=double() ]
for(i in 1:nfiles){
  fileName<-paste0(c("cart_importance_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-readRDS(fileName)
    imp.values<-lapply(imp.values, function(x){
      full<-c(V1 = 0, V2 = 0, V3 = 0, V4 = 0, V5 = 0, V6 = 0, V7 = 0, V8 = 0, V9 = 0)
      full[names(x)]<-x
      full
    })
    imp.values<-as.data.table(t(do.call(cbind, imp.values)))
    set(x=data, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c("ca_V1", "ca_V2", "ca_V3", "ca_V4", "ca_V5", "ca_V6", "ca_V7", "ca_V8", "ca_V9"), value=imp.values)
  }
}

data[,nvars:=model+moderators]
data[model==5, nvars:=nvars-4]

for(nv in as.numeric(names(table(data$nvars)))[-8]){
  set(data, i=which(data$nvars==nv), j=c((min(grep("^ca_V", names(data)))+nv):max(grep("^ca_V", names(data)))), value=NA)
}

#Metacart names
data[,c("cn_V1", "cn_V2", "cn_V3", "cn_V4", "cn_V5", "cn_V6", "cn_V7", "cn_V8", "cn_V9"):=double() ]
for(i in 1:nfiles){
  fileName<-paste0(c("cart_names_importance_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-readRDS(fileName)
    imp.values<-lapply(imp.values, function(x){
      full<-c(V1 = 0, V2 = 0, V3 = 0, V4 = 0, V5 = 0, V6 = 0, V7 = 0, V8 = 0, V9 = 0)
      full[names(x)]<-x
      full
    })
    imp.values<-as.data.table(t(do.call(cbind, imp.values)))
    set(x=data, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c("cn_V1", "cn_V2", "cn_V3", "cn_V4", "cn_V5", "cn_V6", "cn_V7", "cn_V8", "cn_V9"), value=imp.values)
  }
}
for(nv in as.numeric(names(table(data$nvars)))[-8]){
  set(data, i=which(data$nvars==nv), j=c((min(grep("^cn_V", names(data)))+nv):max(grep("^cn_V", names(data)))), value=NA)
}

###
#Check which design factors are relevant for a correlated moderator
#
yvars<-list("mf_V1", "ca_V1")

anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, "(es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2", sep="~")
  thisaov<-aov(as.formula(form), data=data)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})
etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs
names(coef.table)<-yvars

table.names<-row.names(coef.table)
table.names<-gsub("_studies", "", table.names)
table.names<-gsub("moderators", "M", table.names)
table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
table.names<-gsub("es", "$\\beta$", table.names)
table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names
names(coef.table)<-c("MetaForest", "CA importance")

study2.coef.table<-coef.table
sink("output.txt", append = TRUE)
study2.coef.table
sink()
library(xtable)

print(xtable(study2.coef.table), file="table_study2_importance.tex",sanitize.text.function=function(x){x})


#Standardize importances
plotdat<-data

for(var in c("mf", "ca")){
  plotdat[, "imp.sum" := rowSums(abs(.SD), na.rm=TRUE), .SDcols=grep(paste0("^", var, "_V"), names(plotdat), value = TRUE)]
  
  for(thisvar in grep(paste0("^", var, "_V"), names(plotdat), value = TRUE)){
    plotdat[!(imp.sum == 0), (thisvar) := (.SD/imp.sum)*100, .SDcols = thisvar]
    #set(plotdat, i = !(imp.sum == 0), j=var, value=(plotdat[[var]]/plotdat[["imp.sum"]])*100)
  }
}

plotdat[, c("nvars", "imp.sum") := NULL]

plotdat<-melt(plotdat, id.vars = names(plotdat)[c(1:7)])

plotdat[, Algorithm:=substr(variable, 1,2)]

plotdat[, variable:=substr(variable, 5,5)]

library(ggplot2)

names(plotdat)[2:9]<-c("k", "n", "es", "tau2", "Moderators", "Model", "Variable", "Importance")

plotdat[, c(2:8, 10) := lapply(.SD, factor), .SDcols=c(2:8, 10)]

plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat$Model <- factor(plotdat$Model, labels = c("a", "b", "c", "d", "e"))

plotdat <- plotdat[Variable %in% c( "1", "2", "3", "4", "5"), ]

plots<-lapply(levels(plotdat$Algorithm), function(x){
  tmp<-plotdat[Algorithm == x, -10]
  ggplot(tmp, aes(x = Variable, y = Importance)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Model + es ~ k, labeller = 
                 labeller(Model = label_both, es = label_parsed, k=label_both)) + 
    theme_bw()+geom_hline(yintercept = 0)
})

for(i in 1:length(plots)){
  ggsave(paste0("Importance_Plot_", levels(plotdat$Algorithm)[i],".pdf"), plot=plots[[i]], device="pdf", width=210, height=297, units="mm")  
}




##########################################
# Tau2s                                  #
##########################################

#
#Merge files
#
data<-data.table(readRDS(list.files(pattern = "summarydata")))

nfiles<-400
chunk_length<-405

###Read files 
#Metaforest
data[,c("tau2"):=double() ]

for(i in 1:nfiles){
  #i=279
  fileName<-paste0(c("res_tau2s_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-unlist(readRDS(fileName))
    set(x=data, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c("tau2"), value=imp.values)
  }
}


anovas <- aov(tau2 ~ (es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2, data=data)
etasqs <- EtaSq(anovas)[ , 2]
etasqs <- formatC(etasqs, 2, format="f")

coef.table<-etasqs
names(coef.table)<-yvars

plotdat <- data
names(plotdat)[2:8]<-c("k", "n", "es", "tau2", "Moderators", "Model", "tau2_est")

plotdat[, c(2:7) := lapply(.SD, factor), .SDcols=c(2:7)]

plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat$Model <- factor(plotdat$Model, labels = c("a", "b", "c", "d", "e"))
plotdat$tau2 <- factor(plotdat$tau2, labels = c("tau^2:~0", "tau^2:~0.04", "tau^2:~0.28"))
plotdat$n <- factor(plotdat$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))

vlinedata <- data.frame(es=unlist(lapply(levels(plotdat$es), rep, 4)), k = rep(levels(plotdat$k), 3))
vlinedata <- data.frame(rbind(vlinedata, vlinedata, vlinedata))
vlinedata$tau2 <- rep(levels(plotdat$tau2), 12)
vlinedata

plots<-lapply(levels(plotdat$Model), function(x){
  x <-levels(plotdat$Model)[1]
  tmp<-plotdat[Model == x, -7]
  ggplot(tmp, aes(x = tau2_est)) +
    geom_density()+
    geom_vline(aes(xintercept = as.numeric(gsub("tau\\^2\\:\\~", "", tau2))))+
    facet_grid(es + k ~ tau2, labeller = 
                 labeller(es = label_parsed, tau2 = label_parsed, k=label_both)) + 
    theme_bw()
})
for(i in 1:length(plots)){
  ggsave(paste0("Tau2_residual_Plot_", levels(plotdat$Model)[i],".pdf"), plot=plots[[i]], device="pdf", width=210, height=297, units="mm")  
}
