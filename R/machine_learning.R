machine_learning <- function(x, y,
                             models = c("enet"
                                        , "lasso"
                                        , "ridge"
                                        # , "rf"
                                        # , "lightgbm"
                             ),
                             cores = 2){
  y = y[!is.na(y)]
  cells = intersect(rownames(x), names(y))
  x = as.matrix(x[cells, ]); y = y[cells]
  x = x[, colSums(x!=0)>5]

  require(caret)
  require(doMC)
  registerDoMC(cores=cores)
  fits = list()
  lambda <- 10^seq(-3, 3, length = 1000)
  # Elastic net regression
  if("enet" %in% models){
    require(glmnet)
    message(format(Sys.time(), "%X  "), "Train the elastic net regression model ...")
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
    fit_en <- train(x, y, method="glmnet", metric="RMSE", tuneLength=20, trControl=control,
                    tuneGrid = expand.grid(alpha = seq(0.1,0.9,0.1), lambda = lambda))
    fits$enet = fit_en
  }

  # Lasso regression
  if("lasso" %in% models){
    require(glmnet)
    message(format(Sys.time(), "%X  "), "Train the lasso regression model ...")
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
    lasso <- train(x, y, method="glmnet", metric="RMSE", tuneLength=20, trControl=control,
                   tuneGrid = expand.grid(alpha = 1, lambda = lambda))
    fits$lasso = lasso
  }

  # Ridge regression
  if("ridge" %in% models){
    require(glmnet)
    message(format(Sys.time(), "%X  "), "Train the ridge regression model ...")
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
    ridge <- train(x, y, method="glmnet", metric="RMSE", tuneLength=20, trControl=control,
                   tuneGrid = expand.grid(alpha = 0, lambda = lambda))
    fits$ridge = ridge
  }

  # random forest
  if("rf" %in% models){
    require(randomForest)
    message(format(Sys.time(), "%X  "), "Train the random forest model ...")
    control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
    fit_rf <- train(x, y, method="rf", metric="RMSE", tuneLength=20, trControl=control)
    fits$rf = fit_rf
  }

  # lightgbm
  if("lightgbm" %in% models){
    require(lightgbm)
    message(format(Sys.time(), "%X  "), "Train the lightgbm model ...")
    dtrain <- lgb.Dataset(x, label = y, categorical_feature = categorical_feature)
    params <- list(objective = "regression", metric = "l2")
    lightgbm_model <- lgb.train(params, dtrain, 30)
    fits$lightgbm = lightgbm_model
  }
  return(fits)
}
