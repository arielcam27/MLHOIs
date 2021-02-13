
times = 2
times = 5
#---------------------------------

machine_knn_HOI <- function(dataset) {
 
  accuracy    <- rep(NA, times)
  kappa       <- rep(NA, times)
  specificity <- rep(NA, times)
  sensitivity <- rep(NA, times)
  
  for(t in 1:times){
    cat("Progress:", t, "\n")
    inTrain <- createDataPartition(y = dataset$isHOI, p = .80, list = FALSE)
    training <- dataset[ inTrain,]
    testing  <- dataset[-inTrain,]
    
    model       <- train(isHOI~., data=training, method="knn", trControl = trainControl(classProbs =  TRUE), preProcess = c("center","scale"), tuneLength = 20)
    predictions <- predict(model, newdata=testing)
    result      <- confusionMatrix(table(testing$isHOI, predictions))
    
    accuracy[t]    <- result$overall[1]
    kappa[t]       <- result$overall[2]
    sensitivity[t] <- result$byClass[[1]]
    specificity[t] <- result$byClass[[2]]
    
  }
  
  result.predicted.prob <- predict(model, testing, type="prob")
  result.roc            <- roc(testing$isHOI, result.predicted.prob$hoi) # Draw ROC curve.
  #plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
  
  return(list("result" = result,
              "fit"    = model,
              "acc"    = accuracy,
              "kap"    = kappa,
              "spec"   = specificity,
              "sens"   = sensitivity,
              "ROC"    = result.roc))
}

#---------------------------------


#---------------------------------

machine_lin_HOI <- function(dataset) {
  
  accuracy    <- rep(NA, times)
  kappa       <- rep(NA, times)
  specificity <- rep(NA, times)
  sensitivity <- rep(NA, times)
  
  for(t in 1:times){
    cat("Progress:", t, "\n")
    inTrain <- createDataPartition(y = dataset$isHOI, p = .80, list = FALSE)
    training <- dataset[ inTrain,]
    testing  <- dataset[-inTrain,]
    
    model       <- train(isHOI~., data=training, method="glm", trControl = trainControl(classProbs =  TRUE), family = "binomial")
    predictions <- predict(model, newdata=testing)
    result      <- confusionMatrix(table(testing$isHOI, predictions))
    
    accuracy[t]    <- result$overall[1]
    kappa[t]       <- result$overall[2]
    sensitivity[t] <- result$byClass[[1]]
    specificity[t] <- result$byClass[[2]]
  }
  
  result.predicted.prob <- predict(model, testing, type="prob")
  result.roc            <- roc(testing$isHOI, result.predicted.prob$hoi) # Draw ROC curve.
  #plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
  
  return(list("result" = result,
              "fit"    = model,
              "acc"    = accuracy,
              "kap"    = kappa,
              "spec"   = specificity,
              "sens"   = sensitivity,
              "ROC"    = result.roc))
}

machine_svm_HOI <- function(dataset) {
  
  accuracy    <- rep(NA, times)
  kappa       <- rep(NA, times)
  specificity <- rep(NA, times)
  sensitivity <- rep(NA, times)
  
  for(t in 1:times){
    cat("Progress:", t, "\n")
    inTrain <- createDataPartition(y = dataset$isHOI, p = .80, list = FALSE)
    training <- dataset[ inTrain,]
    testing  <- dataset[-inTrain,]
    
    model       <- train(isHOI~., data=training, method="svmRadial", trControl = trainControl(classProbs =  TRUE), preProcess = c("center","scale"), tuneLength = 20)
    predictions <- predict(model, newdata=testing)
    result      <- confusionMatrix(table(testing$isHOI, predictions))
    
    accuracy[t]    <- result$overall[1]
    kappa[t]       <- result$overall[2]
    sensitivity[t] <- result$byClass[[1]]
    specificity[t] <- result$byClass[[2]]
  }
  
  result.predicted.prob <- predict(model, newdata=testing, type="prob")
  result.roc            <- roc(testing$isHOI, result.predicted.prob$hoi) # Draw ROC curve.
  #plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
  
  return(list("result" = result,
              "fit"    = model,
              "acc"    = accuracy,
              "kap"    = kappa,
              "spec"   = specificity,
              "sens"   = sensitivity,
              "ROC"    = result.roc))
}
