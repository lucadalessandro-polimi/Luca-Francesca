.KfoldsCrossValidationCtr <- setRefClass("KfoldsCrossValidationObject",
                                         fields = c(K="integer",       # number of folds
                                                      data="data.frame", # raw data
                                                      folds="list",      # k Folds
                                                      num_obs_kFold="integer" # num_obs_kFold[i] = nrow(folds[[i]])
                                                      ),      
                                         methods = list(
                                           get_Kfold = function(j){
                                             folds[[j]]
                                           },
                                           get_data = function(j){
                                             train_data = list()
                                             for(i in 1:K){
                                               if( i == j){
                                                 test_data = folds[[i]]
                                               }else{
                                                 train_data = rbind(train_data, folds[[i]])
                                                 
                                               }
                                             }
                                             return(list(train_data = train_data, test_data = test_data))
                                           }
                                         ))

setGeneric("KfoldsCrossValidation", function(data, seed, K) standardGeneric("KfoldsCrossValidation"))
setMethod("KfoldsCrossValidation", signature=c(data="data.frame", seed="integer", K="integer"),
          function(data,seed=27182L,K=10L){
            K <- as.integer(K)
            seed <- as.integer(seed)
            set.seed(seed)
            folds <- list()
            num_obs_kFold <- vector(mode="integer", length=K)
            data = data[sample(1:nrow(data)), ]
            num_data = round(nrow(data)/K)
            num_obs_kFold[1:(K-1)] <- rep(num_data, times=(K-1))
            
            for(i in 1:(K-1)){
              folds[[i]] = data[(1 + num_data*(i-1)):(num_data*i),]
              
            }
            folds[[K]] = data[(num_data*(K-1) + 1):nrow(data), ]
            num_obs_kFold[K] <- nrow(folds[[K]])
            storage.mode(num_obs_kFold) <- "integer" 
            .KfoldsCrossValidationCtr(K=K, data=data, folds=folds, num_obs_kFold=num_obs_kFold)
          })

setMethod("KfoldsCrossValidation", signature=c(data="data.frame", seed="missing", K="missing"),
          function(data,seed=27182L,K=10L){
            KfoldsCrossValidation(data, 27182L,10L)
          })
