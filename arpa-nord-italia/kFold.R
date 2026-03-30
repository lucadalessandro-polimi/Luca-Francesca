source("kFoldCV.R")
library(Matrix)


#data = matrix(rnorm(100), ncol=2)
#colnames(data) = c("X1", "X2")
#data = as.data.frame(data)
locs <- read.csv("data/2022.10.789_space_locs.csv", row.names=1)
obs  <- read.csv("data/2022.10.789_PM10.csv",  row.names=1)

data <- data.frame(
  V1 = locs$V1,
  V2 = locs$V2,
  y = obs$y
)

kFold = KfoldsCrossValidation(data) # k = 10 default
print(colnames(data))


errors <- numeric(10)

for(k in 1:10){

  cat("Fold ", k, "\n")
  split_data <- kFold$get_data(k)

  train_df <- split_data$train_data
  test_df  <- split_data$test_data

  # Ricostruiamo locs e obs separati
  locs_train <- train_df[, c("V1","V2")]
  obs_train  <- data.frame(y = train_df$y)

  locs_test  <- test_df[, c("V1","V2")]
  obs_test   <- data.frame(y = test_df$y)

  write.csv(locs_train, "kFold_data/train_locs.csv")
  write.csv(obs_train,  "kFold_data/train_obs.csv")

  write.csv(locs_test,  "kFold_data/test_locs.csv")
  write.csv(obs_test,   "kFold_data/test_obs.csv")

  pred_file <- paste0("kFold_data/pred_test_fold", k, ".csv")

  system(paste(
    "./PM10_space_only_north_italy",
    "mesh_6",
    "kFold_data/train_locs.csv",
    "kFold_data/train_obs.csv",
    "kFold_data/test_locs.csv",
    pred_file
  ))


  pred <- as.numeric(readMM(pred_file))

  errors[k] <- mean((pred - obs_test$y)^2)
}

write.csv(errors, "kFold_data/errors_mesh6.csv", row.names=FALSE)
