
# Define the functions
# F1 score
F1 <- function(mat){
  apply(mat, 1, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    2*TP/(2*TP+FP+FN)
  })
}
# Sensitivity score
Sen <- function(mat){
  apply(mat, 1, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    TP/(TP+FN)
  })
}

# specificity score
Spe <- function(mat){
  apply(mat, 1, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    TN/(TN+FP)
  })
}

# accuracy score
Acc <- function(mat){
  apply(mat, 1, function(x){
    TN <- x[1]
    FP <- x[2]
    TP <- x[3]
    FN <- x[4]
    (TP+TN)/(TP+TN+FP+FN)
  })
}

