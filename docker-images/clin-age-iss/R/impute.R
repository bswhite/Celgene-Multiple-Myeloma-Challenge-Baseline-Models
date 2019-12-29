
## Utilities to impute missing clinical variables (age, ISS, and gender)
## With the exception of age in the Hovon (training study), the percent
## of missing values is very low.  Hence, use a simple imputation strategy:
## For ISS and gender: impute missing values as the mode
## For age: impute missing values as the mean
## Impute _within_ each data set.  This is because each data set may have
## a different patient population.  Hence, values missing at random
## within a data set will not be missing at random across all data sets.

## impute NA values by applying a function to non-NA values
impute.by.applying.method.to.non.nas <- function(X, impute.col, data.set.col, method = "mean") {
    imputed.vals <- rep(NA, nrow(X))
    data.sets <- unique(X[, data.set.col])
    data.sets <- na.omit(data.sets)
    
    for(data.set in data.sets) {
        data.set.flag <- !is.na(X[, data.set.col]) & ( X[, data.set.col] == data.set )
        na.flag <- data.set.flag & is.na(X[, impute.col])
        non.na.flag <- data.set.flag & !is.na(X[, impute.col])
        if(any(non.na.flag)) {
            val <- do.call(method, list(X[non.na.flag, impute.col]))
            imputed.vals[non.na.flag] <- X[non.na.flag, impute.col]
            imputed.vals[na.flag] <- val
        }
    }
    imputed.vals
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

impute.using.mode <- function(X, impute.col, data.set.col) {
    impute.by.applying.method.to.non.nas(X, impute.col, data.set.col,
                                         method = "getmode")
}

impute.using.mean <- function(X, impute.col, data.set.col) {
    impute.by.applying.method.to.non.nas(X, impute.col, data.set.col,
                                         method = "mean")
}
