#' A LASSO model selection function
#'
#' This function allows to select by crossvalidation one or a serie of different models allowing for best performance. The function is a wrapper for glmnet package.
#' 
#' @param datax A data.frame of explicative variables to be studied of length equal to explicative variable vector
#' @param dependent A dependent variable vector
#' @param laglim The model allows to automatically convert the dataset in order to study lagged dependecies
#' @param study.t0 A logical value indicating whether the t0 period should be taken into account
#' @param folds Number of folds for crossvalidation
#' @param modlim A maximum number of additional models to be presented in output
#' @param varlim A maximal number of variables to figure in supplementary models
#' 
#' @return A list of selected models estimated with OLS
#' \describe{
#'  \item{optimal}{The model for best lambda value, selected by cross validation}
#'  \item{other}{A list of other models in order of goodness of fit with limited regressors} 
#' }
#' @keywords LASSOanalyse
#' @export
#' @examples
#' LASSOanalyse(data, dependent, 
#'      laglim = 5, study.t0 = TRUE, folds = 10, 
#'      modlim = 1, varlim = 10)
#' @import tidyverse
#' @import glmnet
#' 
LASSOanalyse = function(datax, # Explanatory variables dataframe
    dependent, # Dependent variable
    laglim = 5, # Limit of the lag to explore to
    study.t0 = TRUE,
    folds = 10,
    modlim = 1,
    varlim = 10 # Limit of the variables to inlude into final model
    ) {
# Testing
# df = read.csv("./groupe_residentiel/data/derived_data/df3.csv")
# dependent = df[,8]
# datax = log(df[,c(14:18)])
# rm(df)
# laglim = 3 # Lag limit 
# plim = 0.1 # Probability trashhold for tests
# limit = 5
# study.t0 = FALSE
# varlim = 8
    ###############
    # Load packages
    cat("Stage 1 \n")
    #####################
    # Database recreation 
    # Taking lags into account
    cat("Stage 2 \n")
    # Taking into accoun t0
    if (study.t0 == TRUE) {
        # For defined lag order
        if (laglim != 0) {
            database = datax
            for (i in 1:laglim) {
                cat(i, "\n")
                databaseX = datax %>% 
                    mutate_all(function(x) lag(x, i))
                names(databaseX) = 
                    paste0("l", i, "_", names(datax))
                database = cbind(database, databaseX)
                rm(databaseX)
            }  
            rm(datax) ; gc()
        # Only t0 data study
        } else {
            database = datax
            rm(datax) ; gc()
        }
    # No taking into account t0
    } else {
        # Defined lag order
        if (laglim != 0) {
            # First lag base creation
            cat(1, "\n")
            database = datax %>% 
                    mutate_all(function(x) lag(x, n = 1))
            names(database) = 
                    paste0("l", 1, "_", names(datax))
            # More than 1 lag account
            if (laglim > 1) {
                for (i in 2:laglim) {
                    cat(i, "\n")
                    databaseX = datax %>% 
                        mutate_all(function(x) lag(x, n = i))
                    names(databaseX) = 
                        paste0("l", i, "_", names(datax))
                    database = cbind(database, databaseX)
                    rm(databaseX)
                } 
            }
            rm(datax) ; gc()
        # Extra case for conflicts
        } else {
            cat("No laglim defined, 
                proceeding with t0 period study \n")
            database = datax
            rm(datax) ; gc()
        }
    }
    ################
    # GLM LASSO part
    cat("Stage 3 \n")
    # Data
    # NA treatment 
    nas = which(is.na(database), arr.ind=TRUE)
    nas = unique(nas[,1])
    # Subsetting rows without NAs
    databaseX = as.matrix(database[-nas,])
    dependent = as.matrix(dependent[-nas,])
    # LASSO Cross validation
    res_lasso = cv.glmnet(databaseX, dependent,
        nfolds = folds)
    # Optimal lambda by variables limit
    ## Get validation information
    ns = data.frame(
        n = 1:length(res_lasso$lambda),
        lambda = res_lasso$lambda,
        cvm = res_lasso$cvm,
        cvsd = res_lasso$cvsd,
        df = head(res_lasso$glmnet.fit$df, 
            length(res_lasso$lambda))) %>% 
        filter(df <= varlim) %>%
        arrange(cvm) %>%
        head(modlim) %>%
        select(n) %>% 
        as.integer()
    # Select lambda proposed by selection
    lambda_opt = res_lasso$lambda.min
    ##########
    # MCO part
    cat("Stage 4.1 \n")
    # Evaluate OLS for every model
    if (length(ns) != 0) {
        models = list()
        for (i in 1:length(ns)) {
            # Get variables positions
            ox_vars = which(res_lasso$glmnet.fit$beta[,ns[i]] != 0)
            # Results OLS
            models[[i]] = lm(dependent ~ databaseX[,ox_vars])
        }
    } else {
        models = "No models satisfying the indicated criteria"
    }
    # LASSO for optimal lambda
    cat("Stage 4.2 \n")
    opt_vars = which(glmnet(databaseX, dependent,
        lambda = lambda_opt)$beta != 0)
    # Run MCO
    ols_opt = lm(dependent ~ databaseX[,opt_vars])
    ################
    # Return results
    results = list(optimal = ols_opt,
        other = models) 
    # Output
    return(results)
}
