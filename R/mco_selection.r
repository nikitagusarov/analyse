#' MCO automation function : Parallelized approach
#'
#' This function allows to select by informational criteria and various tests results (Durbin-Watson, Shapiro-Wilks, Breush-Pagan, F-test) one or a serie of different models allowing for best performance in-sample. The function uses parallelised approach with doParallel package.
#' 
#' @param datax A data.frame of explicative variables to be studied of length equal to explicative variable vector
#' @param dependent A dependent variable vector
#' @param laglim The model allows to automatically convert the dataset in order to study lagged dependecies
#' @param plim The critical value for tests results (default value 0.1)
#' @param limit A maximum number of models to be selected in each category
#' @param study.t0 A logical value indicating whether the t0 period should be taken into account
#' @param varlim A maximal number of variables to make combinations of (large number significantly increases computation time, each supplementary variable increases the tyme by 2)
#' @param valid LOGICAL, indicates whether the crossvalidation procedure should be applied during the model selection 
#' @param part The part of the sample to be used for crossvalidation (by default 15%)
#' 
#' @return A list of selected models estimated with OLS with additional performance information
#' \describe{
#'  \item{choices}{A data.frame containing information for selected models testing results}
#'  \item{models}{A list of models in order as they figure in choices dataframe}
#'  \item{univariate}{Results of first stage univariate model testing for selected variables}
#' }
#' @keywords MCOanalyse
#' @export
#' @examples
#' MCOanalyse(data, dependent, 
#'      laglim = 5, plim = 0.1, limit = 2, 
#'      study.t0 = TRUE, folds = 10, varlim = 10)
#' @import tidyverse
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import filehash
#' @import lmtest
#' @import mctest
#' 
MCOanalyse = function(datax, # Dataframe with ALL variables
        # Must start from t = 0 to t = n
        # head - t0, tail - tn
        # If needed to account for lagged dependent, include it
    dependent, # Name of the dependent variable
    laglim = 2, # Lag limit 
    plim = 0.1, # Probability trashhold for tests
    limit = 1, # Number of models selected limit
    study.t0 = FALSE, # If the model should account for t0
    varlim = 10, # Maximum variables number to study
        # Highter is this number, highter is complexity and time : t = t0*2^n
        ## Avec t0 le temps d'execution pour carlim = 1
        # Advised not to use more than 10 
        # (15? ~ 2 hours computing time)
    valid = TRUE, # If the models should be cross validated
    part = 0.15 # Number of crossvalidation folds
    ) { 
# # Testing
# df = read.csv("./groupe_residentiel/data/derived_data/df3.csv")
# dependent = df[,8]
# datax = log(df[,c(14:18)])
# rm(df)
# laglim = 3 # Lag limit 
# plim = 0.1 # Probability trashhold for tests
# limit = 5
# study.t0 = FALSE
# varlim = 8
    #################
    # Parallelisation 
    cat("Stage 1 \n")
    # Get cores 
    ncores = detectCores()
    registerDoParallel(ncores) 
    # clust = makeCluster(ncores)
    # registerDoSNOW(clust)
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
    #####################################
    # Univariate selection over many lags
    # Loop
    tbase = foreach (i = 1:ncol(database),
        .combine = rbind) %dopar% {
        # Data 
        X = database[,i]
        # Model 
        summar = summary(lm(dependent ~ X)) # %>% 
            #summary()
        data.frame(n = i,
                p.tstat = summar$coeff[2,4],
                r2 = summar$r.squared,
                p.fstat = pf(summar$fstat[1],
                    summar$fstat[2],
                    summar$fstat[3],
                    lower.tail=FALSE),
                varname = colnames(database)[i])
    }
    # Selection of variables by varlim
    ## univar_tbase = tbase
    univar_tbase = tbase %>% 
        filter(p.tstat < plim &
            p.fstat < plim) 
    rm(tbase) ; gc()
    univar_tbase = univar_tbase %>% 
        arrange(-r2)  %>% 
        head(varlim)
    vars = as.integer(univar_tbase$n)
    # Complexity reduction achieved
    #######################
    # Database recreation 2
    database = database[,vars]
    ##############
    # Combinations
    cat("Stage 3 \n")
    # Combinations storage
    dbCreate("DataBase")
    Base = dbInit("DataBase")
    Base$combin = vector("list", length = ncol(database))
    # Loop
    Base$combin = foreach (i = 1:ncol(database)) %dopar% {
        # cat(i, "\n")
        combn(length(names(database)), i, 
            simplify = TRUE)
    }
    # Load DB
    dbLoad(Base)
    # Count possible models 
    nmodels = 0
    for (i in 1:length(Base$combin)) {
        nmodels = nmodels + ncol(Base$combin[[i]])
    }
    ########################
    # Models and performance 
    cat("Stage 4 \n")
    model = list()
    # Results frames
    # Model results
    cat("Stage 4.1 \n")
    if (valid == TRUE) {
        resbase = data.frame(n = NA, i = NA, j = NA,
            rse = NA, r2 = NA, mse = NA,
            p.fstat = NA, 
            p.dw = NA, p.sw = NA, p.bp = NA, k2.fg = NA)[-1,]
    } else {
        resbase = data.frame(n = NA, i = NA, j = NA,
            rse = NA, r2 = NA,
            p.fstat = NA, 
            p.dw = NA, p.sw = NA, p.bp = NA, k2.fg = NA)[-1,]
    }
    # Loop modelisation
    cat("Stage 4.2 \n")
    # x = 1 # Counter
    if (valid == TRUE) {
        # Cross-validation
        for (i in 1:length(Base$combin)) {
            # Iteration number
            cat("i = ", i, "\n")
            # Combination preselection
            # Variable combinations
            Vars = Base$combin[[i]] # May cause memory overuse
            Varnames = colnames(database)
            # Limits
            varlength = ncol(Base$combin[[i]])
            rowlength = nrow(database)
            # Select study part
            study = 
            foreach (k = 1:as.integer(rowlength*part)) %do% {
                sample(rowlength, 
                    as.integer(rowlength*(1-part)))
            }
            # Data
            resbaseX = 
            foreach (j = 1:varlength,
                .packages = c("foreach"),
                .combine = rbind,
                .inorder = FALSE #,
                # .options.snow = opts
                ) %dopar% {
                # Crossvalidation
                mods = foreach (k = 1:length(study),
                .packages = c("lmtest", "mctest"),
                .combine = rbind
                ) %do% {
                    # Study 
                    X = data.frame(database[study[[k]],Vars[,j]])
                    colnames(X) = Varnames[Vars[,j]]
                    Y = data.frame(dependent[study[[k]]])
                    colnames(Y) = "y"
                    df = cbind(Y, X)
                    rm(Y, X)
                    # Modelisation
                    model = lm(y ~ ., df)
                    # Test out of sample
                    Xt = as.data.frame(database[-study[[k]],Vars[,j]]) 
                    colnames(Xt) = Varnames[Vars[,j]]
                    Yt = data.frame(O = dependent[-study[[k]]], 
                        P = predict(model, newdata = Xt))
                    Yt = na.omit(Yt)
                    rm(Xt)
                    r2t = mean((Yt$O - Yt$P)^2) 
                    # Summary
                    summar = summary(model)
                    # fgts = ifelse((length(Vars)-1) > 1, 
                    #     omcdiag(df[,-1], df[,1])$odiags[2,1], 0)
                    # Output
                    data.frame(
                        sup_rse = summar$sigma,
                        sup_r2 = summar$r.squared,
                        sup_se = r2t,
                        sup_p.fstat = pf(summar$fstat[1],
                            summar$fstat[2],
                            summar$fstat[3],
                            lower.tail=FALSE),
                        sup_p.dw = dwtest(model, 
                            alternative = "two.sided")$p.value,
                        sup_p.sw = shapiro.test(model$res)$p.value,
                        sup_p.bp = bptest(model)$p.value,
                        sup_k2.fg = ifelse((ncol(df)-1) > 1, 
                            omcdiag(df[,-1], df[,1])$odiags[2,2], 0))
                }
                # Results 
                ## Performance
                data.frame(n = NA, 
                    i = i, j = j,
                    rse = mean(mods$sup_rse),
                    r2 = mean(mods$sup_r2),
                    mse = mean(mods$sup_se),
                    p.fstat = mean(mods$sup_p.fstat),
                    p.dw = mean(mods$sup_p.dw),
                    p.sw = mean(mods$sup_p.sw),
                    p.bp = mean(mods$sup_p.bp),
                    k2.fg = mean(mods$sup_k2.fg))
            }
            resbase = rbind(resbase, resbaseX)
        }
    } else {
        # No crossvalidation
        for (i in 1:length(Base$combin)) {
            # Iteration number
            cat("i = ", i, "\n")
            # Combination preselection
            # Variable combinations
            Vars = Base$combin[[i]] # May cause memory overuse
            # Limits
            varlength = ncol(Base$combin[[i]])
            rowlength = nrow(database)
            Varnames = colnames(database)
            # Data
            resbaseX = 
            foreach (j = 1:varlength,
                .packages = c("lmtest", "mctest"),
                .combine = rbind,
                .inorder = FALSE #,
                # .options.snow = opts
                ) %dopar% {
                # Explicative variables 
                X = data.frame(database[,Vars[,j]])
                colnames(X) = Varnames[Vars[,j]]
                Y = data.frame(dependent)
                colnames(Y) = "y"
                df = cbind(Y, X)
                rm(Y, X)
                # Modelisation
                model = lm(y ~ ., df)
                # Sumary
                summar = summary(model)
                # Results 
                ## Performance
                data.frame(n = NA, 
                    i = i, j = j,
                    rse = summar$sigma,
                    r2 = summar$r.squared,
                    p.fstat = pf(summar$fstat[1],
                        summar$fstat[2],
                        summar$fstat[3],
                        lower.tail=FALSE),
                    p.dw = dwtest(model, 
                        alternative = "two.sided")$p.value,
                    p.sw = shapiro.test(model$res)$p.value,
                    p.bp = bptest(model)$p.value,
                    k2.fg = ifelse((ncol(df)-1) > 1, 
                        omcdiag(df[,-1], df[,1])$odiags[2,2], 0))
            }
            resbase = rbind(resbase, resbaseX)
        }
    }
    #####################
    # Choice and analysis
    cat("Stage 5 \n")
    # Aux
    crit = list()
    types = list()
    # Optimal
    crit[[1]] = "p.fstat < plim &
        p.dw > plim &
        p.sw < plim &
        p.bp > plim &
        k2.fg < plim"
    types[[1]] = "optimal"
    # Autocorrelation
    crit[[2]] = "p.fstat < plim &
        p.dw <= plim & 
        p.sw < plim & 
        p.bp > plim &
        k2.fg < plim"
    types[[2]] = "autocor"
    # Non-normality
    crit[[3]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw < plim & 
        p.bp <= plim &
        k2.fg < plim"
    types[[3]] = "non.norm"
    # Heteroscedastisity
    crit[[4]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw >= plim & 
        p.bp > plim &
        k2.fg < plim"
    types[[4]] = "heterosk"
    # Multicollinearity
    crit[[5]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw < plim & 
        p.bp > plim &
        k2.fg >= plim"
    types[[5]] = "multicol"
    # Autocorrelation and non-normality
    crit[[6]] = "p.fstat < plim &
        p.dw <= plim & 
        p.sw < plim & 
        p.bp <= plim &
        k2.fg < plim"
    types[[6]] = "autocor_non.norm"
    # Autocorrelation and heteroscedasticity
    crit[[7]] = "p.fstat < plim &
        p.dw <= plim & 
        p.sw >= plim & 
        p.bp > plim &
        k2.fg < plim"
    types[[7]] = "autocor_heterosk"
    # Autocorrelation and multicollinearity
    crit[[8]] = "p.fstat < plim &
        p.dw <= plim & 
        p.sw < plim & 
        p.bp > plim &
        k2.fg >= plim"
    types[[8]] = "autocor_multicol"
    # Non-normality and Heterockedasticity
    crit[[9]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw >= plim & 
        p.bp <= plim &
        k2.fg < plim"
    types[[9]] = "non.norm_hetero"
    # Non-normality and Multicollinearity
    crit[[10]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw < plim & 
        p.bp <= plim &
        k2.fg >= plim"
    types[[10]] = "non.norm_multicol"
    # Heteroscedastisity and Multicollinearity
    crit[[11]] = "p.fstat < plim &
        p.dw > plim & 
        p.sw >= plim & 
        p.bp > plim &
        k2.fg >= plim"
    types[[11]] = "hetero_multicol"
    # Eval results
    if (valid == TRUE) {
        results = foreach (s = 1:length(crit),
            .packages = "tidyverse",
            .combine = rbind) %dopar% {
            plim = plim
            resbase %>% 
                dplyr::filter(eval(parse(text = crit[[s]]))) %>%
                arrange(mse) %>%
                head(limit) %>%
                mutate(type = types[[s]])
        } 
        results = results %>%
            arrange(mse) %>%
            mutate(m = 1:nrow(results))
    } else {
        results = foreach (s = 1:length(crit),
            .packages = "tidyverse",
            .combine = rbind) %dopar% {
            plim = plim
            resbase %>% 
                dplyr::filter(eval(parse(text = crit[[s]]))) %>%
                arrange(r2) %>%
                tail(limit) %>%
                mutate(type = types[[s]])
        } 
        results = results %>%
            arrange(-r2) %>%
            mutate(m = 1:nrow(results))
    }
    # Remove aux 
    rm(resbase) ; gc(full = TRUE)
    #########################
    # Selected models storage
    # Loop
    model = foreach (k = 1:nrow(results)) %dopar% {
        # Explicative variables
        vars = Base$combin[[results$i[k]]][,results$j[k]]
        # Explicative variables 
        X = data.frame(database[,vars])
        colnames(X) = Varnames[vars]
        Y = data.frame(dependent)
        colnames(Y) = "y"
        df = cbind(Y, X)
        rm(Y, X)
        # Models
        model[[k]] = lm(y ~ ., df)
    }
    #################
    # Remove database
    rm(Base)
    unlink("DataBase")
    # Close cluster
    stopImplicitCluster() ; gc()
    # stopCluster(clust) ; gc()
    ##################
    # Output reduction
    ret = list(choices = results,
        models = model,
        univariate = univar_tbase)
    return(ret)
}