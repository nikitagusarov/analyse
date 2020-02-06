#' Seasonal data.frames optimisation function
#'
#' This function allows to optimise the time series seasonal data, based on the stationarity tests' results.
#' It is an extention of test.df() function.
#' 
#' @param datax A data.frame of seasonal data to be studied of length equal to explicative variable vector. The first column is advised to be year and the second to be "M" for monthly data, or "T" for trimestrial data.
#' @param season An indicator of seasonality type to be used as present in the data.frame.
#' 
#' @return A data.frame with test results by variable.
#' 
#' @keywords optimise.df test.df
#' @export
#' @examples
#' optimise.df = function(datax, season)
#' 
#' @import tidyverse
#' @import tseries
#' 
optimise.df = function(datax, season) {
    # Aux
    df = 
        data.frame(matrix(NA, nrow = 3, ncol = ncol(datax)))
    # Identity
    colnames(df) = colnames(datax)
    rownames(df) = c("p-ADF", "p-KPSS.L", "p-KPSS.T")
    # Tests
    df[1,] = datax %>% 
        lapply(function(x) tseries::adf.test(na.omit(x))$p.val) %>%
        unlist()
    df[2,] = datax %>% 
        lapply(function(x) tseries::kpss.test(na.omit(x), null = "Level")$p.val) %>%
        unlist()
    df[3,] = datax %>% 
        lapply(function(x) tseries::kpss.test(na.omit(x), null = "Trend")$p.val) %>%
        unlist()
    # Seasonality test
    if (season == "M") {
        m = datax %>% 
            group_by(M) %>%
            summarise_all(function(x) mean(na.omit(x)))
        v = datax %>% 
            group_by(M) %>%
            summarise_all(function(x) var(na.omit(x)))
        GM = datax %>% 
            select(-M) %>%
            summarise_all(function(x) mean(na.omit(x)))
        GM = GM[rep(1, nrow(m)),]
        tm = round((m[,2:ncol(m)]-GM)/v[,2:ncol(v)], 4)
        tm = cbind(tm, M = m$M)
    } else if (season == "T") {
        m = datax %>% 
            group_by(T) %>%
            summarise_all(function(x) mean(na.omit(x)))
        v = datax %>% 
            group_by(T) %>%
            summarise_all(function(x) var(na.omit(x)))
        GM = datax %>% 
            select(-T) %>%
            summarise_all(function(x) mean(na.omit(x)))
        GM = GM[rep(1, nrow(m)),]
        tm = round((m[,2:ncol(m)]-GM)/v[,2:ncol(v)], 4)
        tm = cbind(tm, T = m$T)
    } else {
        return(df) 
        warning("No seasonal component to analyse")
    }
    df = dplyr::union(df, tm)
    # Transform
    tr = datax %>% 
        mutate(t = 1:nrow(datax)) %>%
        na.omit()
    for (i in 3:ncol(df)) {
        if (df[2,i] < 0.1) {
            tr[,i] = tr[,i] - lm(tr[,i] ~ tr$t)$fit
            colnames(tr)[i] = paste0("d_", colnames(tr)[i])
        } else if (df[3,i] < 0.1) {
            tr[,i] = tr[,i] - mean(tr[,i])
            colnames(tr)[i] = paste0("L_", colnames(tr)[i])
        } else if (df[1,i] > 0.1) {
            tr[,i] = tr[,i] - lag(tr[,i])
            colnames(tr)[i] = paste0("T_", colnames(tr)[i])
        }   
    }
    # Results
    res = list(test = df, data = tr)
    return(res)
}