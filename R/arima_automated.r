#' ARIMA model choice automation
#'
#' This function allows to select by informational criteria and performance results one or a serie of different models allowing for best performance in and out of sample. The function is a wrapper for forecast and tseries packages.
#' 
#' @param seriex A univariate time serie to be studied
#' @param laglim The maximum difference limit for I part (it is the same for ordinary and seasonal components)
#' @param crossvalid The point for cross-validation
#' @param limit A maximum number of models to be selected by criteria (3 in sample and 3 out of sample criterias are used)
#' @param method The estimation method to be applied (as in the timeserie package)
#' 
#' @return A list of selected models estimated with SARIMA with additional performance information
#' \describe{
#' \item{crit_study}{A data.frame of main in-sample choice criterias for selected models}
#' \item{crit_test}{A data.frame of main out of sample choice criterias for selected models}
#' \item{plot_study}{Autoplot predictions for best in-sample models}
#' \item{plot_test}{Autoplot predictions for best out of sample models}
#' \item{t_stat}{Test statistics for tests}
#' \item{t_pval}{P-values for tests}
#' \item{models}{The results for selected models}
#' \item{plots}{Grid-arranged Acf and Pacf plots for each difference sequence}
#' }
#' @keywords MCOanalyse
#' @export
#' @examples
#' ARIMAanalyse(seriex, 
#'      laglim = 5, crossvalid = c(2015, 4), limit = 3, 
#'      method = "ML")
#' @import tidyverse
#' @import tseries
#' @import ggfortify
#' @import forecast
#' @import car
#' @import lmtest
#' 
ARIMAanalyse = function(seriex, # Time serie
    laglim = 1, # Difference limit
    crossvalid = c(2015, 4), # Validation division point
    limit = 3, # Number of models selected by criteria
    method = "ML") { # Estimation method
    ###########
    # Libraries 
    ###########
    library(tidyverse)
    library(ggfortify)
    library(forecast)
    library(tseries)
    #######################
    # Combinations creation
    #######################
    bin = c(0:laglim)
    comb = crossing(yar = bin, yi = bin, yma = bin,
            sar = bin, si = bin, sma = bin) 
    comb = comb %>% 
        mutate(n = 1:nrow(comb))
    ##############
    # Modelisation
    ##############
    # Auxilary list
    model = list()
    ## seriex = ts_data[,2]
    serie = window(seriex, end = crossvalid)
    # Cycle
    for (i in 1:nrow(comb)) {
        # Order selection
        yearly = as.numeric(comb[i, 1:3])
        season = as.numeric(comb[i, 4:6])
        # Model 
        model[[i]] = Arima(serie, 
            order = yearly,
            seasonal = season,
            method = method)
    }
    # Results extraction 
    # Auxilary frame
    frame = data.frame(AIC = NA, BIC = NA, 
        loglik = NA, pseudoR2 = NA, sigma2 = NA, n = NA)[-1,]
    # Model description
    for (i in 1:nrow(comb)) {
        frame[i,] = frame[i,] %>% 
            mutate(AIC = model[[i]]$aic,
                BIC = model[[i]]$bic, 
                loglik = model[[i]]$loglik,
                sigma2 = model[[i]]$sigma2,
                pseudoR2 = 0,
                n = i)
    }
    # Loglik normalisation
    lago = crossing(yi = bin, si = bin)
    # Loop
    for (i in 1:nrow(lago)) {
        sell = left_join(lago[i,], comb, by = c("yi", "si")) %>% 
            mutate(ko = i, 
                kk = as.integer(yar == 0 &
                    yma == 0 &
                    sar == 0 &
                    sma == 0)) %>% 
            select(n, ko, kk)
        frame[sell$n[sell$kk == 0],] = 
            frame[sell$n[sell$kk == 0],] %>% 
            mutate(pseudoR2 = 
                1 - loglik/frame$loglik[sell$n[sell$kk == 1]])
        frame[sell$n[sell$kk == 1],] = 
            frame[sell$n[sell$kk == 1],] %>% 
            mutate(pseudoR2 = 0)
    }
    # Reorganisation 
    frameX = frame %>% 
        dplyr::arrange(AIC) %>% 
        mutate(ak = 1:nrow(frame)) %>% 
        dplyr::arrange(BIC) %>% 
        mutate(bk = 1:nrow(frame)) %>%
        dplyr::arrange(pseudoR2) %>% 
        mutate(lk = nrow(frame):1) %>%
        dplyr::arrange(sigma2) %>% 
        mutate(sk = 1:nrow(frame)) %>%
        mutate(CC4 = ak + bk + lk + sk,
            CC3 = ak + bk + lk,
            CC2 = ak + bk)
    ################
    # Optimums study
    ################
    # Criteria
    cc2 = frameX %>%
        dplyr::arrange(CC2) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "AIC+BIC")
    ccl = frameX %>%
        dplyr::arrange(lk) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "Log-Lik")
    ccs = frameX %>%
        dplyr::arrange(sigma2) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "Sigma2")
    ccx = rbind(cc2, ccl) %>% 
        rbind(ccs)
    # Presentation results
    cc = left_join(ccx, comb)
    ##########################
    # Testing models forecasts
    ##########################
    # Aux
    forec_list = list()
    res = data.frame(erm = NA, er2m = NA, 
        sigma2 = NA, somme = NA, n = NA)[-1,]
    # Loop
    for (i in 1:length(model)) {
        # Forecasts
        forec = forecast(model[[i]], 
            h = (length(seriex) - length(serie)), 
            level = c(50, 90))
        # Forecasts organisation
        forec_list[[i]] = data.frame(pred = forec$mean,
                orig = seriex[(length(serie)+1):length(seriex)]) %>% 
            mutate(er = pred - orig,
                er2 = er^2)
        # Results
        res[i, ] = data.frame(erm = sum(forec_list[[i]]$er),
            er2m = sum(forec_list[[i]]$er2),
            sigma2 = model[[i]]$sigma2,
            somme = sum(forec_list[[i]]$er2) + model[[i]]$sigma2,
            n = i)
    }
    ##############
    # New criteria
    ##############
    ## Square error
    rsq = res %>%
        dplyr::arrange(er2m) %>%
        head(limit)
    ## Mean error
    rm = res %>%
        dplyr::arrange(pmax(erm, -erm)) %>%
        head(limit)
    ## Overall square error
    rosq = res %>%
        dplyr::arrange(somme) %>%
        head(limit)
    # Combination
    rtot = rbind(rsq, rm) %>% 
        rbind(rosq)
    # Merge 
    rcc = left_join(rtot, comb)
    #######
    # Plots
    #######
    ################
    # Study criteria
    # Aux
    I = 1:length(unique(cc$n))
    x = 1
    # Ploting frame
    toplot1 = data.frame(matrix(NA, 
        ncol = length(I)+1, 
        nrow = length(window(seriex, start = crossvalid))))
    names(toplot1)[2:(length(I)+1)] = 
        paste0("n", unique(as.integer(cc$n)))
    names(toplot1)[1] = "original"
    # toplot1[,x] = seriex
    toplot1[,x] = window(seriex, start = crossvalid)
    # Loop
    for (i in unique(as.integer(cc$n))) {
        x = x + 1
        # Forecasts
        forec = forecast(model[[i]], 
            h = (length(seriex) - length(serie)), 
            level = c(50, 90))
        # series = append(serie, forec$mean)
        # toplot1[, x] = series
        toplot1[, x] = append(serie[length(serie)], forec$mean)
    }
    # Ts conversion
    toplot1 = ts(toplot1, 
        start = attributes(serie)$tsp[2],
        freq = attributes(serie)$tsp[3])
    ###############
    # Test criteria
    # Aux
    I = 1:length(unique(rcc$n))
    x = 1
    # Ploting frame
    toplot2 = data.frame(matrix(NA, 
        ncol = length(I)+1, 
        nrow = length(window(seriex, start = crossvalid))))
    names(toplot2)[2:(length(I)+1)] = 
        paste0("n", unique(as.integer(rcc$n)))
    names(toplot2)[1] = "original"
    # toplot2[,x] = seriex
    toplot2[,x] = window(seriex, start = crossvalid)
    # Loop
    for (i in unique(as.integer(rcc$n))) {
        x = x + 1
        # Forecasts
        forec = forecast(model[[i]], 
            h = (length(seriex) - length(serie)), 
            level = c(50, 90))
        # series = append(serie, forec$mean)
        # toplot2[, x] = series
        toplot2[, x] = append(serie[length(serie)], forec$mean)
    }
    # Ts conversion
    toplot2 = ts(toplot2, 
        start = attributes(serie)$tsp[2],
        freq = attributes(serie)$tsp[3])
    #######+###################    
    # Tests+# Exploratory plots
    #######+###################
    # Packages
    library(car)
    library(lmtest)
    # Lagorder 
    lago = crossing(Y_lag = bin, S_Lag = bin)
    # Frame creation
    tests = data.frame(matrix(NA, 
        nrow = nrow(lago), 
        ncol = 5))
    names(tests) = c("Y_Lag", "S_Lag",
        "ADF", "KPSS.L", "KPSS.T")
    tests[,1] = lago[,1]
    tests[,2] = lago[,2]
    # Loop
    for(i in 0:(nrow(lago)-1)) {
        if (tests[i+1,1] != 0) {
            ser = diff(seriex, 
                lag = attributes(seriex)$tsp[3], 
                dif = tests[i+1,1])
        } else {ser = seriex}
        if (tests[i+1,2] != 0) {
            ser = diff(ser,
                lag = 1, 
                dif = tests[i+1,2])
        } else {ser = ser}
        tests[i+1, 3] = adf.test(ser)$p.value
        tests[i+1, 4] = kpss.test(ser,null = "Level")$p.value
        tests[i+1, 5] = kpss.test(ser,null = "Trend")$p.value
    }
    # Statistics frame creation
    statistics = data.frame(matrix(NA, 
        nrow = nrow(lago), 
        ncol = 6))
    names(statistics) = c("Y_Lag", "S_Lag",
        "ADF", "KPSS.L", "KPSS.T",
        "DW")
    statistics[,1] = lago[,1]
    statistics[,2] = lago[,2]
    plots = list()
    # Loop
    for(i in 0:(nrow(lago)-1)) {
        if (statistics[i+1,1] != 0) {
            ser = diff(seriex, 
                lag = attributes(seriex)$tsp[3], 
                dif = tests[i+1,1])
        } else {ser = seriex}
        if (statistics[i+1,2] != 0) {
            ser = diff(ser,
                lag = 1, 
                dif = tests[i+1,2])
        } else {ser = ser}
        statistics[i+1, 3] = adf.test(ser)$statistic
        statistics[i+1, 4] = kpss.test(ser,null = "Level")$statistic
        statistics[i+1, 5] = kpss.test(ser,null = "Trend")$statistic
        statistics[i+1, 6] = durbinWatsonTest(as.integer(ser))
        # Plots
        placf = ggAcf(ser) + 
            ggtitle(paste("Differences :", tests[i+1,1], tests[i+1,2], sep = " "))
        plpacf = ggPacf(ser) + 
            ggtitle(paste("Differences :", tests[i+1,1], tests[i+1,2], sep = " "))
        plots[[i+1]] = grid.arrange(placf, plpacf, nrow = 1)
    }
    #########
    # Results
    #########    
    results = list(crit_study = cc, 
        crit_test = rcc,
        plot_study = toplot1, 
        plot_test = toplot2,
        t_stat = statistics,
        t_pval = tests,
        models = model,
        plots = plots)
    # Output
    return(results)
}
