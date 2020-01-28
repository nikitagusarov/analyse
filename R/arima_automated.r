#' ARIMA model choice automation
#'
#' This function allows to select by informational criteria and performance results one or a serie of different models allowing for best performance in and out of sample. The function is a wrapper for forecast and tseries packages.
#' 
#' @param seriex A univariate time serie to be studied
#' @param order The maximum limit for p, d and q for ordinary ARIMA part
#' @param seasonallim The maximum limit for P, D and Q for SARIMA part
#' @param crossvalid The point for cross-validation
#' @param limit A maximum number of models to be selected by criteria (3 in sample and 3 out of sample criterias are used)
#' @param method The estimation method to be applied (as in the timeserie package)
#' @param adf.plim ADF selection probability limit
#' @param partlim Part of autocorrelated lags in lags
#' 
#' @return A list of selected models estimated with SARIMA with additional performance information
#' 
#' @keywords ARIMAanalyse
#' @export
#' @examples
#' ARIMAanalyse = function(seriex, 
#'    orderlim = c(1, 1, 1), seasonallim = c(1, 1, 1), 
#'    crossvalid = c(2015, 4), limit = 3, 
#'    method = "ML", adf.plim = 0.1, partlim = 0.3)
#'    
#' @import tidyverse
#' @import tseries
#' @import ggfortify
#' @import gridExtra
#' @import forecast
#' @import car
#' @import lmtest
#' 
ARIMAanalyse = function(seriex, # Time serie
    orderlim = c(1, 1, 1), # Difference limit
    seasonallim = c(1, 1, 1), # Seasonal part limit
    crossvalid = c(2015, 4), # Validation division point
    limit = 3, # Number of models selected by criteria
    method = "ML", # Estimation method
    adf.plim = 0.1, # ADF selection probability limit
    partlim = 0.3 # Part of autocorrelated lags in lags
    ) { 
    #######+###################    
    # Tests+# Exploratory plots
    #######+###################
    # Lagorder 
    lago = crossing(d = 0:orderlim[2], D = 0:seasonallim[2])
    # Frame creation
    tests = data.frame(matrix(NA, 
        nrow = nrow(lago), 
        ncol = 5))
    names(tests) = c("d", "D",
        "ADF", "KPSS.L", "KPSS.T")
    tests[,1] = lago[,1]
    tests[,2] = lago[,2]
    # Statistics frame creation
    statistics = data.frame(matrix(NA, 
        nrow = nrow(lago), 
        ncol = 6))
    names(statistics) = c("d", "D",
        "ADF", "KPSS.L", "KPSS.T",
        "DW")
    statistics[,1] = lago[,1]
    statistics[,2] = lago[,2]
    # Correlation 
    acorrelation = list()
    pcorrelation = list()
    # Plots
    cfplots = list()
    # Loop
    for (i in 0:(nrow(lago)-1)) {
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
        # Tests P-VALUES dataframe
        tests[i+1, 3] = adf.test(ser)$p.value
        tests[i+1, 4] = kpss.test(ser,null = "Level")$p.value
        tests[i+1, 5] = kpss.test(ser,null = "Trend")$p.value
        statistics[i+1, 3] = adf.test(ser)$statistic
        # Tests STATISTICS dataframe
        statistics[i+1, 4] = kpss.test(ser,null = "Level")$statistic
        statistics[i+1, 5] = kpss.test(ser,null = "Trend")$statistic
        statistics[i+1, 6] = durbinWatsonTest(as.integer(ser))
        # Values
        s_level = 
            qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(ser)))
        acorrelation[[i+1]] = 
            sum(abs(acf(ser, plot = F)$acf)/s_level > 1)/as.integer(10*log10(length(ser)))
        pcorrelation[[i+1]] = 
            sum(abs(pacf(ser, plot = F)$acf)/s_level > 1)/as.integer(10*log10(length(ser)))
        # CFs
        cf = ggAcf(ser, plot = F)
        pcf = ggPacf(ser, plot = F)
        # Plots
        placf = autoplot(cf) + 
            ggtitle(paste("Differences :", tests[i+1,1], tests[i+1,2], sep = " "))
        plpacf = autoplot(pcf) + 
            ggtitle(paste("Differences :", tests[i+1,1], tests[i+1,2], sep = " "))
        cfplots[[i+1]] = grid.arrange(placf, plpacf, nrow = 1)
    }
    ##############################
    # Analyse and preselect models
    ##############################
    # Selection by correlation  
    nn = intersect(which(acorrelation < partlim), 
        which(pcorrelation < partlim))
    # Selection by test results
    an = tests[nn,] %>% 
        filter(ADF < adf.plim & 
            (KPSS.L >= 0.1 | KPSS.T >= 0.1)) %>%
        select(d, D)
    #######################
    # Combinations creation
    #######################
    dbin = as.integer(an$d)
    Dbin = as.integer(an$D)
    comb = crossing(p = c(0:orderlim[1]), 
        d = dbin, 
        q = c(0:orderlim[3]),
        P = c(0:seasonallim[1]), 
        D = Dbin, 
        Q = c(0:seasonallim[3]))
    if (nrow(comb) == 0) {
        results = list(tests = list(pval = tests,
                stat = statistics), 
            plots = cfplots)
        return(results)
        stop("No satisfying difference orders detected")
        }
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
    lago = crossing(d = dbin, D = Dbin)
    # Loop
    for (i in 1:nrow(lago)) {
        sell = left_join(lago[i,], comb, by = c("d", "D")) %>% 
            mutate(ko = i, 
                kk = as.integer(p == 0 &
                    q == 0 &
                    P == 0 &
                    Q == 0)) %>% 
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
    cca = frameX %>%
        dplyr::arrange(AIC) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "AIC")
    ccb = frameX %>%
        dplyr::arrange(BIC) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "BIC")
    # ccl = frameX %>%
    #     dplyr::arrange(lk) %>%
    #     head(limit) %>% 
    #     select(n) %>%
    #     mutate(TYPE = "Log-Lik")
    ccs = frameX %>%
        dplyr::arrange(sigma2) %>%
        head(limit) %>% 
        select(n) %>%
        mutate(TYPE = "Sigma2")
    ccx = rbind(cca, ccb) %>% 
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
    #########
    # Results
    #########   
    crit = list(study = cc, 
        test = rcc) 
    plotx = list(study = toplot1, 
        test = toplot2)
    tests = list(stat = statistics,
        pval = tests)
    results = list(crit = crit,
        predict = plotx,
        tests = tests,
        models = model,
        cfplots = cfplots)
    # Output
    return(results)
}