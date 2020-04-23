featureCalibration = function (ref, comp, params, LOG) {
    initMzTol = params$initMzTol
    sdWidth = params$sdWidth
    
    ref = ref[order(ref$Intensity, decreasing = TRUE), ] ## Sort by intensity
    comp = comp[order(comp$Intensity, decreasing = TRUE), ] ## Sort by intensity
    
    ########################################
    ## Calibration of RT and m/z globally ##
    ########################################
    ## Calculate rtShift and mzShift in a global manner, and calibrate features in the compared sample
    cat("  Global calibration of RT and m/z is being performed\n")
    cat("  Global calibration of RT and m/z is being performed\n", file = LOG)
    res = globalCalibration(ref, comp, initMzTol)
    rtShifts = res$rtShifts ## unit of second
    mzShifts = res$mzShifts ## unit of ppm
    cat(sprintf("    Based on the matched features within +/- %i ppm\n", initMzTol))
    cat(sprintf("    The global RT-shift is %.4f second\n", median(rtShifts)))
    cat(sprintf("    The global m/z-shift is %.4f ppm\n", median(mzShifts)))
    cat("    RT and m/z of the compared sample is calibrated according to the global RT- and m/z-shift\n")
    cat(sprintf("    Based on the matched features within +/- %i ppm\n", initMzTol), file = LOG)
    cat(sprintf("    The global RT-shift is %.4f second\n", median(rtShifts)), file = LOG)
    cat(sprintf("    The global m/z-shift is %.4f ppm\n", median(mzShifts)), file = LOG)
    cat("    RT and m/z of the compared sample is calibrated according to the global RT- and m/z-shift\n", file = LOG)
    comp$RT = comp$RT - median(rtShifts)
    comp$mz = comp$mz / (1 + median(mzShifts)/ 1e6)
    
    ##########################################
    ## Calibration of RT using LOESS curve  ##
    ##########################################
    cat("  Local calibration of RT and m/z is being performed (through LOESS modeling)\n")
    cat("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows\n")
    cat("      RT- and m/z-tolerance =", sdWidth, "x dynamically estimated SD of RT- and m/z-shifts\n")
    cat("    LOESS modeling may take some time. Please be patient ...\n")
    cat("  Local calibration of RT and m/z is being performed (through LOESS modeling)\n", file = LOG)
    cat("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows\n", file = LOG)
    cat("      RT- and m/z-tolerance =", sdWidth, "x dynamically estimated SD of RT- and m/z-shifts\n", file = LOG)
    cat("    LOESS modeling may take some time. Please be patient ...\n", file = LOG)
    ## 1st round of LOESS modeling and calibration of RT
    ## When performing this LOESS modeling,
    ## the matched features are selected based on the global rt- and mz-shift
    rtSd = max(sd(rtShifts), 1e-3)
    mzSd = max(sd(mzShifts), 1e-3)
    res = localCalibration(ref, comp, rtSd, mzSd, params, "RT")
    estRtShift = predict(res$model, data.frame(x = comp$RT))
    comp$RT = comp$RT - estRtShift
    cat("    The 1st round of RT-calibration is done\n")
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n")
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n")
    cat("    The 1st round of RT-calibration is done\n", file = LOG)
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n", file = LOG)
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n", file = LOG)
    
    ## 2nd round of LOESS modeling and calibration of RT
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, params, "RT")
    estRtShift = predict(res$model, data.frame(x = comp$RT))
    comp$RT = comp$RT - estRtShift
    cat("    The 2nd round of RT-calibration is done\n")
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n")
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n")
    cat("    The 2nd round of RT-calibration is done\n", file = LOG)
    cat("      min SD of RT-shifts =", round(min(res$dynRtSd), 4), "second\n", file = LOG)
    cat("      max SD of RT-shifts =", round(max(res$dynRtSd), 4), "second\n", file = LOG)
    
    ##########################################
    ## Calibration of m/z using LOESS curve ##
    ##########################################
    ## 1st round of LOESS modeling and calibration of m/z
    ## When performing this LOESS modeling,
    ## the matched features are selected based on the global rt- and mz-shift
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, params, "m/z")
    estMzShift = predict(res$model, data.frame(x = comp$mz))
    comp$mz = comp$mz / (1 + estMzShift / 1e6)
    cat("    The 1st round of m/z-calibration is done\n")
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n")
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n")
    cat("    The 1st round of m/z-calibration is done\n", file = LOG)
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n", file = LOG)
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n", file = LOG)
    
    ## 2nd round of LOESS modeling and calibration of m/z
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = localCalibration(ref, comp, rtSd, mzSd, params, "m/z")
    estMzShift = predict(res$model, data.frame(x = comp$mz))
    comp$mz = comp$mz / (1 + estMzShift / 1e6)
    cat("    The 2nd round of m/z-calibration is done\n")
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n")
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n")
    cat("    The 2nd round of m/z-calibration is done\n", file = LOG)
    cat("      min SD of m/z-shifts =", round(min(res$dynMzSd), 4), "ppm\n", file = LOG)
    cat("      max SD of m/z-shifts =", round(max(res$dynMzSd), 4), "ppm\n", file = LOG)
    
    ## Output organization
    rtSd = res$dynRtSd
    mzSd = res$dynMzSd
    res = list(calibratedFeature = comp, rtSd = rtSd, mzSd = mzSd)
    return (res)
}

globalCalibration = function(ref, comp, mzTol = 20) {
    nPeaks = round(0.05 * dim(ref)[1]) ## Number of peaks "matched" between a reference and a compared sample
    i = 1 ## row index of peaks in the reference sample. The 1st peak is the strongest peak
    j = 1 ## index of the matched peaks
    subDf = data.frame()
    while (j <= nPeaks) {
        z = ref$z[i]
        mz = ref$mz[i]
        rt = ref$RT[i]
        intensity = ref$Intensity[i]
        if (z == 0) {
            ## For a reference feature with undetermined charge, consider all possible charges
            rowInd = which(comp$mz >= mz - mz * mzTol / 1e6 & comp$mz < mz + mz * mzTol / 1e6)
        } else {
            ## For a reference feature with a specific charge, only consider features with the same charge
            rowInd = which(comp$mz >= mz - mz * mzTol / 1e6 & comp$mz < mz + mz * mzTol / 1e6 & comp$z == z)
        }
        if (length(rowInd) > 0) {
            rowInd = rowInd[1]
            subDf[j, 1] = mz
            subDf[j, 2] = rt
            subDf[j, 3] = intensity
            subDf[j, 4] = comp$mz[rowInd]
            subDf[j, 5] = comp$RT[rowInd]
            subDf[j, 6] = comp$Intensity[rowInd]
            comp = comp[-rowInd, ]
            j = j + 1
        }
        i = i + 1
    }
    colnames(subDf) = c("refMz", "refRT", "refIntensity", "compMz", "compRT", "compIntensity")
    
    ## Calculate the global shift and its standard deviation using
    ## 10% trimmed rt- and mz-shift values
    rtShifts = subDf$"compRT" - subDf$"refRT"
    rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
    rtShifts = rtShifts[rtInd]
    mzShifts = (subDf$"compMz" - subDf$"refMz") / subDf$"compMz" * 1e6
    mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
    mzShifts = mzShifts[mzInd]
    
    res = list(rtShifts = rtShifts, mzShifts = mzShifts)
    return (res)
}

localCalibration = function(ref, comp, rtSd, mzSd, params, cal) {
    ## For each feature the reference sample, look for the corresponding feature
    ## to be aligned in the compared samples
    subDf = data.frame()
    n = dim(ref)[1]
    if (length(rtSd) == 1) {
        rtSd = rep(rtSd, n)
    }
    if (length(mzSd) == 1) {
        mzSd = rep(mzSd, n)
    }
    sdWidth = params$sdWidth
    rtTol = sdWidth * rtSd
    mzTol = sdWidth * mzSd
    j = 1
    for (i in 1:n) {
        z = ref$z[i]
        mz = ref$mz[i]
        rt = ref$RT[i]
        intensity = ref$Intensity[i]
        rtErr = comp$RT - rt
        mzErr = (comp$mz - mz) / comp$mz * 1e6 ## unit of ppm
        if (z == 0) {
            ## For a reference feature with undetermined charge, consider all possible charges
            rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i])
        } else {
            ## For a reference feature with a specific charge, only consider features with the same charge
            rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i] & comp$z == z)
        }
        if (length(rowInd) > 0) {
            ## When multiple peaks correspond to a peak in the reference sample
            ## choose one (with the highest intensity)
            ## Note that the data.frame(comp) is already sorted according to its intensity values
            rowInd = rowInd[1]
            subDf[j, 1] = mz
            subDf[j, 2] = rt
            subDf[j, 3] = intensity
            subDf[j, 4] = comp$mz[rowInd]
            subDf[j, 5] = comp$RT[rowInd]
            subDf[j, 6] = comp$Intensity[rowInd]
            j = j + 1
        }
    }
    colnames(subDf) = c("refMz", "refRT", "refIntensity", "compMz", "compRT", "compIntensity")
    
    if (cal == "RT") {
        ## Build a LOESS model for the RT-shifts between the reference and compared samples
        compRt = subDf$"compRT"
        refRt = subDf$"refRT"
        rtShifts = compRt - refRt
        if (sum(rtShifts == 0) == length(rtShifts)) {
            rtShifts = 1e-6 * rnorm(length(rtShifts))
        }
        mod = loess.as(compRt, rtShifts, degree = 1, criterion = "aicc",
                       control = loess.control(surface = "direct"))
        
        ## Calculate a new rt-tolerance using trimming
        compRt = compRt - mod$fitted ## Calibration according to the LOESS model
        rtShifts = (compRt - refRt)
        rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
        modRtSd = loess.as(compRt[rtInd], rtShifts[rtInd] ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynRtSd = sqrt(pmax(0, predict(modRtSd, data.frame(x = ref$RT))))
        statRtSd = sd(rtShifts[rtInd])
        
        ## Calculate a new mz-tolerance using trimming
        ## Sometimes, the variation of m/z-shifts cannot be captured when trimming is applied
        ## So, the trimming is not used for m/z-shifts
        # mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
        mzShifts = (subDf$"compMz" - subDf$"refMz") / subDf$"compMz" * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        modMzSd = loess.as(subDf$"compMz", mzShifts ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynMzSd = sqrt(pmax(0, predict(modMzSd, data.frame(x = ref$mz))))
        statMzSd = sd(mzShifts)
    } else if (cal == "m/z") {
        ## Build a LOESS model for the mz-shifts between the reference and compared samples
        compMz = subDf$"compMz"
        refMz = subDf$"refMz"
        mzShifts = (compMz - refMz) / compMz * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        mod = loess.as(compMz, mzShifts, degree = 1, criterion = "aicc",
                       control = loess.control(surface = "direct"))
        
        ## Calculate a new mz-tolerance using trimming
        ## Sometimes, the variation of m/z-shifts cannot be captured when trimming is applied
        ## So, the trimming is not used for m/z-shifts
        # mzInd = which(mzShifts >= quantile(mzShifts, 0.1) & mzShifts <= quantile(mzShifts, 0.9))
        compMz = compMz * (1 + mod$fitted / 1e6) ## Calibration according to the LOESS model
        mzShifts = (compMz - refMz) / compMz * 1e6
        if (sum(mzShifts == 0) == length(mzShifts)) {
            mzShifts = 1e-6 * rnorm(length(mzShifts))
        }
        modMzSd = loess.as(subDf$"compMz", mzShifts ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynMzSd = sqrt(pmax(0, predict(modMzSd, data.frame(x = ref$mz))))
        statMzSd = sd(mzShifts)
        
        ## Calculate a new rt-tolerance using trimming
        rtShifts = (subDf$"compRT" - subDf$"refRT")
        if (sum(rtShifts == 0) == length(rtShifts)) {
            rtShifts = 1e-6 * rnorm(length(rtShifts))
        }
        rtInd = which(rtShifts >= quantile(rtShifts, 0.1) & rtShifts <= quantile(rtShifts, 0.9))
        modRtSd = loess.as(subDf$"compRT"[rtInd], rtShifts[rtInd] ^ 2, degree = 1, criterion = "aicc",
                           control = loess.control(surface = "direct"))
        dynRtSd = sqrt(pmax(0, predict(modRtSd, data.frame(x = ref$RT))))
        statRtSd = sd(rtShifts[rtInd])
    }
    
    res = list(model = mod, dynRtSd = dynRtSd, statRtSd = statRtSd,
               dynMzSd = dynMzSd, statMzSd = statMzSd)
    return (res)
}

loess.as = function(x, y, degree=1, criterion=c("aicc", "gcv"), 
                    family = c("gaussian", "symmetric"), user.span=NULL, plot=FALSE, ...) {
    criterion <- match.arg(criterion)
    family <- match.arg(family)
    x <- as.matrix(x)
    
    if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
    if (!is.numeric(x)) stop("argument 'x' must be numeric!")
    if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
        stop("argument 'user.span' must be a numerical number!")
    if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
    if(length(y) < 3) stop("not enough observations!")
    
    data.bind <- data.frame(x=x, y=y)
    if (ncol(x) == 1) {
        names(data.bind) <- c("x", "y")
    } else { names(data.bind) <- c("x1", "x2", "y") }
    
    opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
        as.crit <- function (x) {
            span <- x$pars$span
            traceL <- x$trace.hat
            sigma2 <- sum(x$residuals^2 ) / (x$n-1)
            aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
            gcv <- x$n*sigma2 / (x$n-traceL)^2
            result <- list(span=span, aicc=aicc, gcv=gcv)
            return(result)
        }
        criterion <- match.arg(criterion)
        fn <- function(span) {
            mod <- update(model, span=span)
            as.crit(mod)[[criterion]]
        }
        result <- optimize(fn, span.range)
        return(list(span=result$minimum, criterion=result$objective))
    }
    
    if (ncol(x)==1) {
        if (is.null(user.span)) {
            fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
            span1 <- opt.span(fit0, criterion=criterion)$span
        } else {
            span1 <- user.span
        }		
        fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
    } else {
        if (is.null(user.span)) {
            fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
            span1 <- opt.span(fit0, criterion=criterion)$span
        } else {
            span1 <- user.span
        }		
        fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
    }
    if (plot){
        if (ncol(x)==1) {
            m <- 100
            x.new <- seq(min(x), max(x), length.out=m)
            fit.new <- predict(fit, data.frame(x = x.new))
            plot(x, y, col="lightgrey", xlab="x", ylab="m(x)", ...)
            lines(x.new,fit.new, lwd=1.5, ...)
        } else {
            m <- 50
            x1 <- seq(min(data.bind$x1), max(data.bind$x1), len=m) 
            x2 <- seq(min(data.bind$x2), max(data.bind$x2), len=m) 
            x.new <- expand.grid(x1=x1, x2=x2) 
            fit.new <- matrix(predict(fit, x.new), m, m) 
            persp(x1, x2, fit.new, theta=40, phi=30, ticktype="detailed", xlab="x1", ylab="x2", zlab="y", col="lightblue", expand=0.6)
        }		
    }
    return(fit)
}
