rescueFeatures = function (ref, comp, dfAligned, rtSd, mzSd, rtTolScale, rtTol, mzTolScale, mzTol) {
  nUnaligned = dim(comp)[1]
  cat(paste0("    There are ", nUnaligned, " unaligned features\n"))
  cat("    The features satisfying the following conditions are going to be reviewed for further alignment\n")
  cat(paste0("    There are ", nUnaligned, " unaligned features\n"), file = LOG)
  cat("    The features satisfying the following conditions are going to be reviewed for further alignment\n", file = LOG)
  
  ## 1. Reduce the "reference" run by selecting "unaligned" features
  ##    rtTol and mzTol should be also reduced
  nRescue = 0
  indUnaligned = setdiff(ref$index, dfAligned$"ref_index")
  ind = which(ref$index %in% indUnaligned)
  ref = ref[ind, ]
  rtSd = rtSd[ind]
  mzSd = mzSd[ind]
  
  ## 2. Apply the selection criteria
  ##    - Absolute intensity level: hard-coded (grand median intensity of aligned features)
  ##    - Intensity-ratio between the ref- and comp-runs: hard-coded, within 95% of the ratios of aligned features)
  ##    - RT- and m/z-shifts should be within specified tolerances (e.g. 10SD or 10ppm)
  resIntLevel = median(as.matrix(dfAligned[, grep("Intensity", colnames(dfAligned))]))
  resRatioPct = 95 ## Hard-coded
  intRatioAligned = log(dfAligned$"comp_Intensity", 2) - log(dfAligned$"ref_Intensity", 2)
  lRatio = quantile(intRatioAligned, (1 - resRatioPct / 100) / 2)
  uRatio = quantile(intRatioAligned, 1 - (1 - resRatioPct / 100) / 2)
  cat(paste0("    - Intensity higher than ", resIntLevel, " (median intensity of aligned features)\n"))
  cat(paste0("    - Ratio of intensities within ", resRatioPct, "% of ratios from aligned features\n"))
  cat(paste0("    - Intensity higher than ", resIntLevel, " (median intensity of aligned features)\n"), file = LOG)
  cat(paste0("    - Ratio of intensities within ", resRatioPct, "% of ratios from aligned features\n"), file = LOG)
  if (rtTolScale == 1) { # times of SD unit
    cat(paste0("    - RT-shifts within ", rtTol, " x SD of estimated RT-shifts from aligned features\n"))
    cat(paste0("    - RT-shifts within ", rtTol, " x SD of estimated RT-shifts from aligned features\n"), file = LOG)
    rtTol = rtTol * rtSd
  } else if (rtTolScale == 2) { # second unit
    cat(paste0("    - RT-shifts less than ", rtTol, " seconds\n"))
    cat(paste0("    - RT-shifts less than ", rtTol, " seconds\n"), file = LOG)
    rtTol = rep(rtTol, length(ind))
  } else {
    cat("  - WARNING: check your parameter for RT-tolerance unit. It should be either 1 or 2\n")
    cat("  - Due to the incorrect RT-tolerance unit parameter, the rescue step is skipped\n")
    cat("  - WARNING: check your parameter for RT-tolerance unit. It should be either 1 or 2\n", file = LOG)
    cat("  - Due to the incorrect RT-tolerance unit parameter, the rescue step is skipped\n", file = LOG)
    return (dfAligned)
  }
  if (mzTolScale == 1) { # times of SD unit
    cat(paste0("    - m/z-shifts within ", mzTol, " x SD of estimated m/z-shifts from aligned features\n"))
    cat(paste0("    - m/z-shifts within ", mzTol, " x SD of estimated m/z-shifts from aligned features\n"), file = LOG)
    mzTol = mzTol * mzSd
  } else if (mzTolScale == 2) { # ppm unit
    cat(paste0("    - m/z-shifts less than ", mzTol, " ppm\n"))
    cat(paste0("    - m/z-shifts less than ", mzTol, " ppm\n"), file = LOG)
    mzTol = rep(mzTol, length(ind))
  } else {
    cat("  - WARNING: check your parameter for m/z-tolerance unit. It should be either 1 or 2\n")
    cat("  - Due to the incorrect m/z-tolerance unit parameter, the rescue step is skipped\n")
    cat("  - WARNING: check your parameter for m/z-tolerance unit. It should be either 1 or 2\n", file = LOG)
    cat("  - Due to the incorrect m/z-tolerance unit parameter, the rescue step is skipped\n", file = LOG)
    return (dfAligned)
  }
  
  for (i in 1:dim(ref)[1]) {
    intRatio = log(comp$Intensity, 2) - log(ref$Intensity[i], 2)
    mzShift = (comp$mz - ref$mz[i]) / comp$mz * 1e6
    rtShift = comp$RT - ref$RT[i]
    selInd = which(comp$Intensity > resIntLevel &
                     intRatio >= lRatio & intRatio <= uRatio &
                     abs(mzShift) <= mzTol[i] &
                     abs(rtShift) <= rtTol[i])
    if (length(selInd) > 0) {
      ## Check charge states of the matched features
      charge = ref$z[i]
      if (charge == 0) {
        selInd = selInd[1]
        rescued = cbind(ref[i, ], comp[selInd, ])
        names(rescued) = names(dfAligned)
        dfAligned = rbind(dfAligned, rescued)
        comp = comp[-selInd, ]
        nRescue = nRescue + 1
      } else {
        if (comp$z[selInd[1]] == 0) {
          selInd = selInd[1]
          rescued = cbind(ref[i, ], comp[selInd, ])
          names(rescued) = names(dfAligned)
          dfAligned = rbind(dfAligned, rescued)
          comp = comp[-selInd, ]
          nRescue = nRescue + 1
        } else {
          chargeInd = which(comp$z[selInd] == charge)
          if (length(chargeInd) > 0) {
            selInd = selInd[chargeInd[1]]
            rescued = cbind(ref[i, ], comp[selInd, ])
            names(rescued) = names(dfAligned)
            dfAligned = rbind(dfAligned, rescued)
            comp = comp[-selInd, ]
            nRescue = nRescue + 1
          }
        }
      }
    } ## end if (length(selInd) > 0) 
  } ## end for
  cat(paste0("    Through the rescue procedure, ", nRescue, " features are additionally aligned\n\n"))
  cat(paste0("    Through the rescue procedure, ", nRescue, " features are additionally aligned\n\n"), file = LOG)
  res = list(ref = ref, comp = comp, dfAligned = dfAligned)
  return (res)
}

matchFeatures = function (ref, comp, refName, compName, rtSd, mzSd, params, LOG) {
  sdWidth = params$sdWidth
  rescue = params$rescue
  resSdWidth = params$sdWidth_rescue
  # resIntLevel = params$intensityLevel_rescue
  #resRatioPct = params$pctIntensityRatio_rescue
  
  ## Sort features by their intensities (to consider high intensity features first)
  ref = ref[order(ref$Intensity, decreasing = T), ]
  comp = comp[order(comp$Intensity, decreasing = T), ]
  
  n = dim(ref)[1]
  nc = dim(comp)[1]
  cat(paste0("  ", refName, ":", n, " features (reference run)\n"))
  cat(paste0("  ", compName, ":", nc, " features (compared run)\n"))
  cat(paste0("  ", refName, ":", n, " features (reference run)\n"), file = LOG)
  cat(paste0("  ", compName, ":", nc, " features (compared run)\n"), file = LOG)
  
  dfAligned = data.frame()
  dfUnaligned = data.frame()
  if (length(rtSd) == 1) {
    rtSd = rep(rtSd, n)
  }
  if (length(mzSd) == 1) {
    mzSd = rep(mzSd, n)
  }
  rtSd[rtSd == 0] = min(rtSd[rtSd > 0]) ## Prevent zero-tolerance of RT
  mzSd[mzSd == 0] = min(mzSd[mzSd > 0]) ## Prevent zero-tolerance of m/z
  rtTol = sdWidth * rtSd
  mzTol = sdWidth * mzSd
  
  j = 1
  mzShifts = NULL
  rtShifts = NULL
  for (i in 1:n) {
    mz = ref$mz[i]
    rt = ref$RT[i]
    intensity = ref$Intensity[i]
    charge = ref$z[i]
    rtErr = comp$RT - rt
    mzErr = (comp$mz - mz) / comp$mz * 1e6
    rowInd = which(abs(rtErr) <= rtTol[i] & abs(mzErr) <= mzTol[i])
    
    ## When there is/are matched feature(s)
    if (length(rowInd) > 0) {
      ## Check charge states of the matched features
      ## 0 charge state can be matched to any charge state
      if (charge == 0) {
        rowInd = rowInd[1]
        dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
        comp = comp[-rowInd, ]
      } else {
        if (comp[rowInd[1], ]$z == 0) { ## When the compared features has charge 0, it is okay
          rowInd = rowInd[1]
          dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
          comp = comp[-rowInd, ]
        } else {
          chargeInd = which(comp$z[rowInd] == charge)
          if (length(chargeInd) > 0) {
            rowInd = rowInd[chargeInd[1]]
            dfAligned = rbind(dfAligned, cbind(ref[i, ], comp[rowInd, ]))
            comp = comp[-rowInd, ]
          }
        }
      }
    }
  }
  colnames(dfAligned) = c(paste0("ref_", names(ref)), paste0("comp_", names(comp)))
  nAligned = dim(dfAligned)[1]
  cat(paste0("    ", nAligned, " features are aligned between runs\n"))
  cat(paste0("    ", nAligned, " features are aligned between runs\n"), file = LOG)
  
  ######################################
  ## Rescue some "unaligned" features ##
  ######################################
  if (rescue == 1) {
  	rtTolUnit = params$rtTolUnit
  	mzTolUnit = params$mzTolUnit
  	rtTol = params$rtTolVal
  	mzTol = params$mzTolVal
  	for (i in 1:length(rtTolUnit)) {
  		res = rescueFeatures(ref, comp, dfAligned, rtSd, mzSd, rtTolUnit[i], rtTol[i], mzTolUnit[i], mzTol[i])
  		ref = res$ref
  		comp = res$comp
  		dfAligned = res$dfAligned
  	}
  }
  return (dfAligned)
}
