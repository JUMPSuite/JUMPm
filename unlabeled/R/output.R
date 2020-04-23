findMatchedFeatures = function(matched, calibrated, files, samples) {
  ## Caution
  ## Input arguments files and samples should have the same order as defined in the parameter file
  nFiles = length(files)
  
  refRunNo = NULL
  for (i in 1:length(matched)) {
    if (is.null(dim(matched[[i]]))) {
      refRunNo = i
    }
  }
  
  ## Extract the reference run's indexes of fully- and partially matched features
  fullInd = NULL
  partialInd = NULL
  for (i in 1:length(matched)) {
    if (i == refRunNo) {
      next
    } else {
      if (is.null(fullInd)) {
        fullInd = matched[[i]][, 1]
        partialInd = matched[[i]][, 1]
      } else {
        fullInd = intersect(fullInd, matched[[i]][, 1])
        partialInd = union(partialInd, matched[[i]][, 1])
      }
    }
  }
  fullInd = sort(unique(fullInd))
  partialInd = setdiff(partialInd, fullInd)
  partialInd = sort(unique(partialInd))
  
  ###############################
  ## 1. Fully-matched features ##
  ###############################
  features = list()
  runNo = setdiff(seq(1, length(matched)), refRunNo)[1]
  nCols = ncol(matched[[runNo]]) / 2
  for (i in 1:length(matched)) {
    if (i == refRunNo) {
      tmp = matched[[runNo]][order(matched[[runNo]][, 1]), ]
      tmp = tmp[tmp[, 1] %in% fullInd, ]
      features = c(features, list(tmp[, 1:nCols]))
    } else {
      tmp = matched[[i]][order(matched[[i]][, 1]), ]
      tmp = tmp[tmp[, 1] %in% fullInd, ]
      features = c(features, list(tmp[, (nCols + 1):(2 * nCols)]))
    }
  }
  features = as.data.frame(features, check.names = F)

  ## Manipulation of column names
  header = NULL
  for (i in 1:length(files)) {
    filename = basename(files[i])
    filename = gsub(".\\d+.feature", "", filename)
    header = c(header, rep(filename, nCols))
  }
  colnames(features) = gsub("ref_", "", colnames(features))
  colnames(features) = gsub("comp_", "", colnames(features))
  colnames(features) = paste0(header, "_", colnames(features))
  fullFeatures = features
  
  ####################################
  ## 2. Partially-matched features  ##
  ####################################
  features = data.frame(matrix(NA, nrow = length(partialInd), ncol = nCols * nFiles))
  for (i in 1:length(partialInd)) {
    for (j in 1:length(matched)) {
      if (j == refRunNo) {
        next
      } else {
        ind = which(matched[[j]][, 1] == partialInd[i])
        if (length(ind) > 0) {
          features[i, ((j - 1) * nCols + 1) : (j * nCols)] = matched[[j]][ind, (nCols + 1) : (2 * nCols)]
          features[i, ((refRunNo - 1) * nCols + 1) : (refRunNo * nCols)] = matched[[j]][ind, 1 : nCols]
        } else if (length(ind) > 1) {
          stop ('Error in partially-matched feature identification')
        }
      }
    }
  }

  ## Manipulation of column names
  header = NULL
  for (i in 1:length(files)) {
    filename = basename(files[i])
    filename = gsub(".\\d+.feature", "", filename)
    header = c(header, rep(filename, nCols))
  }
  colnames(features) = rep(colnames(calibrated[[1]]), nFiles)
  colnames(features) = paste0(header, "_", colnames(features))
  partialFeatures = features
  
  ############################
  ## 3. Un-matched features ##
  ############################
  unmatchedFeatures = list()
  for (i in 1:length(files)) {
    filename = basename(files[i])
    filename = gsub(".\\d+.feature", "", filename)
    if (i == refRunNo) {
      unInd = unique(setdiff(calibrated[[i]]$index, union(fullInd, partialInd)))
      # unmatchedFeatures[[fileName]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
    } else {
      unInd = setdiff(calibrated[[i]]$index, matched[[i]][, (nCols + 1)])
      # unmatchedFeatures[[fileName]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
    }
    unmatchedFeatures[[filename]] = calibrated[[i]][calibrated[[i]]$index %in% unInd, ]
  }
  
  ############
  ## Output ##
  ############
  res = list(fullInd = fullInd, fullFeatures = fullFeatures, 
             partialInd = partialInd, partialFeatures = partialFeatures,
             unmatchedFeatures = unmatchedFeatures)
  return (res)
}
