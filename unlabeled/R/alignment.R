rm(list = ls())

################
## Subroutine ##
################
parseParams = function (paramFile) {
    lines = readLines(paramFile)
    params = NULL
    for (i in 1:length(lines)) {
        lines[i] = gsub("\\s", "", lines[i])
        if (grepl("^#", lines[i]) | !grepl("=", lines[i])) {
            next
        } else {
            line_i = gsub("#.*", "", lines[i])
            key = unlist(strsplit(line_i, "="))[1]
            val = unlist(strsplit(line_i, "="))[2]
        }
        if (key == "output_name") {
            params$output_name = as.character(val)
        } else if (key == "mode") {
            params$mode = as.numeric(val)
        } else if (key == "tol_initial") {
            params$initMzTol = as.numeric(val)
        } else if (key == "sd_width") {
            params$sdWidth = as.numeric(val)
        } else if (key == "rescue") {
            params$rescue = as.numeric(val)
        } else if (key == "sd_width_rescue") {
            params$sdWidth_rescue = as.numeric(val)
        } else if (key == "intensity_level_rescue") {
            params$intensityLevel_rescue = as.numeric(val)
        } else if (key == "pct_intensity_ratio_rescue") {
            params$pctIntensityRatio_rescue = as.numeric(val)
        } else if (key == "pct_full_alignment") {
            params$pctFullyAlignment = as.numeric(val)
        } else if (key == "rt_tolerance_unit") {
            val = gsub(" ", "", val)
            params$rtTolUnit = as.numeric(unlist(strsplit(val, ",")))
        } else if (key == "rt_tolerance_value") {
            val = gsub(" ", "", val)
            params$rtTolVal = as.numeric(unlist(strsplit(val, ",")))
        } else if (key == "mz_tolerance_unit") {
            val = gsub(" ", "", val)
            params$mzTolUnit = as.numeric(unlist(strsplit(val, ",")))
        } else if (key == "mz_tolerance_value") {
            val = gsub(" ", "", val)
            params$mzTolVal = as.numeric(unlist(strsplit(val, ",")))
        } else if (key == "reference_feature") {
            val = gsub(" ", "", val)
            params$referenceFeature = as.character(val)
        } else if (key == "skip_loading_bias_correction") {
            val = gsub(" ", "", val)
            params$skipLoadingBiasCorrection = as.numeric(val)
        }
    }
    return (params)
}

##################
## Main routine ##
##################

## Parameters
##  output_name: prefix of output files containing fully- and partially-aligned features
##  initMzTol: initial m/z-tolerance for the global calibration (default = 20 ppm)
##  sdWidth: SD-width for RT- and m/z-tolerances (when finding "matched/aligned" features, default = 5)
##  rescue: whether or not to rescue unaligned features by loosening RT- and m/z-tolerances (default = 1 (yes))
##  sdWidth_rescue: SD-width for RT- and m/z-tolerances (when rescuing "unmatched/unaligned" features, default = 10)
##  intensityLevel_rescue: absolute intensity threshold of the features to be considered in the rescue step (default = 1e5)
##  pctIntensityRatio_rescue: percentage threshold of intensity-ratios of the aligned features (default = 95)
##                            the unaligned features whose intensity-ratios are within this percentage (of the ratios of aligned ones)
##                            will be considred in the rescue step
##  pctFullyAlignment: percentage of samples are grouped into a feature (e.g. 100% = feature is formed by peaks present in all samples)
##  files: array of input .feature files

############################
## Parse input arguments  ##
############################
startTime = Sys.time()
args = commandArgs(trailingOnly = TRUE)
paramFile = args[1]
filenames = args[2]
logFile = args[3]
LOG = file(logFile, "a")
srcDirectory = args[4]
outDirectory = args[5]
params = parseParams(paramFile)
files = unlist(strsplit(filenames, ","))

# ## For testing in a desktop
# startTime = Sys.time()
# paramFile = "../IROAsamples/jumpm_negative.params"
# # filenames = "../IROAsamples/IROA_c18_target1.1.feature"
# filenames = "../IROAsamples/IROA_IS_NEG_1.1.feature,../IROAsamples/IROA_IS_NEG_2.1.feature,../IROAsamples/IROA_IS_NEG_3.1.feature"
# # filenames = "../IROAsamples/old/IROA_IS_NEG_1.1.feature,../IROAsamples/old/IROA_IS_NEG_2.1.feature,../IROAsamples/old/IROA_IS_NEG_3.1.feature"
# srcDirectory = "U:/Research/Projects/7Metabolomics/JUMPm"
# outDirectory = "."
# logFile = "tmplog"
# LOG = file(logFile, 'a')
# params = parseParams(paramFile)
# files = unlist(strsplit(filenames, ","))

#############################################
## Read feature files generated from JUMPm ##
#############################################
features = list()
for (i in 1:length(files)) {
    df = read.table(files[i], header = T, sep = "\t", row.names = NULL,
                    stringsAsFactors = F, comment.char = "", check.names = F)
    colnames(df) = gsub(" ", "", colnames(df))
    colnames(df) = gsub("/", "", colnames(df))
    features[[i]] = df
}

############################################################
## Alignment and so forth when there are multiple samples ##
############################################################
if (length(files) > 1) {
    cat("\n\n  Feature calibration\n")
    cat("  ===================\n\n")
    cat("\n\n  Feature calibration\n", file = LOG)
    cat("  ===================\n\n", file = LOG)
    
    ################################
    ## Select a reference sample  ##
    ################################
    if (params$referenceFeature == "0") {
        ## A sample with the largest median of top 100 intensities is set to a reference run
        refNo = 1
        refIntensity = 0
        for (i in 1:length(features)) {
            tmpIntensity = median(sort(features[[i]]$Intensity, decreasing = T)[1:100])
            if (tmpIntensity >= refIntensity) {
                refNo = i
                refIntensity = tmpIntensity
            }
        }
    } else {
        refNo = which(files == params$referenceFeature)
    }
    cat(paste0("  ", basename(files[refNo]), " is chosen as the reference run\n"))
    cat(paste0("  ", basename(files[refNo]), " is chosen as the reference run\n"), file = LOG)
    
    ##############################################################
    ## Calibration of features against those in a reference run ##
    ##############################################################
    source(paste0(srcDirectory, "/R/calibration.R"))
    calibratedFeatures = list()
    rtSd = list()
    mzSd = list()
    for (i in 1:length(files)) {
        if (i == refNo) {
            calibratedFeatures[[i]] = features[[i]]
            rtSd[[i]] = NA
            mzSd[[i]] = NA
        } else {
            cat(paste0("\n  ", basename(files[i]), " is being aligned against the reference run (it may take a while)\n"))
            cat(paste0("\n  ", basename(files[i]), " is being aligned against the reference run (it may take a while)\n"), file = LOG)
            res = featureCalibration(features[[refNo]], features[[i]], params, LOG)
            calibratedFeatures[[i]] = res$calibratedFeature
            rtSd[[i]] = res$rtSd
            mzSd[[i]] = res$mzSd
        }
    }
    cat("\n")
    cat("\n", file = LOG)
    
    ########################
    ## Calibration summary  ##
    ########################
    cat("  Calibration summary\n")
    cat("    After calibration, RT- and m/z-shifts of each run (against a reference run) are centered to zero\n")
    cat("    Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows\n")
    cat("    Filename\t\t# features\tSD RT-shift[second]\tSD m/z-shift[ppm]\n")
    cat("  Calibration summary\n", file = LOG)
    cat("    After calibration, RT- and m/z-shifts of each run (against a reference run) are centered to zero\n", file = LOG)
    cat("    Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows\n", file = LOG)
    cat("    Filename\t\t# features\tSD RT-shift[second]\tSD m/z-shift[ppm]\n", file = LOG)
    for (i in 1:length(files)) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        nFeatures = nrow(features[[i]])
        if (i == refNo) {
            meanRtSd = "NA (reference)"
            meanMzSd = "NA (reference)"
        } else {
            meanRtSd = sprintf("%.6f", mean(rtSd[[i]]))
            meanMzSd = sprintf("%.6f", mean(mzSd[[i]]))
        }
        cat(paste0("    ", filename, "\t\t", nFeatures, "\t\t", meanRtSd, "\t", meanMzSd, "\n"))
        cat(paste0("    ", filename, "\t\t", nFeatures, "\t\t", meanRtSd, "\t", meanMzSd, "\n"), file = LOG)
    }
    cat("\n")
    cat("\n", file = LOG)
    
    ################################
    ## Identify aligned features  ##
    ################################
    cat("  Feature alignment\n")
    cat("  =================\n\n")
    cat("  Feature alignment\n", file = LOG)
    cat("  =================\n\n", file = LOG)
    source(paste0(srcDirectory, "/R/matchFeatures.R"))
    matchedFeatures = list()
    for (i in 1:length(files)) {
        if (i == refNo) {
            matchedFeatures[[i]] = NA
        } else {
            matchedFeatures[[i]] = matchFeatures(calibratedFeatures[[refNo]], calibratedFeatures[[i]],
                                                 basename(files[refNo]), basename(files[i]),
                                                 rtSd[[i]], mzSd[[i]], params, LOG)
        }
    }
    
    ##############################################
    ## Summarization and organization of output ##
    ##############################################
    ## Fully-, partially- and un-matched/aligned features
    source(paste0(srcDirectory, "/R/output.R"))
    res = findMatchedFeatures(matchedFeatures, calibratedFeatures, files)
    
    ##############################
    ## Alignment/match summary  ##
    ##############################
    cat("  Alignment/matching summary\n")
    cat("    After alignment/feature matching, fully-, partially- and un-aligned features are as follows\n")
    cat("    Filename\t\tFully-aligned\tPartially-aligned\tun-aligned\n")
    cat("  Alignment/matching summary\n", file = LOG)
    cat("    After alignment/feature matching, fully-, partially- and un-aligned features are as follows\n", file = LOG)
    cat("    Filename\t\tFully-aligned\tPartially-aligned\tun-aligned\n", file = LOG)
    for (i in 1:length(files)) {
        filename = basename(files[i])
        filename = gsub(".\\d+.feature", "", filename)
        nFull = length(res$fullInd)
        nPartial = sum(!is.na(res$partialFeatures[[paste0(filename, "_index")]]))
        nUnmatched = dim(res$unmatchedFeatures[[filename]])[1]
        cat(paste0("    ", filename, "\t\t", nFull, "\t\t", nPartial, "\t\t", nUnmatched, "\n"))
        cat(paste0("    ", filename, "\t\t", nFull, "\t\t", nPartial, "\t\t", nUnmatched, "\n"), file = LOG)
    }
    cat("\n")
    cat(paste0("    There are ", nrow(res$partialFeatures), " partially-aligned features\n"))
    cat("\n", file = LOG)
    cat(paste0("    There are ", nrow(res$partialFeatures), " partially-aligned features\n"), file = LOG)
    colInd = grep("index", colnames(res$partialFeatures))
    nRuns = apply(!is.na(res$partialFeatures[, colInd]), 1, sum)
    for (i in (length(files) - 1):2) {
        cat(paste0("    Partially-aligned features over ", i, " runs = ", sum(nRuns == i), "\n"))
        cat(paste0("    Partially-aligned features over ", i, " runs = ", sum(nRuns == i), "\n"), file = LOG)
    }
    
    ## Depending on the parameter, some "partially-aligned" features are merged to "fully-aligned" ones
    if (!is.null(params$pctFullyAlignment)) {
        pctFullyAlignment = params$pctFullyAlignment
        if (pctFullyAlignment < 100) {
            colInd = grep("index", colnames(res$partialFeatures))
            nRuns = apply(!is.na(res$partialFeatures[, colInd]), 1, sum)
            rowInd = which(nRuns >= ceiling(pctFullyAlignment / 100 * length(files)))	  
            if (length(rowInd) > 0) {
                res$fullFeatures = rbind(res$fullFeatures, res$partialFeatures[rowInd, ])
                res$fullInd = c(res$fullInd, res$partialInd[rowInd])
                res$partialFeatures = res$partialFeatures[-rowInd, ]
                res$partialInd = res$partialInd[-rowInd]
                cat("\n")
                cat("    According to the parameter setting,", length(rowInd), "partially-aligned features are regarded as fully-aligned\n")
                cat("\n", file = LOG)
                cat("    According to the parameter setting,", length(rowInd), "partially-aligned features are regarded as fully-aligned\n", file = LOG)
            } else {
                cat("\n")
                cat("    According to the parameter setting, no feature is added to the set of fully-aligned ones\n")
                cat("\n", file = LOG)
                cat("    According to the parameter setting, no feature is added to the set of fully-aligned ones\n", file = LOG)
            }
        }
    }
    cat("\n")
    cat("\n", file = LOG)
    
    ##############################
    ## Processing quantity data ##
    ##############################
    cat("  Quantity information\n")
    cat("  ====================\n")
    cat("  Quantity information\n", file = LOG)
    cat("  ====================\n", file = LOG)
    fullFeatures = res$fullFeatures
    expr = fullFeatures[, grep("Intensity", colnames(fullFeatures))]
    
    ## Calculation and printout of loading-biases (filter out the most variable features +/- 10%)
    cat("  Loading-bias summary\n")
    cat("  Loading-bias summary\n", file = LOG)
    colInd = grep("Intensity", colnames(expr))
    sampleName = colnames(expr)[colInd]
    sampleName = gsub("_Intensity", "", sampleName)
    lexpr = expr[, colInd] ## Feature intensity table for loading-bias calculation
    nFeatures = dim(lexpr)[1]
    nSamples = dim(lexpr)[2]
    sampleMeans = rowMeans(lexpr, na.rm = T)
    lexpr = log2(lexpr / sampleMeans)
    rowInd = seq(1, nFeatures, by = 1)
    for (i in 1:nSamples) {
        rowInd = intersect(rowInd, which(lexpr[, i] < quantile(lexpr[, i], 0.9, na.rm = T) & lexpr[, i] > quantile(lexpr[, i], 0.1, na.rm = T)))
    }
    meanIntensity = round(as.numeric(2 ^ colMeans(lexpr[rowInd, ]) * 100), 2)
    sdVal = as.numeric(apply(lexpr[rowInd, ], 2, sd))
    sdIntensity = round(as.numeric(((2 ^ sdVal - 1) + (1 - 2 ^ (-sdVal))) / 2 * 100), 2)
    semIntensity = round(as.numeric(sdIntensity / sqrt(length(rowInd))), 2)
    cat("    Samplename\tMean[%]\tSD[%]\tSEM[%]\t#features\n")
    cat("    Samplename\tMean[%]\tSD[%]\tSEM[%]\t#features\n", file = LOG)
    for (i in 1:nSamples) {
        cat("   ", sampleName[i], "\t", meanIntensity[i], "\t", sdIntensity[i], "\t", semIntensity[i], "\t", length(rowInd), "\n")
        cat("   ", sampleName[i], "\t", meanIntensity[i], "\t", sdIntensity[i], "\t", semIntensity[i], "\t", length(rowInd), "\n", file = LOG)
    }
    
    ## Normalization (i.e. loading-bias correction) based on trimmed-mean intensities
    expr = log2(expr)
    rowSel = NULL
    if (params$skipLoadingBiasCorrection == 0) {
        ## Parameters for normalization
        intensityThreshold_normalization = quantile(as.matrix(expr), 0.2, na.rm = T) ## 10% percentile of overall intensities
        trimPct_normalization = 0.1 ## 10% lowest and 10% highest intensities are to be trimmed
        rowSel = expr > intensityThreshold_normalization ## Pre-filtering based on the intensity level
        for (i in 1:ncol(expr)) {
            rowSel[, i] = (expr[, i] > quantile(expr[, i], trimPct_normalization, na.rm = T) &
                               expr[, i] < quantile(expr[, i], 1 - trimPct_normalization, na.rm = T))
        }
        rowInd = which(apply(rowSel, 1, sum) == ncol(rowSel))
        meanIntensity = colMeans(expr[rowInd, ])
        normFactor = meanIntensity - mean(meanIntensity)
        expr = expr - rep(normFactor, each = nrow(expr))
    }
    
    ## Replace NAs in expr with the half of grand minimum intensity
    rowInd = which(apply(is.na(expr), 1, sum) > 0)
    if (length(rowInd) > 0) {
        grandMin = min(expr, na.rm = T)
        for (i in 1:length(rowInd)) {
            colInd = which(is.na(expr[rowInd[i], ]))
            rMean = mean(as.numeric(expr[rowInd[i], ]), na.rm = T)
            expr[rowInd[i], colInd] = (grandMin - 1) + runif(length(colInd)) * 1e-3 ## For numerical stability, small random number is added
        }    
    }
    
    colInd = grep("Intensity", colnames(res$fullFeatures))
    for (i in 1:length(files)) {
        ## Replace intensity values in res$fullFeatures with the normalized intensity values
        res$fullFeatures[, colInd[i]] = 2 ^ expr[, i]
    }
    
    ############################
    ## Write output to files  ##
    ############################
    ## Old fully-matched features
    fullFeatures = res$fullFeatures
    fMz = round(apply(fullFeatures[, grep("mz", colnames(fullFeatures))], 1, mean, na.rm = T), 5) # Feature m/z (averaged over runs)
    meanMz = round(apply(fullFeatures[, grep("mz", colnames(fullFeatures))], 1, mean, na.rm = T), 5) # Feature m/z (averaged over runs)
    fullFeatures = cbind("meanMz" = fMz, fullFeatures)
    fullFeatures = fullFeatures[order(fullFeatures$meanMz), ]
    outputFile = paste0(outDirectory, "/.", params$output_name, "_fully_aligned.feature") ## Hidden by adding "." (dot) to the file name
    write.table(fullFeatures, outputFile, sep = "\t", row.names = F, quote = F)
    
    ## Partially-matched features
    outputFile = paste0(outDirectory, "/", params$output_name, "_partially_aligned.feature")
    write.table(res$partialFeatures, outputFile, sep = "\t", row.names = F, quote = F)
    
    ## Unmatched features
    for (i in 1:length(res$unmatchedFeatures)) {
        outputFile = paste0(outDirectory, "/", params$output_name, "_", basename(files[i]))
        outputFile = gsub(".\\d+.feature", "_unaligned.feature", outputFile)
        write.table(res$unmatchedFeatures[[i]], outputFile, sep = "\t", row.names = F, quote = F)
    }
    ## End of the alignment and generation of necessary files
    
    ########################################################
    ## Generation of output file for single-file analysis ##
    ########################################################
} else {
    cat("  Since a single file is used, the feature alignment step is skipped\n")
    cat("  Since a single file is used, the feature alignment step is skipped\n", file = LOG)
    
    ## Old fully-matched features
    fullFeatures = features[[1]]
    meanMz = fullFeatures$mz
    filename = basename(files[1])
    filename = gsub(".\\d+.feature", "", filename)
    colnames(fullFeatures) = paste0(filename, "_", colnames(fullFeatures))
    fullFeatures = cbind(meanMz, fullFeatures)
    fullFeatures = fullFeatures[order(fullFeatures$meanMz), ]
    outputFile = paste0(outDirectory, "/.", params$output_name, "_fully_aligned.feature")
    nFeatures = dim(fullFeatures)[1]
    cat(" ", nFeatures, "features are detected\n")
    cat(" ", nFeatures, "features are detected\n", file = LOG)
    write.table(fullFeatures, outputFile, sep = "\t", row.names = F, quote = F)
}

## New fully-matched features (formatting according to the discussion on 2019/11/19)
if (length(files) > 1) {
    fullFeatures = res$fullFeatures
    fMz = round(apply(fullFeatures[, grep("mz", colnames(fullFeatures))], 1, mean, na.rm = T), 5) # Feature m/z (averaged over runs)
    refColName = sub(".\\d.feature", "", basename(files[refNo]))
    refColName = paste0(refColName, "_")
} else {
    fullFeatures = features[[1]]
    fMz = fullFeatures[, grep("mz", colnames(fullFeatures))] # Feature m/z
    refColName = ""
}
fNum = c(1:dim(fullFeatures)[1])
fIon = NULL
for (i in 1:dim(fullFeatures)[1]) {
    ## Determination of the charge state
    ##   Heuristic rules
    ##   1. Prefer the charge state other than 0
    ##   2. When multiple charge states have the same frequency over samples, choose the lower one
    charges = fullFeatures[i, grep("_z", colnames(fullFeatures))]
    charges = charges[!is.na(charges)]
    charges = setdiff(as.numeric(charges), 0)
    if (length(charges) == 0) {
        charge = 1
    } else {
        temp = table(as.numeric(charges))
        charge = as.numeric(names(temp)[temp == max(temp)])
        if (length(charge) > 1) {
            charge = charge[1]
        }
    }
    if (params$mode == -1) {
        if (charge > 1) {
            fIon[i] = paste0("[M-", charge, "H]", charge, "-")
        } else {
            fIon[i] = "[M-H]-"
        }
    } else {
        if (charge > 1) {
            fIon[i] = paste0("[M+", charge, "H]", charge, "+")
        } else {
            fIon[i] = "[M+H]+"
        }
    }
}
fRT = round(fullFeatures[, which(colnames(fullFeatures) == paste0(refColName, "RT"))] / 60, 2) # Feature RT (from reference run), unit of minute, 2 decimal
fMaxRT = fullFeatures[, which(colnames(fullFeatures) == paste0(refColName, "maxRT"))] # Feature max RT (from reference run)
fMinRT = fullFeatures[, which(colnames(fullFeatures) == paste0(refColName, "minRT"))] # Feature min RT (from reference run)
fWidth = round((fMaxRT - fMinRT) / 60, 2) # Feature width (maxRT - minRT of feature), unit of minute, 2 decimal
fSNratio = fullFeatures[, which(colnames(fullFeatures) == paste0(refColName, "SN"))] # Feature S/N ratio (from reference run)
if (length(files) > 1) {
    fIntensities = fullFeatures[, grep("Intensity", colnames(fullFeatures))]
} else {
    fIntensities = data.frame("Feature_intensity" = fullFeatures[, grep("Intensity", colnames(fullFeatures))])
}
outDf = cbind("Feature_ion" = fIon, "Feature_m/z" = fMz, "Feature_RT"= fRT, "Feature_width" = fWidth, "Feature_SNratio" = fSNratio, fIntensities)
outDf = outDf[order(outDf$`Feature_m/z`), ]
fNum = c(1:dim(outDf)[1])
outDf = cbind("Feature_num" = fNum, outDf)
outputFile = paste0(outDirectory, "/", params$output_name, "_fully_aligned.feature")
write.table(outDf, outputFile, sep = "\t", row.names = F, quote = F)

endTime = Sys.time()
elapsedTime = endTime - startTime
elapsedTime = round(as.numeric(elapsedTime, units = "mins"), 4)
cat("\n")
cat(paste0("  Feature alignment takes ", elapsedTime, " minutes\n"))
cat("\n", file = LOG)
cat(paste0("  Feature alignment takes ", elapsedTime, " minutes\n"), file = LOG)
close (LOG)
