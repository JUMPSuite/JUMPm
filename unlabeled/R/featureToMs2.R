rm(list = ls())
library(mzR)

##################
## Subroutines  ##
##################
parseParams = function (paramFile) {
    lines = readLines(paramFile)
    params = NULL
    params$files = NULL
    params$samples = NULL
    for (i in 1:length(lines)) {
        lines[i] = gsub("\\s", "", lines[i])
        if (grepl("^#", lines[i]) | lines[i] == "") {
            next
        } else {
            line_i = gsub("#.*", "", lines[i])
            key = unlist(strsplit(line_i, "="))[1]
            val = unlist(strsplit(line_i, "="))[2]
        }
        if (key == "feature_file") {
            params$featureFile = as.character(val)
        } else if (key == "mzXML") {
            params$mzXMLs = as.character(val)
        } else if (key == "isolation_window") {
            params$tolIsolation = as.numeric(val)
        } else if (key == "tol_precursor") {
            params$tolPrecursor = as.numeric(val)
        } else if (key == "tol_intra_ms2_consolidation") {
            params$tolIntraMs2Consolidation = as.numeric(val)
        } else if (key == "tol_inter_ms2_consolidation") {
            params$tolInterMs2Consolidation = as.numeric(val)
        } else if (key == "mode") {
            params$mode = as.numeric(val)
        }
    }
    return (params)
}

ms2Consolidation = function(input, tol, type) {
    if (type == "intra") {
        msData = input$msData
        h = input$h
        scans = input$scans ## scans = MS2 scan numbers within a run for a feature
        scans = scans[order(h$totIonCurrent[scans], decreasing = T)]
        spec = peaks(msData, scans[1])
    } else if (type == "inter") {
        scans = names(input) ## scans = MS2 spectra from different runs for a feature
        tic = NULL
        for (i in 1:length(scans)) {
            tic = c(tic, sum(input[[scans[i]]][, 2]))
        }
        scans = scans[order(tic, decreasing = T)]
        spec = input[[scans[1]]]
    } else {
        stop("MS2 consolidation type should be either 'intra' or 'inter'")
    }
    if (length(scans) > 1) {
        for (i in 2:length(scans)) {
            if (type == "intra") {
                p = peaks(msData, scans[i])
            } else if (type == "inter") {
                p = input[[scans[i]]]
            } else {
                stop("MS2 consolidation type should be either 'intra' or 'inter'")
            }
            for (j in 1:nrow(spec)) {
                if (is.null(p)) {
                    break
                }
                mz = spec[j, 1]
                intensity = spec[j, 2]
                lL = mz - mz * tol / 1e6
                uL = mz + mz * tol / 1e6
                if (class(p) == "numeric") {
                    ind = which(p[1] >= lL & p[1] <= uL)
                    if (length(ind) == 0) {
                        next
                    } else {
                        ## Note that m/z (in "spec") is not changed
                        intensity = intensity + p[2]
                        spec[j, 1] = mz
                        spec[j, 2] = intensity
                        p = NULL
                    }
                } else {
                    ind = which(p[, 1] >= lL & p[, 1] <= uL)
                    if (length(ind) == 0) {
                        next
                    } else {
                        ind = ind[order(p[ind, 2], decreasing = T)[1]] ## Choose the strongest peak
                        intensity = intensity + p[ind, 2]
                        spec[j, 1] = mz
                        spec[j, 2] = intensity
                        p = p[-ind, ]
                    }
                }
            }
            spec = rbind(spec, p)
        }
    }
    return (spec)
}

##########
## Main ##
##########

## Parameters for the identification and consolidation of MS2 spectra
##  tolIsolation = 1.25 ## Half of isolation window (instrument setting)
##  tolPrecursor = 10 ## PPM tolerance for finding MS1 peaks corresponding to a feature
##  tolIntraMs2Consolidation = 10 ## PPM tolerance for merging MS2 spectra within a run
##  tolInterMs2Consolidation = 20 ## PPM tolerance for merging MS2 spectra between runs
##  ppiThreshold = "max" (hard-coded) PPI (Percentage of precursor ion for selecting a feature responsible for a MS2 scan)
##  pctTfThreshold = 50 (hard-coded) Threshold of percentage of TF (true feature)

startTime = Sys.time()
args = commandArgs(T)
paramFile = args[1]
featureFile = args[2]
mzXMLs = unlist(strsplit(args[3], ","))
outDirectory = args[4]
logFile = args[5]
LOG = file(logFile, "a")

## For testing in a desktop
# args = commandArgs(T)
# paramFile = "../IROAsamples/jumpm_negative.params"
# featureFile = "../IROAsamples/.IROA_IS_NEG_fully_aligned.feature"
# mzXMLs = "../IROAsamples/IROA_IS_NEG_1.mzXML,../IROAsamples/IROA_IS_NEG_2.mzXML,../IROAsamples/IROA_IS_NEG_3.mzXML"
# mzXMLs = unlist(strsplit(mzXMLs, ","))
# outDirectory = "../IROAsamples/MS2"
# logFile = "tmplog"
# LOG = file(logFile, "a")

params = parseParams(paramFile)
cat("\n  Identification and consolidation of MS2 spectra for the features\n")
cat("  ================================================================\n\n")
cat("\n  Identification and consolidation of MS2 spectra for the features\n", file = LOG)
cat("  ================================================================\n\n", file = LOG)

## Parameters
ppiThreshold = "max" ## Hard-coded
pctTfThreshold = 50 ## Hard-coded
tolIsolation = params$tolIsolation
tolPrecursor = params$tolPrecursor
tolIntraMs2Consolidation = params$tolIntraMs2Consolidation
tolInterMs2Consolidation = params$tolInterMs2Consolidation

## Read a feature file
features = read.table(featureFile, header = T, sep = "\t", row.names = NULL,
                      stringsAsFactors = F, comment.char = "", check.names = F)
colnames(features) = gsub(" ", "", colnames(features))
colnames(features) = gsub("/", "", colnames(features))

featureToSpectra = vector("list", nrow(features))
featureToScanNum = data.frame(matrix(NA, nrow = nrow(features), ncol = length(mzXMLs)))
colnames(featureToScanNum) = c(basename(mzXMLs))

for (m in 1:length(mzXMLs)) {
    cat(paste0("  Reading ", basename(mzXMLs[m]), " (it may take a while)\n"))
    cat(paste0("  Reading ", basename(mzXMLs[m]), " (it may take a while)\n"), file = LOG)
    msData = openMSfile(mzXMLs[m])
    h = header(msData)
    fileName = sub(".mzXML", "", basename(mzXMLs[m]))
    df = features[, grep(fileName, colnames(features))]
    colnames(df) = gsub(paste0(fileName, "_"), "", colnames(df))
    
    ##############################################################
    ## Identification of MS2 scans responsible for each feature ##
    ##############################################################
    pb = txtProgressBar(1, nrow(h), style = 3)
    cat("    Identifying MS2 scans responsible for the features\n")
    cat("    Identifying MS2 scans responsible for the features\n", file = LOG)
    for (i in 1:nrow(h)) {
        setTxtProgressBar(pb, i)
        if (h$msLevel[i] == 1) {
            ms1ScanNum = i
        } else if (h$msLevel[i] == 2) {
            ind = which(i > df$`minMS1Scan#` & i < df$`maxMS1Scan#` &
                            df$PercentageofTF < pctTfThreshold &
                            df$mz >= (h$precursorMZ[i] - tolIsolation) &
                            df$mz <= (h$precursorMZ[i] + tolIsolation))
            if (length(ind) == 0) {
                next
            }
            ## Check the intensities of candidate features at the very preceding MS1 scan
            ## For example, let's assume that candidate features are as follows for a MS2 scan #140
            ## index  mz        z MS1scan#  RT   minRT maxRT Intensity
            ## 211    218.1498  0 136       136  1     951   37544
            ## 169    220.0705  0 126       126  1     1446  91709
            ## 113    218.8597  6 18        18   1     764   91745
            ## 3      220.1052  0 1         1    1     1248  355843
            ## Also, suppose that the very preceding MS1 scan is at scan#140
            ## Then, we need to check the intensities of those candidate features at scan#140
            ## example) For the 1st feature whose representative m/z = 218.1498,
            ##          1. Get MS1 spectrum of scan#140; p = peaks(msData, 140)
            ##          2. Open a m/z window with a tolerance; [218.1498 - 10ppm, 218.1498 + 10ppm]
            ##          3. Look for the MS1 peak with the highest intensity within the window, and record the intensity
            ##          4. If the intensity is the higest among candidate features, choose it for the MS2 scan
            ##          5. Otherwise, check the next candidate feature
            ##          6. Instead of choosing one feature with the highest intensity,
            ##             PPI can be used and multiple features may be chosen for each MS2
            p = peaks(msData, ms1ScanNum)
            ppi = rep(0, length(ind))
            for (j in 1:length(ind)) {
                mz = df$mz[ind[j]]
                lL = mz - mz * tolPrecursor / 1e6
                uL = mz + mz * tolPrecursor / 1e6
                ppi[j] = max(0, p[(p[, 1] >= lL & p[, 1] <= uL), 2])
            }
            ## When there's no reasonable candidate feature for the MS2 scan
            if (sum(ppi) == 0) {
                next
            }
            ppi = ppi / sum(ppi) * 100
            ## Selection of feature(s)
            if (ppiThreshold == "max") {
                ind = ind[which.max(ppi)]
            } else {
                ind = ind[ppi > ppiThreshold]
            }
            ## When there's no reasonable candidate feature of which PPI is greater than a threshold
            if (length(ind) == 0) {
                next
            }
            for (j in 1:length(ind)) {
                rowInd = ind[j]
                if (is.na(featureToScanNum[rowInd, m])) {
                    featureToScanNum[rowInd, m] = i
                } else {
                    featureToScanNum[rowInd, m] = paste(c(featureToScanNum[rowInd, m], i), collapse = ";")
                }
            }
        } else {
            next ## skip if msLevel is neither 1 nor 2
        }
    }
    cat("\n")
    
    ################################################
    ## Consolidation of MS2 spectra within a run  ##
    ################################################
    pb = txtProgressBar(1, dim(featureToScanNum)[1], style = 3)
    cat("    Merging MS2 spectra within each feature\n")
    cat("    Merging MS2 spectra within each feature\n", file = LOG)
    for (i in 1:dim(featureToScanNum)[1]) {
        setTxtProgressBar(pb, i)
        if (is.na(featureToScanNum[i, m])) {
            next
        } else {
            input = list(msData = msData, h = h, scans = as.numeric(unlist(strsplit(featureToScanNum[i, m], ";"))))
            spec = ms2Consolidation(input, tolIntraMs2Consolidation, "intra")
        }
        spec = spec[order(spec[, 1]), ]
        featureToSpectra[[i]][[fileName]] = spec
    }
    cat("\n")
    cat("\n", file = LOG)
}

################################################
## Consolidation of MS2 spectra between runs  ##
################################################
cat("    Merging MS2 spectra between runs for each feature\n")
cat("    Merging MS2 spectra between runs for each feature\n", file = LOG)
pb = txtProgressBar(1, length(featureToSpectra), style = 3)
dir.create(outDirectory, showWarnings = F)
for (i in 1:length(featureToSpectra)) {
    setTxtProgressBar(pb, i)
    if (is.null(featureToSpectra[[i]])) {
        next
    } else {
        ## Determination of the charge state
        ##   Heuristic rules
        ##   1. Prefer the charge state other than 0
        ##   2. When multiple charge states have the same frequency over samples, choose the lower one
        charges = features[i, grep("_z", colnames(features))]
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
        
        ## Generation of .dta files
        ## Limit the number of peaks in each .dta file to 100
        ## Divide m/z-range into 10 bins (e.g. 0~100, 100~200, etc.) and retain the 10 largest peaks in each bin
        spec = ms2Consolidation(featureToSpectra[[i]], tolInterMs2Consolidation, "inter")
        if (dim(spec)[1] > 100) {
            nBins = 10
            bins = seq(min(spec[, 1]), max(spec[, 1]), length.out = (nBins + 1))
            filteredSpec = NULL
            for (j in 1:nBins) {
                rowInd = which((spec[, 1] >= bins[j] & spec[, 1] < bins[j + 1]))
                if (length(rowInd) > 0) {
                    peaksInBin = matrix(spec[rowInd, ], nrow = length(rowInd), ncol = 2)
                    if (dim(peaksInBin)[1] > 1) {
                        peaksInBin = peaksInBin[order(peaksInBin[, 2], decreasing = T), ]
                    }
                    if (dim(peaksInBin)[1] > 10) {
                        peaksInBin = peaksInBin[1:10, ]
                    }
                    filteredSpec = rbind(filteredSpec, peaksInBin)
                }
            }
            spec = filteredSpec
        }
        ## In .MS2 file, the header includes mass and charge information
        ## Mass information should be (M+H)+ for positive mode or (M-H)- for negative mode, where M is a neutral mass
        ## Charge information is the same as one in .feature file
        Hmass = 1.007276466812;
        spec = spec[order(spec[, 1]), ]
        mz = features[i, 1]
        if (params$mode == 1) {
            neutralMass = (mz - Hmass) * charge
            mh = neutralMass + Hmass
        } else if (params$mode == -1) {
            neutralMass = (mz + Hmass) * charge
            mh = neutralMass - Hmass
        }
        spec = rbind(c(mh, charge), spec)
        filename = paste0(outDirectory, "/f", i, ".MS2")
        write.table(spec, file = filename, col.names = F, row.names = F, sep = "\t", quote = F)
    }
}

## Write featureToScanNum to a file
featureToScanNum = cbind(rownames(featureToScanNum), featureToScanNum)
rownames(featureToScanNum) = NULL
colnames(featureToScanNum)[1] = "featureNumber"
filename = paste0(outDirectory, "/featureToScanNum.txt")
write.table(featureToScanNum, file = filename, row.names = F, sep = "\t", quote = F)

endTime = Sys.time()
elapsedTime = endTime - startTime
elapsedTime = round(as.numeric(elapsedTime, units = "mins"), 4)
cat("\n")
cat(paste0("  Processing of MS2 spectra takes ", elapsedTime, " minutes\n"))
cat("\n", file = LOG)
cat(paste0("  Processing of MS2 spectra takes ", elapsedTime, " minutes\n"), file = LOG)
close (LOG)
