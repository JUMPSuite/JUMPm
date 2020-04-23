limmaTest = function(df, comparison) {
  ## Retrieve comparison group information
  samplesInGroups = unlist(strsplit(comparison, ":"))
  nGroups = length(samplesInGroups)
  groups = list()
  nSamples = 0
  usedSamples = NULL
  for (g in 1:nGroups) {
    groups[[g]] = unlist(strsplit(samplesInGroups[g], ","))
    nSamples = nSamples + length(groups[[g]])
    usedSamples = c(usedSamples, groups[[g]])
  }
  
  ## Generate a design matrix (which contains the information of comparison)
  colnames(df) = gsub("_Intensity", "", colnames(df))
  colNames = colnames(df)[which(colnames(df) %in% usedSamples)]
  design = matrix(0, nrow = nSamples, ncol = nGroups)
  for (g in 1:nGroups) {
    design[which(colNames %in% groups[[g]]), g] = 1
  }
  colnames(design) = paste("group", seq(1, nGroups), sep = "")
  
  ## Generate a contrast matrix and new column names for the LIMMA result table
  contVec = NULL
  combMatrix = combn(seq(1, nGroups), 2)
  for (j in 1:ncol(combMatrix)) {
    contVec = c(contVec, paste(paste("group", combMatrix[1, j], sep = ""), paste("group", combMatrix[2, j], sep = ""), sep = "-"))
  }
  contMatrix = makeContrasts(contrasts = contVec, levels = design)
  
  ## Perform LIMMA
  data = df[, which(colnames(df) %in% colNames)]
  fit = lmFit(data, design)
  fit = contrasts.fit(fit, contMatrix)
  fit = eBayes(fit)
  res = topTable(fit, n = nrow(df), adjust = "BH", sort = "none")
  return (res)
}
