rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
# args[1] = "../alignment_residual.txt"
# args[2] = "../measured_rt_error.txt"
modX = read.table(args[1], sep = "\t", row.names = NULL,
                  stringsAsFactors = F, comment.char = "", check.names = F)
x = read.table(args[2], sep = "\t", row.names = NULL,
               stringsAsFactors = F, comment.char = "", check.names = F)
f = ecdf(as.numeric(unlist(modX)))
p = f(as.numeric(x[, 3])) ## similar to the concept of p-value (when error is small, y becomes small)
x[, 3] = p
write.table(file = "calculated_rt_pvalue.txt", x, sep = "\t", row.names = F, col.names = F, quote = F)