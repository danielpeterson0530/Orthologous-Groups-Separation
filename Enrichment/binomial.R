#/usr/bin/env Rscript
#
# R script to conduct exact binomial comparison
# ** Can be run independant, however can be used for multiple files using binomial_run.py

args = commandArgs(trailingOnly=TRUE)
if (length(args)<6) {
  stop("improper usage: Rscript tool.R ./filename portal colint_xobs colint_yobs x_total y_total", call.=FALSE)
}

filename = args[1]
portal = args[2]
xint = as.integer(args[3])
yint = as.integer(args[4])
x_total = as.integer(args[5])
y_total = as.integer(args[6])
data <- read.table(filename, header=FALSE, sep = "\t")

term <- vector()
x_obs <- vector()
y_obs <- vector()
pval <- vector()
x_totals <- vector()
y_totals <- vector()

for (i in 1:nrow(data)) {
  term_tmp <- paste(data[i,2])
  x_tmp <- as.integer(data[i,xint])
  y_tmp <- as.integer(data[i,yint])

  yprop_tmp <- (y_tmp / y_total)

  result_tmp <- binom.test(x=x_tmp, n=x_total, p=yprop_tmp, alternative = "greater", conf.level = 0.95)
  pval_tmp <- result_tmp[["p.value"]]

  term <- c(term, term_tmp)
  x_obs <- c(x_obs, x_tmp)
  y_obs <- c(y_obs, y_tmp)
  x_totals <- c(x_totals, x_total)
  y_totals <- c(y_totals, y_total)
  pval <- c(pval, pval_tmp)
}

all_results <- data.frame(term, x_obs, x_totals, y_obs, y_totals, pval)
all_results$adjpval <- p.adjust(all_results$pval,method="BH")

for (i in 1:nrow(all_results)) {
#   if (all_results[i,"adjpval"] < 0.05 & all_results[i,"y_obs"] > 0) {
    if (all_results[i,"y_obs"] > 0 & all_results[i,"x_obs"] > 0) {
      cat( paste( toString(portal),
                  toString(gsub(" ", "_", all_results[i,"term"])),
                  toString(all_results[i,"x_obs"]),
                  toString(all_results[i,"x_totals"]),
                  toString(all_results[i,"y_obs"]),
                  toString(all_results[i,"y_totals"]),
                  toString(all_results[i,"pval"]),
                  toString(all_results[i,"adjpval"]),
                  "\n"), sep="\t")
   }
}



#outtable = all_results[order(all_results$acc_pval_adjust),]
#printtable= outtable[ which(outtable$acc_pval_adjust < 0.05), ]

#printtable= outtable[ which(outtable$pcore_pval_adjust < 0.05), ]
#printtable= outtable[ which(outtable$npcore_pval_adjust < 0.05), ]

#print(printtable, quote=FALSE, row.names=FALSE)

#outfile = paste(filename, "_goenrichment_results.txt", sep="")
#outtable = all_results[order(all_results$pval_adjust),]
#printtable= outtable[ which(outtable$pval_adjust < 0.05), ]
#write.table(printtable, outfile, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)
