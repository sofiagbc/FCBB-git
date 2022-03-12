args <- commandArgs(trailingOnly = TRUE)
print(args)

input <- args[2]
output <- args[4]

input <- "p05pvalues.tsv"
output <- "pvalues.corrected.tsv"

data <- read.table(input, header=TRUE, sep= "\t")
data$pvalue <- sort (data$pvalue, decreasing = FALSE)

fdr=0.05
m=length(data$index)
data$qvalue <- rep(0,m)
flag <- 0
for (i in 1:m) {
  data$qvalue[i] <- data$pvalue[i]*m/i 
  if (data$qvalue[i]>fdr & flag==0) {
    index <- i
    flag <- 1
  }
    
}

print(paste("Number of discoveries at fdr<=0.05 is", index))

cat("index\tpvalue\tqvalue\n",file=output,append=FALSE)

for (i in 1:m) {
  cat("",data$index[i],"\t",round(data$pvalue[i],6),"\t",round(data$qvalue[i],6),"\n",file=output,append=TRUE)
}

# round values from beginning or only in the output file?
# check rounding, not correct
    
    

