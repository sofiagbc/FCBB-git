library('getopt')

args <- commandArgs(TRUE)
#print(args)

spec = matrix(c("input", "i", 1, "character","output", "o", 1, "character"), byrow=TRUE, ncol=4)
opt=getopt(spec)
#print(opt)

# get flags here with getopt

input <- opt$input
output <- opt$output

data_input <- read.table(input, header=TRUE, sep= "\t")
m=length(data_input$index)

data_input$pvalue <- as.numeric(as.character(data_input$pvalue))
data_sorted <- order(data_input$pvalue)


data <- data_input
data$index <- rep(0,m)
data$pvalue <- rep(0,m)
data$qvalue <- rep(0,m)
for (i in 1:m) {
  num=data_sorted[i]
  data$index[i]=num
  data$pvalue[i]=data_input$pvalue[num]
}

fdr=0.05
data$qvalue <- rep(0,m)
flag <- 0
for (i in 1:m) {
  data$qvalue[i] <- data$pvalue[i]*m/i 
  if (data$qvalue[i]>fdr & flag==0) {
    index <- i-1
    flag <- 1
  }
    
}

print(paste("Number of discoveries at fdr<=0.05 is", index))

data$pvalue <- round(data$pvalue,6)
data$qvalue <- round(data$qvalue,6)
write.table(data,output,quote=FALSE, sep="\t",row.names=FALSE)



