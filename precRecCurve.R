library(ROCR)
args = commandArgs(trailingOnly = TRUE)
print(args[1])
df = read.delim(args[1])
pred = prediction(df$score,df$label)
# precision/recall plot
perf = performance(pred, "prec", "rec")
coords = as.matrix(data.frame(fpr=perf@x.values, fnr=perf@y.values))
write.table(coords, file=args[2], row.names=FALSE, col.names=FALSE, sep=" ")

# ROC plot
perf = performance(pred, "tpr", "fpr")
coords = as.matrix(data.frame(fpr=perf@x.values, fnr=perf@y.values))
write.table(coords, file=args[3], row.names=FALSE, col.names=FALSE, sep=" ")
