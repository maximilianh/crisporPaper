require("cocor")

dataFrame = read.delim("out/signFig5.tab")
for(i in 1:nrow(dataFrame)) {
        row = dataFrame[i, ]
        cocor.dep.groups.overlap(row$, r.jh=0.299678784731, r.kh=0.229522146089, n=163, alternative="less", alpha=0.05, conf.level=0.95, null.value=0)
}
