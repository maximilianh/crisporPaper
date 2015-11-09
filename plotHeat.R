# large parts part copied from http://sebastianraschka.com/Articles/heatmaps_in_r.html
require("gplots")
require("RColorBrewer")

system("python heatFilter.py")

data <- read.delim("out/heatData.tsv", comment.char="#", check.names=FALSE)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 

# creates a 5 x 5 inch image
pdf("out/heatMap.pdf",    # create PNG for the heat map        
  width = 5,        # 5 x 300 pixels
  height = 5,
  #res = 300,            # 300 pixels per inch
  pointsize = 7)        # smaller font size

#lmat = rbind(c(0,3),c(2,1),c(0,4))
#lwid = c(1.5,4)
#lhei = c(1.5,4,1)

heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  #hclustfun = function(x) hclust(x),
  #hclustfun = function(x) hclust(x,method = 'ward.D'),
  #hclustfun = function(x) hclust(x,method = 'ward.D2'),
  #hclustfun = function(x) hclust(x,method = 'centroid'),
  #hclustfun = function(x) hclust(x,method = 'median'),
  hclustfun = function(x) hclust(x,method = 'average'),
  #hclustfun = function(x) hclust(x,method = 'single'),
  #hclustfun = function(x) hclust(x,method = 'complete'),
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(14,16),     # widens margins around plot
  keysize = 0.8,
  cexRow=1.1,
  cexCol=1.1,
  #lmat = lmat, 
  #lwid = lwid, lhei = lhei,
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device

system("convert -density 300 -quality 100 out/heatMap.pdf out/heatMap.png")
print("output written to out/heatMap.png and .pdf")
