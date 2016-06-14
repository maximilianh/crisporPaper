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
  pointsize = 7)        # smaller font size

heatmap.2(mat_data,
  #distfun=function(c) dist(c, method="euclidian"),
  cellnote = mat_data,  # same data set for cell labels
  #hclustfun = function(x) hclust(x),
  #hclustfun = function(x) hclust(x,method = 'ward.D'),
  #hclustfun = function(x) hclust(x,method = 'ward.D2'),
  #hclustfun = function(x) hclust(x,method = 'centroid'),
  #hclustfun = function(x) hclust(x,method = 'median'),
  #hclustfun = function(x) hclust(x,method = 'average'),
  #hclustfun = function(x) hclust(x,method = 'mcquitty'),
  #hclustfun = function(x) hclust(x,method = 'complete'),
  #hclustfun = function(x) hclust(x,method = 'single'),
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(14,16),     # widens margins around plot
  keysize = 0.8,
  cexRow=1.1,
  cexCol=1.1,
  #lmat = lmat, 
  #lwid = lwid, lhei = lhei,
  #dendrogram="row",     # only draw a row dendrogram
  dendrogram='none', Rowv=FALSE, # remove the row reordering and dendrogram
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device

system("convert -density 300 -quality 100 out/heatMap.pdf out/heatMap.png")
print("output written to out/heatMap.png and .pdf")

#  ----- 2ND HEAT MAP -----

system("python heatFilter.py part2")

data <- read.delim("out/heatData.tsv", comment.char="#", check.names=FALSE)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 

# creates a 5 x 5 inch image
pdf("out/heatMap-part2.pdf",    # create PNG for the heat map        
  width = 5,        # 5 x 300 pixels
  height = 3,
  pointsize = 7)        # smaller font size

heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  dendrogram='none', Rowv=FALSE, # remove the row reordering and dendrogram
  hclustfun = function(x) hclust(x,method = 'single'),
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(14,14),     # widens margins around plot
  #keysize = 0.8,
  cexRow=0.8,
  cexCol=0.8,
  #lmat = lmat, 
  #lwid = lwid, lhei = lhei,
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device

system("convert -density 300 -quality 100 out/heatMap-part2.pdf out/heatMap-part2.png")
print("output written to out/heatMap-part2.png and .pdf")

# ----- PART 3 ---- for suppl fig 

system("python heatFilter.py part3")

data <- read.delim("out/heatData.tsv", comment.char="#", check.names=FALSE)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 

# creates a 5 x 5 inch image
pdf("out/heatMap-part3.pdf",    # create PNG for the heat map        
  width = 7,        # 5 x 300 pixels
  height = 6,
  pointsize = 7)        # smaller font size

heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  dendrogram='none', Rowv=FALSE, # remove the row reordering and dendrogram
  hclustfun = function(x) hclust(x,method = 'single'),
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(14,14),     # widens margins around plot
  #keysize = 0.8,
  cexRow=0.8,
  cexCol=0.8,
  key = FALSE,
  #lmat = lmat, 
  #lwid = lwid, lhei = lhei,
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device

system("convert -density 300 -quality 100 out/heatMap-part3.pdf out/heatMap-part3.png")
print("output written to out/heatMap-part3.png and .pdf")
