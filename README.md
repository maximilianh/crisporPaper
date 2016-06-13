Required software
=================

* python2.7 with pip 
  * brew install python
  * apt-get install python2.7
* matplotlib and scipy 
  * pip install matplotlib scipy
* R 
  * brew tap homebrew/science; brew install r
* gplots, colorbrewer, ROCR and e1071
  * Rscript -e 'install.packages(c("gplots", "RColorBrewer","ROCR","e1071"),  repos="http://cran.rstudio.com/")'
* limma
  * Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("limma")'

The "crisporWebsite" github repo has to be located in ../crisporWebsite. This because the scripts in here need the 
crispor python libraries for the scoring.

Offtarget-prediction
====================

offtargets.tsv is a big table of crispr offtargets collected from six studies.

The data conversion is done in subdirectories of origData:
each has a convert.py file which does the conversion.
Or alternatively a file log.txt that documents how the conversion was done.
The output is always a file convert.tab which were concatenated into offtargets.tsv
with this command:
    cat origData/*/convert.tab | sort -k1,1 -k3,3gr | cat origData/headers.txt - > offtargets.tsv

The score column is the fraction of reads that fall on the off-target, if reads were determined, otherwise the number of lentiviral insertion sites relative to all sites for the off-target.

The off-target studies are:

Cho2014
Frock2015
Hsu2013
Kim2015
Ran2015
Tsai2015
Wang2015

annotateOffs.py is the central library with helper functions for all scripts.

scripts
-------

Each script creates a tab-sep file and/or a plot in PNG and PDF format.

The scripts have to be run in this order:

filtAnnotateOfftargets.py: creates SUPPL FILE 1 out/annotOfftargets.tsv from offtargets.tsv
   Creates out/annotOfftargets.tsv, a version of offtargets.tsv with more 
   annotations, like spec. score, mismatch count, etc

   Creates out/offtargetsFilt.tsv, a filtered version of offtargets.tsv
   offtargetsFilt.tsv has all the two GC-rich guides removed and all 
   off-targets < 0.01% also removed.
   This file does contain guides that have no off-targets > 0.01%.

   The script also creates out/annotFiltOfftargets.tsv with 413 filtered off-targets.
   The file has guides without any off-targets removed.

  697 offtargets, 36 guides
  wrote out/annotOfftargets.tsv, 697 rows
  removed these guides: Tsai_HEK293_sgRNA4,Tsai_VEGFA_site2
  kept 34 guide names (->tests of identical guides in diff. cells)
  kept 28 different guide sequences
  kept 447 off-targets
  output written to out/offtargetsFilt.tsv
  413 offtargets, 34 guides
  wrote out/annotFiltOfftargets.tsv, 413 rows

plotGcOtRatio.py: scatter plot of on/off-target ratio vs GC content. 
        Creates SUPPL FIGURE 1.
        Output: 

plotGcSpecScore.py: scatter plot GC content vs specificity score.
        Output: out/gcSpecScores.pdf

plotVenn.py: Creates out/venn.pdf
        Input are offtargets.tsv
        = SUPPL FIGURE S2: venn diagram of EMX1 VEGFA off-target overlaps.
        Also creates SUPPL FILES 2 (EMX) and 3 (VEGFA) with the details of the overlaps.
        (tabCat commands in log.txt)

plotSpecScoreChange.py: Creates out/specScoreMMComp.pdf 
        = SUPPL FIGURE S3. 
        Input are out/annotFiltOfftargets.tsv and the files in crisporOfftargets/,
        cropitOfftargets/ and mitOfftargets/.

plotMismatchFraction.py : creates a plot off-target strength versus mismatch count
        = FIGURE 1 
        Input: out/offtargetsFilt.tsv and the files in crisporOfftargets/
        Output: out/out/mismatchFraction-all.pdf
        With the argument "supp", shows only the two Tsai outliers in the plot.

plotRoc.py: create ROC plot 
        = FIGURE 2
        Input: the files in crisporOfftargets, mitOfftargets and cropitOfftargets
        Output: out/roc.pdf
        With the argument "supp", it adds a dataset "with two outliers" to the ROC plot
        (question by referee)

compareOfftargetTools.py: 
        Creates a table with MIT versus CRispor versus CassOffFinder.
        = Supplemental Table 2
        Input: the files in crisporOfftargets/, casOffOfftargets/ and mitOfftargets/
        Output: out/mitCrisporSensDiff.tsv

compSpecScoreVsOtCount_split.py:
        Creates a scatterplot that shows specificity score histogram, off-target frequency
        and off-target count, 
        = SUPPL FIG S4 and FIGURE 3
        Output: out/specScoreVsOtCount-CRISPOR.pdf and out/specScoreVsOtCount-MIT.pdf
        The version of the script without _split.py is an older version that
        overlays both diagrams into a single figure which didn't pass 
        the coauthors/colleagues.

Efficiency scoring
================

This part requires 3rd party software packages:
- an R package for the Wang/Lander score, to install it:
    Rscript -e 'install.packages("e1071", repos="http://cran.rstudio.com/")'
- usually not needed: the python svmlight library for the fast calculation of Chari scores:
    pip install svmlight

effData/:
Contains the data from the efficiency studies.
Each subdirectory contains supplementary files or converted excel tables and a conversion script
for one paper.
convert.py or log.txt that converts these to a file. the output is written to a file <author><year>.tab with always the three same columns.  (guideName, sequence, cleavageFrequency). 

WangSabatiniSvm/ - contains the R code of the Wang et al 2013 paper SVM for
efficiency scoring

To run all code, we need the library crisporEffScores.py in the directory ../crispor.
To git clone the github repo to ../crispor:
        git clone https://github.com/maximilianh/crisporWebsite ../crispor/

corrCellType.py: creates SUPPL TABLE 4
        Creates a table with inter-cell variation between efficiency experiments.
        Input: various files in effData/
        Output: out/corrCellType.tsv

effDataAddContext.py: extends 34mers or 20mers of files in the effData directory
        to the 100mers required for score calculation. For most datasets, 34mers 
        are available from an older version of this script, they were mapped with 
        BLAT from the 23mers provided by the papers. 
        The 34mers are mapped with gfClient, the 23mers with local BLAT.
        Input is effData/*.ext.tab (34mers) or effData/*.guides.tab (20mers)
        and Output is effData/*.context.tab
        All newer datasets are in effData/*.guides.tab files.

effDataAddScores.py: 
        runs the scoring models over effData/*.context.tab and outputs
        effData/*.scores.tab
        USED TO CREATE SUPPL FILE 2, see effData/log.txt

compEffScores.py
        Main script. Compares predicted scores from effData/*.scores.tab with
        obtained efficiency from the same files. 
        Input: effData/*.scores.tab
        Output: out/effScoreComp.tsv and out/compEffScores-train.pdf and
        out/compEffScores-valid.pdf

plotHeat.R: 
        Creates the heatmap of correlations 
        Input: out/effScoreComp.tsv
        Output: out/heatData.tsv and out/heatMap.pdf and .png
               out/heatMap-2.pdf and .png
        = FIGURE 4 and FIGURE 5
        This is an R script!
        Has to be run with "Rscript plotHeat.R" unlike all other scripts.

scoreCutoffs.py
        obtain the 75 percentiles from the scores and print to stdout.
        Has been manually copied already into binClass.py

binClass.py 
        discretizes the predictions scores and outcomes to generate binary
        classification scenario: the top 25% of the outcomes are considered
        "positives" and the top 25% percent of the scores considered
        positive predictions. 
        Input: effData/*.scores
        Output: out/binClassMetrics.tsv, with precision, recall and f1 
        values for every combination of score and dataset.

plotPrecRecall.py:
        plot precision recall for SUPPLEMENTAL FIGURE 5
        Input: out/binClassMetrics.tsv
        Output: 

monteCarloAlena.py
monteCarloSchoenig.py
        Calculate P-Values to obtain a certain result from the Schoenig 
        or Alena (="Shkumatava") datasets. 
        Input: effData/alenaAll.scores.tab or effData/schoenig.scores.tab
        Output: To stdout

compRelativeScores.py: pick out pairs of guides on the same gene and compare efficiency scores for them
simulate.py: pick 2,3,4 guides and determine how if one gave a better result than X percent in the experiment.
        repeat 1000 times and give percentage. Also try using only guides with higher than 0.1-1.0 efficiency
        scores.

compareTools.py: for each guide, get recall, sens, spec for CRISPOR and MIT
compareMitCrisporSpecScore.py: compare the MIT with the CRISPOR specificity
        score. CRISPOR uses the same formula but a more sensitive aligner, so
        the scores aren't always the same.
calcSvmScores.py: calculate all possible SVM scores and write them to a cache file, so all other scripts
        don't have to startup R again. SVM scores are relatively slow to calculate.
calcChariScores.py: same, for Chari et al SVM scores.

plotBigRocPrecRecall.py: create ROC and precision-recall plots for all big datasets (based on cell cultures)
        -> SUPPLEMENTAL FIGURE 7

compareOfftargetTools.py: compare MIT, CRISPOR and CassOffFinder to figure out which
        off-targets MIT is missing and how they are distributed over the mismatches.

corrActivityOof.py: determine correlation between KO activity and the Bae et al OOF score for all datasets
        in the effData directory. 

compHartParameters.py: obtain correlations of all Hart 2016 cell-culture+timepoint against all scores

The efficiency studies are:

xu2015Train.tab - The Wang/Sabatini/Lander Science 2014 data, as used to train the Xu SVM model.
The data was also used to train the Wang et al SVM model, which they provided to us

chari2015Train.tab
chari2015Valid + seven different cell lines
doench2014-Hs.tab
doench2014-Mm.tab

farboud2015.tab
gagnon2014.tab

ren2015.tab
liu2016.tab
varshney2015.tab

xu2015.tab

museumIC50.tab
museumT7.tab

Add another efficiency score file
=================================

* create a file `effData/dataset_genome.guides.tab` with the columns guide, seq, modFreq
* make sure you have /gbdb/genome/genome.2bit and /gbdb/genome/genome.sizes (fetchChromSizes from UCSC if needed or via rsync from hgdownload.soe.ucsc.edu)
* run python effDataAddContext.py, will create `effData/dataset_genome.context.tab`
* run python effDataAddScores.py, will create `effData/dataset_genome.scores.tab`
* run python compEffScores.py, will create out/compEffScores-*.pdf and .png and out/effScoreComp.tsv
* run Rscript plotHeat.R, will create out/heatMap*.pdf and .png
