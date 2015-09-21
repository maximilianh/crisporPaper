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
   Creates out/offtargetsFilt.tsv, a filtered version of offtargets.tsv
   Creates out/annotFiltOfftargets.tsv with 413 filtered off-targets.
   Prints a few basic stats to stdout.

plotGcOtCount.py: scatter plot of off-target count vs GC content. Creates FIGURE 1.

plotVenn.py: Creates out/venn.pdf
        Input are offtargets.tsv
        = SUPPL FIGURE 1: venn diagram of EMX1 VEGFA off-target overlaps.
        Also creates SUPPL FILES 2 (EMX) and 3 (VEGFA) with the details of the overlaps.

plotGuideMismatchCount.py: Creates out/specScoreMMComp.pdf 
        = SUPPL FIGURE 3. 
        Input are out/offtargetsFilt.tsv and the files in crisporOfftargets/,
        cropitOfftargets/ and mitOfftargets/.

plotMismatchFraction.py : creates a plot of
        off-target strength versus mismatch count, out/mismatchFraction-all.pdf
        = FIGURE 2 
        Input: out/offtargetsFilt.tsv and the files in crisporOfftargets/

plotRoc.py: create ROC plot = FIGURE 3
        Input: the files in crisporOfftargets, mitOfftargets and cropitOfftargets

compareOfftargetTools.py: 
        Creates a table with MIT versus CRispor versus CassOffFinder.
        Input: the files in crisporOfftargets/, casOffOfftargets/ and mitOfftargets/
        Output: out/mitCrisporSensDiff.tsv

Efficacy scoring
================

effData/:
Contains the data from the efficacy studies.
Each subdirectory contains supplementary files or converted excel tables and a conversion script
for one paper.
convert.py or log.txt that converts these to a file. the output is written to a file <author><year>.tab with always the three same columns.  (guideName, sequence, cleavageFrequency). 

WangSabatiniSvm/ - contains the R code of the Wang et al 2013 paper SVM for
efficiency scoring

To run all code, we need the library crisporEffScores.py in the directory ../crispor.
To git clone the github repo to ../crispor:
        git clone https://github.com/maximilianh/crisporWebsite ../crispor/

corrCellType.py: creates SUPPL TABLE 4
        Creates a table with inter-cell variation between efficacy experiments.
        Input: various files in effData/
        Output: out/corrCellType.tsv

effDataAddContext.py: extends 34mers or 20mers of files in the effData directory
        to the 100mers required for score calculation. For most datasets, 34mers 
        are available from an older version of this script, they were mapped with 
        BLAT from the 23mers provided by the papers. 
        The 34mers are mapped with gfClient, the 23mers with local BLAT.
        Input is effData/*.ext.tab (34mers) of effData/*.guides.tab (20mers)
        and Output is effData/*.context.tab

effDataAddScores.py: 
        runs the scoring models over effData/*.context.tab and outputs
        effData/*.scores.tab
        USED TO CREATE SUPPL FILE 2, see effData/log.txt

compEffScores.py
        Main script. Compares predicted scores from effData/*.scores.tab with
        obtained efficacy from the same files. Creates 

scoreCutoffs.py:
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


compareOfftargetTools.py: compare MIT, CRISPOR and CassOffFinder to figure out which
        off-targets MIT is missing and how they are distributed over the mismatches.

corrActivityOof.py: determine correlation between KO activity and the Bae et al OOF score for all datasets
        in the effData directory. 

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
varshney2015.tab
xu2015.tab
museumIC50.tab
museumT7.tab

