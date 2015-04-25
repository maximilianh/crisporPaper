offtargets.tsv is a big table of crispr offtargets collected from six studies.

The data conversion is done in subdirectories of origData:
each has a convert.py file which does the conversion.
Or alternatively a file log.txt that documents how the conversion was done.
The output is always a file convert.tab which are concatenated into offtargets.tsv
column headers are in origData/headers.txt

The score column is usually the fraction of reads that fall on the off-target.

The studies are:

Cho2014
Frock2015
Hsu2013
Kim2015
Ran2015 - not included currently. Unsure where the data is.
Tsai2015
Wang2015
