download:
	wget http://hgwdev.soe.ucsc.edu/~max/crispor/analysisData/offtargets2.tgz -O - | tar xvz 

pack:
	tar cvfz ~/public_html/crispor/bigData2.tgz crisporOfftargets/*.tsv wholeGenome/specScores*.tab cropitOfftargets/ mitOfftargets/
