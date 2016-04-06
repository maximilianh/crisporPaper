download:
	wget http://hgwdev.soe.ucsc.edu/~max/crispor/analysisData/bigData.tgz -O - | tar xvz 

pack:
	tar cvfz ~/public_html/crispor/bigData.tgz crisporOfftargets/*.tsv wholeGenome/specScores*.tab cropitOfftargets/ mitOfftargets/
