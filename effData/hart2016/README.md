Received TKOFoldChange.tgz from Trevor Hart by email

convert.py creates ../hart* files. 
python effDataAddContext / effDataAddScores was run to map these to the genome and score them.
the resulting ../hart*.scores.tab files were then moved to ./mapped/

python filterScores.py was run to filter these files down to only the essential genes, the output are
../hart*.scores.tab files

