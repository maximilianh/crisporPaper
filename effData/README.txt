The data in museum*.tab was manually curated by JP Concordet from experiments
at the Paris Museum. Concordet2 is another set from him.
schoenigHs.tab and schoenigMm.tab are from Kai Schoenig, University of Cologne.
eschstruth.tab is from Alexis Eschstruth, by email.

every .tab file is created by a directory of the same name with scripts or
commands to convert the original data to this tab file.
The commands to create the .tab are either called log.txt or convert.py in this
directory.  The .tab files are only the first step of the data import.

The .ext.tab files include the genomic position of the sequences and were
created from the .tab files with BLAT or from the original data in the
supplemental files. The .ext.tab files include manual selection of the best
match for a 23mer and should be considered the input data.

The .context.tab files include the +- 50bp context around the PAM, based on the
.ext.tab files.  They were created by the script ../effDataAddContext.py

The .scores.tab files include all the former fields + all efficiency scores.
They were created with the script ../effDataAddScores.py

CHANGES:

One sequence in eschstruth.context was manually changed from 
CCTCAGTTGACACTTTTGAGCGG = the sequence in the reference genome to 
CCTCAGTTGACAGGTTTCAGCGG = the sequence in this strain 
So if you regenerate this .context.tab file from guides.tab, 
you would have to introduce this mutation again.



