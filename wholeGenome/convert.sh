# used table S5 from Chari 2015 to get the 20mers
cut -f2 nmeth.3473-S5.txt | grep -v ^$ | grep -v Seq > allGuides.txt
