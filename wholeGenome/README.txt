scripts to get all scores for all sp cas9 guides in the hg19 genome

  # get sequences of extended exons
  hgsql hg19 -NB -e 'select * from refGene' | cut -f2- > refGene.gp
  genePredToBed refGene.gp stdout | grep -v hap | grep -v chrUn | bedToExons stdin stdout -cdsOnly | awk '{$2=$2-17; $3=$3+17; $6="+"; print}' > refGene.cdsExons.bed
  twoBitToFa -bed=refGene.cdsExons.bed /gbdb/hg19/hg19.2bit extExonSeqs.fa -bedPos &
  # search for the guides
  python findGuides.py extExonSeqs.fa guides.bed > guides.fa
  cat guides.fa | grep -v \> > allGuides.txt

  # filter the guides and split for parasol
  python splitGuides.py
  # create the joblist
  python makeJoblist.py
  
  mv jobFiles jobList ../../crispor/
  cd  ../../crispor/
  mkdir out/guides out/offs -p
  ssh ku
  para create jobList
  para push

  tabCat out/guides/ | grep ^[0-9]*fw | cut -f2,3 | grep -v ^specScore > specScores.tab
  # concat all spec scores into one file
  cat specScores.*.tab > specScores.tab

