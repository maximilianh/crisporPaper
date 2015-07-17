scripts to get all scores for all sp cas9 guides in the hg19 genome
Based on the guide file provided by Chari 2015


  # source log.txt - not needed anymore. .gz file is provided now
  python splitGuides.py
  python makeJoblist.py
  
  mv jobFiles jobList ../../crispor/
  cd  ../../crispor/
  mkdir out/guides out/offs -p
  ssh ku
  para create jobList
  para push
