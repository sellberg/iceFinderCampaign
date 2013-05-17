#!/usr/bin/awk -f

BEGIN {
  default=-1
  header=1    #print header
  run="NORUN?"#run number
  nxtc=-1     #Number of XTC files
  sizextc=-1  #Total size of XTC files
  darksub=-1  #Dark Cal Subtraction
  nfr=-1      #Number of Frames Processed
  nhitf=-1    #number of hits
  hitf=-1     #hitfinder
  hitfADC=-1  #hitfinderADC
  hitfNAT=-1  #hitfinderNAT
  nicef=-1    #icefinder
  icefADC=-1  #icefinderADC
  icefNAT=-1  #icefinderNAT
  hitfr=-1    #hitrate
  icefr=-1    #hitrate ice
  type0=-1    #sorted water hits
  type0_below50=-1    #sorted water hits below 50 ADUs
  type0_below100=-1    #sorted water hits below 100 ADUs but above 50 ADUs
  type1=-1    #sorted ice hits
  type1_below50=-1    #sorted ice hits below 50 ADUs
  type1_below100=-1    #sorted ice hits below 100 ADUs but above 50 ADUs
  type2=-1    #sorted weak hits
  type3=-1    #sorted hits out-of-focus

  if (length(comment) == 0){comment "no comment"}
} 
{
  if ($0~/r0/)            {                         run=$1     }
  if ($0~/Number of XTC/) {                         nxtc=$5    }
  if ($0~/XTC total size/){                         sizextc=$5 }
  if ($0~/useDarkcalSub/) {gsub(/useDarkcal.*=/,"");darksub=$0 }
  if ($0~/Frames pro/)    {                         nfr=$3     }
  if ($0~/hitfinder=/)    {gsub(/hitfinder=/,"")   ;hitf=$0    }
  if ($0~/hitfinderADC=/) {gsub(/hitfinderADC=/,"");hitfADC=$0 }
  if ($0~/hitfinderNAT=/) {gsub(/hitfinderNAT=/,"");hitfNAT=$0 }
  if ($0~/Number of hits: /) {                      nhitf=$4   }
  if ($0~/Aver.*hit r/)   {                         hitfr=$4   }
  if ($0~/nFr.*ice powd/) {                         nicef=$6   }
  if ($0~/icefinder=/)    {gsub(/icefinder=/,"")   ;icef=$0    }
  if ($0~/icefinderADC=/) {gsub(/icefinderADC=/,"");icefADC=$0 }
  if ($0~/icefinderNAT=/) {gsub(/icefinderNAT=/,"");icefNAT=$0 }
  if ($0~/Num.* sorted water hits: /) {             type0=$6   }
  if ($0~/Num.* sorted water.*50 ADUs:/) {          type0_below50=$9   }
  if ($0~/Num.* sorted water.*100 ADUs:/) {         type0_below100=$9   }
  if ($0~/Num.* sorted ice hits: /) {               type1=$6   }
  if ($0~/Num.* sorted ice.*50 ADUs:/) {            type1_below50=$9   }
  if ($0~/Num.* sorted ice.*100 ADUs:/) {           type1_below100=$9   }
  if ($0~/Num.* sorted weak hits: /) {              type2=$6   }
  if ($0~/Num.* sorted out.* hits: /) {             type3=$6   }
}
END {
  if (nfr > 0 ) {  icefr=(nicef/nfr)*100 }
  format1="%8s ,"
  format1e="%8s"
  format2="%8.2f ,"
  format3="%8d ,"
  format3e="%8d"
  format4="%s ,"
  format4b="\"%s\" ,"

  runnum = substr(run,2,4)
#  print runnum

  printf(format1,"runnum");
  printf(format1,"nxtc");
  printf(format1,"sizextc");
  printf(format1,"nfr");
  printf(format1,"darksub")
  printf(format1,"hitf")
  printf(format1,"nhitf")
  printf(format1,"hitfADC")
  printf(format1,"hitfNAT")
  printf(format1,"hitfr")
  printf(format1,"icef ")
  printf(format1,"nicef ")
  printf(format1,"icefADC")
  printf(format1,"icefNAT")
  printf(format1,"icefr")
  printf(format4,"comment")
  printf(format1,"type0")
  printf(format1,"type0 < 50ADUs")
  printf(format1,"50ADUs < type0 < 100ADUs")
  printf(format1,"type1")
  printf(format1,"type1 < 50ADUs")
  printf(format1,"50ADUs < type1 < 100ADUs")
  printf(format1,"type2")
  printf(format1e,"type3")

  printf("\n")

#  printf(format1,substr(FILENAME,0,5));
  printf(format3,runnum);
  printf(format3,nxtc);
  printf(format3,sizextc)
  printf(format3,nfr);
  printf(format3,darksub)
  printf(format3,hitf)
  printf(format3,nhitf)
  printf(format3,hitfADC)
  printf(format3,hitfNAT)
  printf(format2,hitfr)
  printf(format3,icef )
  printf(format3,nicef )
  printf(format3,icefADC)
  printf(format3,icefNAT)
  printf(format2,icefr)
  printf(format4b,comment)
  printf(format3,type0)
  printf(format3,type0_below50)
  printf(format3,type0_below100)
  printf(format3,type1)
  printf(format3,type1_below50)
  printf(format3,type1_below100)
  printf(format3,type2)
  printf(format3e,type3)

  printf("\n")
}

