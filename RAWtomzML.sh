while read f;do mono ThermoRawFileParser.exe -i="$f" -o=/hwfssz5/ST_CANCER/CGR/USER/zhuyafeng/MS_Data/LiverCancer/ -f=1 -m=1; done < liver.cancer.filelist.2.txt
#liver.cancer.filelist.2.txt 原始数据存放目录列表


