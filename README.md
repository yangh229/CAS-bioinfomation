# CAS-bioinfomation
keep studying, record the way of doing reseach

##this is my first day to make reserch dairy

过滤数据，某个染色体
bcftools view -r chr6 Sample_NA18525.gvcf.gz -o chr6.vcf.gz -O z

vcf文件，比较位点，相同差异
vcf-compare GPU.testM.vcf.gz PGG.testM.vcf.gz

查看进程
 ps aux | grep yanghong
 
 加索引 
tabix

统计行数目
wc -l filename

一.PCA图
1.plink  将所有plink二进制文件合并在一起
输入文件 bed  bim  fam  
输出文件 bed  bim  fam  log nosex
命令行
plink --bfile data1 --bmerge data2.bed data2.bim data2.fam --make-bed --out merge
plink2 --bfile $filename --export vcf      --out $filename




2.flashpc 对合并好的文件进行pca分析
输入文件 bed文件
输出文件 pcs.txt（详细文件） pve.txt（每个主成分的重要性，降序） eigvalues.txt eigenvectors.txt
命令行
flashpca --bfile data_pruned

#flashpca --bfile  data_pruned --check --outvec eigenvectors.txt --outval eigenvalues.txt

3.R包ggplot2 处理一下pcs.txt数据格式后进行画图（未做完）
library(ggplot2)
setwd("d:r.job")
#pcs.out，pro_color_shape.txt ##pcs.txt.guolv.province,pro_color_shape_guolv.txt 
##pcs.guolv.002.txt.province,pro_color_shape_guolv.002.txt
##005.nosh.pcs.txt
pcs <- read.table(file = '005.nosh.pcs.txt.p',header = T, stringsAsFactors = F)
info<- read.table(file = '005.pro_color_shape.txt',header = T, stringsAsFactors = F)
colnames(info)<-c('province','color','shape')
Data<-factor(pcs$province,levels = info$province)

p<- ggplot(pcs)+geom_point(aes(x=PC1,y=PC2,color=Data,shape=Data),size=4)

p<-p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor= element_blank())
p<-p+labs(x='PC1(1.30%)',y='PC2(0.83%)')

##1.28%   1.00%     #0.808%  0.698%

p<-p+scale_color_manual(values = info$color)+scale_shape_manual(values = info$shape)
p


做出PCA图
筛选掉 PCA比较离谱的点
返回源文件进行筛个体：Plink神舟:
plink2  --vcf data_pruned.vcf --keep  003.pcs.guolv --make-bed --out 003.data_pruned_guolv 
###--double id
二.计算FST

1.利用plink 将所有bed文件转成vcf格式
plink2 --bfile data_pruned --export vcf --out data_pruned

2.对vcf格式进行过滤(过滤掉 indel 以及多肽位点)
somepy/guolv.py

2.两两对亚群体进行fst计算
somepy/a_b_fst.sh

3.计算后进行排序sort  fst_guolv.out_sort.sh
sort -r -k 2 -t \t  $file -o  $filename

4.排序后取前1%突变的位点
somepy/sort_1%.sh

5.利用vep对位点进行注释  #物种（未解决）
a.把需要的数据提取出来重做一遍 003.vcf.guolv.fst.py
b.所有.1%文件合成一个文件   由于 文件很小  使用  cat  *.1%> han_taiwan.vcf.guolv.fst.sort.1%.together
c.转成师兄的数据格式han_taiwan.vcf.guolv.fst.sort.1%.together.vep.py
d.vep运行命令行
/picb/humpopg/zhangxiaoxi/software/anaconda/bin/vep -i AX-AM_AX-AT.fst_guolv.out_sort.1%  -o AX-AM_AX-AT.fst_guolv.out_sort.1%.vep --species human --fasta /picb/humpopg/zhangxiaoxi/Han.taiwan/human_g1k_v37.fasta

#
infile='han_taiwan.vcf.guolv.fst.sort.1%.together.vep'
outfile='han_taiwan.vcf.guolv.fst.sort.1%.together.vep.result'

/picb/humpopg/zhangxiaoxi/software/anaconda/bin/vep -cache --cache_version 92 --dir_cache /picb/humpopg/zhangxiaoxi/vep.cache -i ${infile} --port 3337 --force_overwrite --sift b --canonical --symbol -o ${outfile}  --vcf
#
###fst的概率密度图
library(ggplot2)
data=read.table('d:r.job/1000fst',header = F,col.names = 'fst')
p<-ggplot(data,aes(x=fst,alpha = 1/10000))
p<-p+  geom_density()+xlim(-0.0015,0.0025)
p


三.进行F3  test  软件qp3pop
命令行：/picb/humpopg/zhangxiaoxi/software/anaconda/bin/qp3Pop -p parfile >logfile
对合并好的han_taiwan.vcf.guolv.fst.sort.1%.together文件进行格式转换
parfile: Name of parameter file
DESCRIPTION OF EACH PARAMETER in parfile:
genotypename:   input genotype file (in eigenstrat format)
snpname:   input snp file      (in eigenstrat format)
indivname:   input indiv file    (in eigenstrat format)
popfilename:  list_qp3test (contains 3 populations on each line <Source1 (A)> <Source2 (B)> < Target (C)>

首先需要3个输入文件genotype,snp，indiv
snp文件：/somepy/qp3pop.snp.sh


genetype文件：
对.vcf.guolv文件处理  vcf.guolv.qp3pop.genetype.py  and  .sh
together处理后的genetype文件  


indiv文件
vcf.guolv.sq3pop.indiv1.py
vcf.guolv.sq3pop.indiv1.sh

f3test绘图
#########qp3pop test 绘图
#install.packages("Hmisc")
library('Hmisc')
f3dat = read.table("d:/r.job/qp3Pop.out.1",header = F ,col.names=c("result","PopA", "PopB", "PopC", "F3", "StdErr", "Z", "SNPs"))
s = f3dat[f3dat$PopC == "TW-HB",] ###target
head(s[order(-s$F3),])
sOrdered = s[order(s$F3),]

##1:21   21 ge popA popB
errbar(1:21, sOrdered$F3[1:21],
       (sOrdered$F3+sOrdered$StdErr)[1:21],
       (sOrdered$F3-sOrdered$StdErr)[1:21], pch=20, las=2, cex.axis=0.8, xaxt='n',
       xlab="Target_TW-HB", ylab="F3")
a<-sOrdered$PopA[1:21]
b<-sOrdered$PopB[1:21]
c<-paste(a,b,sep='_')
axis(1, at=1:21, labels=c, las=2, cex.axis=0.8,pos = 0.005)
##pos   side


TIPS：
外网登陆 ssh wangchenghu@gate2.picb.ac.cn
登录上自己服务器 ssh wangchenghu@10.10.118.135




 


