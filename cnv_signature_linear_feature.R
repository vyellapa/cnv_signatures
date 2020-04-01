library(dplyr)
library(hdp)
library(stringr)
require(plyr)
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library(stringr)
options(digits=5) 
library(RColorBrewer)
library(ConsensusClusterPlus)
library(survival)
library(survminer)

#################################################################################
###
### create color for plotting
###
#################################################################################

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n)

cnv_colors<- cols[c(1,4,12, 18,23,29 )]
cnv_colors_final<- c(rep(cnv_colors[1], 3), rep(cnv_colors[2], 8), rep(cnv_colors[3], 7), 
                     rep(cnv_colors[4], 5), rep(cnv_colors[5], 3), rep(cnv_colors[6], 10))

barplot(1:30,col=cnv_colors_final, names.arg = 1:30, las=2)

########################################################################################################
###
### upload the references for chromosome and band sizes
###
#########################################################################################################

####### ref chromosomes

setwd("~/Desktop/project/UCSC_ref_files/")
snps = read.table("hg19.chrom_sizes.txt", header=F, sep="\t")
levels(snps$V2) <- c(levels(snps$V2), "chr23","chr24")
snps$V2[snps$V2 == "chrX"] <- "chr23"
snps$V2[snps$V2 == "chrY"] <- "chr24"
chr_list <- as.character(1:24) 

cyto <- list()
for  (i in (1:24))
{
  chr_num = i
  chr_filter <- paste("chr",chr_num,sep="")
  limit <- snps[snps$V2 == chr_filter,]
  post <- limit$V3
  x <- seq(from = 0, to = (post - 1000000), by = 1000000)
  p <- seq(from = 1000000, to = post, by = 1000000)
  alfa<- length(p)
  tail <-  c(p[alfa],post)
  y <- data.frame(x,p)
  y <- rbind(y,tail)
  y$chr <- rep(chr_filter,nrow(y))
  colnames(y)<-c("start","end","chr")  
  cyto[[chr_filter]] <- y
}
caryo <- do.call("rbind", cyto) # merge results from all chromosomes  
caryo$chr<-gsub("chr", "",caryo$chr)

gr1 = with(caryo, GRanges(chr, IRanges(start=start, end=end)))

### ref chromosome 10 Mb

cyto_10mb <- list()
for  (i in c(1:22,"X"))
{
  chr_filter = i
  limit <- snps[snps$V1 == chr_filter,]
  post <- limit$V3
  x <- seq(from = 0, to = (post - 10000000), by = 10000000)
  p <- seq(from = 10000000, to = post, by = 10000000)
  alfa<- length(p)
  tail <-  c(p[alfa],post)
  y <- data.frame(x,p)
  y <- rbind(y,tail)
  y$chr <- rep(chr_filter,nrow(y))
  colnames(y)<-c("start","end","chr")  
  cyto_10mb[[chr_filter]] <- y
}
cyto_10mb_all <- do.call("rbind", cyto_10mb)
gr1_10mb = with(cyto_10mb_all, GRanges(chr, IRanges(start=start, end=end)))

##########################################################################################################################################
###
### upload CoMMpass CNV
###
#############################################################################################################################################


### commpass orgnal calls (for a original mistake colum "major" = total alleles)
cnv<- read.delim("~/Desktop/project/SV_project/2020_rebuttal/data/commpass_cnv_new_2020_merged_old.txt", stringsAsFactors = F)
cnv_clean<- cnv[cnv$major %in% c(0,1.5, 2, 2.5, 3:10),  ]
cna_all_all<- cnv_clean[,c(1:6)]
cna_mmrf<- cna_all_all[order(cna_all_all$IDA, cna_all_all$seqnames, cna_all_all$startA),]


### before collapsing segmnent with the same copy number you need to adjust some values
cna_mmrf$minor[is.na(cna_mmrf$minor)]<-0
cna_mmrf$minor[(cna_mmrf$minor)==0.5]<-0 ### transform the NA in 0

cna_mmrf$major<- ceiling(cna_mmrf$major) #### transform the 0.5 in 1
colnames(cna_mmrf)[1:4]<-c("sample","Chrom","start","end")
cnv_mmrf_sel<-cna_mmrf

#### remove IGH, IGH, IGL to reduce VDJ-CSR noise

igh_cnv <- cnv_mmrf_sel[cnv_mmrf_sel$Chrom == 14 & cnv_mmrf_sel$start >106032614 &  cnv_mmrf_sel$start< 108288051 | 
                      cnv_mmrf_sel$Chrom == 22 & cnv_mmrf_sel$start >22080474. &  cnv_mmrf_sel$end< 24065085 | 
                      cnv_mmrf_sel$Chrom == 2 & cnv_mmrf_sel$start >88090568 &  cnv_mmrf_sel$end< 91274235,]

cnv_mmrf_no_igh<- cnv_mmrf[- as.numeric(rownames(igh_cnv)),]

#### remove chromosome X to avoid overestimation of deletion

cnv_mmrf_no_gih_no_x<- cnv_mmrf_no_igh[cnv_mmrf_no_igh$Chrom!="X",]

###############################################################################################################
##
## small event problem - not clear how they handle this in the paper. 
## FFPE generates a lot of artifacts and small gains/losses, but in the Nat Gen paper this is not reported
## 
################################################################################################################

### if you want to run all segments
cnv_mmrf_final<- cnv_mmrf_no_gih_no_x
plot(density(cnv_mmrf_no_gih_no_x$end-cnv_mmrf_no_gih_no_x$start))

### if you want to run segments >100Kb - probably what you want to do
cnv_mmrf_final<- cnv_mmrf_no_gih_no_x[(cnv_mmrf_no_gih_no_x$end-cnv_mmrf_no_gih_no_x$start)>100000,]
plot(density(cnv_mmrf_final$end-cnv_mmrf_final$start))


###############################################################################################################
##
## collapse near segmenets with same copy number
## 
################################################################################################################

cnv_mmrf2<- list()
sample_list<- unique(cnv_mmrf_final$sample)
for(j in (1:length(sample_list)))
{
  cna_mmrf_single<- cnv_mmrf_final[cnv_mmrf_final$sample == sample_list[j],]
  sam_cnv_list<- list()
  chr_list<- unique(cna_mmrf_single$Chrom)
  for(i in (1:length(chr_list)))
  {
    cna_mmrf_single_chr<- cna_mmrf_single[cna_mmrf_single$Chrom == chr_list[i],]

    list_chr<- list()
    # vec<- rle((paste(cna_mmrf_single_chr$major, cna_mmrf_single_chr$minor)))$length
    
    ### just use the major to avoi possibile LOH issue
    ### LOH infact was not considered in the original paper
    
    vec<- rle((paste(cna_mmrf_single_chr$major)))$length 
    for(w in (1:length(vec)))
    {
      if(w==1){
        int<- cna_mmrf_single_chr[1:vec[w],]
        cna_mmrf_single_row<- c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
      }else{
        int<- cna_mmrf_single_chr[(sum(vec[1:(w-1)])+1):sum(vec[1:(w)]),]
        cna_mmrf_single_row<- c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
      }
      list_chr[[w]]<- cna_mmrf_single_row
    }
    list_chr2<- do.call("rbind",list_chr)
    sam_cnv_list[[i]]<- list_chr2
  }
  sam_cnv_list2<- do.call("rbind", sam_cnv_list)
  cnv_mmrf2[[j]]<-  sam_cnv_list2

}
cnv_mmrf<- do.call("rbind", cnv_mmrf2)
cnv_mmrf<- as.data.frame.matrix(cnv_mmrf)
colnames(cnv_mmrf)<-c("sample","Chrom","start","end", "major","minor")

cnv_mmrf$sample<- as.character(as.character(cnv_mmrf$sample))
cnv_mmrf$Chrom<- as.character(as.character(cnv_mmrf$Chrom))
cnv_mmrf[,3:ncol(cnv_mmrf)]<- apply(cnv_mmrf[,3:ncol(cnv_mmrf)], 2, function(x){as.numeric(as.character(x))})



###################################################################################################################
###################################################################################################################
###################################################################################################################
#####
##### COPY NUMBER SIGNATURES matrix generation
#####
###################################################################################################################
###################################################################################################################
###################################################################################################################

### list of good samples from SV paper (latest update March 2020)

sv_good<- read.delim("~/Desktop/project/SV_project/2020_rebuttal/sv_final/200323_sv_sample_list.txt", stringsAsFactors = F)
sample_list<- unique(cnv_mmrf$sample)[unique(cnv_mmrf$sample) %in% sv_good$sample]


#### this is to create an identical matrix to the one that they used. I think, it will be nice to add 
#### all these columns despite some of them will be 0 in all cases. This can be used for a biological explanation
#### right now the code works without 
# mat_sig_cnv_final<- matrix(c(paste(c(rep("1", 3), rep("2", 8), rep("3",7), rep("4",5), rep("6",3), rep("7", 10)),
#                                    c("1" , "2" , "3" , "1" , "2" , "3" , "4" , "5" , "6" , "7" , "8" , "1" , "2" , "3" , "4" ,
#                                      "5" , "6" , "7" , "1" , "2" , "3" , "4" , "5" , "1" , "2" , "3" , "1" , "2" , "3" , "4" ,
#                                      "5" , "6" , "7" , "8" , "9" , "10"), sep="_")),
#                            ncol=1)
# colnames(mat_sig_cnv_final)[1]<-"code"
# mat_sig_cnv_final<- as.data.frame(mat_sig_cnv_final)


mat_sig_cnv_final<- list()

max_10mb<- 100
max_copy_number<- max(cnv_mmrf$major)

### list for each feature

count_10mb_all<-list()
size_all<- list()
count_jump_all<- list()
count_cnv_all<- list()
band_rate_all<-list()
osci_all<- list()

## start loop for each sample

for(w in (1:length(sample_list)))
{

cnv_mmrf2<- cnv_mmrf[cnv_mmrf$sample ==sample_list[w],]

##################################################################################################################
### 
### segment size
###
##################################################################################################################

cnv_mmrf2$seg_size<- (cnv_mmrf2$end - cnv_mmrf2$start) 
size_all[[w]]<-cnv_mmrf2[,c("sample","seg_size")]

##################################################################################################################
###
### numbe rof break in 10 Mb - in this loop 10 Mb without CNV breaks will be kepts as 0
### All zero values will be removed later
###
##################################################################################################################

cnv_temp_brk<- cnv_mmrf2[,c(1,2,3,5,6,7)]
cnv_mmrf2_second<- cnv_mmrf2[,c(1,2,4,5,6,7)]
# colnames(cnv_mmrf2_first)<- colnames(cnv_mmrf2_second)
# cnv_temp_brk<- rbind.data.frame(cnv_mmrf2_second, cnv_mmrf2_first)

### remove diploid whole chromosome regions
int_dipl<- as.data.frame.matrix(table(cnv_mmrf2_second$Chrom, cnv_mmrf2_second$major))
diploid_chr<- rownames(int_dipl[int_dipl$`2`==1 & rowSums(int_dipl)==1,])
cnv_temp_brk<- cnv_mmrf2_second[! cnv_mmrf2_second$Chrom %in% diploid_chr,]
###

# cnv_temp_brk<- cnv_mmrf2_second ### just keep the start to not count two times the same breakpoint

gr_cna_comm = with(cnv_temp_brk, GRanges(Chrom, IRanges(start=end, end=end)))
values(gr_cna_comm) <- DataFrame(sample = cnv_temp_brk$sample, major = cnv_temp_brk$major, minor= cnv_temp_brk$minor, seg_size= cnv_temp_brk$seg_size)
range_10mb <- merge(as.data.frame(gr_cna_comm),as.data.frame(gr1_10mb),by="seqnames",suffixes=c("A","B"))
range_dri_10mb <- range_10mb[with(range_10mb, startB <= startA & endB >= endA),]
count_brk_10mb<- as.data.frame(table(paste(range_dri_10mb$seqnames, range_dri_10mb$startB)))
#take the ref for 10 Mb without CNV breaks
cyto_10mb_all$Var1<- paste(cyto_10mb_all$chr, cyto_10mb_all$start)
count_10mb_file <-join(cyto_10mb_all, count_brk_10mb, by="Var1")
count_10mb_file[is.na(count_10mb_file)]<-0
count_10mb_file$sample<- sample_list[w]
count_10mb_df<- count_10mb_file[,c("sample","Freq")]
colnames(count_10mb_df)[2]<-"count"
count_10mb_all[[w]]<-count_10mb_df 


#############################################################################################
###
### change in copy number adjacent- correct for chromosome - not total 
### just consider either the fiirst or the second breaks in order to not count two time the same break
###
#############################################################################################


cnv_mmrf2_second$jump<- NA
summary_jump<- as.data.frame(table(cnv_mmrf2_second$Chrom))

### keep also whole chromosome gain/loss
# summary_jump_whole_chrom<- summary_jump[summary_jump$Freq<2,]# we don't want to count for whole chromosome gain- loss or diploid
# cnv_mmrf2_second_jump<- cnv_mmrf2_second[! cnv_mmrf2_second$Chrom %in% summary_jump_whole_chrom$Var1,]

cnv_mmrf2_second_jump<- cnv_mmrf2_second

### function from the Nat Gen paper - this is just to check what they did
# allcp<-c()
# for(c in c(1:22))
# {
#   currseg<-as.numeric(cnv_mmrf2_second_jump[cnv_mmrf2_second_jump$Chrom==c,"major"])
#   allcp<-c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))
# }

if(length(unique(cnv_mmrf2_second_jump$Chrom))!=0){
chr_list_jump<- unique(cnv_mmrf2_second_jump$Chrom)
cnv_mmrf2_second_jump$jump<- NA
all_chr_jump<- list()
for(jj in chr_list_jump){
  
  cnv_mmrf2_second_jump_int<- cnv_mmrf2_second_jump[cnv_mmrf2_second_jump$Chrom == jj,]
  for(z in (1:nrow(cnv_mmrf2_second_jump_int)))
  {
    if(z==1){
      cnv_mmrf2_second_jump_int$jump[1]=NA
    }else{
      
      cnv_mmrf2_second_jump_int$jump[z]<- abs((cnv_mmrf2_second_jump_int$major[z]) - 
                                             (cnv_mmrf2_second_jump_int$major[z-1]))
    }
    
    
  }
  
  all_chr_jump[[jj]]<-cnv_mmrf2_second_jump_int
}
all_chr_jump2<- do.call("rbind", all_chr_jump)

}else{
  all_chr_jump2<- cnv_mmrf2_second[1,]
  all_chr_jump2$jump<-0
}

### from the function that they used it seems they did not count the first segment

all_chr_jump2<- all_chr_jump2[! is.na(all_chr_jump2$jump),]
temp_jump<- all_chr_jump2[,c("sample","jump")]
colnames(temp_jump)<- c("sample","count")
count_jump_all[[w]]<-temp_jump 


#############################################################################################
###
### Segment copy number: the observed absolute copy number state of each segment;
###
#############################################################################################

count_cnv_final_df<- cnv_mmrf2[,c("sample","major")]
colnames(count_cnv_final_df)[2]<-"count"
count_cnv_all[[w]]<- count_cnv_final_df[,c("sample","count")]

#############################################################################################
###
### Breakpoint count per chromosome arm: the number of breaks occurring per chromosome arm;
###
#############################################################################################

chrom_arms<- read.delim("~/Desktop/project/UCSC_ref_files/CentromerePosition_hg19.txt")
chrom_arms$chrom<- gsub("chr", "",chrom_arms$chrom)

### just keep the start to not count two times the same breakpoint

### if you have one diploid it mean that you don't have breaks - this would explain the high density peak in 0 in the suppl. figure 5
### create a file where the diploid whole chromosomes Are kept. They will be later removed
### just keep the start to not count two times the same breakpoint

cnv_temp_brk_arm <- cnv_mmrf2_second

if(nrow(cnv_temp_brk_arm)!=0){
gr_cna_comm = with(cnv_temp_brk_arm, GRanges(Chrom, IRanges(start=end, end=end)))
values(gr_cna_comm) <- DataFrame(sample = cnv_temp_brk_arm$sample, 
                                 major = cnv_temp_brk_arm$major, minor= cnv_temp_brk_arm$minor, seg_size= cnv_temp_brk_arm$seg_size)

gr_band = with(chrom_arms, GRanges(chrom, IRanges(start=chromStart, end=chromEnd)))

range_arm <- merge(as.data.frame(gr_cna_comm),as.data.frame(gr_band),by="seqnames",suffixes=c("A","B"))
range_arm$arm<- NA
range_arm$arm[range_arm$startA> range_arm$endB]<-"q_arm"
range_arm$arm[range_arm$startA< range_arm$startB]<-"p_arm"
range_10mb$arm[range_10mb$startB <= range_10mb$startA & range_10mb$endB >= range_10mb$startA]<-"centro"

table(paste(range_arm$seqnames, range_arm$arm))
db_arm_counts<- as.data.frame(table(paste(range_arm$seqnames, range_arm$arm)))
}else{
  db_arm_counts<- matrix(c("13 q_arm", 0), nrow=1)
  db_arm_counts<- as.data.frame(db_arm_counts)
  colnames(db_arm_counts)<- c("Var1","Freq")
}

file_int_band<- as.data.frame(c(paste0(c(1:22), (" p_arm")), paste0(c(1:22), (" q_arm"))))
colnames(file_int_band)<-"Var1"
band_rate<- join(file_int_band, db_arm_counts, by="Var1")
band_rate[is.na(band_rate)]<-0
band_rate$sample<- sample_list[w]
colnames(band_rate)[2]<-"count"
band_rate_all[[w]]<-band_rate[,c("sample","count")]


#############################################################################################
###
### oscillating copy number length
###
#############################################################################################

out<-c()
chrs<-unique(cnv_mmrf2$Chrom)
cnv_mmrf2$tot<- cnv_mmrf2$major 

oscCounts<-c()
for(c in chrs)
{
  currseg<-cnv_mmrf2[cnv_mmrf2$Chrom==c,"tot"]
  currseg<-round(as.numeric(currseg))
  if(length(currseg)>3)
  {
    # print(c)
    prevval<-currseg[1]
    count=0
    for(j in 3:length(currseg))
    {
      if(j==length(currseg)){
        # print("FALSE")
        oscCounts<-rbind(oscCounts,c(c,count))
        count=0
      }else{
            if(currseg[j]==prevval & currseg[j]!=currseg[j-1])
              {
              count<-count+1
              }else{
                oscCounts<-rbind(oscCounts,c(c,count))
                count=0
             }
        }
    
      prevval<-currseg[j-1]
    }
  }else{
    oscCounts<- rbind(oscCounts,c(c,0))
  }
}

oscCounts_df<- as.data.frame(oscCounts)
oscCounts_df$sample<- sample_list[w]
osci_all[[w]]<-oscCounts_df

### END OF THE LOOP FOR FEATURE COUNT EXTRACTION
}

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
###
###
### Fra code for definition of component
###
###
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

### i use a simpler mixed effect model (mclust) to define the main group for each feature.

###### size
size_all22<-do.call("rbind", size_all)
size_all22<- as.data.frame(size_all22)
size_all22$seg_size<- as.numeric(as.character(size_all22$seg_size))
size_all22_alt<- size_all22
myMclust_size <- Mclust(size_all22_alt$seg_size,G=2:10,verbose=FALSE)
# in case you want to plot the density
# pod5 = densityMclust(size_all22_alt$seg_size)
# plot(pod5, what = "density", type = "image", col = "dodgerblue3")
size_all22_alt$size_code <- myMclust_size$classification
boxplot(size_all22_alt$seg_size~size_all22_alt$size_code)
table(size_all22_alt$size_code)

###### count cnv
count_cnv_all2<-do.call("rbind", count_cnv_all)
count_cnv_all2<- as.data.frame(count_cnv_all2)
count_cnv_all2$count<- as.numeric(as.character(count_cnv_all2$count))
count_cnv_all2_alt<- count_cnv_all2
myMclust_count_cnv <- Mclust(count_cnv_all2_alt$count,G=2:5,verbose=FALSE)
# in case you want to plot the density
# pod5_cnv = densityMclust(count_cnv_all2_alt$count,G=2:5)
plot(pod5_cnv, what = "density", type = "image", col = "dodgerblue3")
count_cnv_all2_alt$cnv_count_code <- myMclust_count_cnv$classification
boxplot(count_cnv_all2_alt$count~count_cnv_all2_alt$cnv_count_code)
table(count_cnv_all2_alt$cnv_count_code)

###### count cnv per 10mb
count_10mb_all2<-do.call("rbind", count_10mb_all)
count_10mb_all2<- as.data.frame(count_10mb_all2)
count_10mb_all2$count<- as.numeric(as.character(count_10mb_all2$count))
count_10mb_all2_alt<- count_10mb_all2[count_10mb_all2$count!=0,]
myMclust_10mb <- Mclust(count_10mb_all2_alt$count,G=2:3,verbose=FALSE)
# in case you want to plot the density
# pod5_cnv = densityMclust(count_10mb_all2_alt$count,G=2:2)
# plot(pod5_cnv, what = "density", type = "image", col = "dodgerblue3")
count_10mb_all2_alt$mb_code <- myMclust_10mb$classification
boxplot(count_10mb_all2_alt$count~count_10mb_all2_alt$mb_code)
table(count_10mb_all2_alt$mb_code)

###### count jump 
count_jump_all2<-do.call("rbind", count_jump_all)
count_jump_all2<- as.data.frame(count_jump_all2)
count_jump_all2$count<- as.numeric(as.character(count_jump_all2$count))
count_jump_all2_alt<- count_jump_all2
myMclust_jump<- Mclust(count_jump_all2_alt$count,G=2:3,verbose=FALSE)
# in case you want to plot the density
# pod5_cnv = densityMclust(count_jump_all2_alt$count,G=2:2)
# plot(pod5_cnv, what = "density", type = "image", col = "dodgerblue3")
count_jump_all2_alt$cnv_count_code <- myMclust_jump$classification
boxplot(count_jump_all2_alt$count~count_jump_all2_alt$cnv_count_code)
table(count_jump_all2_alt$cnv_count_code)

###### band rate
band_rate_all2<-do.call("rbind", band_rate_all)
band_rate_all2<- as.data.frame(band_rate_all2)
band_rate_all2$Freq<- as.numeric(as.character(band_rate_all2$count))
band_rate_all2_alt<- band_rate_all2[band_rate_all2$count!=0,]
myMclust_band<- Mclust(band_rate_all2_alt$count,G=5,verbose=FALSE)
# in case you want to plot the density
# pod5_cnv = densityMclust(band_rate_all2_alt$count,G=2:3)
# plot(pod5_cnv, what = "density", type = "image", col = "dodgerblue3")
band_rate_all2_alt$band_code <- myMclust_band$classification
boxplot(band_rate_all2_alt$count~band_rate_all2_alt$band_code)
table(band_rate_all2_alt$band_code)

###### oscillation
osci_all2<-do.call("rbind", osci_all)
osci_all2<- as.data.frame(osci_all2[,c(3,2)])
colnames(osci_all2)[2]<-"count"
osci_all2$count<- as.numeric(as.character(osci_all2$count))
osci_all2_alt<- osci_all2
myMclust_osci<- Mclust(osci_all2_alt$count,G=3:3,verbose=FALSE)
# in case you want to plot the density
# pod5_cnv = densityMclust(band_rate_all2_alt$count,G=2:3)
# plot(pod5_cnv, what = "density", type = "image", col = "dodgerblue3")
osci_all2_alt$osci_code <- myMclust_osci$classification
boxplot(osci_all2_alt$count~osci_all2_alt$osci_code)


#### create matrix using mclust cluster output to define channel and counts

osci_tab<- as.data.frame.matrix(table(osci_all2_alt$sample, osci_all2_alt$osci_code))
colnames(osci_tab)<- paste0("osci_", colnames(osci_tab))
osci_tab$sample<- rownames(osci_tab)

band_tab<- as.data.frame.matrix(table(band_rate_all2_alt$sample, band_rate_all2_alt$band_code))
colnames(band_tab)<- paste0("band_", colnames(band_tab))
band_tab$sample<- rownames(band_tab)

jump_tab<- as.data.frame.matrix(table(count_jump_all2_alt$sample, count_jump_all2_alt$cnv_count_code))
colnames(jump_tab)<- paste0("jump_", colnames(jump_tab))
jump_tab$sample<- rownames(jump_tab)

mb_10_tab<- as.data.frame.matrix(table(count_10mb_all2_alt$sample, count_10mb_all2_alt$mb_code))
colnames(mb_10_tab)<- paste0("mb_10_", colnames(mb_10_tab))
mb_10_tab$sample<- rownames(mb_10_tab)

count_tab<- as.data.frame.matrix(table(count_cnv_all2_alt$sample, count_cnv_all2_alt$cnv_count_code))
colnames(count_tab)<- paste0("count_cnv_", colnames(count_tab))
count_tab$sample<- rownames(count_tab)

size_tab<- as.data.frame.matrix(table(size_all22_alt$sample, size_all22_alt$size_code))
colnames(size_tab)<- paste0("size_cnv_", colnames(size_tab))
size_tab$sample<- rownames(size_tab)

### create final matrix to be run on hdp

hdp_final<-Reduce(merge, list(size_tab, count_tab, band_tab, jump_tab, mb_10_tab,osci_tab))

########################################################################################################################
###
### run hdp for cnv extraction - 
###
#######################################################################################################################

rownames(hdp_final)<- hdp_final$sample
hdp_final2<- hdp_final[,-1]

## plot heatmap to have a visualisation of different entity and distribution
pheatmap(hdp_final2)

### set the right order following the Nat Gen paper
mat_final<- hdp_final2[,c( 23:25, 11:15, 20:22, 16:18,26:27,1:10)]

### select the right number of colors considering that here we don;t have some of the Nat Gen channels
cnv_colors_final<- c(rep(cnv_colors[1], 3), rep(cnv_colors[2], 5), rep(cnv_colors[3], 3), 
                     rep(cnv_colors[4], 3), rep(cnv_colors[5], 2), rep(cnv_colors[6], 10))

names_id<- colnames(hdp_final2)

### plot entire cohort
par(mar=c(9,4,4,4))
barplot(colSums(mat_final), las=2,
        col=cnv_colors_final, names.arg = colnames(mat_final))

### plot profile of each patient
pdf("~/Desktop/project/SV_signatures/cnv/Commpass/PLOT_PROFILE.PDF")
for(i in (1:nrow(mat_final))){
  barplot(as.numeric(mat_final[i,]), las=2,
          col=cnv_colors_final, names.arg = colnames(mat_final))
  
}
dev.off()

#################################################################################################################################
###
### run hdp
###
#################################################################################################################################

id_names<- colnames(mat_final)
genomicData<- (mat_final) ### HDP doesn't accept decimals
names_id<- colnames(genomicData)
all<- genomicData
n<- ncol(genomicData)
shape<- 1
invscale<- 1
hdp<- hdp_init(ppindex=0, #index of the parent DP for initial DP
               cpindex=1, #index of alphaa and alphab for initial DP
               hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
               alphaa=shape,
               alphab=invscale)

hdp<- hdp_adddp(hdp,
                numdp=nrow(genomicData),
                pp=1,
                cp=1)

hdp<- hdp_setdata(hdp= hdp,dpindex=1:nrow(genomicData)+1,data=genomicData)
hdp<- dp_activate(hdp,1:(nrow(genomicData)+1),10)

# Run four independent posterior sampling chains
# Takes ~30 minutes - Load pre-made copy once run the first time.

chlist <- vector("list", 4)
for (i in 1:4){
  chlist[[i]] <- hdp_posterior(hdp,
                               burnin=20000,
                               n=100,
                               space=50,
                               cpiter=3,
                               seed=i*1e4)
}

mut_example_multi <- hdp_multi_chain(chlist)

### QC HDP plots
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=500)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")


### extract signature with high cosine similarities
mut_example_multi_0.8_10_all <- hdp_extract_components(mut_example_multi, cos.merge = 0.9, min.sample = 10) #### 5 component
saveRDS(mut_example_multi_0.8_10_all, "~/Desktop/project/SV_signatures/cnv/Commpass/mut_example_multi_0.8_10_all.RDS")

### extract signature with lower cosine similarities
mut_example_multi_0.85_10_all <- hdp_extract_components(mut_example_multi, cos.merge = 0.85, min.sample = 10) #### 5 component
saveRDS(mut_example_multi_0.85_10_all, "~/Desktop/project/SV_signatures/cnv/Commpass/mut_example_multi_0.85_10_all.RDS")


#################################################################################################################################\
###
### run with cos.merge = 0.9, min.sample = 10 solution
###
#################################################################################################################################\

mut_example_multi<- mut_example_multi_0.8_10_all
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")


group_factor <- factor(names_id,  levels = names_id)
posteriorMeans<- t(comp_categ_distn(mut_example_multi)[[1]])

### print number of clusters extracted - remember that "0" can be a real signature or not. 
ncol(posteriorMeans)


rownames(posteriorMeans)<-  names_id

### create pdf with all signature plot

pdf("~/Desktop/project/SV_signatures/cnv/Commpass/cnv_signatures__no_low_0.90_10.pdf", width = 10, height = 5)
plot_comp_distn(mut_example_multi,  col=cnv_colors_final, col_nonsig = "grey", 
                cat_names=as.character(names_id))
dev.off()

x<-((mut_example_multi@comp_dp_distn))
kk<- x[["mean"]]

### plot first 10 cases contribution to see heterogenity
par(mar=c(10,4,4,4))
barplot(as.matrix(t(kk[1:10,])), col=brewer.pal(8, "Set1"), las=2, 
        names.arg = c("offset",sample_list)[1:10])
legend("topright", legend=c(0:8), col=brewer.pal(8, "Set1"), pch=15)

### plot un-annotated heatmap to see heterogenity
pheatmap(kk)

###################################################################################
###
### all possible linear regressiion combination of signatures
###
#####################################################################################

## here we test both lm cnv signatureeach combination

combinations<- combn(colnames(kk), m=2)
all_lm<- list()
for(i in (1:ncol(combinations))){
  
  fit<- lm(kk[,combinations[1,i]]~ kk[,combinations[2,i]])
  slope<- fit$coefficients[2]
  slope_st_error<- summary(fit)$coefficients[2,2]
  r_square <- summary(fit)$r.squared
  p_value<- anova(fit)$'Pr(>F)'[1]
  all_lm[[i]]<- c(combinations[,i],p_value,r_square, slope, slope_st_error )
}
all_lm2<- do.call("rbind", all_lm)
all_lm2<- as.data.frame(all_lm2)
colnames(all_lm2)<- c("first_sig","second_sig","p_value","R_square","slope", "sd_err")
all_lm2$p_value<- as.numeric(as.character(all_lm2$p_value))
all_lm2$R_square<- as.numeric(as.character(all_lm2$R_square))
comb_siig<- all_lm2[all_lm2$p_value<0.05 & all_lm2$R_square>0.1,]

#####################################################################################
###
### clinical tests - combine signature data with driver and clinical data frame 
###
#####################################################################################

mat_heat_all<- x[["mean"]]
mat_heat2<- mat_heat_all[-1,]
rownames(mat_heat2)<- rownames(mat_final)
colnames(mat_heat2)<- paste0("cnv_sig_", colnames(mat_heat2))
mat_heat2<- as.data.frame(mat_heat2)
mat_heat2$sample<- row.names(mat_heat2)

clin<- read.delim("~/Desktop/project/AMoritz/prediction_dec_2019.txt", stringsAsFactors = F)
clin_mmrf<- clin[clin$study=="MMRF",]
clin_mmrf$sample<- paste0(clin_mmrf$sample, "_1_BM")
cllust_clin<- merge(mat_heat2, clin_mmrf,by="sample" )

options("scipen"=-10, "digits"=4) 


########################################################################################
####
#### start correlation analysis -  (using both lm and wilcoxon)
####
########################################################################################

corr_file<- NULL
for(j in (2:10)){
  
  for(w in 33:ncol(cllust_clin)){
  x<- cllust_clin[,j]
  y<- cllust_clin[,w]
  fit<- lm(x~ y)
  slope<- fit$coefficients[2]
  slope_st_error<- summary(fit)$coefficients[2,2]
  r_square <- summary(fit)$r.squared
  p_value<- anova(fit)$'Pr(>F)'[1]
  if(length(unique(y))>2){p_wilcox=1}else{
    p_wilcox<- wilcox.test(x~y)$p.value}
  corr_file<- rbind(corr_file, c(colnames(cllust_clin)[j], colnames(cllust_clin)[w], 
                      p_value, r_square, slope, slope_st_error, p_wilcox))
  }
}
corr_file2<- as.data.frame.matrix(corr_file)
colnames(corr_file2)<-c("cnv_sig","genomic","p_value_lm","R_square",
                        "slope","slope_st","p_wilcox")
corr_file2[,3:ncol(corr_file2)]<- apply(corr_file2[,3:ncol(corr_file2)], 
                                        2, function(x){as.numeric(as.character(x))})

corr_file2$p_lm_adj<- p.adjust(corr_file2$p_value_lm, method = "fdr")
corr_file2$p_wilco_adj<- p.adjust(corr_file2$p_wilcox, method = "fdr")

options("scipen"=-100, "digits"=4)
corr_file2[corr_file2$p_lm_adj<0.1 & corr_file2$p_wilco_adj<0.1,]
corr_file2[corr_file2$p_value_lm<0.05 & corr_file2$R_square>0.1,]
corr_file2[corr_file2$p_value_lm<0.05 & corr_file2$p_wilcox<0.05,]

corr_file2[corr_file2$p_value_lm<0.05,]

### check double hit

boxplot(cllust_clin$cnv_sig_4~(cllust_clin$TP53+cllust_clin$delTP53))
pairwise.wilcox.test(cllust_clin$cnv_sig_4, (cllust_clin$TP53+cllust_clin$delTP53))
boxplot(cllust_clin$cnv_sig_4~(cllust_clin$gain1q21))
pairwise.wilcox.test(cllust_clin$cnv_sig_4, (cllust_clin$gain1q21))


####################################################################################################################################
####
#### check chromothripsis and other SV classes
####
####################################################################################################################################
sv_all<- read.delim("~/Desktop/project/SV_project/2020_rebuttal/sv_final/200324_sv_final.txt", stringsAsFactors = F)
missing_sample_no_sv<- sv_good$sample[! sv_good$sample %in% sv_all$sample]## samples without SVs
sv_classes<- as.data.frame.matrix(table(sv_all$sample, sv_all$PCAWG_class))
int_missing<- matrix(rep(0,ncol(sv_classes)*length(missing_sample_no_sv)), 
                                                     nrow=length(missing_sample_no_sv))
int_missing<- as.data.frame(int_missing)
colnames(int_missing)<- colnames(sv_classes)
rownames(int_missing)<- missing_sample_no_sv

sv_classes_final<- rbind.data.frame(sv_classes, int_missing)
sv_classes_final$sample<- rownames(sv_classes_final)

cllust_sv<- merge(sv_classes_final, cllust_clin,by="sample" )

corr_sv<- NULL
for(j in (2:12)){
  
  for(w in 13:21){
    x<- cllust_sv[,j]
    y<- cllust_sv[,w]
    fit<- lm(x~ y)
    slope<- fit$coefficients[2]
    slope_st_error<- summary(fit)$coefficients[2,2]
    r_square <- summary(fit)$r.squared
    p_value<- anova(fit)$'Pr(>F)'[1]
    corr_sv<- rbind(corr_sv, c(colnames(cllust_sv)[j], colnames(cllust_sv)[w], 
                                   p_value, r_square, slope, slope_st_error))
  }
}
corr_sv2<- as.data.frame.matrix(corr_sv)
colnames(corr_sv2)<-c("sv_class","cnv_sig","p_value_lm","R_square",
                        "slope","slope_st")
corr_sv2[,3:ncol(corr_sv2)]<- apply(corr_sv2[,3:ncol(corr_sv2)], 
                                        2, function(x){as.numeric(as.character(x))})

corr_sv2$p_lm_adj<- p.adjust(corr_sv2$p_value_lm, method = "fdr")
options("scipen"=-100, "digits"=6)
corr_sv2[corr_sv2$p_lm_adj<0.1,]

plot(cllust_sv$chromothripsis, cllust_sv$cnv_sig_4)

### take signature from commpass to use apobec

signature_file<- read.delim("~/Desktop/project/clock/fm6_clock_data/signatures_SNV_COMMPASS_cosine_no_mm1.txt", stringsAsFactors = F)
signature_file$apobec<- signature_file$Signature.Subs.02 + signature_file$Signature.Subs.13
signature_file2<- signature_file[,c("sampleID","apobec")]
colnames(signature_file2)[1]<-"sample"
cllust_sv<- merge(cllust_sv, signature_file2, by="sample")
cllust_sv$code_apobec<- cllust_sv$apobec
cllust_sv$code_apobec[cllust_sv$code_apobec<0.15]<-0
cllust_sv$code_apobec[cllust_sv$code_apobec>=0.15]<-1

####################################################################################################################################
####
#### annotated heatmap - this is cool
####
####################################################################################################################################

mat_heat2<- mat_heat[-1,]
kk_clust$cluster<- as.factor(kk_clust$cluster)

annot_df<- cllust_sv[,c("double_hit","chromothripsis_code","complex_code",
             "gain1q21","t_CCND1","t_MMSET","t_MAF","HRD","code_apobec")]
rownames(annot_df)<- cllust_sv$sample
mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
pheatmap(t(mat_heat2),  annotation_col = annot_df)


####################################################################################################################################
####
#### start clinical analysis
####
####################################################################################################################################
cllust_sv$clust_bad<- cllust_sv$cnv_sig_4 + cllust_sv$cnv_sig_8 + cllust_sv$cnv_sig_7
cllust_sv$double_hit<- cllust_sv$delTP53 + cllust_sv$TP53
cllust_sv$chromothripsis_code<- cllust_sv$chromothripsis
cllust_sv$chromothripsis_code[cllust_sv$chromothripsis_code>1]<-1
cllust_sv$complex_code<- cllust_sv$complex
cllust_sv$complex_code[cllust_sv$complex_code>1]<-1
options(scipen = 999)

### test some knonw variables

### chromothripsis (yes - no)
fit<- survfit(Surv(pfs_time, pfs_code) ~ chromothripsis_code, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### TP53 double hit
fit<- survfit(Surv(pfs_time, pfs_code) ~ double_hit, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### ISS
fit<- survfit(Surv(pfs_time, pfs_code) ~ ISS, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### apobec
fit<- survfit(Surv(pfs_time, pfs_code) ~ code_apobec, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### gain1q21
fit<- survfit(Surv(pfs_time, pfs_code) ~ gain1q21, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### t_MMSET
fit<- survfit(Surv(pfs_time, pfs_code) ~ t_MMSET, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### p value for cnv signatures in univariate

surv_object_os_nos <-Surv(time = cllust_sv$pfs_time, event = cllust_sv$pfs_code)

coxph(surv_object_os_nos ~ clust_bad, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_0, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_1, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_2, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_3, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_4, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_5, data = cllust_sv) ## not with 0.85 and 10
coxph(surv_object_os_nos ~ cnv_sig_6, data = cllust_sv) ## not with 0.85 and 10
coxph(surv_object_os_nos ~ cnv_sig_7, data = cllust_sv) ## not with 0.85 and 10


### with 0.9 and 10
# c2_os<- paste("cnv_sig_4","cnv_sig_6","ISS","age","ecog","TP53","t_MMSET","delTP53",
#               "gain1q21",sep="+")
# 
# c2_os<- paste("cnv_sig_1","cnv_sig_2","cnv_sig_3","cnv_sig_5","cnv_sig_6",
#               "cnv_sig_4","cnv_sig_7",
#               "ISS","age","ecog","TP53","t_MMSET","delTP53",
#               "gain1q21",sep="+")
# 
# c2_os<- paste("cnv_sig_1","cnv_sig_2","cnv_sig_3","cnv_sig_5","cnv_sig_6",
#               "cnv_sig_4","cnv_sig_7",sep="+")

c2_os<- paste("cnv_sig_4","chromothripsis_code","complex_code",
              "ISS","age","ecog","TP53","t_MMSET","delTP53",
              "gain1q21",sep="+")


c2_os<- paste("clust_bad","chromothripsis_code","complex_code",
              "ISS","age","ecog","TP53","t_MMSET","delTP53",
              "gain1q21","code_apobec",sep="+")

### removing chromothripsis and complex clust_bad is indipendent from apobec

f1_os <- as.formula(paste("Surv(pfs_time, pfs_code) ~ ",
                          c2_os, collapse= "+"))
f1_nos_os_single<- coxph(f1_os, data=cllust_sv)
ggforest(f1_nos_os_single, fontsize =1, data = cllust_sv)



#########################################################################################
###
### kaplan meier for clusters based on quartile
###
########################################################################################
quantile(cllust_sv$clust_bad)
library(gtools)
id_quant<- quantcut(cllust_sv$clust_bad, q=3, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)


quantile(cllust_sv$cnv_sig_2)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_2, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)


quantile(cllust_sv$cnv_sig_1)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_1, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)

quantile(cllust_sv$cnv_sig_3)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_3, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)

quantile(cllust_sv$cnv_sig_4)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_4, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)


quantile(cllust_sv$clust_bad)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_4, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)
# 
# 
# 
# #####################################################################################
# ###
# ### low CI for hdp - too strict probably
# ###
# #####################################################################################
# 
# ## select low CI values from hdp
# x<-((mut_example_multi@comp_dp_distn))
# kk_low<- x[["cred.int"]]
# kk2<- do.call("rbind", kk_low)
# kk2<- as.data.frame(kk2)
# dim(kk2)
# low<- as.data.frame(kk2[seq(1,nrow(kk2), by=2),])
# low[,1:ncol(low)]<- apply(low[,1:ncol(low)],2, function(x){as.numeric(as.character(x))})
# jj<- apply(low[,1:ncol(low)],2, function(x){(max(x))})
# low_high<-low
# 
# barplot(as.matrix(t(low_high[1:10,])), col=brewer.pal(8, "Set1"), las=2, 
#         names.arg = c("offset",sample_list)[1:10])
#####################################################################################
###
### create clustering based on signature - not sure is the best approach
###
#####################################################################################


# #### matrix selection
# kk<- x[["mean"]] ## mean
# # kk<- low ## low CI
# 
# mat_heat<- kk
# pheatmap(kk)
# rownames(kk)<- c("offset",rownames(mat_final))
# title=tempdir()
# d<- data.matrix(t(kk))
# d[,1:ncol(d)]<- apply(d[,1:ncol(d)],2, function(x){as.numeric(as.character(x))})
# # d = sweep(d,1, apply(d,1,median,na.rm=T))
# results = ConsensusClusterPlus(d,maxK=8,
#                                pFeature=1,
#                                title="test", 
#                                clusterAlg="hc",
#                                innerLinkage="ward.D2", 
#                                finalLinkage="ward.D2", 
#                                distance="euclidean", 
#                                seed=123456789)
# kk_clust<- as.data.frame((results[[5]]$consensusClass)) ##### 4 significant cluster
# kk_clust$sample<- c("offset",rownames(mat_final))
# colnames(kk_clust)[1]<- "cluster"
# table(kk_clust$cluster)
# kk_clust<- kk_clust[-which(kk_clust$sample == "offset"),]
# kk_clust<- kk_clust[order(kk_clust$cluster, kk_clust$sample),]

# 
# names_file<- as.data.frame(kk_clust[,-2])
# rownames(names_file)<- kk_clust$sample
# colnames(names_file)[1]<-"clusters"
# names_file$clusters<- paste0("cnv_sig_", names_file$clusters)
# mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
# mat_heat22<- mat_heat2[match(rownames(mat_heat2), kk_clust$sample),]
# 
# colnames(mat_heat22)<-paste0("cnv_sig_", colnames(mat_heat22))
# mat_heat22<- as.data.frame(mat_heat22)
# mat_heat22$sample<- rownames(mat_heat22)
# kk_clust_mat<- merge(kk_clust, mat_heat22, by="sample")
# rownames(kk_clust_mat)<- kk_clust_mat$sample
# kk_clust_mat<- kk_clust_mat[order(kk_clust_mat$cluster),]
# pheatmap(t(as.matrix(kk_clust_mat[,-c(1,2)])), annotation_col = names_file, cluster_cols = F,
#          annotation_colors = list(clusters = c("cnv_sig_1"= mycol_plus[1], 
#                                                "cnv_sig_2" = mycol_plus[2], 
#                                                "cnv_sig_3" = mycol_plus[3], 
#                                                "cnv_sig_4" = mycol_plus[4],
#                                                "cnv_sig_5" = mycol_plus[5])))




################################################################################################################################################################
################################################################################################################################################
###
###
### try the entire workflow with more collapsed signature solutions
###
###
################################################################################################################################################################
################################################################################################################################################

mut_example_multi<- mut_example_multi_0.85_10_all
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")


group_factor <- factor(names_id,  levels = names_id)
posteriorMeans<- t(comp_categ_distn(mut_example_multi)[[1]])

rownames(posteriorMeans)<-  names_id

pdf("~/Desktop/project/SV_signatures/cnv/Commpass/cnv_signatures__no_low_0.85_10.pdf", width = 10, height = 5)
plot_comp_distn(mut_example_multi,  col=cnv_colors_final, col_nonsig = "grey", 
                cat_names=as.character(names_id))
dev.off()



x<-((mut_example_multi@comp_dp_distn))
kk<- x[["mean"]]
par(mar=c(10,4,4,4))
barplot(as.matrix(t(kk[1:10,])), col=brewer.pal(8, "Set1"), las=2, 
        names.arg = c("offset",sample_list)[1:10])
legend("topright", legend=c(0:8), col=brewer.pal(8, "Set1"), pch=15)


###################################################################################
###
### all possible linear regressiion combination of signatures
###
#####################################################################################

combinations<- combn(colnames(kk), m=2)
all_lm<- list()
for(i in (1:ncol(combinations))){
  
  fit<- lm(kk[,combinations[1,i]]~ kk[,combinations[2,i]])
  slope<- fit$coefficients[2]
  slope_st_error<- summary(fit)$coefficients[2,2]
  r_square <- summary(fit)$r.squared
  p_value<- anova(fit)$'Pr(>F)'[1]
  all_lm[[i]]<- c(combinations[,i],p_value,r_square, slope, slope_st_error )
}
all_lm2<- do.call("rbind", all_lm)
all_lm2<- as.data.frame(all_lm2)
colnames(all_lm2)<- c("first_sig","second_sig","p_value","R_square","slope", "sd_err")
all_lm2$p_value<- as.numeric(as.character(all_lm2$p_value))
all_lm2$R_square<- as.numeric(as.character(all_lm2$R_square))
comb_siig<- all_lm2[all_lm2$p_value<0.05 & all_lm2$R_square>0.1,]

# plot(kk[,5], kk[,5])
# abline(lm(kk[,5]~kk[,5]))
# fit_test<- summary(lm(kk[,5]~(log(kk[,5]+1))))
# fit_test$coefficients[2]
# fit_test2<- summary(lm(kk[,3]~kk[,5]))
# fit_test2$coefficients[2]

#####################################################################################
###
### clinical tests - correlation tests
###
#####################################################################################

mat_heat_all<- x[["mean"]]
mat_heat2<- mat_heat_all[-1,]
rownames(mat_heat2)<- rownames(mat_final)
colnames(mat_heat2)<- paste0("cnv_sig_", colnames(mat_heat2))
mat_heat2<- as.data.frame(mat_heat2)
mat_heat2$sample<- row.names(mat_heat2)

clin<- read.delim("~/Desktop/project/AMoritz/prediction_dec_2019.txt", stringsAsFactors = F)
clin_mmrf<- clin[clin$study=="MMRF",]
clin_mmrf$sample<- paste0(clin_mmrf$sample, "_1_BM")
cllust_clin<- merge(mat_heat2, clin_mmrf,by="sample" )

options("scipen"=-10, "digits"=4)


########################################################################################
####
#### start correlation
####
########################################################################################

corr_file<- NULL
for(j in (2:6)){
  
  for(w in 29:ncol(cllust_clin)){
    x<- cllust_clin[,j]
    y<- cllust_clin[,w]
    fit<- lm(x~ y)
    slope<- fit$coefficients[2]
    slope_st_error<- summary(fit)$coefficients[2,2]
    r_square <- summary(fit)$r.squared
    p_value<- anova(fit)$'Pr(>F)'[1]
    if(length(unique(y))>2){p_wilcox=1}else{
      p_wilcox<- wilcox.test(x~y)$p.value}
    corr_file<- rbind(corr_file, c(colnames(cllust_clin)[j], colnames(cllust_clin)[w], 
                                   p_value, r_square, slope, slope_st_error, p_wilcox))
  }
}
corr_file2<- as.data.frame.matrix(corr_file)
colnames(corr_file2)<-c("cnv_sig","genomic","p_value_lm","R_square",
                        "slope","slope_st","p_wilcox")
corr_file2[,3:ncol(corr_file2)]<- apply(corr_file2[,3:ncol(corr_file2)], 
                                        2, function(x){as.numeric(as.character(x))})

corr_file2$p_lm_adj<- p.adjust(corr_file2$p_value_lm, method = "fdr")
corr_file2$p_wilco_adj<- p.adjust(corr_file2$p_wilcox, method = "fdr")

options("scipen"=-100, "digits"=4)
corr_file2[corr_file2$p_lm_adj<0.1 & corr_file2$p_wilco_adj<0.1,]
corr_file2[corr_file2$p_value_lm<0.05 & corr_file2$R_square>0.1,]
corr_file2[corr_file2$p_value_lm<0.05 & corr_file2$p_wilcox<0.05,]

corr_file2[corr_file2$p_value_lm<0.01,]
corr_file2[corr_file2$p_lm_adj<0.1 ,]

### check double hit
options("scipen"=-999)
boxplot(cllust_clin$cnv_sig_2~(cllust_clin$TP53+cllust_clin$delTP53))
pairwise.wilcox.test(cllust_clin$cnv_sig_2, (cllust_clin$TP53+cllust_clin$delTP53))
boxplot(cllust_clin$cnv_sig_2~(cllust_clin$gain1q21))
pairwise.wilcox.test(cllust_clin$cnv_sig_2, (cllust_clin$gain1q21))


####################################################################################################################################
####
#### check chromothripsis
####
####################################################################################################################################
sv_all<- read.delim("~/Desktop/project/SV_project/2020_rebuttal/sv_final/200324_sv_final.txt", stringsAsFactors = F)
missing_sample_no_sv<- sv_good$sample[! sv_good$sample %in% sv_all$sample]## samples without SVs
sv_classes<- as.data.frame.matrix(table(sv_all$sample, sv_all$PCAWG_class))
int_missing<- matrix(rep(0,ncol(sv_classes)*length(missing_sample_no_sv)), 
                     nrow=length(missing_sample_no_sv))
int_missing<- as.data.frame(int_missing)
colnames(int_missing)<- colnames(sv_classes)
rownames(int_missing)<- missing_sample_no_sv

sv_classes_final<- rbind.data.frame(sv_classes, int_missing)
sv_classes_final$sample<- rownames(sv_classes_final)

cllust_sv<- merge(sv_classes_final, cllust_clin,by="sample" )

corr_sv<- NULL
for(j in (2:12)){
  
  for(w in 13:17){
    x<- cllust_sv[,j]
    y<- cllust_sv[,w]
    fit<- lm(x~ y)
    slope<- fit$coefficients[2]
    slope_st_error<- summary(fit)$coefficients[2,2]
    r_square <- summary(fit)$r.squared
    p_value<- anova(fit)$'Pr(>F)'[1]
    corr_sv<- rbind(corr_sv, c(colnames(cllust_sv)[j], colnames(cllust_sv)[w], 
                               p_value, r_square, slope, slope_st_error))
  }
}
corr_sv2<- as.data.frame.matrix(corr_sv)
colnames(corr_sv2)<-c("sv_class","cnv_sig","p_value_lm","R_square",
                      "slope","slope_st")
corr_sv2[,3:ncol(corr_sv2)]<- apply(corr_sv2[,3:ncol(corr_sv2)], 
                                    2, function(x){as.numeric(as.character(x))})

corr_sv2$p_lm_adj<- p.adjust(corr_sv2$p_value_lm, method = "fdr")
options("scipen"=-100, "digits"=6)
corr_sv2[corr_sv2$p_lm_adj<0.1,]

plot(cllust_sv$chromothripsis, cllust_sv$cnv_sig_2)

### take signature from commpass to use apobec

signature_file<- read.delim("~/Desktop/project/clock/fm6_clock_data/signatures_SNV_COMMPASS_cosine_no_mm1.txt", stringsAsFactors = F)
signature_file$apobec<- signature_file$Signature.Subs.02 + signature_file$Signature.Subs.13
signature_file2<- signature_file[,c("sampleID","apobec")]
colnames(signature_file2)[1]<-"sample"
cllust_sv<- merge(cllust_sv, signature_file2, by="sample")
cllust_sv$code_apobec<- cllust_sv$apobec
cllust_sv$code_apobec[cllust_sv$code_apobec<0.15]<-0
cllust_sv$code_apobec[cllust_sv$code_apobec>=0.15]<-1

### create annotation for binary tests
options("scipen"=-1, "digits"=0)
cllust_sv$clust_bad<- cllust_sv$cnv_sig_2
cllust_sv$double_hit<- cllust_sv$delTP53 + cllust_sv$TP53
cllust_sv$chromothripsis_code<- cllust_sv$chromothripsis
cllust_sv$chromothripsis_code[cllust_sv$chromothripsis_code>1]<-1
cllust_sv$complex_code<- cllust_sv$complex
cllust_sv$complex_code[cllust_sv$complex_code>1]<-1
####################################################################################################################################
####
#### annotated heatmap
####
####################################################################################################################################
mat_heat<- cllust_sv[,c(13:17)]
rownames(mat_heat)<- cllust_sv$sample
annot_df<- cllust_sv[,c("double_hit","chromothripsis_code","complex_code",
                        "gain1q21","t_CCND1","t_MMSET","t_MAF","HRD","code_apobec")]
rownames(annot_df)<- cllust_sv$sample
mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
pheatmap(t(as.matrix(mat_heat)),  annotation_col = annot_df)


####################################################################################################################################
####
#### start clinical analysis
####
####################################################################################################################################

options(scipen = 999)

### test some knonw variables

### chromothripsis (yes - no)
fit<- survfit(Surv(pfs_time, pfs_code) ~ chromothripsis_code, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### TP53 double hit
fit<- survfit(Surv(pfs_time, pfs_code) ~ double_hit, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### ISS
fit<- survfit(Surv(pfs_time, pfs_code) ~ ISS, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### apobec
fit<- survfit(Surv(pfs_time, pfs_code) ~ code_apobec, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### gain1q21
fit<- survfit(Surv(pfs_time, pfs_code) ~ gain1q21, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### t_MMSET
fit<- survfit(Surv(pfs_time, pfs_code) ~ t_MMSET, data = cllust_sv)
ggsurvplot(fit, data = cllust_sv, pval = TRUE)

### p value for cnv signatures in univariate

surv_object_os_nos <-Surv(time = cllust_sv$pfs_time, event = cllust_sv$pfs_code)

coxph(surv_object_os_nos ~ clust_bad, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_0, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_1, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_2, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_3, data = cllust_sv)
coxph(surv_object_os_nos ~ cnv_sig_4, data = cllust_sv)
# coxph(surv_object_os_nos ~ cnv_sig_5, data = cllust_sv) ## not with 0.85 and 10
# coxph(surv_object_os_nos ~ cnv_sig_6, data = cllust_sv) ## not with 0.85 and 10
# coxph(surv_object_os_nos ~ cnv_sig_7, data = cllust_sv) ## not with 0.85 and 10


### with 0.9 and 10
# c2_os<- paste("cnv_sig_4","cnv_sig_6","ISS","age","ecog","TP53","t_MMSET","delTP53",
#               "gain1q21",sep="+")
# 
# c2_os<- paste("cnv_sig_1","cnv_sig_2","cnv_sig_3","cnv_sig_5","cnv_sig_6",
#               "cnv_sig_4","cnv_sig_7",
#               "ISS","age","ecog","TP53","t_MMSET","delTP53",
#               "gain1q21",sep="+")
# 
# c2_os<- paste("cnv_sig_1","cnv_sig_2","cnv_sig_3","cnv_sig_5","cnv_sig_6",
#               "cnv_sig_4","cnv_sig_7",sep="+")

c2_os<- paste("cnv_sig_1","cnv_sig_2","chromothripsis_code","complex_code",
              "ISS","age","ecog","TP53","t_MMSET","delTP53",
              "gain1q21",sep="+")


c2_os<- paste("cnv_sig_2","chromothripsis_code","complex_code",
              "ISS","age","ecog","TP53","t_MMSET","delTP53",
              "gain1q21","code_apobec",sep="+")

### removing chromothripsis and complex clust_bad is indipendent from apobec

f1_os <- as.formula(paste("Surv(pfs_time, pfs_code) ~ ",
                          c2_os, collapse= "+"))
f1_nos_os_single<- coxph(f1_os, data=cllust_sv)
ggforest(f1_nos_os_single, fontsize =1, data = cllust_sv)



#########################################################################################
###
### kaplan meier for clusters based on quartile
###
########################################################################################


quantile(cllust_sv$cnv_sig_2)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_2, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)


quantile(cllust_sv$cnv_sig_1)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_1, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)

quantile(cllust_sv$cnv_sig_3)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_3, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)

quantile(cllust_sv$cnv_sig_4)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_4, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)


quantile(cllust_sv$clust_bad)
library(gtools)
id_quant<- quantcut(cllust_sv$cnv_sig_4, q=4, na.rm=TRUE)
fit2_os <- survfit(surv_object_os_nos ~ id_quant, data = cllust_sv)
ggsurvplot(fit2_os, data = cllust_sv, pval = TRUE)
# 
# 
##############################################################################################################################
###
###
### Nat Gen paper data - the statistiical method that they used generate wired clsuters... hard to use and to make sense
###
###
##############################################################################################################################
# Following copy-number feature extraction we applied mixture modelling to breakdown each 
# feature distribution into mixtures of Gaussian or mixtures of Poisson distributions using the flexmix package.

## create list of all features (the second column is the count of the considered feature)
CN_features<- list()
CN_features[["segsize"]]<- size_all22

CN_features[["changepoint"]]<-count_jump_all2

CN_features[["bp10MB"]]<- count_10mb_all2

CN_features[["bpchrarm"]]<- band_rate_all2[,1:2]

CN_features[["osCN"]] <- osci_all2

CN_features[["copynumber"]]<-  count_cnv_all2


library(flexmix)
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

#### 10bp breakpooint counts
par(xpd=F)
bp10MB_mm<-fitComponent(CN_features[["bp10MB"]][,2],dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)
plot(hist(CN_features[["bp10MB"]][,2]))
bp10MB_mm
parameters(bp10MB_mm)
abline(v=parameters(bp10MB_mm), col="red")


#### segment size
segsize_mm<-fitComponent(CN_features[["segsize"]][,2],seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=10,max_comp=10)
plot(hist(CN_features[["segsize"]][,2]))
parameters(segsize_mm)
segsize_mm
abline(v=parameters(segsize_mm)[1,], col="red")

#### oscillating copy number
osCN_mm<-fitComponent(CN_features[["osCN"]][,2],dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)
plot(hist(CN_features[["osCN"]][,2],breaks =36))
parameters(osCN_mm)
osCN_mm
abline(v=parameters(osCN_mm), col="red")

#### break for arm
bpchrarm_mm<-fitComponent(CN_features[["bpchrarm"]][,2],dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=3,min_comp=5,max_comp=5)
plot(hist(CN_features[["bpchrarm"]][,2],breaks =39))
parameters(bpchrarm_mm)
bpchrarm_mm
abline(v=parameters(bpchrarm_mm), col="red")


#### changepoint 

### if you run "dist  = "BIC" you got  
# *Error in FLXfit(model = model, concomitant = concomitant, control = control,  : 7 Log-likelihood: NaN
changepoint_mm<-fitComponent(CN_features[["changepoint"]][,2],seed=seed,dist="pois",model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=7,max_comp=7)
plot(hist(CN_features[["changepoint"]][,2]))
parameters(changepoint_mm)
changepoint_mm
abline(v=parameters(changepoint_mm), col="red")

#### total copy number
copynumber_mm<-fitComponent(CN_features[["copynumber"]][,2],seed=seed,dist="pois",model_selection=model_selection,
                            nrep=nrep,min_comp=10,max_comp=10,min_prior=0.005,niter=2000)
plot(hist(CN_features[["copynumber"]][,2]))
parameters(copynumber_mm)
copynumber_mm
abline(v=parameters(copynumber_mm), col="red")


# 
# #### example of how to plot parameters that define different classes
# parameters(osCN_mm)
# c1 <- parameters(osCN_mm, component=1)
# c2 <- parameters(osCN_mm, component=2)
# c3 <- parameters(osCN_mm, component=3)
# 
# #### plot example of how to plot parameters that define different classes
# ggplot(CN_features[["osCN"]]) +
#   geom_histogram(aes(count, ..density..),colour = "black", fill = "white") +
#   geom_vline(xintercept = c1, col = "red", size = 2) + 
#   geom_vline(xintercept = c2, col = "blue", size = 2) +
#   geom_vline(xintercept = c3, col = "blue", size = 2)


### combine all the extraction
CN_components<-list(segsize=segsize_mm,
                    bp10MB=bp10MB_mm,
                    osCN=osCN_mm,
                    changepoint=changepoint_mm,
                    copynumber=copynumber_mm,
                    bpchrarm=bpchrarm_mm)



britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features,CN_components,cores=1,subcores=2)
NMF::aheatmap(britroc_sample_component_matrix,fontsize = 7,Rowv=FALSE,Colv=FALSE,legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")


pheatmap((britroc_sample_component_matrix))

######################################################################################################
###
### NNMF run - long run - check how to run it on the cluster
###
#######################################################################################################
# num_cores<- c(1:10)
# nmfalg<-"brunet"
# seed<-77777
# estim.r <- NMF::nmfEstimateRank(t(britroc_sample_component_matrix), 3:12,seed = seed,nrun=10,
#                                 verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
# V.random <- randomize(t(britroc_sample_component_matrix))
# 
# estim.r.random <- NMF::nmfEstimateRank(V.random, 3:12, seed =seed,nrun=10,
#                                        verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
# p<-plot(estim.r,estim.r.random, 
#         what = c("cophenetic", "dispersion","sparseness", "silhouette"),xname="Observed",yname="Randomised",main="")+
#   theme(axis.text=element_text(size=5),axis.title=element_text(size=5),
#         strip.text.x = element_text(size = 5),
#         strip.text.y = element_text(size = 5),
#         legend.text = element_text(size = 5),
#         legend.title = element_text(size = 7))
# g<-ggplotGrob(p)
# g[["grobs"]][[2]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[3]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[4]]$children[[4]]$size[[1]]<-0.5
# g[["grobs"]][[5]]$children[[4]]$size[[1]]<-0.5
# grid::grid.newpage()
# grid::grid.draw(g)
# 
# signature_profile<- basis(sigs)[,]
# 
# signature_profile<- signature_profile[c(11:13, 24:29, 17:23, 30:33, 14:16, 1:10),]
# 
# cnv_colors_final<- c(rep(cnv_colors[1], 3), rep(cnv_colors[2], 6), rep(cnv_colors[3], 7), 
#                      rep(cnv_colors[4], 4), rep(cnv_colors[5], 3), rep(cnv_colors[6], 10))
# pdf("~/Desktop/project/SV_signatures/cnv/TEST.pdf")
# for(i in (1:5))
# {
# x<- barplot(signature_profile[,i], las=2, xaxt="n", border = F,
#             col= cnv_colors_final,
#             space=c(rep(0, 3), 0.5, rep(0, 5), 0.5, rep(0,6), 0.5, rep(0,3), 0.5, rep(0,2),0.5, rep(0, 9)))
# 
# axis(1, at=x, labels  = rownames(signature_profile), line=0)
# abline(v=3.25, lwd=1, lty=2)
# abline(v=10.75, lwd=1, lty=2)
# abline(v=19.25, lwd=1, lty=2)
# abline(v=24.75, lwd=1, lty=2)
# abline(v=28.25, lwd=1, lty=2)
# }
# dev.off()

########################################################################################################################
###
### run hdp for cnv extraction - 
###
#######################################################################################################################

pheatmap(round(britroc_sample_component_matrix))

mat_final<- britroc_sample_component_matrix[,c(11:13, 24:29, 17:23, 30:33, 14:16, 1:10)]

names_id<- c(paste0("brk_10Mb_",1:3),
             paste0("cnv_",1:6),
             paste0("jump",1:7),
             paste0("brk_amr",1:4),
             paste0("oscill_",1:3),
             paste0("sizd",1:10))

### plot entire cohort
par(mar=c(9,4,4,4))
barplot(colSums(mat_final), las=2,
        col=cnv_colors_final, names.arg = colnames(mat_final))

# genomicData[,1:ncol(genomicData)]<- apply(genomicData[,1:ncol(genomicData)],2, function(x){as.numeric(as.character(x))})
### run hdp

genomicData<- round(mat_final[,-c(1:3)]) ### HDP doesn't accept decimals
names_id<- colnames(genomicData)
barplot(colSums(genomicData), las=2,
        col=cnv_colors_final[-c(1:3)], names.arg = colnames(genomicData))
all<- genomicData
n<- ncol(genomicData)
shape<- 1
invscale<- 1
hdp<- hdp_init(ppindex=0, #index of the parent DP for initial DP
               cpindex=1, #index of alphaa and alphab for initial DP
               hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
               alphaa=shape,
               alphab=invscale)

hdp<- hdp_adddp(hdp,
                numdp=nrow(genomicData),
                pp=1,
                cp=1)

hdp<- hdp_setdata(hdp= hdp,dpindex=1:nrow(genomicData)+1,data=genomicData)
hdp<- dp_activate(hdp,1:(nrow(genomicData)+1),10)

# Run four independent posterior sampling chains
# Takes ~50 minutes - don't run. Load pre-made copy.
chlist <- vector("list", 4)
for (i in 1:4){
  chlist[[i]] <- hdp_posterior(hdp,
                               burnin=20000,
                               n=100,
                               space=50,
                               cpiter=3,
                               seed=i*1e4)
}

mut_example_multi <- hdp_multi_chain(chlist)

par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=500)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

mut_example_multi_0.8_10_all <- hdp_extract_components(mut_example_multi, cos.merge = 0.9, min.sample = 10) #### 5 component

mut_example_multi<- mut_example_multi_0.8_10_all
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")

group_factor <- factor(colnames(names_id),  levels = colnames(names_id))
posteriorMeans<- t(comp_categ_distn(mut_example_multi)[[1]])

rownames(posteriorMeans)<-  names_id

pdf("~/Desktop/project/SV_signatures/cnv/Commpass/cnv_signatures__no_low_0.9_10.pdf", width = 10, height = 5)
plot_comp_distn(mut_example_multi,  col=cnv_colors_final[-c(1:3)], col_nonsig = "grey", 
                cat_names=as.character(names_id))
dev.off()



x<-((mut_example_multi@comp_dp_distn))
kk<- x[["mean"]]
par(mar=c(10,4,4,4))
barplot(as.matrix(t(kk[1:10,])), col=brewer.pal(8, "Set1"), las=2, 
        names.arg = c("offset",sample_list)[1:10])
legend("topright", legend=c(0:8), col=brewer.pal(8, "Set1"), pch=15)

#####################################################################################
###
### create clustering based on signature
###
#####################################################################################

### select low CI values from hdp

x<-((mut_example_multi@comp_dp_distn))
kk_low<- x[["cred.int"]]
kk2<- do.call("rbind", kk_low)
kk2<- as.data.frame(kk2)
dim(kk2)
low<- as.data.frame(kk2[seq(1,nrow(kk2), by=2),])
low[,1:ncol(low)]<- apply(low[,1:ncol(low)],2, function(x){as.numeric(as.character(x))})
jj<- apply(low[,1:ncol(low)],2, function(x){(max(x))})
low_high<-low


#### matrix selection
kk<- x[["mean"]] ## mean
# kk<- low ## low CI

mat_heat<- kk
pheatmap(kk)
rownames(kk)<- c("offset",sample_list)
title=tempdir()
d<- data.matrix(t(kk))
d[,1:ncol(d)]<- apply(d[,1:ncol(d)],2, function(x){as.numeric(as.character(x))})
# d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=8,
                               pFeature=1,
                               title="test", 
                               clusterAlg="hc",
                               innerLinkage="ward.D2", 
                               finalLinkage="ward.D2", 
                               distance="euclidean", 
                               seed=123456789)
kk_clust<- as.data.frame((results[[4]]$consensusClass)) ##### 4 significant cluster
kk_clust$sample<- c("offset",sample_list)
colnames(kk_clust)[1]<- "cluster"
table(kk_clust$cluster)
kk_clust<- kk_clust[-which(kk_clust$sample == "offset"),]


mat_heat2<- mat_heat[-1,]
rownames(mat_heat2)<- sample_list
kk_clust$cluster<- as.factor(kk_clust$cluster)
names_file<- as.data.frame(kk_clust[,-2])
rownames(names_file)<- sample_list
colnames(names_file)[1]<-"clusters"
mycol_plus<- c(brewer.pal(11,"Paired"),brewer.pal(6,"Set2"))
pheatmap(t(mat_heat2), annotation_col = names_file,
                  annotation_colors = list(clusters = c("1"= mycol_plus[1], 
                                                      "2" = mycol_plus[2], 
                                                      "3" = mycol_plus[3], 
                                                      "4" = mycol_plus[4],
                                                      "5" = mycol_plus[5])))
                                           
#####################################################################################
###
### clinical tests - correlation tests
###
#####################################################################################
clin<- read.delim("~/Desktop/project/AMoritz/prediction_dec_2019.txt", stringsAsFactors = F)
clin_mmrf<- clin[clin$study=="MMRF",]
clin_mmrf$sample<- paste0(clin_mmrf$sample, "_1_BM")
cllust_clin<- merge(kk_clust, clin_mmrf,by="sample" )
cllust_clin$cluster<- as.factor(cllust_clin$cluster)
surv_object_os_nos <-Surv(time = cllust_clin$pfs_time, event = cllust_clin$pfs_code)
fit2_os <- survfit(surv_object_os_nos ~ cluster, data = cllust_clin)

ggsurvplot(fit2_os, data = cllust_clin, pval = TRUE)

c2_os<- paste("cluster","ISS","age","ecog",
              "gain1q21",sep="+")

c2_os<- paste("cluster",sep="+")

f1_os <- as.formula(paste("Surv(pfs_time, pfs_code) ~ ",
                              c2_os, collapse= "+"))
                        
f1_nos_os_single<- coxph(f1_os, data=cllust_clin)

ggforest(f1_nos_os_single, fontsize =1, data = cllust_clin)

