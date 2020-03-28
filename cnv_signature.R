library(dplyr)
library(hdp)
library(stringr)
require(plyr)
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library(stringr)
options(digits=20) 
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
### upload the references
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

cna_mmrf[cna_mmrf$IDA == "MMRF_1031_1_BM",]

### before collapsing segmnent with the same copy number you need to adjust some values
cna_mmrf$minor[is.na(cna_mmrf$minor)]<-0
cna_mmrf$minor[(cna_mmrf$minor)==0.5]<-0 ### transform the NA in 0

cna_mmrf$major<- ceiling(cna_mmrf$major) #### transform the 0.5 in 1
colnames(cna_mmrf)[1:4]<-c("sample","Chrom","start","end")

#### remove IGH, IGH, IGL 

igh_cnv <- cnv_mmrf[cnv_mmrf$Chrom == 14 & cnv_mmrf$start >106032614 &  cnv_mmrf$start< 108288051 | 
                      cnv_mmrf$Chrom == 22 & cnv_mmrf$start >22080474. &  cnv_mmrf$end< 24065085 | 
                      cnv_mmrf$Chrom == 2 & cnv_mmrf$start >88090568 &  cnv_mmrf$end< 91274235,]

cnv_mmrf_no_igh<- cnv_mmrf[- as.numeric(rownames(igh_cnv)),]

cnv_mmrf_no_gih_no_x<- cnv_mmrf_no_gih[cnv_mmrf_no_gih$Chrom!="X",]


### collapse near segments with identical

cnv_mmrf2<- list()
sample_list<- unique(cnv_mmrf_no_gih_no_x$sample)
for(j in (1:length(sample_list)))
{
  cna_mmrf_single<- cnv_mmrf_no_gih_no_x[cnv_mmrf_no_gih_no_x$sample == sample_list[j],]
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


## small event problem - not clear how they handle this in the paper. 
## FFPE generates a lot of artifacts and small gains/losses
## in the


jj<-cnv_mmrf_final[cnv_mmrf_final == "MMRF_1016_1_BM",]
jj$diff<- jj$end- jj$start
jj

### if you want to run all
cnv_mmrf_final<- cnv_mmrf
plot(density(cnv_mmrf_final$end-cnv_mmrf_final$start))

### if you want to run segments >100Kb
cnv_mmrf_final<- cnv_mmrf[(cnv_mmrf$end-cnv_mmrf$start)>100000,]
plot(density(cnv_mmrf_final$end-cnv_mmrf_final$start))
###################################################################################################################
###################################################################################################################
###################################################################################################################
#####
##### COPY NUMBER SIGNATURES
#####
###################################################################################################################
###################################################################################################################
###################################################################################################################

### het list of good samples from SV paper (latest update March 2020)

sv_good<- read.delim("~/Desktop/project/SV_project/2020_rebuttal/sv_final/200323_sv_sample_list.txt", stringsAsFactors = F)
sample_list<- unique(cnv_mmrf$sample)[unique(cnv_mmrf$sample) %in% sv_good$sample]


####
# mat_sig_cnv_final<- matrix(c(paste(c(rep("1", 3), rep("2", 8), rep("3",7), rep("4",5), rep("6",3), rep("7", 10)),
#                                    c("1" , "2" , "3" , "1" , "2" , "3" , "4" , "5" , "6" , "7" , "8" , "1" , "2" , "3" , "4" ,
#                                      "5" , "6" , "7" , "1" , "2" , "3" , "4" , "5" , "1" , "2" , "3" , "1" , "2" , "3" , "4" ,
#                                      "5" , "6" , "7" , "8" , "9" , "10"), sep="_")),
#                            ncol=1)
# colnames(mat_sig_cnv_final)[1]<-"code"
# mat_sig_cnv_final<- as.data.frame(mat_sig_cnv_final)
which(sample_list=="MMRF_1210_1_BM")


mat_sig_cnv_final<- list()
## start loop for each sample
for(w in (1:length(sample_list)))
{
  
# print(w) 
cnv_mmrf2<- cnv_mmrf_final[cnv_mmrf_final$sample ==sample_list[w],]

##################################################################################################################
### 
### segment size
###
##################################################################################################################

cnv_mmrf2$seg_size<- (cnv_mmrf2$end - cnv_mmrf2$start)/100000 # corrected 10^7

##################################################################################################################
###
### numbe rof break in 10 Mb
###
##################################################################################################################

cnv_mmrf2_first<- cnv_mmrf2[,c(1,2,3,5,6,7)]
cnv_mmrf2_second<- cnv_mmrf2[,c(1,2,4,5,6,7)]
# colnames(cnv_mmrf2_first)<- colnames(cnv_mmrf2_second)
# cnv_temp_brk<- rbind.data.frame(cnv_mmrf2_second, cnv_mmrf2_first)

cnv_temp_brk<- cnv_mmrf2_second ### just keep the start to not count two times the same breakpoint

gr_cna_comm = with(cnv_temp_brk, GRanges(Chrom, IRanges(start=end, end=end)))
values(gr_cna_comm) <- DataFrame(sample = cnv_temp_brk$sample, major = cnv_temp_brk$major, minor= cnv_temp_brk$minor, seg_size= cnv_temp_brk$seg_size)

range_10mb <- merge(as.data.frame(gr_cna_comm),as.data.frame(gr1_10mb),by="seqnames",suffixes=c("A","B"))
range_dri_10mb <- range_10mb[with(range_10mb, startB <= startA & endB >= endA),]
count_brk_10mb<- as.data.frame(table(paste(range_dri_10mb$seqnames, range_dri_10mb$startB)))

#############################################################################################
###
### change in copy number adjacent- correct for chromosome - not total
###
#############################################################################################

# again we don't want to count multiple times the same copy number change
# we don't want to count for whole chromosome gain- loss or diploid

cnv_mmrf2_second$jump<- NA
summary_jump<- as.data.frame(table(cnv_mmrf2_second$Chrom))

summary_jump_whole_chrom<- summary_jump[summary_jump$Freq<2,]
cnv_mmrf2_second_jump<- cnv_mmrf2_second[! cnv_mmrf2_second$Chrom %in% summary_jump_whole_chrom$Var1,]
# cnv_mmrf2$jump<- NA

# for(i in (2:nrow(cnv_mmrf2_second_jump)))
# {
#   if(i==1){
#     cnv_mmrf2_second_jump$jump[1]=0
#   }else{
#     cnv_mmrf2_second_jump$jump[i]<- abs((cnv_mmrf2_second_jump$major[i]) - (cnv_mmrf2_second_jump$major[i+1]))
#   }
# }
# cnv_mmrf2_second_jump$jump[is.na(cnv_mmrf2_second_jump$jump)]<-0
if(length(unique(cnv_mmrf2_second_jump$Chrom))!=0){
chr_list_jump<- unique(cnv_mmrf2_second_jump$Chrom)
cnv_mmrf2_second_jump$jump<- NA
all_chr_jump<- list()
for(jj in chr_list_jump){
  
  cnv_mmrf2_second_jump_int<- cnv_mmrf2_second_jump[cnv_mmrf2_second_jump$Chrom == jj,]
  for(z in (1:nrow(cnv_mmrf2_second_jump_int)))
  {
    if(z==1){
      cnv_mmrf2_second_jump_int$jump[1]=0
    }else{
      
      cnv_mmrf2_second_jump_int$jump[z]<- abs((cnv_mmrf2_second_jump_int$major[z]) - 
                                             (cnv_mmrf2_second_jump_int$major[z-1]))
    }
    
    
  }
  
  all_chr_jump[[jj]]<-cnv_mmrf2_second_jump_int
}
all_chr_jump2<- do.call("rbind", all_chr_jump)
all_chr_jump2$jump[is.na(all_chr_jump2$jump)]<-0
}else{
  all_chr_jump2<- cnv_mmrf2_second[1,]
  all_chr_jump2$jump<-0
}

### if you want to look for chromosome - i'm not really sure that this is correct
int_jump<- as.data.frame.matrix(table(all_chr_jump2$Chrom, all_chr_jump2$jump)) 

#############################################################################################
###
### Segment copy number: the observed absolute copy number state of each segment;
###
#############################################################################################

count_segments<- as.data.frame(table(cnv_mmrf2$major))

#############################################################################################
###
### Breakpoint count per chromosome arm: the number of breaks occurring per chromosome arm;
###
#############################################################################################

chrom_arms<- read.delim("~/Desktop/project/UCSC_ref_files/CentromerePosition_hg19.txt")
chrom_arms$chrom<- gsub("chr", "",chrom_arms$chrom)
# colnames(snps)<-c("chr","chrom","end")
# arm<- merge(snps, chrom_arms, by="chrom")

### just keep the start to not count two times the same breakpoint

### if you have one diploid it mean that you don't have breaks - this would explain the high density peak in 0 in the suppl. figure 5
### create a file where the diploid whole chromosomes are removed
### just keep the start to not count two times the same breakpoint

filter_diploid<- as.data.frame.matrix(table(cnv_mmrf2$Chrom, cnv_mmrf2$major))
filter_diploid$sum<- rowSums(filter_diploid)
chr_diploid_all<- filter_diploid[filter_diploid$sum ==1 & filter_diploid$`2`==1,]


cnv_temp_brk <- cnv_mmrf2_second
cnv_temp_brk_arm<- cnv_temp_brk[! cnv_temp_brk$Chrom %in% rownames(chr_diploid_all), ]

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

#############################################################################################
###
### oscillating copy number length
###
#############################################################################################

# osCN<-getOscilation(CN_data,chrlen)
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

# lengthg_oscillating_cnv<- length(oscCounts_df$V2[! oscCounts_df$V2 %in% 0]) ## this is the number of total oscillating event

############################################################################################
####
#### create final counts
####
############################################################################################

### Breakpoint count per 10 Mb
# plot(hist(count_brk_10mb$Freq))

brrk_10vec<-  count_brk_10mb$Freq
brrk_10vec[brrk_10vec>3]<-3
brrk_10_final<- (table(brrk_10vec)) ### final count
missing_10<- as.character(1:3)[! as.character(1:3) %in% names(brrk_10_final)  ]
missing_10_2<- rep(0,length(missing_10))
names(missing_10_2)<- missing_10
brrk_10vec_seg<- c(brrk_10_final, missing_10_2)[order(as.numeric(names(c(brrk_10_final, missing_10_2))))] ### final count

### Copy number - copy number segments
# plot(density(cnv_mmrf2$major))

cnv_seg_all<- cnv_mmrf2$major
cnv_seg_all[cnv_seg_all>8]<-8
# cnv_seg_all[cnv_seg_all==0]<- 1
cnv_seg_all<- cnv_seg_all[cnv_seg_all!=0]
missing<- as.character(1:8)[! as.character(1:8) %in% names(table(cnv_seg_all))  ]
missing_2<- rep(0,length(missing))
names(missing_2)<- missing
cnv_count_seg<- c(table(cnv_seg_all), missing_2)[order(as.numeric(names(c(table(cnv_seg_all), missing_2))))] ### final count

### Copy number change point - adiacent segment
# plot(hist(cnv_mmrf2$jump[cnv_mmrf2$jump !=0]), xlim=c(0,10))

all_chr_jump2$jump[all_chr_jump2$jump>7]<-7
cnv_jump <- table(all_chr_jump2$jump[all_chr_jump2$jump !=0])
missing_j<- as.character(1:7)[! as.character(1:7) %in% names(cnv_jump)  ]
missing_j_2<- rep(0,length(missing_j))
names(missing_j_2)<- missing_j
cnv_count_jump<- c(cnv_jump, missing_j_2)[order(as.numeric(names(c(cnv_jump, missing_j_2))))] ### final count

###Breakpoint count per chromosome arm
# barplot(table(db_arm_counts$Freq))

arm_count<- as.numeric(db_arm_counts$Freq)
arm_count[arm_count>5]<- 5
arm_count_tab<- table(arm_count)
missing_arm<- as.character(1:5)[! as.character(1:5) %in% names(arm_count_tab)  ]
missing_arm_j_2<- rep(0,length(missing_arm))
names(missing_arm_j_2)<- missing_arm
cnv_count_arm<- c(arm_count_tab, missing_arm_j_2)[order(as.numeric(names(c(arm_count_tab, missing_arm_j_2))))] ### final count


### Oscillating copy number length
# lengthg_oscillating_cnv
# plot(hist(as.numeric(oscCounts_df$V2)))

osci_cnv<- as.numeric(as.character(oscCounts_df$V2))
osci_cnv[osci_cnv>3]<-3
osci_cnv<- osci_cnv[osci_cnv !=0]

osci_cnv_tab<- table(osci_cnv)
missing_osci<- as.character(1:3)[! as.character(1:3) %in% names(osci_cnv_tab)  ]
missing_osci_j_2<- rep(0,length(missing_osci))
names(missing_osci_j_2)<- missing_osci
cnv_count_osci<- c(osci_cnv_tab, missing_osci_j_2)[order(as.numeric(names(c(osci_cnv_tab, missing_osci_j_2))))] ### final count

### Segment size

# plot(hist(cnv_mmrf2$seg_size))

## remove diploid reegions
cnv_mmrf2_size<- cnv_mmrf2[cnv_mmrf2$major!=2,]

seg_zie_count<- round(cnv_mmrf2_size$seg_size)
seg_zie_count[seg_zie_count>10]<-10
seg_zie_count<- seg_zie_count[seg_zie_count!=0]
size_cnv_tab<- table(seg_zie_count)
missing_size<- as.character(1:10)[! as.character(1:10) %in% names(size_cnv_tab)  ]
missing_size_j_2<- rep(0,length(missing_size))
names(missing_size_j_2)<- missing_size
cnv_count_size<- c(size_cnv_tab, missing_size_j_2)[order(as.numeric(names(c(size_cnv_tab, missing_size_j_2))))] ### final count

### final vector

vec_cnv_Sig<- c(brrk_10vec_seg,cnv_count_seg, cnv_count_jump,  cnv_count_arm,cnv_count_osci,cnv_count_size)
if(length(vec_cnv_Sig)==36){
mat_sig_cnv_final[[w]]<- vec_cnv_Sig

mat_barplot<- as.data.frame(cbind(as.numeric(vec_cnv_Sig), 
                                  c(rep("1", 3), rep("2", 8), rep("3",7), rep("4",5), rep("6",3), rep("7", 10))))
colnames(mat_barplot)<-c("num","id")
mat_barplot$num<- as.numeric(as.character(mat_barplot$num))
mat_barplot$ref<- 1:nrow(mat_barplot)


pdf(sprintf("%s_cnv_sig_profile.pdf",paste0("~/Desktop/project/SV_signatures/cnv/Commpass/PLOT_COMMPASS/",sample_list[w])), 
    height=6, width=15)
x<- barplot(vec_cnv_Sig, las=2, xaxt="n", border = F,
            col= cnv_colors_final, main=sample_list[w],
            space=c(rep(0, 3), 0.5, rep(0, 7), 0.5, rep(0,6), 0.5, rep(0,4), 0.5, rep(0,2),0.5, rep(0, 9)))
axis(1, at=x, labels  = names(vec_cnv_Sig), line=0)
abline(v=3.25, lwd=1, lty=2)
abline(v=11.75, lwd=1, lty=2)
abline(v=19.25, lwd=1, lty=2)
abline(v=24.75, lwd=1, lty=2)
abline(v=28.25, lwd=1, lty=2)
dev.off()
}else{print (w)
  break}
}

mat_sig_cnv_final2<- do.call("rbind", mat_sig_cnv_final)
mat_sig_cnv_final2<- as.data.frame(mat_sig_cnv_final2)
rownames(mat_sig_cnv_final2)<- sample_list
colnames(mat_sig_cnv_final2)<- 1:length(mat_sig_cnv_final2)
mat_sig_cnv_final2[,1:length(mat_sig_cnv_final2)]<- apply(mat_sig_cnv_final2[,1:length(mat_sig_cnv_final2)], 2, 
                                                          function(x){as.numeric(as.character(x))})

########################################################################################################################
###
### run hdp for cnv extraction
###
#######################################################################################################################
names_id<- c(paste0("brk_10Mb_",1:3),
             paste0("cnv_",1:8),
             paste0("jump",1:7),
             paste0("brk_amr",1:5),
             paste0("oscill_",1:3),
             paste0("sizd",1:10 ))

### plot entire cohort
par(mar=c(9,4,4,4))
barplot(colSums(mat_sig_cnv_final2), las=2,
        col=cnv_colors_final, names.arg = names_id)
colnames(mat_sig_cnv_final2)<-names_id
pheatmap(as.matrix(t(mat_sig_cnv_final2)))


### run hdp
genomicData<- mat_sig_cnv_final2
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
                               n=50,
                               space=50,
                               cpiter=3,
                               seed=i*1e4)
}

mut_example_multi <- hdp_multi_chain(chlist)

par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=500)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

mut_example_multi_0.8_10_all <- hdp_extract_components(mut_example_multi, cos.merge = 0.85, min.sample = 10) #### 5 component

mut_example_multi<- mut_example_multi_0.8_10_all
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")

group_factor <- factor(colnames(names_id),  levels = colnames(names_id))
posteriorMeans<- t(comp_categ_distn(mut_example_multi)[[1]])

rownames(posteriorMeans)<-  names_id

pdf("~/Desktop/project/SV_signatures/cnv/Commpass/cnv_signatures__no_low_0.9_10.pdf", width = 10, height = 5)
plot_comp_distn(mut_example_multi,  col=cnv_colors_final, col_nonsig = "grey", 
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
kk_clust<- as.data.frame((results[[5]]$consensusClass)) ##### 4 significant cluster
kk_clust$sample<- c("offset",sample_list)
colnames(kk_clust)[1]<- "cluster"
table(kk_clust$cluster)
kk_clust<- kk_clust[-which(kk_clust$sample == "offset"),]

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

c2_os<- paste("cluster","ISS","age","ecog",sep="+")

f1_os <- as.formula(paste("Surv(pfs_time, pfs_code) ~ ",
                              c2_os, collapse= "+"))
                        
f1_nos_os_single<- coxph(f1_os, data=cllust_clin)

ggforest(f1_nos_os_single, fontsize =1, data = cllust_clin)

