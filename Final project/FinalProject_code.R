# FCBB - Final project Spring 2022
# Adam Kenet, Sophia Nehs and Sofía García del Barrio
# 05/02/2022

# Associations Between Mutational Signatures and Immune Cell Infiltration in Lung 
# Adenocarcinoma (LUAD) and Lung Squamous Cell Carcinoma (LUSC)

# ========================= 1. Get data for LUAD and LUSC =======================
library('maftools')

# Load the TCGA cohort
luad.tcga = tcgaLoad(study = "LUAD")
lusc.tcga = tcgaLoad(study = "LUSC")

# Finding mutational signatures

library('BSgenome') # Requires BSgenome object - human genome reference: hg38
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
ref_gen19 <- BSgenome.Hsapiens.UCSC.hg19

luad.tnm = trinucleotideMatrix(maf = luad.tcga, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
lusc.tnm = trinucleotideMatrix(maf = lusc.tcga, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

# get cophenetic plot to know number of signatures
library('NMF')
#luad.sign = estimateSignatures(mat = luad.tnm, nTry = 6,pConstant = 1e-9) # drop in n=2
#lusc.sign = estimateSignatures(mat = lusc.tnm, nTry = 6) # drop in n=2

# get contribution to each signature. De novo analysis
luad.sig = extractSignatures(mat = luad.tnm, n = 2)
lusc.sig = extractSignatures(mat = lusc.tnm, n = 2)

# calculate cosine similarity to COSMIC
luad.cosm = compareSignatures(nmfRes = luad.sig, sig_db = "legacy") # COSMIC version 2, the one we use
lusc.cosm = compareSignatures(nmfRes = lusc.sig, sig_db = "legacy")
#luad.cosm2 = compareSignatures(nmfRes = luad.sig, sig_db = "SBS") # COSMIC version 3, not used
#lusc.cosm2 = compareSignatures(nmfRes = lusc.sig, sig_db = "SBS")

# =============== 2. Show results from TCGA data ========================

# TCGA ============================
# Shows sample summary
getSampleSummary(luad.tcga)
getSampleSummary(lusc.tcga)
# Shows gene summary
getGeneSummary(luad.tcga)
getGeneSummary(lusc.tcga)
# shows clinical data associated with samples
getClinicalData(luad.tcga)
getClinicalData(lusc.tcga)
# Shows all fields in MAF
getFields(luad.tcga)
getFields(lusc.tcga)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = luad.tcga, basename = 'luad')
write.mafSummary(maf = lusc.tcga, basename = 'lusc')

# plot all summaries
plotmafSummary(maf = luad.tcga, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = lusc.tcga, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# plot TITV summary: mutation distribution
luad.titv = titv(maf = luad.tcga, plot = FALSE, useSyn = TRUE)
plotTiTv(res = luad.titv)
lusc.titv = titv(maf = lusc.tcga, plot = FALSE, useSyn = TRUE)
plotTiTv(res = lusc.titv)

# Signatures ================
# show cosine similarity heatmaps
library('pheatmap')
pheatmap::pheatmap(mat = luad.cosm$cosine_similarities, cluster_rows = FALSE, main = "LUAD: cosine similarity against validated COSMIC signatures")
pheatmap::pheatmap(mat = lusc.cosm$cosine_similarities, cluster_rows = FALSE, main = "LUSC: cosine similarity against validated COSMIC signatures")

 # mutation distribution per predominant signature
maftools::plotSignatures(nmfRes = luad.sig, title_size = 1.2, sig_db = "legacy")
maftools::plotSignatures(nmfRes = lusc.sig, title_size = 1.2, sig_db = "legacy")

# ========================= 3. Get ESTIMATE scores =============================

# read table obtained from ESTIMATE website. Meaning: 
  # stromal score: presence of stroma in tumor tissue,
  # immune score: infiltration of immune cells in tumor tissue
  # estimate score: infers tumor purity
luad_estimate <- read.table("luad_estimate.txt", header=TRUE, sep= "\t")
lusc_estimate <- read.table("lusc_estimate.txt", header=TRUE, sep= "\t")

# shorten the ID to match TCGA
for (i in 1:nrow(lusc_estimate)) {
  long_ID <- lusc_estimate$ID[i]
  chop_ID <- substr(long_ID,1,12)
  lusc_estimate$ID[i] <- chop_ID
  
}
for (i in 1:nrow(luad_estimate)) {
  long_ID <- luad_estimate$ID[i]
  chop_ID <- substr(long_ID,1,12)
  luad_estimate$ID[i] <- chop_ID
  
}

# ===================== 4. Create datafame for patient summary ============================
# has ID, Contribution to main signatures, Immune Score and Signature group

# patient ID and contributions to each signature (de novo)
lusc_contr <- t(lusc.sig[["contributions"]])
lusc_contr <- as.data.frame(lusc_contr)
lusc_contr <- cbind(rownames(lusc_contr), lusc_contr)
names(lusc_contr) <- c('ID', 'COSMIC_4', 'COSMIC_7') # sig1=Cosmic4, sig2=Cosmic7

luad_contr <- t(luad.sig[["contributions"]])
luad_contr <- as.data.frame(luad_contr)
luad_contr <- cbind(rownames(luad_contr), luad_contr)
names(luad_contr) <- c('ID', 'COSMIC_2', 'COSMIC_4') # sig1=Cosmic2, sig2=Cosmic4

# Shorten the ID tag name to match the naming in the immune score data
for (i in 1:nrow(lusc_contr)) {
  long_ID <- lusc_contr$ID[i]
  chop_ID <- substr(long_ID,1,12)
  lusc_contr$ID[i] <- chop_ID
  
}
for (i in 1:nrow(luad_contr)) {
  long_ID <- luad_contr$ID[i]
  chop_ID <- substr(long_ID,1,12)
  luad_contr$ID[i] <- chop_ID
  
}

# Only pulling out the patients that are in both data sets
lusc_contr <-lusc_contr[lusc_contr$ID %in% lusc_estimate$ID, ]
luad_contr <-luad_contr[luad_contr$ID %in% luad_estimate$ID, ]
lusc_estimate <-lusc_estimate[lusc_estimate$ID %in% lusc_contr$ID, ]
luad_estimate <-luad_estimate[luad_estimate$ID %in% luad_contr$ID, ]

# Sort by ID to match all data (just in case)
lusc_contr=lusc_contr[order(lusc_contr$ID),]; 
luad_contr=luad_contr[order(luad_contr$ID),]; 

# Add the immune scores for the patients to the contr dataframe
lusc_contr$ImmuneScore <- lusc_estimate[,4]
luad_contr$ImmuneScore <- luad_estimate[,4]

# ===================== 5. Get tobacco history from TCGA ============================

# creating data frame for patient ID and number of pack-years smoked (column 36)
lusc_tobacco <- as.data.frame(lusc.tcga@clinical.data[,36])
ID <- as.data.frame(lusc.tcga@clinical.data[["Tumor_Sample_Barcode"]])
lusc_tobacco <- cbind(ID, lusc_tobacco)
names(lusc_tobacco)[1] <- c('ID') # sig1=Cosmic4, sig2=Cosmic7

luad_tobacco <- as.data.frame(luad.tcga@clinical.data[,36])
ID <- as.data.frame(luad.tcga@clinical.data[["Tumor_Sample_Barcode"]])
luad_tobacco <- cbind(ID, luad_tobacco)
names(luad_tobacco)[1] <- c('ID') # sig1=Cosmic4, sig2=Cosmic7

# Shorten the ID tag name to match the naming in the immune score data
for (i in 1:nrow(lusc_tobacco)) {
  long_ID <- lusc_tobacco$ID[i]
  chop_ID <- substr(long_ID,1,12)
  lusc_tobacco$ID[i] <- chop_ID
  
}
for (i in 1:nrow(luad_tobacco)) {
  long_ID <- luad_tobacco$ID[i]
  chop_ID <- substr(long_ID,1,12)
  luad_tobacco$ID[i] <- chop_ID
  
}

# Only pulling out the patients that are in both data sets
lusc_tobacco <-lusc_tobacco[lusc_tobacco$ID %in% lusc_estimate$ID, ]
luad_tobacco <-luad_tobacco[luad_tobacco$ID %in% luad_estimate$ID, ]

# Sort by ID to match all data (just in case)
lusc_tobacco=lusc_tobacco[order(lusc_tobacco$ID),]; 
luad_tobacco=luad_tobacco[order(luad_tobacco$ID),]; 

# ============ 6. Get contributions to each COSMIC and divide into groups ===================

# Create dataframe for count of substitutions and APOBEC per patient
lusc_count <- as.data.frame(lusc.tnm[["nmf_matrix"]])
lusc_count <- cbind(rownames(lusc_count), lusc_count)
names(lusc_count)[1] <- c('ID') # name of first column 

luad_count <- as.data.frame(luad.tnm[["nmf_matrix"]])
luad_count <- cbind(rownames(luad_count), luad_count)
names(luad_count)[1] <- c('ID') 

# Shorten the ID tag name to match the naming in the count data
for (i in 1:nrow(lusc_count)) {
  long_ID <- lusc_count$ID[i]
  chop_ID <- substr(long_ID,1,12)
  lusc_count$ID[i] <- chop_ID
  
}
for (i in 1:nrow(luad_count)) {
  long_ID <- luad_contr$ID[i]
  chop_ID <- substr(long_ID,1,12)
  luad_count$ID[i] <- chop_ID
  
}

# Only pulling out the patients that are in both data sets
lusc_count <-lusc_count[lusc_count$ID %in% lusc_estimate$ID, ]
luad_count <-luad_count[luad_count$ID %in% luad_estimate$ID, ]

# Sort by ID to match all data (just in case)
lusc_count=lusc_count[order(lusc_count$ID),]; 
luad_count=luad_count[order(luad_count$ID),];  

# calculate contributions for each COSMIC signature
library(sigminer)
lusc_allsig=as.data.frame(sig_fit(t(lusc_count[,2:97]), sig_index = 1:30,return_class = "data.table"))
luad_allsig=as.data.frame(sig_fit(t(luad_count[,2:97]), sig_index = 1:30,return_class = "data.table"))

# get COSMIC which has the maximum contribution
lusc_allsig$max=apply(lusc_allsig[,2:31],1,which.max)
luad_allsig$max=apply(luad_allsig[,2:31],1,which.max)

# Shorten the ID tag name to match the naming in the immune score data
for (i in 1:nrow(lusc_allsig)) {
  long_ID <- lusc_allsig$sample[i]
  chop_ID <- substr(long_ID,1,12)
  lusc_allsig$sample[i] <- chop_ID
  
}
for (i in 1:nrow(luad_allsig)) {
  long_ID <- luad_allsig$sample[i]
  chop_ID <- substr(long_ID,1,12)
  luad_allsig$sample[i] <- chop_ID
  
}

# Only pulling out the patients that are in both data sets
lusc_contr <-lusc_contr[lusc_contr$ID %in% lusc_allsig$sample, ]
luad_contr <-luad_contr[luad_contr$ID %in% luad_allsig$sample, ]
lusc_allsig <-lusc_allsig[lusc_allsig$sample %in% lusc_contr$ID, ]
luad_allsig <-luad_allsig[luad_allsig$sample %in% luad_contr$ID, ]

# write absolute contribution to predominant signatures into the Contr dataframe
lusc_contr$COSMIC_4 <- lusc_allsig$COSMIC_4
lusc_contr$COSMIC_7 <- lusc_allsig$COSMIC_7

luad_contr$COSMIC_2 <- luad_allsig$COSMIC_2
luad_contr$COSMIC_4 <- luad_allsig$COSMIC_4

# Divide into groups depending on maximum signature
lusc_contr$Group <- lusc_allsig$max
luad_contr$Group <- luad_allsig$max

# ========================== 7. Plot signature group histograms =======================

h <- hist(lusc_contr$Group,breaks=30,
          main="LUSC Histogram of signature groups",
          xlab="signature group", ylab="# patients")+
  geom_boxplot()

h <- hist(luad_contr$Group,breaks=30,
          main="LUAD Histogram of signature groups",
          xlab="signature group", ylab="# patients")+
  geom_boxplot()

# ===================== 8. ANOVA one-way test with 8 groups ============================

# LUAD ===========
# transform into numbers and factors
luad_contr2=luad_contr # new dataframe to do factor
luad_contr2$Group=factor(luad_contr$Group);
luad_levels=levels(luad_contr2$Group);
luad_contr2$Group <- ordered(luad_contr2$Group,levels=luad_levels)

# get summary of number of patients, mean and std per group
library(dplyr)
group_by(luad_contr2,Group) %>%
  summarise(
    count = n(),
    mean = mean(ImmuneScore, na.rm = TRUE),
    sd = sd(ImmuneScore, na.rm = TRUE)
  )

# manually filter the groups with >=5 patients
luad_analysis=filter(luad_contr2,(Group==1 | Group==2 | Group==3 | Group==4 | Group==6 | Group==13 | Group==18 | Group==24))

# ANOVA test - H0: all groups have same mean for IS
luad_anova <- aov(ImmuneScore ~ Group, data = luad_analysis)

# check conditions for ANOVA
# Levene's test - H0: variance is homogeneous
library(car)
leveneTest(ImmuneScore ~ Group, data = luad_analysis) # pvalue=0.9381

# Shapiro-Wilk test- H0: distributions are normal
luad_residuals <- residuals(object = luad_anova)
shapiro.test(x = luad_residuals ) # pvalue: 0.1271

# both conditions satisfied, so we can perform ANOVA
luad_anova <- aov(ImmuneScore ~ Group, data = luad_analysis) # pvalue=0.152
summary(luad_anova)

# pvalue>0.05 --> no significant difference

# LUSC ===================
# transform into numbers and factors
lusc_contr2=lusc_contr # new dataframe as factor
lusc_contr2$Group=factor(lusc_contr2$Group)
lusc_levels=levels(lusc_contr2$Group);
lusc_contr2$Group <- ordered(lusc_contr2$Group,levels=lusc_levels)

# get summary of number of patients, mean and std per group
group_by(lusc_contr2,Group) %>%
  summarise(
    count = n(),
    mean = mean(ImmuneScore, na.rm = TRUE),
    sd = sd(ImmuneScore, na.rm = TRUE)
  )

# manually filter the groups with >=5 patients
lusc_analysis=filter(lusc_contr2,(Group==1 | Group==2 | Group==3 | Group==4 | Group==5 | Group==13 | Group==16 | Group==24))

# ANOVA test - H0: all groups have same mean for IS
lusc_anova <- aov(ImmuneScore ~ Group, data = lusc_analysis)

# check conditions for ANOVA
# Levene's test - H0: variance is homogeneous
leveneTest(ImmuneScore ~ Group, data = lusc_analysis) # pvalue=0.34

# Shapiro-Wilk test- H0: distributions are normal
lusc_residuals <- residuals(object = lusc_anova )
shapiro.test(x = lusc_residuals ) # pvalue: 0.004569

# distributions are not normal, so we can't perform ANOVA. 
# Kruskal-Wallis (non parametric) - H0: all groups have same mean for IS
kruskal.test(ImmuneScore ~ Group, data = lusc_analysis) # pvalue: 0.01503

# as pvalue<0.05, do pairwise
# Wilcoxon rank sum test - H0: mean is the same in pairs of groups
pairwise.wilcox.test(lusc_analysis$ImmuneScore, lusc_analysis$Group,p.adjust.method = "BH")

# plot mean of 8 groups
library(ggplot2)
ggplot(lusc_analysis, aes(x = Group, y = ImmuneScore, fill = Group),
       ylab = "Immune Score", xlab = "Group",
       main="LUSC: Signature distribution") +  # Change filling color
  geom_boxplot()


# ========================== 8. ANOVA one-way test with 3 groups =========================

# set the non-predominant groups equal to 0
for (i in 1:nrow(lusc_contr)) {
  if (lusc_contr$Group[i]!=4 && lusc_contr$Group[i]!=7) {
    lusc_contr$Group[i]=0
  }
  
}

for (i in 1:nrow(luad_contr)) {
  if (luad_contr$Group[i]!=2 & luad_contr$Group[i]!=4) {
    luad_contr$Group[i]=0
  }
  
}

# LUAD ===========
# transform into numbers and factors
luad_contr$Group=factor(luad_contr$Group);
luad_levels=levels(luad_contr$Group);
luad_contr$Group <- ordered(luad_contr$Group,levels=luad_levels)

# get summary of number of patients, mean and std per group
group_by(luad_contr,Group) %>%
  summarise(
    count = n(),
    mean = mean(ImmuneScore, na.rm = TRUE),
    sd = sd(ImmuneScore, na.rm = TRUE)
  )

# ANOVA test - H0: all groups have same mean for IS
luad_anova <- aov(ImmuneScore ~ Group, data = luad_contr)

# check conditions for ANOVA
# Levene's test - H0: variance is homogeneous
leveneTest(ImmuneScore ~ Group, data = luad_contr) # pvalue=0.386

# Shapiro-Wilk test- H0: distributions are normal
luad_residuals <- residuals(object = luad_anova )
shapiro.test(x = luad_residuals ) # pvalue: 0.1436

# both conditions satisfied, so we can perform ANOVA
luad_anova <- aov(ImmuneScore ~ Group, data = luad_contr) # pvalue=0.0328
summary(luad_anova)

# pvalue<0.05 --> significant difference. do pairwise
TukeyHSD(luad_anova) # pairwise ANOVA

# plot mean of 3 groups
ggplot(luad_contr, aes(x = Group, y = ImmuneScore, fill = Group),
       ylab = "Immune Score", xlab = "Signatures",
       main="LUAD: Signature distribution") +  # Change filling color
  geom_boxplot()

# LUSC =============
# transform into numbers and factors
lusc_contr$Group=factor(lusc_contr$Group)
lusc_levels=levels(lusc_contr$Group);
lusc_contr$Group <- ordered(lusc_contr$Group,levels=lusc_levels)

# get summary of number of patients, mean and std per group
lusc_summary=group_by(lusc_contr,Group) %>%
  summarise(
    count = n(),
    mean = mean(ImmuneScore, na.rm = TRUE),
    sd = sd(ImmuneScore, na.rm = TRUE)
  )

# ANOVA test - H0: all groups have same mean for IS
lusc_anova <- aov(ImmuneScore ~ Group, data = lusc_contr)

# check conditions for ANOVA
# Levene's test - H0: variance is homogeneous
leveneTest(ImmuneScore ~ Group, data = lusc_contr) # pvalue=0.2035

# Shapiro-Wilk test- H0: distributions are normal
lusc_residuals <- residuals(object = lusc_anova )
shapiro.test(x = lusc_residuals ) # pvalue: 0.001396

# distributions are not normal, so we can't perform ANOVA. 
# Kruskal-Wallis (non parametric) - H0: all groups have same mean for IS
kruskal.test(ImmuneScore ~ Group, data = lusc_contr) # pvalue: 0.7574

# as pvalue>0.05 --> no significant difference

# ====================== 9. Tobacco analysis: correlation for IS and COSMIC_4 ========================

# match IDs
luad_contr <-luad_contr[luad_contr$ID %in% luad_tobacco$ID, ] 
luad_tobacco <-luad_tobacco[luad_tobacco$ID %in% luad_contr$ID, ]

lusc_contr <-lusc_contr[lusc_contr$ID %in% lusc_tobacco$ID, ] 
lusc_tobacco <-lusc_tobacco[lusc_tobacco$ID %in% lusc_contr$ID, ]

# Sort by ID
lusc_contr=lusc_contr[order(lusc_contr$ID),]; 
lusc_tobacco=lusc_tobacco[order(lusc_tobacco$ID),];  

luad_contr=luad_contr[order(luad_contr$ID),]; 
luad_tobacco=luad_tobacco[order(luad_tobacco$ID),];  

# turn character data from TCGA into numeric
luad_tobacco$number_pack_years_smoked=as.numeric(luad_tobacco$number_pack_years_smoked)
lusc_tobacco$number_pack_years_smoked=as.numeric(lusc_tobacco$number_pack_years_smoked)

# turn NA into 0
luad_tobacco$number_pack_years_smoked[is.na(luad_tobacco$number_pack_years_smoked)] = 0
lusc_tobacco$number_pack_years_smoked[is.na(lusc_tobacco$number_pack_years_smoked)] = 0

# plot histograms for number of pack-years smoked
h <- hist(lusc_tobacco$number_pack_years_smoked,breaks=150,
          main="LUSC Histogram of # pack-years smoked",
          xlab="# pack-years smoked", ylab="# patients")+
  geom_boxplot()

h <- hist(luad_tobacco$number_pack_years_smoked,breaks=150,
          main="LUAD Histogram of # pack-years smoked",
          xlab="# pack-years smoked", ylab="# patients")

# Correlation test - H0: variables are not correlated
# IS
cor.test(luad_contr$ImmuneScore, luad_tobacco$number_pack_years_smoked) # pvalue=0.09686
cor.test(lusc_contr$ImmuneScore, lusc_tobacco$number_pack_years_smoked) # pvalue=0.8289

# COSMIC 4 exposure
cor.test(luad_contr$COSMIC_4, luad_tobacco$number_pack_years_smoked) # pvalue=0
cor.test(lusc_contr$COSMIC_4, lusc_tobacco$number_pack_years_smoked) # pvalue=0.4042

# ======================= 10. Tobacco analysis: t-test smokers vs non-smokers =================

# filter dataframes between on-smokers (pack-years=0) and smokers
luad_nonsmoke=filter(luad_tobacco,number_pack_years_smoked==0)
luad_smoke=filter(luad_tobacco,number_pack_years_smoked>0)
luad_nonsmoke2 <-luad_contr[luad_contr$ID %in% luad_nonsmoke$ID, ] 
luad_smoke2 <-luad_contr[luad_contr$ID %in% luad_smoke$ID, ] 

lusc_nonsmoke=filter(lusc_tobacco,number_pack_years_smoked==0)
lusc_smoke=filter(lusc_tobacco,number_pack_years_smoked>0)
lusc_nonsmoke2 <-lusc_contr[lusc_contr$ID %in% lusc_nonsmoke$ID, ] 
lusc_smoke2 <-lusc_contr[lusc_contr$ID %in% lusc_smoke$ID, ] 

# Shapiro-Wilk test- H0: distributions are normal
shapiro.test(luad_nonsmoke2$ImmuneScore) # pvalue: 0.02481
shapiro.test(luad_smoke2$ImmuneScore) # pvalue: 0.06134
shapiro.test(lusc_nonsmoke2$ImmuneScore) # pvalue: 0.0241
shapiro.test(lusc_smoke2$ImmuneScore) # pvalue: 0.02681
shapiro.test(luad_nonsmoke2$COSMIC_4) # pvalue: 0
shapiro.test(luad_smoke2$COSMIC_4) # pvalue: 0
shapiro.test(lusc_nonsmoke2$COSMIC_4) # pvalue: 0
shapiro.test(lusc_smoke2$COSMIC_4) # pvalue: 0

# they are not normal, so can't do t-test

# Wilcoxon signed-rank test (non parametric) - H0: mean is the same in both groups
# IS
wilcox.test(luad_nonsmoke2$ImmuneScore, luad_smoke2$ImmuneScore, alternative = "two.sided") # pvalue=0.01768
wilcox.test(lusc_nonsmoke2$ImmuneScore, lusc_smoke2$ImmuneScore, alternative = "two.sided") # pvalue=0.5904

# COSMIC 4 exposure
wilcox.test(luad_nonsmoke2$COSMIC_4, luad_smoke2$COSMIC_4, alternative = "two.sided") # pvalue=0
wilcox.test(lusc_nonsmoke2$COSMIC_4, lusc_smoke2$COSMIC_4, alternative = "two.sided") # pvalue=0.8217

# only LUAD shows significant difference in both IS and COSMIC 4. Plot for LUAD

# assign 0 to non-smokers, 1 to smokers to divide into groups for plotting
luad_nonsmoke2$smoke=0
luad_smoke2$smoke=1

# merge into one dataframe
luad_a=merge(luad_nonsmoke2,luad_smoke2,all='TRUE')

# plot mean of smokers vs non-smokers
library("ggpubr")
# IS
ggboxplot(luad_a, x = "smoke", y = "ImmuneScore", 
          color = "smoke", palette = c("#00AFBB", "#E7B800"),
          ylab = "Immune Score", xlab = "Non smoker / Smoker",
          main='LUAD: IS in smokers vs non-smokers')
# COSMIC 4
ggboxplot(luad_a, x = "smoke", y = "COSMIC_4", 
          color = "smoke", palette = c("#00AFBB", "#E7B800"),
          ylab = "COSMIC 4", xlab = "Non smoker / Smoker",
          main='LUAD: COSMIC 4 exposure in smokers vs non-smokers')
