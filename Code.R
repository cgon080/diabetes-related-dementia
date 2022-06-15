setwd("~/DRD")

dir()

dat1 = read.csv("CMDHBML.csv", header = TRUE, sep = ",", dec = ".")
summary(dat1)
names(dat1)

dat2 = dat1[,c(2,6:7,29,129:133,166,168,169,171,52,57,38:40,5)]
names(dat2)

for (i in 1:19) {
  dat2[,i] <- replace(dat2[,i],dat2[,i]=='NI',NA)
}

dat2$date_referral <- as.Date(dat2$date_referral)
dat2$hba1c_before_date <- as.Date(dat2$hba1c_before_date)
dat2$hba1c_after_date <- as.Date(dat2$hba1c_after_date)
dat2$imaging_before_date <- as.Date(dat2$imaging_before_date)
dat2$imaging_after_date <- as.Date(dat2$imaging_after_date)

dat2$dif_ref_bef = dat2$date_referral-dat2$hba1c_before_date
dat2$dif_ref_aft = dat2$hba1c_after_date-dat2$date_referral

dat2$dif_ref_bef[is.na(dat2$dif_ref_bef)] <- 9999
dat2$dif_ref_aft[is.na(dat2$dif_ref_aft)] <- 9999

dat2$hba1c = rep(0, 3949)
for (i in 1:3949){
  if (dat2$dif_ref_bef[i]<dat2$dif_ref_aft[i]){
    dat2$hba1c[i] = dat2$hba1c_before[i]
  } else if (dat2$dif_ref_bef[i]>dat2$dif_ref_aft[i]){
    dat2$hba1c[i] = dat2$hba1c_after[i]
  } else {
    dat2$hba1c[i] = dat2$hba1c_before[i]
  } 
}

dat2$hba1c[3938] <- 23
dat2$hba1c <- as.numeric(dat2$hba1c)
summary(dat2$hba1c)

# 1st filter people with hba1c>50

dat3 = subset(dat2, hba1c>=50)
dat3$rudas_memory <- as.numeric(dat3$rudas_memory)
dat3$rudas_total <- as.numeric(dat3$rudas_total)
dat3$RUDAS_mem_sc = dat3$rudas_memory/8

dat3$ace_total <- as.numeric(dat3$ace_total)
dat3$ace_memory <- as.numeric(dat3$ace_memory)
dat3$ace_attention <- as.numeric(dat3$ace_attention)
dat3$ACE_mem_sc = dat3$ace_memory/26

dat3$mem_att_rat = dat3$ace_memory/dat3$ace_attention

# With results about CT-Scans

# 1) Select the closest result to the Referral 
dat3$dif_bef_ima = dat3$date_referral-dat3$imaging_before_date
dat3$dif_aft_ima = dat3$imaging_after_date-dat3$date_referral

dat3$dif_bef_ima[is.na(dat3$dif_bef_ima)] <- 9999
dat3$dif_aft_ima[is.na(dat3$dif_aft_ima)] <- 9999

dat3$CTScan_Result = rep(0, 895)
for (i in 1:895){
  if (dat3$dif_bef_ima[i]<dat3$dif_aft_ima[i]){
    dat3$CTScan_Result[i] = dat3$imaging_before_result[i]
  } else if (dat3$dif_bef_ima[i]>dat3$dif_aft_ima[i]){
    dat3$CTScan_Result[i] = dat3$imaging_after_result[i]
  } else {
    dat3$CTScan_Result[i] = dat3$imaging_before_result[i]
  } 
}

library(dplyr)
library(stringr)

dat3$CTScan_Result = tolower(dat3$CTScan_Result)

SCScan_Code = dat3 %>% mutate(I1 = +str_detect(CTScan_Result,str_c(c('medial atrophy','temporal atrophy','hippocampal atrophy'), collapse = '|')),
                I2 = +str_detect(CTScan_Result,str_c(c('vascular lesions','vascular lesion','small vessel disease','ischaemic changes','ischaemic change'), collapse = '|')),
                I3 = +str_detect(CTScan_Result,str_c(c('generalised atrophy'), collapse = '|')),
                I5 = +str_detect(CTScan_Result,str_c(c('the patient did not present'),collapse = '|')))

dat4 = subset(SCScan_Code, !is.na(CTScan_Result) & dementia_diagnosis==1)
dat4$I4 = rep(0, nrow(dat4))
for (i in 1:nrow(dat4)) {
  if (rowSums(dat4[i,29:32])==0){
    dat4$I4[i]=1
  } else {
    dat4$I4[i]=0
  }
}

dat4$I2[dat4$I2==1]<-2
dat4$I3[dat4$I3==1]<-3
dat4$I4[dat4$I4==1]<-4
dat4$I5[dat4$I5==1]<-5

dat4$IndexCT = rowSums(dat4[,c(31,33)])

dat5 = subset(dat4, IndexCT!=0 & I2==0)
dat6 = subset(dat5, !is.na(rudas_total)|!is.na(rudas_memory)|!is.na(ace_total)|!is.na(ace_attention)|!is.na(ace_memory))

dat7 = subset(dat6, !(is.na(RUDAS_mem_sc) & is.na(ACE_mem_sc)))
dat7$Cognitive = rowSums(dat7[,23:24],na.rm = T)
dat6['1028', 'Cognitive']<-0.30769231
dat6['1061', 'Cognitive']<-0.65384615
dat6['3423', 'Cognitive']<-0.03846154

dat8 = subset(dat7, Cognitive>0.5)

#Table with results

dat8 %>% group_by(ethnicity1) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))
dat8 %>% group_by(gender) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))

#------------------
#For control group
#------------------

dat9 = subset(dat2, hba1c<50)

dat9$rudas_memory <- as.numeric(dat9$rudas_memory)
dat9$rudas_total <- as.numeric(dat9$rudas_total)
dat9$RUDAS_mem_sc = dat9$rudas_memory/8

dat9$ace_total <- as.numeric(dat9$ace_total)
dat9$ace_memory <- as.numeric(dat9$ace_memory)
dat9$ace_attention <- as.numeric(dat9$ace_attention)
dat9$ACE_mem_sc = dat9$ace_memory/26

dat9$mem_att_rat = dat9$ace_memory/dat9$ace_attention

# With results about CT-Scans

# 1) Select the closest result to the Referral 
dat9$dif_bef_ima = dat9$date_referral-dat9$imaging_before_date
dat9$dif_aft_ima = dat9$imaging_after_date-dat9$date_referral

dat9$dif_bef_ima[is.na(dat9$dif_bef_ima)] <- 9999
dat9$dif_aft_ima[is.na(dat9$dif_aft_ima)] <- 9999

dat9$CTScan_Result = rep(0, 3035)
for (i in 1:3035){
  if (dat9$dif_bef_ima[i]<dat9$dif_aft_ima[i]){
    dat9$CTScan_Result[i] = dat9$imaging_before_result[i]
  } else if (dat9$dif_bef_ima[i]>dat9$dif_aft_ima[i]){
    dat9$CTScan_Result[i] = dat9$imaging_after_result[i]
  } else {
    dat9$CTScan_Result[i] = dat9$imaging_before_result[i]
  } 
}

library(dplyr)
library(stringr)

dat9$CTScan_Result = tolower(dat9$CTScan_Result)

dat10 = dat9 %>% mutate(I1 = +str_detect(CTScan_Result,str_c(c('medial atrophy','temporal atrophy','hippocampal atrophy'), collapse = '|')),
                              I2 = +str_detect(CTScan_Result,str_c(c('vascular lesions','vascular lesion','small vessel disease','ischaemic changes','ischaemic change'), collapse = '|')),
                              I3 = +str_detect(CTScan_Result,str_c(c('generalised atrophy'), collapse = '|')),
                              I5 = +str_detect(CTScan_Result,str_c(c('the patient did not present'),collapse = '|')))

dat11 = subset(dat10, !is.na(CTScan_Result) & dementia_diagnosis==1)
dat11$I4 = rep(0, nrow(dat11))
for (i in 1:nrow(dat11)) {
  if (rowSums(dat11[i,29:32])==0){
    dat11$I4[i]=1
  } else {
    dat11$I4[i]=0
  }
}

dat11$I2[dat11$I2==1]<-2
dat11$I3[dat11$I3==1]<-3
dat11$I4[dat11$I4==1]<-4
dat11$I5[dat11$I5==1]<-5

dat11$IndexCT = rowSums(dat11[,c(31,33)])

dat12 = subset(dat11, IndexCT!=0 & I2==0)
dat13 = subset(dat12, !is.na(rudas_total)|!is.na(rudas_memory)|!is.na(ace_total)|!is.na(ace_attention)|!is.na(ace_memory))

dat14 = subset(dat13, !(is.na(RUDAS_mem_sc) & is.na(ACE_mem_sc)))
dat14.1 = subset(dat13, !is.na(RUDAS_mem_sc) & !is.na(ACE_mem_sc))
dat14$Cognitive = rowSums(dat14[,23:24],na.rm = T)
dat14['173', 'Cognitive']<-0.1153846
dat14['322', 'Cognitive']<-0.4615385
dat14['505', 'Cognitive']<-0.5000000
dat14['1985', 'Cognitive']<-0.7307692
dat14['2037', 'Cognitive']<-0.6153846
dat14['2088', 'Cognitive']<-0.3076923
dat14['2580', 'Cognitive']<-0.3461538
dat14['2629', 'Cognitive']<-0.5000000
dat14['3072', 'Cognitive']<-0.2692308

dat15 = subset(dat14, Cognitive>0.5)

dat15 %>% group_by(ethnicity1) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))
dat15 %>% group_by(gender) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))

set.seed(111)
dat15$Rand = runif(86)

dat16 = dat15[order(dat15$Rand),]
dat16.1 = dat16[1:34,]

dat16.1 %>% group_by(ethnicity1) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))
dat16.1 %>% group_by(gender) %>% summarise(mean(age_at_referral),mean(hba1c),mean(rudas_total, na.rm=T),mean(rudas_memory, na.rm=T),mean(ace_total, na.rm=T),mean(ace_memory, na.rm=T))

#------------------
#Exercise
#------------------

dat8$Group = 1
dat15$Group = 0
dat15.1 = dat15[,-36]

dat17 = rbind(dat8, dat15.1)

fit1 = glm(as.factor(Group) ~ as.factor(ethnicity1) + as.factor(gender) + age_at_referral, 
           family = binomial(link = "logit"), data = dat17)
summary(fit1)

#------------------
#Text analysis
#------------------
library(tm)
library(dplyr)
library(radiant.data)
library(wordcloud)
library(SnowballC)
library(RColorBrewer)

docs <- Corpus(VectorSource(dat8$CTScan_Result))
inspect(docs)

# 1. Stripping any extra white space:
docs <- tm_map(docs, stripWhitespace)

# 2. Transforming everything to lowercase
docs <- tm_map(docs, content_transformer(tolower))

# 3. Removing numbers 
docs <- tm_map(docs, removeNumbers)

# 4. Removing punctuation
docs <- tm_map(docs, removePunctuation)

# 5. Removing stop words
docs <- tm_map(docs, removeWords, stopwords("english"))
docs <- tm_map(docs, removeWords, c("findings", "indication", "technique", "conclusion", 
                                    "radiologists", "radiologist", "consultant"))

dtm <- DocumentTermMatrix(docs)

sums <- as.data.frame(colSums(as.matrix(dtm)))
sums <- rownames_to_column(sums) 
colnames(sums) <- c("term", "count")
sums <- arrange(sums, desc(count))
head <- sums[1:100,]

wordcloud(words = head$term, freq = head$count, min.freq = 10,
          max.words=100, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

library(wordcloud2)
wordcloud2(head)

#----------------------
#ACM
#----------------------

dat18 = dat17
dat18$Age = cut(dat17$age_at_referral, 4)
dat18$gender = as.factor(dat18$gender)
dat18$ethnicity1 = as.factor(dat18$ethnicity1)
dat18$Group = as.factor(dat18$Group)

dat19 = dat18[,c(2,3,36,37)]

library(FactoClass)

acm1 = FactoClass(dat19, dudi.acm)
plotFactoClass(acm1)
acm1$carac.cate
