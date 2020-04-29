library(plyr)
library(dplyr)	# MAKE SURE THIS LOADS AFTER plyr, or errors can arise
library(grid)
library(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)
library(gridExtra)
library(gtools)
library(ggthemes)
library(gsubfn)
library(sjmisc)
library(timeDate)
library(data.table)

## --------------------------------- ##
## read in data
## --------------------------------- ##
setwd("../data/raw data/")
# read in and tidy data for original conditions involving item repetition during training
d1 <- read.csv("exp1-TS-CNTRL-run1-Nov4_2015-48subj.csv", sep = "\t")
d2 <-read.csv("exp1-TS-CNTRL-run2-Nov9_2015-52subj.csv", sep = "\t")
d3 <- read.csv("exp1-MT-ST-Nov11_2015-80subj.csv", sep = "\t")
d4 <- read.csv("exp1-balance-Dec1_2015-62subj.csv", sep = "\t")
d5 <- read.csv("exp1-talker37-March9_2016-40subj.csv", sep = ",")
d6 <- read.csv("new-talker-balance-3-10-2016-22subj.csv", sep = "\t")
d7 <- read.csv("feb2017-new-talker-balance.csv")


d.rep <- rbind(d1, d2, d3, d4, d5, d6, d7)
d.rep <- within(d.rep, {
  ItemRepetition <- "yes"
})

# read in and tidy data for no repetition conditions
d.norep.1 <- read.csv("exp2-NOREP-MT-ST-Dec1_2015-80subj.csv", sep = "\t")
d.norep.2 <- read.csv("exp2-NOREP-MT-ST-CNTRL-March9_2016-64subj.csv", sep = ",")

d.norep <- rbind.fill(d.norep.1, d.norep.2)
d.norep <- within(d.norep, {
  Answer.Submit = NULL #Not really sure why this column appeared
  ItemRepetition <- "no"
})

d.new.exp2 <- rbind.fill(read.csv("replication/replicate-april12.csv", sep = "\t"),
                read.csv("replication/replicate-fillup.csv", sep = "\t"),
                read.csv("replication/replicate-fillup2.csv", sep = "\t"), 
                read.csv("replication/replicate-fillup3.csv", sep = "\t"),
                read.csv("replication/replicate-april15-fillup.csv", sep = "\t"),
                read.csv("replication/replicate-MT2.csv", sep = "\t"),
                read.csv("replication/replicate-april16.csv", sep = "\t"))

d.new.exp2 <- within(d.new.exp2, {
  ItemRepetition <- "Exp2Rep"
})


# bind all data together
d <- rbind.fill(d.rep, d.norep, d.new.exp2)

# read in sentence keyword info
key.df <- read.delim("../../exp1_BB08_replication/lists/usable_sentences_96.txt", header = T)

key.df <- within(key.df, {
	SentenceID <- gsub("ENG_", "", SentenceID)
	Sentence <- as.character(Sentence)
})

## --------------------------------- ##
## helper functions
## --------------------------------- ##

#Get nth smallest values
getNFirst <- function(v, n) {
  sort(unique(v))[1:n]
}

# Multi-gsub function: takes a list of patterns and an equal-length list of replacement characters
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

# To code accuracy based on the number of correctly transcribed keywords
num.keywords.correct <- function(transcription, keywords){
	# remove punctuation and ensure all characters are lower case
	keys <- gsub("[[:punct:]]", "", tolower(as.character(keywords)))
	trans <- gsub("[[:punct:]]", "", tolower(as.character(transcription)))
	
	# split keyword list and transcription into word arrays
	keys.list <- unlist(strsplit(keys, " "))	
	trans.list <- unlist(strsplit(trans, " "))
	
	# count number of keywords that occur in the transcription array
	return(sum(!is.na(pmatch(keys.list, trans.list))))
}


## --------------------------------- ##
## variables for later
## --------------------------------- ##

chars <- c("<", ">", '"', ",", ";", "|", "'")
char_codes <- c("&lt&", "&gt&", "&quot&", "&comma&", "&semicol&", "&pipe&", "&apos&")

Talker_map <- list('1'='ALL_021_M_CMN', 
                '2'='ALL_016_M_CMN', 
                '3'='ALL_032_M_CMN',
                '4'='ALL_043_M_CMN',
                '5'='ALL_035_M_CMN',
                '6'='ALL_037_M_CMN')

####################################
#Survey
#####################################
d.survey = d %>%
  dplyr::select(contains("Answer"), workerid, assignmentsubmittime, ItemRepetition, -c(contains("Resp")))

d.survey <- within(d.survey, {
  TemporaryID = rleid(workerid)
  WorkerID = workerid
  workerid = NULL
  Answer.Talker = NULL
  Answer.userAgent = NULL
  Answer.set = NULL
})

####################################
#Test Block
#####################################
d.test = ldply(strsplit(unlist(strsplit(as.character(d$Answer.testResp), ";")), "\\|"))
names(d.test) = c("PartOfExp", "Trial", "ListPosition", "Filename", "Sentence", "Transcription", "Condition", "TrainingTestSet", "ListNum", "TestTalkerNum", "WorkerID")
d.test <- within(d.test, {
  Trial = as.integer(Trial)
  TemporaryID = rleid(WorkerID)
  ListPosition = as.integer(ListPosition)
  TestTalkerID = as.integer(ifelse(Condition %in% c("5AmEng_1Acc", "5Acc_Diff1Acc", 
                                        "NoRep_5AmEng_1Acc", "1Acc_Same1Acc"), 
                                    TestTalkerNum, gsub(".*,", "", TestTalkerNum)))
  TestTalkerName = as.character(Talker_map[TestTalkerID])
  TestTalkerNameShort = as.integer(gsub("[^0-9]", "", TestTalkerName))
  TestTalkerNum = NULL
  PresentationBlock = NA
})

#####################################
#Training Block
#####################################
d.training = ldply(strsplit(unlist(strsplit(as.character(d$Answer.trainingResp), ";")), "\\|"))
names(d.training) = c("PartOfExp", "Trial", "ListPosition", "Filename", "Sentence", "Transcription", "Condition", "TrainingTestSet", "ListNum", "TestTalkerNum", "WorkerID")

d.training <- within(d.training, {
  Trial = as.integer(Trial)
  ListPosition = as.integer(ListPosition)
  TemporaryID = rleid(WorkerID)
  TestTalkerID = as.integer(ifelse(Condition %in% c("5AmEng_1Acc", "5Acc_Diff1Acc", 
                                                     "NoRep_5AmEng_1Acc", "1Acc_Same1Acc"), 
                                    TestTalkerNum, gsub(".*,", "", TestTalkerNum)))
  TestTalkerName = as.character(Talker_map[TestTalkerID])
  TestTalkerNameShort = as.integer(gsub("[^0-9]", "", TestTalkerName))
  TestTalkerNum = NULL
  PresentationBlock <- ifelse(Trial %in% c(0:15), 1, 
                              ifelse(Trial %in% c(16:31), 2,
                                     ifelse(Trial %in% c(32:47), 3,
                                            ifelse(Trial %in% c(48:63), 4,
                                                   ifelse(Trial %in% (64:79), 5, "")))))
  
})

#Join together training and test
d.both <- rbind(d.training, d.test)

#Join together test/training and survey
d.all <- merge(d.both, d.survey, by = c("WorkerID", "TemporaryID"))

d.all <- d.all %>%
  arrange(assignmentsubmittime) %>%
  mutate(
    OrderOfSubmission = as.integer(as.factor(rank(assignmentsubmittime, ties.method = "min")))
)

## The following code executes much faster when using dplyr (as abov)
# d.test$OrderOfSubmission = as.integer(as.factor(rank(d.test$AssignmentSubmitTime,ties.method = "min")))

#Here's a sanity check -- the dates should be in order
unique(d.all[with(d.all, order(OrderOfSubmission)),]$assignmentsubmittime)


d.all.m <- merge(d.all, subset(key.df, select = -c(SentenceSet)), by = c("Sentence"), sort = FALSE)

d.all.m <- within(d.all.m, {
	SentenceID <- gsub(".*(HT[1-2].*)", "\\1", Filename)
	CurrentTalkerID <- gsub("talker_(.*)/.*", "\\1", Filename)
  CurrentTalkerIDShort <- gsub("[^0-9]", "", CurrentTalkerID)
	# adjust trial numbering to start at 1 rather than 0
	Trial <- Trial + 1

	AudioType = factor(Answer.audio_type, levels = c("", "computer Talkers", "external", "in-ear", "over-ear"))
	AudioType2_WoreHeadphones = as.factor(ifelse(AudioType %in% c("in-ear", "over-ear"), "headphones", "other"))
	Answer.audio_type = NULL
	
  AudioQual = factor(Answer.audio_qual, levels = c("", "okay", "good", "excellent", "professional"))
  Answer.audio_qual = NULL
  
	AccentFamFreq <- factor(Answer.accent_familiarity_frequency, levels = c("", "all_time", "day", "week", "month", "year", "never"))
	Answer.accent_familiarity_frequency <- NULL
	
  Condition2 <- as.character(Condition)
    Condition2 [Condition == "1Acc_Same1Acc"] <- "Talker-specific"
    Condition2 [Condition %in% c("1Acc_Diff1Acc", "NoRep_1Acc_Diff1Acc")] <- "Single talker"
    Condition2 [Condition %in% c("5Acc_Diff1Acc", "NoRep_5Acc_Diff1Acc")] <- "Multi-talker"
  	Condition2 [Condition %in% c("5AmEng_1Acc", "NoRep_5AmEng_1Acc")] <- "Control"
	Condition2 <- factor(Condition2, levels = c("Talker-specific", "Single talker", "Multi-talker", "Control"))
  
	Condition.long <- NA
		Condition.long [Condition == "1Acc_Same1Acc"] <- "Talker-specific\ntraining"
	  Condition.long [Condition %in% c("1Acc_Diff1Acc", "NoRep_1Acc_Diff1Acc")] <- "Single talker\ntraining"
	  Condition.long [Condition %in% c("5Acc_Diff1Acc", "NoRep_5Acc_Diff1Acc")] <- "Multi-talker\ntraining"
		Condition.long [Condition %in% c("5AmEng_1Acc", "NoRep_5AmEng_1Acc")] <- "Control\n(training = 5 AmEng)"
	Condition.long <- as.factor(Condition.long)
  Condition.long <- factor(Condition.long, levels = as.character(levels(Condition.long)[c(1,3,2,4)]))

  ItemRepetition <- factor(ItemRepetition, levels = c("yes", "no", "Exp2Rep"))
  
	# remove capitalization and convert special punctuation characters (char_codes) to actual punctuation
	Transcription <- tolower(mgsub(char_codes, chars, Transcription))

	NumKeywordsCorrect <- mapply(num.keywords.correct, Transcription, Keywords)
	
	PropKeywordsCorrect <- NumKeywordsCorrect / NumKeywords	
})

###############################################################################################################

## --------------------------------- ##
# Effect of audio variables on accuracy (independent of condition)
## --------------------------------- ##

#xtabs(~AudioType + AudioQual, subset(d.all.m, !duplicated(WorkerID)))
#xtabs(~ AccentFamFreq, subset(d.all.m, !duplicated(WorkerID)))

# p_audiotype <- d.all.m %>%
#   mutate(
#     PartOfExp = factor(PartOfExp, levels = c("training", "test")),
#     AudioType = sjmisc::word_wrap(as.character(AudioType), wrap = 6)
#     ) %>%
#   group_by(WorkerID, AudioType, PartOfExp) %>%
#   dplyr::summarise(PropKeywordsCorrect = mean(PropKeywordsCorrect)) %>%
#   ggplot(aes(x = AudioType, y = PropKeywordsCorrect)) +
#     geom_dotplot(binaxis = "y", stackdir = "center", alpha = .2, dotsize = 0.3) +
#     stat_summary(fun.y = "mean", geom = "point", size = 3) + 
#     stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .25) +
#     facet_wrap(~ PartOfExp) +
#     ggtitle("Effect of Audio Type on Accuracy\n(independent of condition)") +
#     theme_bw() +
#     theme(axis.text.x = element_text(size = 9))
#ggsave(p_audiotype, file = "../figures/accuracy_by_audiotype.pdf")



# p_audioqual <- d.all.m %>%
#   mutate(
#     PartOfExp = factor(PartOfExp, levels = c("training", "test"))
#   ) %>%
#   filter(AudioType %in% c("in-ear", "over-ear")) %>%
#   group_by(WorkerID, AudioQual, PartOfExp) %>%
#   dplyr::summarise(PropKeywordsCorrect = mean(PropKeywordsCorrect)) %>%
#   ggplot(aes(x = AudioQual, y = PropKeywordsCorrect)) +
#     geom_dotplot(binaxis = "y", stackdir = "center", alpha = .2, dotsize = 0.3) +
#     stat_summary(fun.y = "mean", geom = "point", size = 3) + 
#     stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .25) +
#     facet_wrap(~ PartOfExp) +
#     ggtitle("Effect of Audio Quality on Accuracy\n(among participants using in-ear or over-ear headphones,\nindependent of condition)") +
#     theme_bw()
#ggsave(p_audioqual, file = "../figures/accuracy_by_audioquality.pdf", width = 11)



# p_accentfam <- d.all.m %>%
#   mutate(
#     PartOfExp = factor(PartOfExp, levels = c("training", "test"))
#   ) %>%
#   group_by(WorkerID, AccentFamFreq, PartOfExp) %>%
#   dplyr::summarise(PropKeywordsCorrect = mean(PropKeywordsCorrect)) %>%
#   ggplot(aes(x = AccentFamFreq, y = PropKeywordsCorrect)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", alpha = .2, dotsize = 0.3) +
#   stat_summary(fun.y = "mean", geom = "point", size = 3) + 
#   stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .25) +
#   facet_wrap(~ PartOfExp) +
#   xlab("Accent familiarity\n(frequency of exposure to accent)")
#   ggtitle("Effect of Accent Familiarity on Accuracy\n(independent of condition)") +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 9))
#ggsave(p_audiotype, file = "../figures/accuracy_by_accentfam.pdf")


## --------------------------------- ##
#Exclusion Criteria: Puts rejected worker ids in data frame called rejectedWorkers
## --------------------------------- ##

rejectedWorkers = data.frame()
#Repeat Worker
workers = data.frame(table(d$workerid))
repeatWorkers = workers[workers$Freq > 1,]

#TODO add repeated workers to rejectWorkers (none now so it doesn't matter)

#Exclude people who speak an East Asian Language
d.lgbackground = data.frame(cbind(as.character(d$Answer.language_background),
                                  as.character(d$Answer.language_background_free), 
                                  as.character(d$Answer.accent_familiarity_frequency),
                                  as.character(d$Answer.accent_familiarity_place),
                                  as.character(d$Answer.audio_type),
                                  as.character(d$Answer.audio_qual),
                                  as.character(d$workerid)))
names(d.lgbackground) = c("LgBack", "LgBackFree","AccentFreq", "AccentPlace", "AudioType", "AudioQual", "WorkerId")

#Check if they checked "I Talker Mandarin or another east asian"
if (nrow(subset(d.lgbackground, LgBack == "chinese"))) { #Note: chinese refers to any east asian lg here
  rejectedWorkers = rbind(rejectedWorkers,
                          data.frame(cbind(as.character(subset(d.lgbackground, LgBack == "chinese")$WorkerId), "SpeaksEastAsian")))
}

#Some people's free responses don't agree with their multi. Maybe we need to make this clearer in the questions.
if (nrow(subset(d.lgbackground, grepl("Japan|Chinese|China|Mandarin|Cantonese|Korean", 
                                      d.lgbackground$LgBackFree,
                                      ignore.case=TRUE)))) { 
  rejectedWorkers = rbind(rejectedWorkers,
                          data.frame(cbind(as.character(
                            subset(d.lgbackground, 
                                   grepl("Japan|Chinese|China|Mandarin|Cantonese|Korean", 
                                         d.lgbackground$LgBackFree,
                                         ignore.case=TRUE))$WorkerId), "SpeaksEastAsianFree")))
}
#You might have noticed that the last talker you heard had an accent. 
#How often do you hear talkers with a similar accent, whether in person or in movies/TV shows etc.?
if (nrow(subset(d.lgbackground, AccentFreq %in% c("all_time")))) {
  rejectedWorkers = rbind(rejectedWorkers,
                          data.frame(cbind(as.character(subset(d.lgbackground, AccentFreq %in% c("all_time"))$WorkerId), "HearsAccentAllTime")))
}
#If you somewhat regularly hear listeners with a similar accent, 
#please tell us in what context you encounter these talkers (select all that apply)
if (nrow(subset(d.lgbackground, grepl("in_my_family|among_my_close_friends", AccentPlace)))) {
  rejectedWorkers = rbind(rejectedWorkers,
                          data.frame(cbind(as.character(subset(d.lgbackground, grepl("in_my_family|among_my_close_friends", AccentPlace))$WorkerId), "FamilyCloseFriendsHaveAccent")))
}
#What kind of audio equipment did you use for the experiment?
# if (nrow(subset(d.lgbackground, !(AudioType %in% c("in-ear", "over-ear"))))) {
#   rejectedWorkers = rbind(rejectedWorkers,
#                           data.frame(cbind(as.character(subset(d.lgbackground, !(AudioType %in% c("in-ear", "over-ear")))$WorkerId), "NotWearingHeadphones")))
# }


#Add in condition and list info
names(rejectedWorkers) = c("WorkerID", "ReasonRejected")
rejectedWorkers = aggregate(rejectedWorkers$ReasonRejected, by=list(WorkerID=rejectedWorkers$WorkerID), FUN=paste)

rejectedWorkers$ReasonRejected = as.character(rejectedWorkers$x)
rejectedWorkers$x = NULL
d.all.m$RejectWorker = ifelse(d.all.m$WorkerID %in% rejectedWorkers$WorkerID, TRUE, FALSE)
d.all.m = left_join(d.all.m, rejectedWorkers)

# Identify the test and training talkers that each worker heard
d.training <- subset(d.all.m, PartOfExp == "training")
training_talkers_byWorker <- d.training %>%
  dplyr::mutate(cb = paste(WorkerID, CurrentTalkerID, sep = "-")) %>%
  filter(!duplicated(cb)) %>%
  dplyr::select(WorkerID, CurrentTalkerID) %>%
  arrange(WorkerID)# %>%  
 # group_by(WorkerId)%>%
#  summarise(TrainingTalkers = paste(TalkerID, collapse = ", "))
training_talkers_byWorker = aggregate(training_talkers_byWorker$CurrentTalkerID, list(training_talkers_byWorker$WorkerID), paste, collapse=", ")
names(training_talkers_byWorker) <- c("WorkerID", "TrainingTalkerID")

# Add test and training talker IDs to rejectedWorkers data frame
d.all.m <- left_join(d.all.m, training_talkers_byWorker)

# how many workers were excluded from each condition and for each reason
addmargins(xtabs(~ ReasonRejected + Condition, subset(d.all.m, !duplicated(WorkerID))))

# Data frame of usable subjects
d.test.m.usable = subset(d.all.m, !WorkerID %in% rejectedWorkers$WorkerID & PartOfExp == "test")

## --------------------------------- ##
#Identify workers for each experiment
## --------------------------------- ##

# Experiment 1: Talker-specific vs. Control. 
# Get only the first 8 participants in the talker-specific and control conditions for each of the 6 test talkers


exp1_workers <- d.test.m.usable %>%
  dplyr::select(WorkerID, Condition2, TestTalkerID, ItemRepetition, OrderOfSubmission) %>%
  distinct(WorkerID, .keep_all = TRUE) %>%
  filter(Condition2 %in% c("Talker-specific", "Control"),
         ItemRepetition == "yes") %>%
  group_by(Condition2, TestTalkerID) %>%
  arrange(OrderOfSubmission) %>%
  do(head(., n = 8))
  
xtabs(~ Condition2 + TestTalkerID , exp1_workers)  


# Experiment 2: All conditions vs. Control (with item repetition during training)
exp2_workers <- d.test.m.usable %>%
  dplyr::select(WorkerID, Condition2, TestTalkerID, ItemRepetition, OrderOfSubmission, ListNum, TrainingTestSet) %>%
  distinct(WorkerID, .keep_all = TRUE) %>%
  filter(TestTalkerID%in% c(3,4,5,6),
         ItemRepetition == "yes")
  
xtabs(~ TrainingTestSet + ListNum + TestTalkerID + Condition2, exp2_workers)

# Experiment2 Replication 
exp_replication_workers <- d.test.m.usable %>%
  filter(ItemRepetition == "Exp2Rep") %>%
  #filter(ItemRepetition == "Exp2Rep" & AudioType2_WoreHeadphones != "other" & Condition2 %in% c("Control", "Talker-specific")) %>%
  #filter(ItemRepetition == "Exp2Rep" & PartOfExp == "test" & AudioType2_WoreHeadphones == "other" & Condition2 == "Control") %>%  
  filter(ItemRepetition == "Exp2Rep" & PartOfExp == "test") %>%
  dplyr::group_by(WorkerID, Condition2, TestTalkerID, TrainingTalkerID, ItemRepetition, TrainingTestSet, ListNum, AudioType2_WoreHeadphones) %>%
  dplyr::summarise(MeanProp = mean(PropKeywordsCorrect)) %>%
  distinct(WorkerID, .keep_all = TRUE)

xtabs(~ Condition2 + AudioType2_WoreHeadphones, exp_replication_workers)


xtabs(~ TrainingTestSet + ListNum + Condition2 + TestTalkerID +  TrainingTalkerID, exp_replication_workers)


# Experiment 3: Cross-talker generalization conditions (single talker and multi-talker) vs. Control, with no item repetition during training
exp3_workers <- d.test.m.usable %>%
  dplyr::select(WorkerID, Condition2, TestTalkerNameShort, ItemRepetition) %>%
  distinct(WorkerID, .keep_all = TRUE) %>%
  filter(TestTalkerNameShort %in% c(32, 43, 37, 35),
         ItemRepetition == "no")

xtabs(~ Condition2 + TestTalkerNameShort + ItemRepetition, exp3_workers)

d.all.m %>%
  filter(ItemRepetition == "Exp2Rep" & RejectWorker == TRUE & PartOfExp == "test") %>%
  distinct(Condition2, WorkerID) %>%
  group_by(Condition2, ReasonRejected) %>%
  dplyr::summarise(a = n())

## --------------------------------- ##
## save data frames to RData file
## --------------------------------- ##

d.all.m.final <- within(d.all.m, {
  Exp1_subj <- ifelse(WorkerID %in% exp1_workers$WorkerID, "yes", "no")
  Exp2_subj <- ifelse(WorkerID %in% exp2_workers$WorkerID, "yes", "no")
  Exp3_subj <- ifelse(WorkerID %in% exp3_workers$WorkerID, "yes", "no")
  ExpRep_subj <- ifelse(WorkerID %in% exp_replication_workers$WorkerID, "yes", "no")
})

library(digest)
#Remove unnecessary columns
d.all.m.final$TemporaryID = NULL
d.all.m.final$Answer.Submit = NULL
# TFJ: I commented out these lines so that order of submission remains available information
# d.all.m.final$OrderOfSubmission = NULL
# d.all.m.final$assignmentsubmittime = NULL
d.all.m.final$WorkerID = sapply(d.all.m.final$WorkerID, digest, algo="md5")
saveRDS(d.all.m.final, file = "../replication_data_all-February-2019.RDS")
#save(d.training.m.usable, file = "exp1-2-3_data_training-March10_2016.RData")
#save(d.test.m.usable, file = "exp1-2-3_data_test-March10_2016.RData")
