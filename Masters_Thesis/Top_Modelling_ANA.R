
##### TOPIC MODEL #######
# this code describes the topic model used to adress RQ2

library(ggplot2)
library(data.table)
library(slam)         #needed for perplexity calc
library(quanteda)     #pre-process text and create corpus
library(topicmodels)
library(ldatuning) #calculate griffiths metric for number of topics
library(tidytext)
library(cowplot)
library(tidyverse)
library(kableExtra)
library(gt)
library(gridExtra)
library(grid)
library(janitor) #clean names


#load df with just retweets
load("month_RT_tm.Rda")


##### PRE-PROCESSING CORPUS #######
Tweet_text_df <- month_RT_tm %>% filter(sourcetweet_type == "retweeted") %>%  
  select(sourcetweet_id, full_text, month) %>% distinct
# I did a lot of the preprocessing manually to deal with twitter specific words, symbols and syntax

# remove RT: infront of retweeted username
Tweet_text_df$full_text <- gsub("(RT|via)((?:\\b\\W*@\\w+)+)", "", Tweet_text_df$full_text)

#remove mentions
Tweet_text_df$full_text = gsub("@\\w+", "", Tweet_text_df$full_text)

#remove Hashtags
Tweet_text_df$full_text = gsub("#\\w+","", Tweet_text_df$full_text)

#remove urls
Tweet_text_df$full_text = gsub("http.+ |http.+$", " ", Tweet_text_df$full_text)

#when processing in a previous test, its discovered that quanteda doesnt seem to handle contractions very well. 
#Lematize doesn't appear to help it's something to with the apostrophe. 
Tweet_text_df$full_text = gsub("don’t", "do not", Tweet_text_df$full_text)
Tweet_text_df$full_text = gsub("can’t", "cannot", Tweet_text_df$full_text)
Tweet_text_df$full_text = gsub("isn’t", "is not", Tweet_text_df$full_text)
Tweet_text_df$full_text = gsub("it’s", "it is", Tweet_text_df$full_text)
Tweet_text_df$full_text = gsub("you’re", "you are", Tweet_text_df$full_text) 
Tweet_text_df$full_text = gsub("I’m", "I am", Tweet_text_df$full_text) 

#other pre-processing
Tweet_text_df$full_text = gsub("[[:punct:]]", " ", Tweet_text_df$full_text)  # Remove punctuation
Tweet_text_df$full_text = gsub("[[ |\t]]{2,}", " ", Tweet_text_df$full_text)  # Remove tabs
Tweet_text_df$full_text = gsub("\\b[a-z]\\b", "", Tweet_text_df$full_text) #remove single letters
Tweet_text_df$full_text = gsub("^ ", "", Tweet_text_df$full_text)  # Leading blanks
Tweet_text_df$full_text= gsub(" $", "", Tweet_text_df$full_text)  # Lagging blanks
Tweet_text_df$full_text = gsub(" +", " ", Tweet_text_df$full_text) # General spaces
Tweet_text_df$full_text = gsub("[[:digit:]]", "", Tweet_text_df$full_text) #numbers
Tweet_text_df$full_text = tolower(Tweet_text_df$full_text) #lowercase
Tweet_text_df$full_text = gsub("\n", " ", Tweet_text_df$full_text) # remove \n for new line

# previous run showed a tendency of users to repeat phrases multiple times that was biasing the data
Tweet_text_df$full_text = gsub("\\b(\\w+)\\b(?=.*\\b\\1\\b)", "", Tweet_text_df$full_text, perl = TRUE) 

#remove duplicates 
Tweet_text_df <- Tweet_text_df[!duplicated(Tweet_text_df$full_text),]


##### MAKE CORPUS #####

#corpus with metadata
RT_corpus <- corpus(x = Tweet_text_df,
                        text_field = "full_text",
                        meta = list( "month"),
                        docid_field = "sourcetweet_id")

# tokenize
RT_tokens <- quanteda::tokens(x = RT_corpus) 

# Remove stopwords
RT_tokens <- quanteda::tokens_remove(x = RT_tokens, stopwords("en"), padding = FALSE)

# lemmatize
RT_tokens <- quanteda::tokens_replace(RT_tokens, 
                                          pattern = lexicon::hash_lemmas$token, 
                                          replacement = lexicon::hash_lemmas$lemma)

#limit to works over 3 characters to remove noise and meaningless words
RT_tokens <- tokens_select(RT_tokens, min_nchar=3L, max_nchar=79L)

####create document term matrix
dtm <- dfm(x = RT_tokens)
#Document-feature matrix of: 7,119 documents, 5,498 features (99.86% sparse) and 1 docvars.

#trim doc matrix according to Maier's 99% and 0.05%
ndocs <- length(RT_corpus)

# ignore overly sparse terms (appearing in less than 1% of the documents)
minDocFreq <- ndocs * 0.005
# ignore overly common terms (appearing in more than 80% of the documents)
maxDocFreq <- ndocs * 0.99

dtm <- dfm_trim(dtm, min_termfreq = minDocFreq, max_termfreq = maxDocFreq)
#Document-feature matrix of: 7,129 documents, 340 features (98.29% sparse) and 33 docvars.

raw_sum <- apply(dtm,1,FUN=sum) #sum by raw each raw of the table
dtm_new <- dtm[raw_sum > 0, ] # remove empties
#Document-feature matrix of: 6,827 documents, 354 features (98.12% sparse) and 33 docvars.

# find ideal number of topics
# first try griffiths2004 metric

edtops_griff <- FindTopicsNumber(dtm_new, topics = seq(from = 2, to = 60, by = 1), metrics = c("Griffiths2004"), method = "Gibbs", control = list(seed = 22),
                            verbose = TRUE)

# plot with ggplot
ggplot(edtops_griff, aes(topics, Griffiths2004)) +
  geom_point(color = "#66A182") +
  geom_line(color = "#66A182") +
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,1)) +
  theme(axis.text = element_text(size = 5)) +
  ggtitle("Griffiths 2004 Metric for Number of Topics")

#find ideal K with perplexity plot 

#make test and training set
rowids <- 1:nrow(dtm_new)
trids <- sample(x = rowids, size = 0.8*length(rowids), replace = F)
tstids <- rowids[!rowids %in% trids]
tr_dt<- dtm_new[trids,]
tst_dt <- dtm_new[-trids,]
dim(tst_dt)
#test Ks
ks <- c(3,10, 20,25, 26, 27,30,40,41,50,60,100)
perp <- c()
for(i in 1:length(ks)){
  
  # Train
  current_tm <- LDA(x = tr_dt,
                        k = ks[i], 
                        method="Gibbs",
                        control=list(alpha = 0.01,
                                     iter = 1000, 
                                     seed = 2222, 
                                     verbose = FALSE))
  # Compute perplexity
  perp[i] <- perplexity(object = current_tm, 
                            newdata = as.simple_triplet_matrix(tst_dt))
  print(perp)
}

#make df for ggplot
perp_plot <- data.frame(ks = ks,
                        perp = perp)

# perplexity scores
ggplot(perp_plot, aes(x=ks, y = perp)) +
  geom_point(color= "#4682B4") +
  geom_line(color = "#4682B4") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10) ) +
  ggtitle("Perplexity Scores") +
  ylab("perplexity score") +
  xlab("number of topics")

### LDA MODEL #####

lda_26 <-  LDA(x = dtm_new,
                   k = 26, 
                   method="Gibbs",
                   control=list(alpha = 0.01,
                                iter = 1000, 
                                seed = 7777, 
                                verbose = FALSE))

#betas 0.001 and 0.01 were tested using delta paramter in other iterations. 
#Didn't effect results so left at default 0.01
#alphas 0.05, 0.001, 0.1 were test with little effect on the results

#get top 15 terms for each topic
print(topicmodels:: get_terms(lda_26, 15))

#get term-topic distribution of terms
lda_26_posterior <- topicmodels::posterior(object = lda_26)
#create column for terms
top_dist_over_words_26 <- lda_26_posterior$terms
#make dataframe using datatable
top_dist_over_words_dt_26 <- data.table(topic=1:26, 
                                            top_dist_over_words_26)

#reshape data to columns topic, term and beta value
top_dist_over_words_dt_26 <- melt.data.table(top_dist_over_words_dt_26,
                                                 id.vars = 'topic')
#arrange in descending order
top_dist_over_words_dt_26 <- top_dist_over_words_dt_26[order(value,decreasing = T)]
#subset top 15
top15per_topic_26 <- top_dist_over_words_dt_26 %>% 
  group_by(topic) %>% 
  slice_max(order_by = value, n = 15)

#extract topic proportions from each document
doc_top_prop_26 <- lda_26_posterior$topic

#create data table with columns sourcetweet_id to identify the tweet and its topic proportions
doc_top_prop_26 <- data.table(sourcetweet_id = lda_26@documents,
                                  doc_top_prop_26)

#rename columns by their topic number
setnames(x = doc_top_prop_26, 
         old = 2:ncol(doc_top_prop_26),
         new = paste0("Topic", 1:26))

#plot term-topic probabilties

#separate into groups for better viewing
#topics 1 to 13 in facet plot
top15per_topic_26_plot %>%
  filter(topic %in% 1:13) %>% 
  group_by(topic) %>%
  slice_max(value, n = 15) %>% 
  ungroup() %>%
  arrange(topic, -value) %>% 
  mutate(variable = reorder_within(variable, value, topic)) %>%
  ggplot(aes(value, variable, fill = factor(label)))+
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free_y") +
  scale_y_reordered() + 
  theme_minimal() +
  xlab("beta") +
  ylab("terms") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme(plot.margin = unit(c(2,1.8,2,1.8), "cm"))


#Topics are given names based on their terms and corresponding tweets. I used this code, going 
#through each topic individually, here's topic 1 as an example:

#make empty column for topic labels
top15per_topic_26$label <- NA
#get top 10 tweets with highest topic proportions
t1_ids <- doc_top_prop_26[order(Topic1, decreasing = T)][,.(sourcetweet_id, Topic1)][1:10,]

#inspect top tweets
Tweet_text_df %>% filter(sourcetweet_id %in% t1_ids$sourcetweet_id) %>% select(sourcetweet_id, full_text)
#discomfort/restriction/food/drink

#assign topic label 
top15per_topic_26 <- top15per_topic_26 %>%
  mutate(label = ifelse(topic == 1, "Discomfort/Restriction/Water", label))

#add rest of labels
top15per_topic_26 <- top15per_topic_26 %>% 
  mutate(label = case_when(topic == 1 ~ "Discomfort/Restriction/Water",
                          topic == 2 ~ "Aspiration/Clothing",
                           topic == 3 ~ "Sweetspo/Goals",
                           topic == 4 ~ "Burning Calories/Diets",
                           topic == 5 ~ "Group Chat/General",
                           topic == 6 ~ "Community Maintenance",
                           topic == 7 ~ "Meanspo",
                           topic == 8 ~ "Proana Coach/Support",
                           topic == 9 ~ "Undefined/Mixed 2",
                           topic == 10 ~ "Meanspo Severe",
                           topic == 11 ~ "New User/Edtwt",
                           topic == 12 ~ "Binge/Restriction/Tips",
                           topic == 13 ~ "Thinspo/Clothing",
                           topic == 14 ~"Aspiration/Body Comparison",
                           topic == 15 ~"New User/Thinspo/Meanspo",
                           topic == 16 ~ "Thoughts/Emotions",
                           topic == 17 ~ "New User/Friends/Mutuals",
                           topic == 18 ~ "Thinspo/Representation",
                           topic == 19 ~ "Feel/Hunger/Fasting",
                           topic == 20 ~ "Food/Low Calorie",
                           topic == 21 ~ "Meanspo Harsh",
                           topic == 22 ~ "Group Chat/Pro-ED Topics",
                           topic == 23 ~ "Undefined/Mixed Topic",
                           topic == 24 ~ "Weight loss",
                           topic == 25 ~ "Thinspo/Motivation",
                           topic == 26 ~ "Aspiration/Body Parts"))

# factor into themes
top15per_topic_26_fact <- top15per_topic_26 %>%
  mutate(label = factor(label, levels = c(
    "Thoughts/Emotions",
    "Sweetspo/Goals",
    "Feel/Hunger/Fasting",
    "Discomfort/Restriction/Water",
    "Binge/Restriction/Tips",
    "Burning Calories/Diets",
    "Food/Low Calorie",
    "Weight loss",
    "Aspiration/Clothing",
    "Aspiration/Body Comparison",
    "Aspiration/Body Parts",
    "Thinspo/Clothing",
    "Thinspo/Motivation",
    "Thinspo/Representation",
    "Community Maintenance",
    "Proana Coach/Support",
    "Group Chat/General",
    "Group Chat/Pro-ED Topics",
    "New User/Edtwt",
    "New User/Thinspo/Meanspo",
    "New User/Friends/Mutuals",
    "Meanspo",
    "Meanspo Harsh",
    "Meanspo Severe",
    "Undefined/Mixed Topic",
    "Undefined/Mixed 2")))

top15per_topic_26_fact %>% 
  filter(label %in% head(levels(label), 13)) %>% 
  group_by(label) %>%
  slice_max(value, n = 15) %>% 
  ungroup() %>%
  arrange(label, -value) %>% 
  mutate(variable = reorder_within(variable, value, label)) %>%
  ggplot(aes(value, variable, fill = factor(label)))+
  geom_col(show.legend = FALSE) +
  facet_wrap(~ label, scales = "free_y") +
  scale_y_reordered() + 
  theme_minimal() +
  xlab("beta") +
  ylab("terms") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme(plot.margin = unit(c(2,1.8,2,1.8), "cm"))

top15per_topic_26_fact %>% 
  filter(label %in% head(levels(label), 13)) %>% 
  group_by(label) %>%
  slice_max(value, n = 15) %>% 
  ungroup() %>%
  arrange(label, -value) %>% 
  mutate(variable = reorder_within(variable, value, label)) %>%
  ggplot(aes(value, variable, fill = factor(label)))+
  geom_col(show.legend = FALSE) +
  facet_wrap(~ label, scales = "free_y") +
  scale_y_reordered() + 
  theme_minimal() +
  xlab("beta") +
  ylab("terms") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme(plot.margin = unit(c(2,1.8,2,1.8), "cm"))

#### Mean topic proportions ######
### Make a dataframe of topic proportions for each tweet

#make vector of label names
lab_vec <- top15per_topic_26 %>% distinct(label)

#make vector sourcetweet_id and each label that will be columns in  as columns
lab_vec <- c("sourcetweet_id", lab_vec$label)

# change column names to topic names
colnames(doc_top_prop_26) <- lab_vec

#calculate topic means
topic_means <- colMeans(doc_top_prop_26[, c(2:26)])
#descending order
topic_means <- topic_means[order(topic_means, decreasing = T)]
#put into a table
topic_means_dt <- data.table(topic = names(topic_means),
                             mean = topic_means)

#plot means for source tweets
ggplot(topic_means_dt, aes(x = mean, y = reorder(topic,mean))) +
  geom_point(size =3, color = "#C2b280") +
  ylab("topics") + 
  theme(plot.margin = unit(c(2,1,2,1), "cm"))

#expand proportions to whole Retweet dataframe
Full_RT_Topics_df <- left_join(month_RT_tm, doc_top_prop_26, by = "sourcetweet_id")
RT_doc_top_prop <- Full_RT_Topics_df[,36:61]
RT_doc_top_prop <- RT_doc_top_prop %>% drop_na()
topic_means_RT <- colMeans(RT_doc_top_prop[, c(1:26)])
topic_means_RT <- topic_means[order(topic_means_RT, decreasing = T)]
topic_means_dt_RT <- data.table(topic = names(topic_means_RT),
                             mean = topic_means_RT)

#plot
ggplot(topic_means_dt, aes(x = mean, y = reorder(topic,mean))) +
  geom_point(size =3, color = "#8d6b94") + 
  ylab("topics")+
  theme(plot.margin = unit(c(2,1,2,1), "cm"))

month_props <- Full_RT_Topics_df %>% select(35:61) %>% drop_na()

#### Measure Topic Proportions by month ####
# aggregate by months
topic_proportion_per_month <- aggregate(month_props[,-1], by = list(month = month_props$month), mean)

#reshape
month_props_plot_df <- melt(topic_proportion_per_month, id.vars = "month")


#colors for better identification
my_colors <- c(
  "#C2b280", "#d9d9ce",
           "#5b3256","violetred4", "orchid3","#8d6b94","#D8BFD8",
           "#FFA07A",
           "#66A182","#6B8E23","#c3cfa0",
           "deepskyblue4","#4682B4","#87CEEB",
           "#c08081","#EAC9c1", "darkred", "firebrick2",
           "#C71585", "#FF1493", "palevioletred",
           "goldenrod", "gold4", "darkgoldenrod1",
           "slategray", "grey80")
           

#factor so colors are in the right order
           month_props_plot_df <- month_props_plot_df %>% 
             mutate(variable = factor(variable, c("Thoughts/Emotions",
                                                  "Sweetspo/Goals",
                                                  "Feel/Hunger/Fasting",
                                                  "Discomfort/Restriction/Water",
                                                  "Binge/Restriction/Tips",
                                                  "Burning Calories/Diets",
                                                  "Food/Low Calorie",
                                                  "Weight loss",
                                                  "Aspiration/Clothing",
                                                  "Aspiration/Body Comparison",
                                                  "Aspiration/Body Parts",
                                                  "Thinspo/Clothing",
                                                  "Thinspo/Motivation",
                                                  "Thinspo/Representation",
                                                  "Community Maintenance",
                                                  "Proana Coach/Support",
                                                  "Group Chat/General",
                                                  "Group Chat/Pro-ED Topics",
                                                  "New User/Edtwt",
                                                  "New User/Thinspo/Meanspo",
                                                  "New User/Friends/Mutuals",
                                                  "Meanspo",
                                                  "Meanspo Harsh",
                                                  "Meanspo Severe",
                                                  "Undefined/Mixed Topic",
                                                  "Undefined/Mixed 2")))

#plot
 ggplot(month_props_plot_df, aes(x=month, y=value, fill=variable)) + 
             geom_bar(stat = "identity") + ylab("proportion") + 
             scale_fill_manual(values = my_colors, name = "topics") + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1),
                   legend.key.size = unit(.5, 'cm'), #change legend key size
                   legend.key.height = unit(.5, 'cm'), #change legend key height
                   legend.key.width = unit(.5, 'cm'), #change legend key width
                   legend.title = element_text(size=8), #change legend title font size
                   legend.text = element_text(size=5)) 
 
 ## Top 50 proportions 
 
 #subset top 50
 #load top degree user list
 load("top_50_users.Rdata")
 #these are top 20
 
 #break full RT df with topic proportions into month dfs
 months_tm_list <- Full_RT_Topics_df %>% as.data.frame() %>% 
   group_split(month) 
 
 #subset top 50 in rt_username
 top_50_tm_full_df_list <- list()
 for(i in seq_along(months_tm_list)){
   top_50_tm_full_df_list[[i]] <- months_tm_list[[i]] %>% 
     filter(rt_username %in% top_50_users[[i]]$name)
 }
 
 #bind list of dataframes into one dataframe 
 top_50_tm_full_df <- do.call(rbind, top_50_tm_full_df_list)

 #select month and topics columns
 month_props_50 <- top_50_tm_full_df %>% select(35:61) %>% drop_na()
 
 #aggregate by month
 topic_proportion_per_month_50 <- aggregate(month_props_50[,-1], by = list(month = month_props_50$month), mean)
 
 #reshape
 month_props_plot_df_50 <- melt(topic_proportion_per_month_50, id.vars = "month")

 #factor again for this df so colors are right
  month_props_plot_df_50 <- month_props_plot_df_50 %>% 
   mutate(variable = factor(variable, c("Thoughts/Emotions",
                                        "Sweetspo/Goals",
                                        "Feel/Hunger/Fasting",
                                        "Discomfort/Restriction/Water",
                                        "Binge/Restriction/Tips",
                                        "Burning Calories/Diets",
                                        "Food/Low Calorie",
                                        "Weight loss",
                                        "Aspiration/Clothing",
                                        "Aspiration/Body Comparison",
                                        "Aspiration/Body Parts",
                                        "Thinspo/Clothing",
                                        "Thinspo/Motivation",
                                        "Thinspo/Representation",
                                        "Community Maintenance",
                                        "Proana Coach/Support",
                                        "Group Chat/General",
                                        "Group Chat/Pro-ED Topics",
                                        "New User/Edtwt",
                                        "New User/Thinspo/Meanspo",
                                        "New User/Friends/Mutuals",
                                        "Meanspo",
                                        "Meanspo Harsh",
                                        "Meanspo Severe",
                                        "Undefined/Mixed Topic",
                                        "Undefined/Mixed 2")))
 
 #plot
 ggplot(month_props_plot_df_50, aes(x=month, y=value, fill=variable)) + 
   geom_bar(stat = "identity") + ylab("proportion") +
   scale_fill_manual(values = my_colors, name = "topics") + 
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         legend.key.size = unit(.5, 'cm'), #change legend key size
         legend.key.height = unit(.5, 'cm'), #change legend key height
         legend.key.width = unit(.5, 'cm'), #change legend key width
         legend.title = element_text(size=8), #change legend title font size
         legend.text = element_text(size=5)) 
 
### Raw tweet count for selected topics 
 
 #create threshold for proportion above which it is counted as being mostly about
 #that topic
 threshold <- 0.5
 
 
 # new column indicating if the tweet is associated with the topic
 columns_to_process <- 2:27
 
 # Iterate over each column
 for (col_index in columns_to_process) {
   col_name <- colnames(retweet_props)[col_index]
   
   # calcultate of tweets are mostly about that topic assigning them 1 is they are 0 if they aren't
   retweet_props[[paste0("Topic",col_name)]] <- ifelse(retweet_props[[col_name]] > threshold, 1, 0)
 }
 
 #do same procedure for top 50 nodes 
 rt_props_50 <- top_50_tm_full_df %>% select(35:61) %>% drop_na()
 for (col_index in columns_to_process) {
   col_name <- colnames(rt_props_50)[col_index]
   
   rt_props_50[[paste0("Topic",col_name, "50")]] <- ifelse(rt_props_50[[col_name]] > threshold, 1, 0)
 }
 
 # plot weight 
 
 ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_weight_loss, color = "Overall", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_weight_loss50, color = "Top 50", linetype = "Top 50"), group = 1) +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   scale_color_manual(name = "Lines", values = c("Overall" = "#FFA07A", "Top 50" = "#FFA07A")) + ggtitle("Weight loss") +
   ylab("tweet count") +
   xlab("month") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 # plot thinspo themes 
 
 hin_topic_trad <-  ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_thinspo_clothing, color = "Thinspo/Clothing", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_thinspo_clothing50, color = "Thinspo/Clothing", linetype = "Top 50"), group = 1) +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   scale_color_manual(name = "Topic", values =c("Thinspo/Clothing" = "midnightblue",
                                                "Thinspo/Representation" = "#87CEEB")) + 
   ggtitle("Monthly Count of Tweets with > 50% proportion of Thinspo Themes") +
   geom_line(aes(y = topic_thinspo_representation, color = "Thinspo/Representation", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_thinspo_representation50, color = "Thinspo/Representation", linetype = "Top 50"), group = 1)+
   ylab("tweet count") +
   xlab("month") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text = element_text(size = 5),
         legend.key.size = unit(.2, 'cm'), 
         legend.key.height = unit(.2, 'cm'), 
         legend.key.width = unit(.2, 'cm'), 
         legend.title = element_text(size=5), 
         legend.text = element_text(size=3),
         plot.title = element_text(size = 6),
         axis.title = element_text(size = 7))
 
 #plot aspiration/body parts
 
 ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_aspiration_body_parts, color = "Aspiration/Body Parts", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_aspiration_body_parts50, color = "Aspiration/Body Parts", linetype = "Top 50"), group = 1) +
   ggtitle("Aspirational/Body Parts") +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   scale_color_manual(name = "Topic", values =c(
     "Aspiration/Body Parts" = "olivedrab3")) +
   ylab("tweet count") +
   xlab("month") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text = element_text(size = 5),
         legend.key.size = unit(.2, 'cm'), 
         legend.key.height = unit(.2, 'cm'), 
         legend.key.width = unit(.2, 'cm'), 
         legend.title = element_text(size=5), 
         legend.text = element_text(size=3),
         plot.title = element_text(size = 6),
         axis.title = element_text(size = 7))
 
 #plot group chat themese
 ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_group_chat_general, color = "Group Chat/General", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_group_chat_general50, color = "Group Chat/General", linetype = "Top 50"), group = 1) +
   geom_line(aes(y = topic_group_chat_pro_ed_topics, color = "Group Chat/Pro-ED Topics", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_group_chat_pro_ed_topics50, color = "Group Chat/Pro-ED Topics", linetype = "Top 50"), group = 1) +
   ggtitle("Group Chat Themes") +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   scale_color_manual(name = "Topic", values =c("Group Chat/General" = "darkred",
                                                "Group Chat/Pro-ED Topics" = "firebrick2")) +
   ylab("tweet count") +
   xlab("month") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text = element_text(size = 5),
         legend.key.size = unit(.2, 'cm'), 
         legend.key.height = unit(.2, 'cm'), 
         legend.key.width = unit(.2, 'cm'), 
         legend.title = element_text(size=5), 
         legend.text = element_text(size=3),
         plot.title = element_text(size = 6),
         axis.title = element_text(size = 7))

 
 #New User Themes
 
 ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_new_user_edtwt, color = "New User/Edtwt", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_new_user_edtwt50, color = "New User/Edtwt", linetype = "Top 50"), group = 1) +
   geom_line(aes(y = topic_new_user_thinspo_meanspo, color = "New User/Thinspo/Meanspo", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_new_user_thinspo_meanspo50, color = "New User/Thinspo/Meanspo", linetype = "Top 50"), group = 1) +
   ggtitle("New Users Themes") +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   scale_color_manual(name = "Topic", values =c("New User/Edtwt" = "maroon4",
                                                "New User/Thinspo/Meanspo" = "#FF1493")) +
   ylab("tweet count") +
   xlab("month") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text = element_text(size = 5),
         legend.key.size = unit(.2, 'cm'), #change legend key size
         legend.key.height = unit(.2, 'cm'), #change legend key height
         legend.key.width = unit(.2, 'cm'), #change legend key width
         legend.title = element_text(size=5), #change legend title font size
         legend.text = element_text(size=3),
         plot.title = element_text(size = 6),
         axis.title = element_text(size = 7))
 
 ggplot(all_topics_counts, aes(x = month)) +
   geom_line(aes(y = topic_new_user_friends_mutuals, color = "New User/Friends/Mutuals", linetype = "Overall"), group = 1) +
   geom_line(aes(y = topic_new_user_friends_mutuals50, color = "New User/Friends/Mutuals", linetype = "Top 50"), group = 1) +
   ylab("tweet count") +
   xlab("month") +
   scale_linetype_manual(name = "Lines", values = c("Overall" = "solid", "Top 50" = "dotted")) +
   ggtitle("New User/Friends/Mutuals Topic") +
   scale_color_manual(name = "Topic", values =c("New User/Friends/Mutuals" = "lightpink2"))+
   scale_y_continuous(limits = c(0,600), breaks = seq(0,600,200)) +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text = element_text(size = 5),
         legend.key.size = unit(.2, 'cm'), #change legend key size
         legend.key.height = unit(.2, 'cm'), #change legend key height
         legend.key.width = unit(.2, 'cm'), #change legend key width
         legend.title = element_text(size=5), #change legend title font size
         legend.text = element_text(size=3),
         plot.title = element_text(size = 6),
         axis.title = element_text(size = 7))