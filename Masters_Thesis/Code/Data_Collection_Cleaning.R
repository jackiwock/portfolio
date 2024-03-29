
########################################### DATA COLLECTION AND CLEANING #####################################################
# The following code was executed as part of my master's thesis influencers in confinement, in which I examined the influence of 
# covid lockdown conditions on pro-eating disorder twitter's leadership structure. In this portion I collected the data from 
# Twitter and cleaned and tranformed it for further analysis in the following parts:

#           # Part 1: Data collection
#            - connecting to twitter's API
#            - collect tweets with desired hashtags and timeframe
#            - transform tweets into dataframe 
#           # Part 2: Data Cleaning and Tranformation
#            - filter out non-english tweets and bot
#            - new variables: retweeted user, full text, reformat date, hashtags
#           # Part 3: Hashtag Co-occurrence network for further filtering
#            - create co-occurrence network
#            - louvain community detection
#            - filter unrelated communities
#           # Part 4: Remove duplicates



# Required Packages

library(tidyverse)         # data cleaning and transformation   
library(academictwitteR)   # access api
library(igraph)            # hashtag co-occurrence network
library(reshape2)          # transform hashtag data from long to wide format



                  ########################## Part 1: Data Collection ###########################

#### connect to twitter's API using bearer token 

set_bearer()
bearer_token <- Sys.getenv("bearer_token_here")
headers <- c(`Authorization` = sprintf('Bearer %s', bearer_token))


#### collect tweets with desired hashtags and timeframe 
#"thinspo", "edtwt", and "proana" between 2019-11-01 and 2020-09-30

Tweets <-
  get_all_tweets(
    query = c("#thinspo", "#edtw", "#proana"),
    start_tweets = "2019-11-01T00:00:00Z",
    end_tweets = "2020-09-30T00:00:00Z",
    file = "Tweets",
    data_path = "Tweets/",
    n = 1000000
  )


### Transform tweets into dataframe

Tweets_df <- bind_tweets(data_path = "Tweets/", output_format = "tidy")



                  ########################## Part 2: Data Cleaning and Tranformation ###########################
 

### filter out non-english tweets and bot

Tweets_df <- Tweets_df %>% 
  filter(lang == "en",
        user_username != "Ana___Winter")


### new variables: retweeted user, full text, reformat date, hashtags

# There are three kinds of tweets in the dataframe: retweets, orginal tweets and quote tweets. 
# Retweets begin with the syntax RT:@username, which is the user that is being retweeted. This user is extracted and place in a new column "rt_username". 
# The text column contains the text of the tweet, but retweets's text is cut off after a certain number of characters. The full text is in the sourcetweet text, but original tweets' sourcetweet text are NA. 
# To have the full text of all the tweets, a new column full_text is created, using the what is in the text column for original tweets and what is in the sourcetweet text for retweets.

Tweets_df <- Tweets_df %>%
  mutate(rt_username =  str_extract(text, "RT @\\w+"),                          #extract RT:@username and move to new column
         rt_username = str_replace_all(rt_username, "RT @", ""),                #remove RT:@ before username
         full_text =  if_else(is.na(sourcetext_text), text, sourcetweet_text))  #new column with full text of both retweets and original tweets



# create column from "created_at" with yyyy/mm/dd date format

Ana_Tweets_df$Date <- as.Date(Ana_Tweets_df$created_at, format = "%Y-%m-%d")



# extract hashtags from full_text variable 

Tweets_df <- Tweets_df %>% 
  mutate(hashtags = str_extract_all(full_text, "#(\\w+)(?=(?:\\s|,|\\.|$|#\\w))")) %>%  # regex accounts for hts followed by a space, another #, commas or periods or occurs at the end of the character string. 
  rowwise() %>% 
  mutate(hashtags = paste(hashtags, collapse = " ")) %>% 
  mutate(hashtags =tolower(hashtags)) %>%  # make lowercase
  hashtags = str_replace_all(hashtags, "#", "") #remove "#" symbol



                  ################ Part 3: Hashtag Co-occurrence Network for Further Filtering ##################

# edtwt picks up some education and tech tweets unrelated to proana. To try to tease some of these out of the data, I created
# a hashtag co-occurrence network to make connections between hashtags that appear together in tweets, then used louvain community
# detection to separate clusters of related hashtags to create hashtag topics. Topics unrelated to eating disorders were then filtered
# out of the data.


###  create co-occurrence network

# unnest the hashtags column, separating each hashtag into a new row

hashnet_all <- Tweets_df %>% separate_rows(hashtags, sep = " ")


# create dataframe with tweet_id, username and hashtags

edge_hash <- hashnet_all %>% select(tweet_id, hashtags)]


# create adjacency matrix with tweet_id as row and hashtags as columns

edgelist <- acast(edge_hash, formula = edge_hash$tweet_id ~ edge_hash$hashtags, length, value.var = "hashtags")


# remove rare hashtags that only appear once in data

edgelist <- edgelist[,apply(edgelist, MARGIN = 2, FUN = sum, na.rm = TRUE) > 1]


# transpose and multiply edgelist by itself for co-occurences

edgelist <- t(edgelist) %*% edgelist


# change diagnals to 0 to remove self loops

diag(edgelist_1) <- 0


# remove zeros from rows and columns, hashtags with no relationships to other hashtags

edgelist <- edgelist [,apply(edgelist, MARGIN = 2, FUN = sum, na.rm = TRUE)>0]
edgelist <- edgelist [apply(edgelist, MARGIN = 1, FUN = sum, na.rm = TRUE)>0,]


# make undirected weighted network

whole_cc <- graph.adjacency(edgelist, mode="undirected", weighted=TRUE)

### louvain community detection

set.seed(1254)
cc_lou <- cluster_louvain(whole_cc)


# assign community membership as node attribute

V(whole_cc)$membership <- cc_lou$membership


# calculate node degree and add as node attribute

V(whole_cc)$degree <- degree(whole_cc)


# convert network to dataframe for better viewing 

whole_df <- as_data_frame(whole_cc, what = "vertices")


### filter unrelated communities 


# create separate dataframes for each communtity

whole_mem <- whole_df %>% group_split(membership)


#take a look at each ordered by highest degree to get a sense of topics 

lapply(whole_mem, function(x){
  x %>% arrange(desc(degree))
})

#19 communties of hashtags are detected 
# community 1 has body related hashtags, some seem more pornographic in nature and may not be related to EDs 


#transform community 1 into a vector of hashtags and surround each term with \\b to search exact terms in hashtag column

comm_1_bound <- sapply(whole_mem[[1]]$name, function(word) paste0("\\b", word, "\\b"))

                       
#look at tweets with hashtags from community 1
                     
Tweets_df %>% filter(grepl(paste(comm_1_bound, collapse = "|"), hashtags)) %>% 
  select(user_username, user_description, full_text, hashtags)

                       
# most are ED related, but a handful attached to thinspo are pornographic/fetish or pertaining to erectile dyfunction
# and don't appear to be related to proed
                       
# make into vector and add \\b
                       
comm_1_fet <- c("dick","boner", "teendick", "veinydick", "youngdick", "menshealth")                      
fet_bound <- sapply(comm_1_fet, function(word) paste0("\\b", word, "\\b"))

                    
# find and capture tweet ids with these hashtags and put into vector
                    
fet_tweets <- Tweets_df %>% 
  filter(grepl(paste(fet_bound, collapse = "|"), hashtags)) %>%
 select(tweet_id, user_username, user_description, full_text, hashtags)

                    
fet_tweet_ids <- fet_tweets$tweet_id

# 10  other communities have hashtags seemingly unrelated to pro-ed twitter
# put hashtags into a list of vectors
                    
comm_to_check <- list(whole_mem[[4]]$name, whole_mem[[7]]$name, whole_mem[[11]]$name, whole_mem[[12]]$name, whole_mem[[13]]$name,
                      whole_mem[[14]]$name, whole_mem[[16]]$name, whole_mem[[17]]$name, whole_mem[[19]]$name)

                    
# //b for exact terms
                    
comm_to_check_bound <- list()
for(i in seq_along(comm_to_check)){
comm_to_check_bound[[i]] <- sapply(comm_to_check[[i]], function(word) paste0("\\b", word, "\\b"))
}

                                   
# find tweets that contain hashtag lists and check content
                                   
check_ht_dfs <- list()
for(i in seq_along(comm_to_check_bound)){
check_ht_dfs[[i]] <- Tweets_df %>% filter(grepl(paste(comm_to_check_bound[[i]], collapse = "|"), hashtags)) %>% 
  select(user_username, user_description, full_text, hashtags)
}
                            
# most of these are swept up in #edtwt and have no identifiable proed tweets in them:
# Community 7: religious 
#          11: autism research
#          12: foreign language 
#          13: bilingual education 
#          14: manchester org 
#          16: education 
#          17: doctors
#          19: nursing 

# create a vector of tweet_ids of the tweets where the unrelated hashtags were found                                
                                   
non_proana_ht <- do.call(c, lapply(check_ht_dfs[2:10], function(df) df$tweet_ids))

                                   
# community 4 is education and technology related, but some proed users get swept up with hashtags like edchat and selfcare 
# most are accompanied by known proed hashtags
                                   
# make vector of words that would exclude proed users from this dataframe
                                   
remove_school <- c("edchat", "emotions", "fitnessgoals", "selfcare", "covid", "thelockdown",
                   "thinspo", "proana", "meanspo", "anorexia", "abwtbs", "bulimia", "eatingdisorder")
                                 
rem_school_bound <- sapply(remove_school, function(word) paste0("\\b", word, "\\b"))

                           
# subset data into new dataframe
                           
School_tweets <- Tweets_df %>% 
  filter(grepl(paste(comm_to_check_bound, collapse = "|"), hashtags)) %>%
  filter(!grepl(paste(rem_school_bound, collapse = "|"), hashtags)) %>% 
  select(tweet_id, user_username, user_description, full_text, hashtags)

                           
# make vector of all tweet_ids to be removed
                           
rem_tweet_ids <- c(School_tweets$tweet_id, non_proana_ht, fet_tweet_ids)
#remove from df
                           
clean_tweets <- Tweets_df %>% filter(!tweet_id %in% rem_tweet_ids)

                           
                           
                  ############################ Part 4: Remove Duplicates ##################################
                           
# found duplicates in the dataframe. Seems to be disparity around the like_count and the user_followers count. Remove duplicates and
# keep tweet with higher user_followers_count

final_ana <- clean_tweets %>%
  group_by(tweet_id) %>%
  slice(which.max(user_followers_count)) %>%
  ungroup()

# Original count 73,896. Now 71,001

                           

