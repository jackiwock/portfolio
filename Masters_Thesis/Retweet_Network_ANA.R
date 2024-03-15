
###### RETWEET NETWORK #######

# This code handles the first three sections in methods concerning RQ1: 
  # forming retweet network
  # fitting powerlaw lines to degree distributions
  # User Turnover 

library(tidyverse)
library(igraph)
library(ggplot2)
library(modelsummary) 
library(poweRlaw)
library(kableExtra)
library(Hmisc) #spearman's correlation


# load data 
load("final_ana.Rdata")
# create variable that classes tweets by month
final_ana <- final_ana %>% mutate(month = format(Date, "%Y-%m"))

####### PREPARE DATA FOR RT NTW #########

# Separate data into months
months_df_list <- final_ana %>% as_tibble() %>% 
  group_split(month) 

# make df with just retweets 
month_RT_ntw_df <- lapply(months_df_list, function(x){
  x %>% filter(sourcetweet_type == "retweeted")
})

#drop self loops
month_RT_ntw_df <- lapply(month_RT_ntw_df, function(x){
  x %>% filter(user_username != rt_username)
})

######## MAKE RT NTW #########
# Created edgelist
ALL_RT <- lapply(month_RT_ntw_df, function(x){
  x %>% 
    select(user_username, rt_username)
})
#make weighted
RT_edgelist <- lapply(ALL_RT, function(x){
  x %>% 
    group_by(user_username, rt_username) %>% 
    summarise(weight = sum(n()))
})

#create networks from edgelist
RT_ntws <- lapply(RT_edgelist, graph_from_data_frame)

#remove self loops (just to be sure)
RT_ntws <- lapply(RT_ntws, function(x){
  igraph::simplify(x, remove.loop = T)
})
# remove isolates 
RT_ntws <- lapply(RT_ntws, function(x){
  x <- delete.vertices(x, which(degree(x) == 0))
})

########## CALCULATE NTW STATS FUNCTION ########
calculate_network_stats <- function(network) {
  nodes <- vcount(network)
  ties <- ecount(network)
  density <- graph.density(network)
  reciprocity <- reciprocity(network)
  transitivity_global <- transitivity(network, type = "global")
  clustering_coef <- transitivity(network, type = "average")
  ave_path_length <- mean_distance(network)
  
  # Create a data frame with the calculated statistics
  stats_df <- data.frame(
    Nodes = nodes,
    Ties = ties,
    Density = density,
    Reciprocity = reciprocity,
    Transitivity = transitivity_global,
    Clustering_Coef = clustering_coef,
    Mean_Distance = ave_path_length
  )
  
  return(stats_df)
}


##### CALCULATE NTW STATS MONTHLY #######
network_stats <- lapply(RT_ntws, calculate_network_stats)

##### CALCULATE NTW STATS OVERALL #######
average_stats <- do.call(rbind,network_stats)
#save(average_stats, file = "average_net_stats.Rdata")
#summarize data
Months <- c("2019-11", "2019-12", "2020-01", "2020-02", "2020-03", "2020-04", "2020-05", "2020-06", "2020-07", "2020-08", "2020-09")

Months_names <- c("November", "December", "January", "February", "March", "April", "May", "June", "July", "August", "September")

#Add month do each df 

#convert networks in to dataframes
RT_add_month_df <- lapply(RT_ntws, as_data_frame, what = "vertices")

#add month variable
for(i in seq_along(RT_ntws)){
  RT_add_month_df[[i]] <- RT_add_month_df[[i]] %>% mutate(month = Months[i])
}

##### CREATE FULL RETWEET NETWORK DF ########
RT_ntw_full <- do.call(rbind, RT_add_month_df)

##### ADD MONTH JOINED ####### 
#Add month_joined variable by grouping usernames and assigning variable based on first value that appears in df
RT_ntw_full <- RT_ntw_full %>% 
  arrange(month) %>%  # Make sure the dataframe is sorted by date.
  group_by(name) %>% 
  mutate(month_joined = format(first(month))) %>%
  ungroup()


# Add as a node attribute to the network
for (i in seq_along(RT_ntws)) {
  # Extract the node names in the current monthly network
  network_nodes <- V(RT_ntws[[i]])$name
  
  # Match node names from the full data frame with the node names in the network
  matched_indices <- match(network_nodes, RT_ntw_full$name)
  
  # Filter out NA indices
  matched_indices <- matched_indices[!is.na(matched_indices)]
  
  # Extract the month_joined values from the full data frame
  attributes_to_add <- RT_ntw_full$month_joined[matched_indices]
  
  # Set the month_joined attribute in the igraph network
  V(RT_ntws[[i]])$month_joined <- attributes_to_add
}

##### ANONYMIZE USER ID ##########

#make a list of all users in networks
unique_users <- unique(RT_ntw_full$name)

#make a function that generates 4 random letters
random_letters <- replicate(length(unique_users), paste(sample(letters, 4), collapse = ''))

# Create a user_id column from combining month_joined with random letters
RT_ntw_full$user_id <- paste(RT_ntw_full$month_joined, random_letters[match(RT_ntw_full$name, unique_users)], sep = '-')
# Rename the columns if needed

#add as node attribute in similar way to month_joined
for (i in seq_along(RT_ntws)) {
  #Extract the node names in the current monthly network
  network_nodes <- V(RT_ntws[[i]])$name
  
  # Match node names from the full data frame with the node names in the network
  matched_indices <- match(network_nodes, RT_ntw_full$name)
  
  # Filter out NA indices
  matched_indices <- matched_indices[!is.na(matched_indices)]
  
  # Extract the month_joined values from the full data frame
  node_id_to_add <- RT_ntw_full$user_id[matched_indices]
  
  # Set the month_joined attribute in the igraph network
  V(RT_ntws[[i]])$user_id <- node_id_to_add
}

###### DEGREE CENTALITY ########
degrees <- lapply(RT_ntws, degree, mode = "in")
for(i in seq_along(RT_ntws)){
  V(RT_ntws[[i]])$degree <- degrees[[i]]
}
####### MONTHLY RT NTW DATAFRAMES WITH ALL VARIABLES #########
RT_ntws_df <-  lapply(RT_ntws, as_data_frame, what = "vertices")


degree_stats <- lapply(seq_along(RT_ntws_df), function(i) {
  data.frame(
    Month = Months_names[i],
    Mean = mean(RT_ntws_df[[i]]$degree),
    SD = sd(RT_ntws_df[[i]]$degree),
    Min = min(RT_ntws_df[[i]]$degree),
    Max = max(RT_ntws_df[[i]]$degree)
  )
})
degree_stats <- do.call(rbind, degree_stats)

######## FITTING POWERLAW LINES ############

est_pl_RT <- list()
pl_results_list <- list()
pl_RT <- list()
for (i in seq_along(RT_ntws_df)) {
  # Create a power-law object for the degree distribution
  pl_RT[[i]] <- displ$new(RT_ntws_df[[i]]$degree)
  
  # Estimate xmin
  est_pl_RT[[i]] <- estimate_xmin(pl_RT[[i]])
  pl_RT[[i]]$setXmin(est_pl_RT[[i]])
  
  # Add the power-law object to the results list
  pl_results_list[[i]] <- pl_RT[[i]]
}


# make df for each month with results 
est_pl_RT_df <- Map(function(x, month) {
  data.frame(
    month = month,  # Use the corresponding month value
    gof = x$gof,
    xmin = x$xmin,
    alpha = x$pars
  )
}, est_pl_RT, Months)

#bind them all together 
est_pl_RT_full <- do.call(rbind, est_pl_RT_df)



# plot in ggplot for better display 
pl_plots_gg <- list()
pl_lines <- list()
for (i in seq_along(pl_results_list)) {
  pl_plots_gg[[i]] <- plot(pl_results_list[[i]])
  pl_lines[[i]] <- lines(pl_results_list[[i]], col = 2)
}

###### PLOTS FOR TOP 20 and top 50 DEGREE USERS ##########
top_20_users <- list()
top_50_users <- list()
plots_deg <- list()
for (i in seq_along(RT_ntws_df)) {
  # Convert name to factor with levels sorted in descending order of degree
  RT_ntws_df[[i]]$user_id <- fct_reorder(RT_ntws_df[[i]]$user_id, RT_ntws_df[[i]]$degree, .desc = TRUE)
  
  #save top degree users
  top_20_users[[i]] <- RT_ntws_df[[i]] %>% 
    mutate(month_joined = as.factor(month_joined)) %>% 
    arrange(desc(degree)) %>% 
    slice_head(n = 20)
  #save top 50 for TM
  top_50_users[[i]] <- RT_ntws_df[[i]] %>% 
    mutate(moth_joined = as.factor(month_joined)) %>% 
    arrange(desc(degree)) %>% 
    slice_head(n = 50)
  
  #mark ppl who are present in the last month 
  #create an empty vector to put results in 
  top_20_users[[i]]$retained <- 0
  
  for (i in 1:length(top_20_users)) {
    if (i == 1) {
      next  # Skip the first iteration as information for Nov turnover is not available
    } 
    
    if (!is.data.frame(top_20_users[[i]])) {
      cat("top_20_users[[i]] is not a data frame in iteration", i, "\n") 
      next
    }
    #mark retained from last month 1, otherwise 0
    top_20_users[[i]] <- top_20_users[[i]] %>% 
      mutate(retained = ifelse(top_20_users[[i]]$name %in% top_20_users[[i - 1]]$name, "1", "0"))
  }
  
  #plot 
  plots_deg[[i]] <- ggplot(top_20_users[[i]], aes(x = fct_rev(user_id), y = degree, fill = month_joined)) +
    geom_col() +
    coord_flip() +
    theme(axis.text.y = element_text(size = 10))
  
  for(i in length(plots_deg)){
    
    if (i == 1) {
      next  # Skip the first iteration
    } 
    plots_deg[[i]] <- plots_deg[[i]] +
      # add border to retained users to visually track turnover
      aes(color = retained) +
      scale_color_manual(values = c("lightgrey", "coral"))
  }
}

############# TURNOVER AND NEW NODES FUNCTIONS. #######
#### TURNOVER #####
# function to calculate turnover
calculate_turnover <- function(data_frames) {
  
  # empty vector to store turnover for each month
  turnover_rates <- c()
  #start at 2nd month in data
  for (i in 2:length(data_frames)) {
    prev_month_df <- data_frames[[i - 1]]
    current_month_df <- data_frames[[i]]
    
    # Calculate nodes that left the network
    nodes_left <- setdiff(prev_month_df$name, current_month_df$name)
    
    # Calculate turnover rate for this month
    turnover_rate <- length(nodes_left) / length(prev_month_df$name)
    
    #store turnover rates
    turnover_rates <- c(turnover_rates, turnover_rate)
  }
  
  return(turnover_rates)
}

#calculate overall turnover
turnover_overall <- calculate_turnover(RT_ntws_df)
#high turnover, relatively stable between 60 and 76%, except July, at 54%

#calculate turnover in top 50 
turnover_leaders <- calculate_turnover(RT_ntw_top_50)

### CACULATE NEW NODES ####### 
calculate_new <- function(x){
  new_percents <- c()
  for(i in 2:length(x)){
    #number of members in current month
    current_month_df <- x[[i]]
    previous_month_df <- x[[i - 1]]
    
    new_nodes <- setdiff(current_month_df$name, previous_month_df$name)
    new_percent <- length(new_nodes)/length(current_month_df$name)
    new_percents <- c(new_percents, new_percent)
  }
  return(new_percents)
}

# new nodes overall
new_overall <- calculate_new(RT_ntws_df)
# new nodes in top 50
new_50 <- calculate_new(RT_ntw_top_50)

#Turnover/New Nodes DF
turnover_df <- data.frame(month = month_no_NOV,
                           turnover_all = turnover_overall,
                           turnover_50 = turnover_leaders, 
                           new_overall = new_overall,
                           new_50 = new_50)


