
################################################ RETWEET NETWORK ###############################################

# This code handles the first three sections of my thesis's methods portion concerning RQ1: 

#        How did the structure of Pro-Eating disorder retweet networks evolve as far as the
#        shape of leadership and network structure before, during and after the first wave of
#        Covid-19 lockdowns following March 2020?

# This question is explored in the following parts 

  # Part 1: Constructing a Retweet Network
  # Part 2: Calculating Network Statistics
  # Part 3: Adding month_joined variable and anonymize users
  # Part 4: Fitting Powerlaw Lines to Degree Distributions
  # Part 5: Visualizing Top Users' Degree Distributions
  # Part 6: Calculating User Turnover 


# Required Packages 

library(tidyverse)     # data tranformation
library(igraph)        # create retweet network
library(ggplot2)       # visualization
library(modelsummary)  # network statistics
library(poweRlaw)      # fit powerlaw lines to the empirical data


                ########################## Part 1: Constructing Retweet Networks #########################

### Prepare data for monthly retweet networks

# load data cleaned and processed data 

load("final_ana.Rdata")


# create variable that classes tweets by month

final_ana <- final_ana %>% mutate(month = format(Date, "%Y-%m"))


# Separate data into monthly data frames

months_df_list <- final_ana %>% as_tibble() %>% 
  group_split(month) 


# make df with just retweets 

month_RT_ntw_df <- lapply(months_df_list, function(x){
  x %>% filter(sourcetweet_type == "retweeted")
})


### Make retweet netowork

# Created edgelist with user name and retweeted username

ALL_RT <- lapply(month_RT_ntw_df, function(x){
  x %>% 
    select(user_username, rt_username)
})

# make weighted edgelist

RT_edgelist <- lapply(ALL_RT, function(x){
  x %>% 
    group_by(user_username, rt_username) %>% 
    summarise(weight = sum(n()))
})

# create networks from edgelist

RT_ntws <- lapply(RT_edgelist, graph_from_data_frame)


# remove self loops and isolates

RT_ntws <- lapply(RT_ntws, function(x){
  igraph::simplify(x, remove.loop = T)
  x <- delete.vertices(x, which(degree(x) == 0))
})


                ######################## Part 2: Calculating Network Statistics #########################

# create a function that calculates network statistics and puts then into a dataframe 

calculate_network_stats <- function(network) {
  nodes <- vcount(network)
  ties <- ecount(network)
  density <- graph.density(network)
  reciprocity <- reciprocity(network)
  transitivity_global <- transitivity(network, type = "global")
  clustering_coef <- transitivity(network, type = "average")
  ave_path_length <- mean_distance(network)
  
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


# Calculate network statistics for each month

network_stats <- lapply(RT_ntws, calculate_network_stats)


                ############## Part 3: Add month_joined variable and anonymize users #################



#convert networks in to dataframes
RT_add_month_df <- lapply(RT_ntws, as_data_frame, what = "vertices")


# create a vector of months 

Months <- c("2019-11", "2019-12", "2020-01", "2020-02", "2020-03", "2020-04", "2020-05", "2020-06", "2020-07", "2020-08", "2020-09")


# add month variable to corresponding network 
for(i in seq_along(RT_ntws)){
  RT_add_month_df[[i]] <- RT_add_month_df[[i]] %>% mutate(month = Months[i])
}

# combine network data frames into one 

RT_ntw_full <- do.call(rbind, RT_add_month_df)


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


# anonymize user id 

# make a list of all users in networks
unique_users <- unique(RT_ntw_full$name)

# make a function that generates 4 random letters for the number of unique users in the full network

random_letters <- replicate(length(unique_users), paste(sample(letters, 4), collapse = ''))


# Create a column with user ids combining month_joined with random letters

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


                ############## Part 4: Fitting Powerlaw Lines to Degree Distributions #################

# Calculate in-degree centrality of each node in each network

degrees <- lapply(RT_ntws, degree, mode = "in")
for(i in seq_along(RT_ntws)){
  V(RT_ntws[[i]])$degree <- degrees[[i]]
}

# Convert networks into dataframes with all node attributing including degree distributions

RT_ntws_df <-  lapply(RT_ntws, as_data_frame, what = "vertices")

### Fit Powelaw lines to degree distributions for each network

# theories about patterns of leadership during covid were partially based on retweet
# networks' tendency to be scale-free, which have a degree distribution that follows a power-law
# distribution

# store results of analysis

est_pl_RT <- list()
pl_RT <- list()

for (i in seq_along(RT_ntws_df)) {
  
  # Create a power-law object for the degree distribution for each network
  pl_RT[[i]] <- displ$new(RT_ntws_df[[i]]$degree)
  
  # Estimate xmin: most organically occurring networks only express power-law behavior
  # after a certain value, this calculates that value along with the alpha scaling parameter
  #that determines the rate of decrease of the power-law line and goodness-of-fit of the actual 
  #network data to theoretical powerlaw distribution
  
  est_pl_RT[[i]] <- estimate_xmin(pl_RT[[i]])
  
  #add to powerlaw object
  pl_RT[[i]]$setXmin(est_pl_RT[[i]])
  
}


# make a dataframe for each month with results 

est_pl_RT_df <- Map(function(x, month) {
  data.frame(
    month = month,  # Use the corresponding month value
    gof = x$gof, # goodness of fit of each month's degree distribution to theoretical powerlaw distribution 
    xmin = x$xmin,   # minimun x value where powerlaw behavior begins
    alpha = x$pars  # alpha, the slope of the line
  )
}, est_pl_RT, Months)

# bind them all together 
est_pl_RT_full <- do.call(rbind, est_pl_RT_df)

# visualiztion


pl_plots_gg <- list()   # store the lists of monthly plots 
pl_lines <- list()      # store the powerlaw lines for each month to add to plots


for (i in seq_along(pl_results_list)) {
  pl_plots_gg[[i]] <- plot(pl_results_list[[i]])
  pl_lines[[i]] <- lines(pl_results_list[[i]], col = 2)
}

                ############## Part 5: Visualizing Top Users' Degree Distributions #################

# Create visualizations of top 20 and 50 users 

# Store dataframes of top 20 and top 50 users and plots

top_20_users <- list()
top_50_users <- list()
plots_deg <- list()

# 
for (i in seq_along(RT_ntws_df)) {
  # Convert name to factor with levels sorted in descending order of degree
  RT_ntws_df[[i]]$user_id <- fct_reorder(RT_ntws_df[[i]]$user_id, RT_ntws_df[[i]]$degree, .desc = TRUE)
  
  # save top degree users
  
  top_20_users[[i]] <- RT_ntws_df[[i]] %>% 
    mutate(month_joined = as.factor(month_joined)) %>% 
    arrange(desc(degree)) %>% 
    slice_head(n = 20)
  
  #save top 50 
  
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
  
  # plot results 
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

                ############## Part 6: Calculating User Turnover #################

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




