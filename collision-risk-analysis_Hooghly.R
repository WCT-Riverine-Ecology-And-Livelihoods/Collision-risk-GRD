##The script below is adapted from Martins et al. 2016, - A quantitative framework for investigating risk of deadly collisions between marine wildlife and boats
##To investiagte per km collision risk of Ganges river dolphins with 3 types of boats - small, medium and large boats, in the river Hooghly of India
##Number of boats per km has been counted based on Google earth images and 
##mean length and width of boats were calculated based on 10 measurements taken at random for each of the three boat categories using the 'ruler' measuring tool of Google earth
library(dplyr)
library(reshape2)
library(ggplot2)
options(scipen = 999)

##load data
df_hooghly <- read.csv("Data/Hooghly_boat_data_final.csv")
df_hooghly <- df_hooghly %>% select(-c(google_earth_image_date, person_involved, remarks))
  
#mean encounter rate (lambda_e) for 1 animal-1 boat (Martins et al., 2016), when animal speed is assumed to be fixed
#vb: boat speed; vm: animal speed; rc: critical distance of encounter; S: surface area of the study region 
getLambda <- function(vb, vm, rc, S)
{
  getJ <- Vectorize(function(vb, vm) {
    alpha <- (2*vm*vb)/(vm^2+vb^2) 
    f.theta <- function(theta, alpha) sqrt(1 - alpha * cos(theta))/(2 * pi)
    integrate(f.theta, lower = 0, upper = 2 * pi, alpha = alpha)$value})
  (2*rc/S)*sqrt(vm^2 + vb^2)*getJ(vb, vm)
}

##input variables
distance <- 1000 #distance traveled by a boat in each river segment in m
n.days <- 365
n.dolph <- 0.5 #number of dolphins per km 
v.dolph <- 0.97222 #3.5 km/hr in m/s based on  Sasaki-Yamamoto et al. (2012) 
##add column for area - considered as distance * width of channel
df_hooghly$area <- df_hooghly$river_width_m * distance
##melt the data frame to long format
df_long <- melt(df_hooghly, 
                id.vars = c("km", "river_width_m", "area"), 
                variable.name = "boat_type", 
                value.name = "boat_number") 
df_long <- df_long %>% arrange(km)
df_long$l_boat <- rep(c(5.4, 12.3, 27.8), times = nrow(df_hooghly)) #length of boat (m)
df_long$w_boat <- rep(c(1.5, 4.3, 6.9), times = nrow(df_hooghly)) #width of boat (m)
l_dolph <- 2.25 #length of dolphin (m)
w_dolph <- 0.5 #width of dolphin (m)

##calculating radius of encounter (rc) using areas of boat and animal
df_long$radius_boat <- sqrt((df_long$l_boat * df_long$w_boat/pi)) 
df_long$radius_dolph <- sqrt((l_dolph * w_dolph/pi))
df_long$rc_radius <- df_long$radius_boat + df_long$radius_dolph 

##generate speeds for the boats per km using a normal distribution with the following parameters
##small boat speed: mean of 5 kmph and 3 sd, converted into m/s
##medium boat speed: mean of 10 kmph and 5 sd, converted into m/s
##large boat speed: mean of 15 kmph and 5 sd, converted into m/s
set.seed(121)
list_boatspeed <- vector(mode = 'list', length = nrow(df_long))
for(i in 1:nrow(df_long)){
  if(df_long$boat_number[i] > 0) {
    list_boatspeed[[i]] <- if(df_long$boat_type[i] == "small_boats") {
      abs(rnorm(df_long$boat_number[i], 1.4, 0.8)) ##abs to convert negative to positive values 
    }
    else if (df_long$boat_type[i] == "medium_boats"){
      abs(rnorm(df_long$boat_number[i], 2.8, 1.4)) 
    }
    else {
      abs(rnorm(df_long$boat_number[i], 4.2, 1.4)) 
    }
  }
}

mean_boatspeed <- unlist(lapply(list_boatspeed, FUN = mean))

##########----------------------------##########
##Simulation 1: vary striking depth parameter alone
p_striking_depth <- seq(0.05, 0.3, by = 0.05)

##prepare outputs for simulations
lam_e <- c()
lam_FD <- c()
lam_final <- c()
n_encounters <- c()
n_collisions <- c()
p_avoid <- c()
p_avoid_final <- c()
p_strike <- c()
n_sim <- 10000

##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km
idx <- 1  # index tracker
for(x in 1:length(p_striking_depth)){
  p_strike[x] <- p_striking_depth[x]
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0) {
      #resetting vectors for subsequent iterations of i (nrow(df_trial))
      lam_FD <- numeric(length(list_boatspeed[[i]]))
      p_avoid <- numeric(length(list_boatspeed[[i]]))
      for(j in 1:length(list_boatspeed[[i]])){
        lam_e[j] <- getLambda(vb = list_boatspeed[[i]][j], 
                              vm = v.dolph, 
                              S = df_long$area[i], 
                              rc = df_long$rc_radius[i]) 
        #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
        lam_FD[j] <- lam_e[j] * (distance/list_boatspeed[[i]][j])
        ##Probability to avoid
        if(list_boatspeed[[i]][j] <= 1.5)
        {
          p_avoid[j] <- 0.95
        }
        else if(list_boatspeed[[i]][j] > 1.5 & list_boatspeed[[i]][j] < 3.5){
          p_avoid[j] <- 0.7
        }
        else {
          p_avoid[j] <- 0.5
        }
      }
      #printing results to verify (can do it for a subset of the df)
      print(lam_FD) 
      #encounter rate per km per year
      lam_final[idx] <- sum(lam_FD) * n.days * n.dolph 
      p_avoid_final[idx] <- mean(p_avoid)
      #number of encounters follows a poisson distribution of final encounter rate
      n_encounters[idx] <- mean(rpois(n_sim, lam_final[idx]))
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(n_encounters[idx]), 
                                       p = p_strike[x] * (1 - p_avoid_final[idx])))
    }
    else {
      lam_final[idx] <- 0
      p_avoid_final[idx] <- NA
      n_encounters[idx] <- 0
      n_collisions[idx] <- 0
    }
    idx <- idx + 1 
  }
}

result_df <- data.frame(
  lam_final = lam_final,
  p_avoid_final = p_avoid_final,
  n_encounters = n_encounters,
  n_collisions = n_collisions
)

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = length(p_striking_depth)), ])
mean_boatspeed_fordf <- rep(mean_boatspeed, times = length(p_striking_depth))
p_striking_depth_fordf <- rep(p_striking_depth, each = nrow(df_long))
df_final <- cbind(df_final, mean_boatspeed_fordf, p_striking_depth_fordf, result_df)

df_final$collision_rate <- df_final$n_collisions/df_final$area

graph_variable_pstrikingdepth <- ggplot(df_final, aes(p_avoid_final, collision_rate)) + 
  geom_point() + 
  facet_wrap(~factor(p_striking_depth_fordf))+
  labs(x = "p-avoid", y = "Collision rate (no. of collisions/m2)")

graph_variable_pstrikingdepth2 <- ggplot(df_final, aes(p_avoid_final, collision_rate)) + 
  geom_point() + 
  ylim(0, 0.00025) + ##restrict ylim to remove outliers
  facet_wrap(~factor(p_striking_depth_fordf)) +
  labs(x = "p-avoid", y = "Collision rate (no. of collisions/m2)") 
  
ggsave(path = "Figs",
       filename = "collision_variable_pstrikingdepth2.png",
       plot = graph_variable_pstrikingdepth2, 
       width = 30, height = 20, units = "cm",dpi = 300)
write.csv(df_final, file = paste0("results_variable_pstrikingdepth.csv"))

##########----------------------------##########
##Simulation 2: vary p_avoid for large boats for a seq(0.3, 0.6, by = 0.05)
p_striking_depth <- 0.05
p_avoid_highspeed <- seq(0.3, 0.6, by = 0.05)

##prepare outputs for simulations
lam_e <- c()
lam_FD <- c()
lam_final <- c()
n_encounters <- c()
n_collisions <- c()
p_avoid <- c()
p_avoid_final <- c()
n_sim <- 10000

##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km
idx <- 1  # index tracker
for(x in 1:length(p_avoid_highspeed)){
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0) {
      #resetting vectors for subsequent iterations of i (nrow(df_trial))
      lam_FD <- numeric(length(list_boatspeed[[i]]))
      p_avoid <- numeric(length(list_boatspeed[[i]]))
      for(j in 1:length(list_boatspeed[[i]])){
        lam_e[j] <- getLambda(vb = list_boatspeed[[i]][j], 
                              vm = v.dolph, 
                              S = df_long$area[i], 
                              rc = df_long$rc_radius[i]) 
        #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
        lam_FD[j] <- lam_e[j] * (distance/list_boatspeed[[i]][j])
        ##Probability to avoid
        if(list_boatspeed[[i]][j] <= 1.5)
        {
          p_avoid[j] <- 0.95
        }
        else if(list_boatspeed[[i]][j] > 1.5 & list_boatspeed[[i]][j] < 3.5){
          p_avoid[j] <- 0.7
        }
        else {
          p_avoid[j] <- p_avoid_highspeed[x]
        }
      }
      #printing results to verify (can do it for a subset of the df)
      print(p_avoid_highspeed[x]) 
      #encounter rate per km per year
      lam_final[idx] <- sum(lam_FD) * n.days * n.dolph 
      p_avoid_final[idx] <- mean(p_avoid)
      #number of encounters follows a poisson distribution of final encounter rate
      n_encounters[idx] <- mean(rpois(n_sim, lam_final[idx]))
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(n_encounters[idx]), 
                                       p = p_striking_depth * (1 - p_avoid_final[idx])))
    }
    else {
      lam_final[idx] <- 0
      p_avoid_final[idx] <- NA
      n_encounters[idx] <- 0
      n_collisions[idx] <- 0
    }
    idx <- idx + 1 
  }
}

result_df <- data.frame(
  lam_final = lam_final,
  p_avoid_final = p_avoid_final,
  n_encounters = n_encounters,
  n_collisions = n_collisions
)

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = length(p_avoid_highspeed)), ])
mean_boatspeed_fordf <- rep(mean_boatspeed, times = length(p_avoid_highspeed))
p_avoid_highspeed_fordf <- rep(p_avoid_highspeed, each = nrow(df_long))
df_final <- cbind(df_final, mean_boatspeed_fordf, p_avoid_highspeed_fordf, result_df)

df_final$collision_rate <- df_final$n_collisions/df_final$area

graph_variable_pavoid <- ggplot(df_final, aes(p_avoid_final, collision_rate)) + 
  geom_point() + 
  facet_wrap(~factor(p_avoid_highspeed_fordf))+
  labs(x = "p-avoid", y = "Collision rate (no. of collisions/m2)")

ggsave(path = "Figs",
       filename = "collision_variable_pavoid.png",
       plot = graph_variable_pavoid, 
       width = 30, height = 20, units = "cm",dpi = 300)
write.csv(df_final, file = paste0("results_variable_pavoid_highspeed.csv"))

##########----------------------------##########
##Simulation 3: vary p_striking depth and p_avoid
p_strike <- seq(0.05, 0.3, by = 0.05)
p_avoid_highspeed <- seq(0.3, 0.6, by = 0.05)
new_df <- expand.grid(p_strike, p_avoid_highspeed)
colnames(new_df) <- c("p_striking_depth", "p_avoid_highspeed")
  
##prepare outputs for simulations
lam_e <- c()
lam_FD <- c()
lam_final <- c()
n_encounters <- c()
n_collisions <- c()
p_avoid <- c()
p_avoid_final <- c()
n_sim <- 10000

##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km
idx <- 1  # index tracker
for(x in 1:nrow(new_df)){
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0) {
      #resetting vectors for subsequent iterations of i (nrow(df_trial))
      lam_FD <- numeric(length(list_boatspeed[[i]]))
      p_avoid <- numeric(length(list_boatspeed[[i]]))
      for(j in 1:length(list_boatspeed[[i]])){
        lam_e[j] <- getLambda(vb = list_boatspeed[[i]][j], 
                              vm = v.dolph, 
                              S = df_long$area[i], 
                              rc = df_long$rc_radius[i]) 
        #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
        lam_FD[j] <- lam_e[j] * (distance/list_boatspeed[[i]][j])
        ##Probability to avoid
        if(list_boatspeed[[i]][j] <= 1.5)
        {
          p_avoid[j] <- 0.95
        }
        else if(list_boatspeed[[i]][j] > 1.5 & list_boatspeed[[i]][j] < 3.5){
          p_avoid[j] <- 0.7
        }
        else {
          p_avoid[j] <- new_df$p_avoid_highspeed[x]
        }
      }
      #printing results to verify (can do it for a subset of the df)
      print(new_df$p_avoid_highspeed[x]) 
      #encounter rate per km per year
      lam_final[idx] <- sum(lam_FD) * n.days * n.dolph 
      p_avoid_final[idx] <- mean(p_avoid)
      #number of encounters follows a poisson distribution of final encounter rate
      n_encounters[idx] <- mean(rpois(n_sim, lam_final[idx]))
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(n_encounters[idx]), 
                                       p = new_df$p_striking_depth[x] * (1 - p_avoid_final[idx])))
    }
    else {
      lam_final[idx] <- 0
      p_avoid_final[idx] <- NA
      n_encounters[idx] <- 0
      n_collisions[idx] <- 0
    }
    idx <- idx + 1 
  }
}

result_df <- data.frame(
  lam_final = lam_final,
  p_avoid_final = p_avoid_final,
  n_encounters = n_encounters,
  n_collisions = n_collisions
)

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = nrow(new_df)), ])
mean_boatspeed_fordf <- rep(mean_boatspeed, times = nrow(new_df))
new_df2 <- as.data.frame(new_df[rep(seq_len(nrow(new_df)), each = nrow(df_long)), ])
df_final <- cbind(df_final, mean_boatspeed_fordf, new_df2, result_df)

df_final$collision_rate <- df_final$n_collisions/df_final$area

graph_variable_both <- ggplot(df_final, aes(p_avoid_final, collision_rate)) + 
  geom_point() + 
  facet_grid(factor(p_avoid_highspeed)~factor(p_striking_depth), switch = "y") +
  scale_y_continuous(position = "right") +
  labs(x = "p-avoid", y = "Collision rate (no. of collisions/m2)")

ggsave(path = "Figs",
       filename = "collision_variable_both.png",
       plot = graph_variable_both, 
       width = 40, height = 30, units = "cm",dpi = 300)
write.csv(df_final, file = paste0("results_variable_both.csv"))




