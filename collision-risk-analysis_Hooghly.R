##The script below is adapted from Martins et al. 2016, - A quantitative framework for investigating risk of deadly collisions between marine wildlife and boats
##To investigate per km collision risk of Ganges river dolphins with 3 types of boats - small, medium and large boats, in the river Hooghly of India
##Number of boats per km has been counted based on Google earth images and 
##mean length and width of boats were calculated based on 10 measurements taken at random for each of the three boat categories using the 'ruler' measuring tool of Google earth
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
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
v.dolph <- 1 #based on  Sasaki-Yamamoto et al. (2012) 
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
  else{
    list_boatspeed[[i]] <- NA
  }
}

##prepare dataframe of base encounter rates and fixed distance encounter rates for all boats
km <- vector(mode = 'list', length = nrow(df_long))
area <- vector(mode = 'list', length = nrow(df_long))
boat_speed_categories <- vector(mode = 'list', length = nrow(df_long))
lam_e <- vector(mode = 'list', length = nrow(df_long))
lam_FD <- vector(mode = 'list', length = nrow(df_long))

for(i in 1:nrow(df_long)){
  if(df_long$boat_number[i] > 0) {
    #resetting vector length for each i 
    lam_FD[[i]] <- numeric(length(list_boatspeed[[i]]))
    boat_speed_categories[[i]] <- numeric(length(list_boatspeed[[i]]))
    km[[i]] <- numeric(length(list_boatspeed[[i]]))
    area[[i]] <- numeric(length(list_boatspeed[[i]]))
    for(j in 1:length(list_boatspeed[[i]])){
      lam_e[[i]][j] <- getLambda(vb = list_boatspeed[[i]][j], 
                                 vm = v.dolph, 
                                 S = df_long$area[i], 
                                 rc = df_long$rc_radius[i]) 
      #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
      lam_FD[[i]][j] <-  lam_e[[i]][j] * (distance/list_boatspeed[[i]][j]) 
      boat_speed_categories[[i]][j] <- as.character(df_long$boat_type[i])
      km[[i]][j] <- df_long$km[i]
      area[[i]][j] <- df_long$area[i]
    }
    #printing results to verify (can do it for a subset of the df)
    print(km) 
  }
  else {
    lam_e[[i]] <- NA
    lam_FD[[i]] <- NA
    boat_speed_categories[[i]] <- as.character(df_long$boat_type[i])
    km[[i]] <- df_long$km[i]
    area[[i]] <- df_long$area[i]
  }
}

##checking some values
getLambda(vb = list_boatspeed[[1000]][2], vm = v.dolph, S = df_long$area[1000], rc = df_long$rc_radius[1000]) 
lam_e[[1000]][2]

df_speeds_encounters <- data.frame(km = unlist(km), area = unlist(area), boat_type = unlist(boat_speed_categories), boatspeed = unlist(list_boatspeed), lam_e = unlist(lam_e), lam_FD = unlist(lam_FD)) 
df_speeds_encounters$boat_type <- factor(df_speeds_encounters$boat_type, 
                                         levels = c("small_boats", "medium_boats", "large_boats" ))
df_speeds_encounters$lam_e_area <- df_speeds_encounters$lam_e/df_speeds_encounters$area
df_speeds_encounters$lam_FD_area <- df_speeds_encounters$lam_FD/df_speeds_encounters$area

##graphs
my_theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  strip.text = element_text(size = 10),
  legend.position = "bottom")

##boxplot of boatspeeds per boat type
boxplot_boatspeed <- ggplot(df_speeds_encounters, aes(x = boat_type, y = boatspeed, fill = boat_type)) +
                     geom_boxplot() +
                     labs(x = "Boat type", y = "Boat speed (m/s)", fill = "") +
                     my_theme

##boxplot of fixed distance encounter rate per boat type
boxplot_lamFD_area <- ggplot(df_speeds_encounters, aes(x = boat_type, y = lam_FD_area, fill = boat_type)) +
                geom_boxplot() +
                labs(x = "Boat type", y = "Fixed distance encounter rate per m2", fill = "") +
                my_theme

##graph of fixed distance encounter rate per m2 vs. boatspeed, for all 3 boat types
graph_speed_enc_small <- ggplot(df_speeds_encounters[df_speeds_encounters$boat_type == "small_boats", ], 
                           aes(x = boatspeed, y = lam_FD_area)) +
  facet_wrap(~boat_type) +
  geom_line(color = "#0072b2") +
  geom_point(color = "#0072b2") + 
  labs(x = "Boat speed (m/s)", y = "Fixed distance encounter rate per m2", fill = "") +
  coord_cartesian(ylim = c(0, 0.0000025)) +
  theme(legend.position = "none")

graph_speed_enc_medium <- ggplot(df_speeds_encounters[df_speeds_encounters$boat_type == "medium_boats", ], 
                                 aes(x = boatspeed, y = lam_FD_area)) +
  facet_wrap(~boat_type) +
  geom_line(color = "#f0e442") +
  geom_point(color = "#f0e442") + 
  labs(x = "Boat speed (m/s)", y = "Fixed distance encounter rate per m2", fill = "") +
  coord_cartesian(ylim = c(0, 0.0000025)) +
  theme(legend.position = "none")

graph_speed_enc_large <- ggplot(df_speeds_encounters[df_speeds_encounters$boat_type == "large_boats", ], 
                                aes(x = boatspeed, y = lam_FD_area)) +
  facet_wrap(~boat_type) +
  geom_line(color = "#d55e00") +
  geom_point(color = "#d55e00") + 
  labs(x = "Boat speed (m/s)", y = "Fixed distance encounter rate per m2", fill = "") +
  theme(legend.position = "none")

graph_all <- ggarrange(graph_speed_enc_small, graph_speed_enc_medium, graph_speed_enc_large, ncol = 3, nrow = 1)

ggsave(path = "Figs/For_MS",
       filename = "boxplot_boatspeed.png",
       plot = boxplot_boatspeed, 
       width = 30, height = 20, units = "cm",dpi = 300)

ggsave(path = "Figs/For_MS",
       filename = "graph_encFD.png",
       plot = boxplot_lamFD_area, 
       width = 30, height = 20, units = "cm",dpi = 300)

ggsave(path = "Figs/For_MS",
       filename = "graph_boatspeed_encounters.png",
       plot = graph_all, 
       width = 40, height = 15, units = "cm",dpi = 300)

write.csv(df_speeds_encounters, file = paste0("df_boatspeeds_encounter_rate.csv"))

##########----------------------------##########
##Simulation 1: vary probability of being within striking depth (p_striking_depth)
p_striking_depth <- seq(0.05, 0.3, by = 0.05)

##prepare outputs for simulations
lam_e <- c()
lam_FD <- c()
lam_FD_mean <- c()
lam_final <- c()
n_encounters <- c()
n_collisions <- c()
p_avoid <- c()
p_avoid_final <- c()
p_strike <- c()
n_sim <- 10000

##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km
for(i in 1:nrow(df_long)){
  if(df_long$boat_number[i] > 0) {
    #resetting vector length for each i 
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
    #mean encounter rate per boat
    df_long$lam_FD_mean[i] <- mean(lam_FD)
    #final encounter rate for one year
    df_long$lam_final[i] <- sum(lam_FD) * n.days * n.dolph 
    df_long$p_avoid_final[i] <- mean(p_avoid)
    #number of encounters follows a poisson distribution of final encounter rate
    df_long$n_encounters[i] <- mean(rpois(n_sim, df_long$lam_final[i]))
  }
  else {
    df_long$lam_FD_mean[i] <- NA
    df_long$lam_final[i] <- NA
    df_long$p_avoid_final[i] <- NA
    df_long$n_encounters[i] <- NA
  }
}

##summary of encounter results per boat in each km segment 
df_long$er_perboat <- df_long$n_encounters/df_long$boat_number
er_perboat_sum <- psych::describeBy(df_long$er_perboat, df_long$boat_type)

ggplot(df_long, aes(x = boat_type, y = er_perboat, fill = boat_type)) +
  geom_boxplot() +
  labs(x = "Boat type", y = "Mean expected number of encounters", fill = "") +
  coord_cartesian(ylim = c(0, 100)) +
  my_theme

##estimate number of collisions with all values of p_striking_depth
idx <- 1  # index tracker
for(x in 1:length(p_striking_depth)){
  p_strike[x] <- p_striking_depth[x]
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0){
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(df_long$n_encounters[i]), 
                                       p = p_strike[x] * (1 - df_long$p_avoid_final[i])))
    }
    else {
      n_collisions[idx] <- NA
    }
    idx <- idx + 1 
  }
}

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = length(p_striking_depth)), ])
p_striking_depth_fordf <- rep(p_striking_depth, each = nrow(df_long))
df_final <- cbind(df_final, p_striking_depth_fordf, n_collisions)

##mean number of encounters and collisions per boat per m2 area and summary
df_final$collision_percentage <- df_final$n_collisions/df_final$n_encounters
df_final$cr_perboat <- df_final$n_collisions/df_final$boat_number 

ggplot(df_final, aes(x = boat_type, y = cr_perboat, fill = boat_type)) +
  geom_boxplot() +
  labs(x = "Boat type", y = "Mean number of collisions per boat", fill = "") +
  coord_cartesian(ylim = c(0, 10)) +
  my_theme

df_final_summary <- df_final %>% group_by(boat_type, p_striking_depth_fordf) %>%
                    summarise(total_encounters = sum(n_encounters, na.rm = T), total_collision = sum(n_collisions, na.rm = T))

df_summary <- df_final %>% group_by(boat_type, p_striking_depth_fordf) %>% 
              summarise(mean_CR = mean(collision_rate, na.rm = T), median_CR = median(collision_rate, na.rm = T), min_CR = min(collision_rate, na.rm = T), max_CR = max(collision_rate, na.rm = T))
  
##graphs 
##violinplot of collision rate1 ~ p_striking_depth
graph_variable_pstrikingdepth_CR1 <- ggplot(df_final, aes(x = as.factor(p_striking_depth_fordf), y = collision_rate1, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  labs(x = "Probability of being within striking depth", y = "Collision rate (expected no. of collisions per boat)", fill = "Boat type") +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw() +
  my_theme 

ggsave(path = "Figs",
       filename = "collision_variable_pstrikingdepth_CR1.png",
       plot = graph_variable_pstrikingdepth_CR1, 
       width = 30, height = 20, units = "cm",dpi = 300)

write.csv(df_final, file = paste0("results_variable_pstrikingdepth_05.03.25.csv"))
write.csv(df_summary, file = paste0("summary_variable_pstrikingdepth.csv"))

##########----------------------------##########
##Simulation 2: vary probability of avoidance by the animal (p_avoid) for large boats for a seq(0.3, 0.6, by = 0.05)
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
for(i in 1:nrow(df_long)){
  if(df_long$boat_number[i] > 0) {
    #resetting vector length for each i 
    lam_FD <- numeric(length(list_boatspeed[[i]]))
    for(j in 1:length(list_boatspeed[[i]])){
      lam_e[j] <- getLambda(vb = list_boatspeed[[i]][j], 
                            vm = v.dolph, 
                            S = df_long$area[i], 
                            rc = df_long$rc_radius[i]) 
      #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
      lam_FD[j] <- lam_e[j] * (distance/list_boatspeed[[i]][j])
    }
    #printing results to verify (can do it for a subset of the df)
    print(lam_FD) 
    #mean encounter rate per boat
    df_long$lam_FD_mean[i] <- mean(lam_FD)
    #final encounter rate for one year
    df_long$lam_final[i] <- sum(lam_FD) * n.days * n.dolph 
    #number of encounters follows a poisson distribution of final encounter rate
    df_long$n_encounters[i] <- mean(rpois(n_sim, df_long$lam_final[i]))
  }
  else {
    df_long$lam_FD_mean[i] <- NA
    df_long$lam_final[i] <- NA
    df_long$n_encounters[i] <- NA
  }
}

##estimate number of collisions with all values of p_avoid_highspeed
idx <- 1  # index tracker
for(x in 1:length(p_avoid_highspeed)){
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0) {
      #resetting vector length for each i 
      p_avoid <- numeric(length(list_boatspeed[[i]]))
      for(j in 1:length(list_boatspeed[[i]])){
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
      #mean p_avoid
      p_avoid_final[idx] <- mean(p_avoid)
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(df_long$n_encounters[i]), 
                                       p = p_striking_depth * (1 - p_avoid_final[idx])))
    }
    else {
      p_avoid_final[idx] <- NA
      n_collisions[idx] <- NA
    }
    idx <- idx + 1 
  }
}

result_df <- data.frame(
  p_avoid_final = p_avoid_final,
  n_collisions = n_collisions
)

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = length(p_avoid_highspeed)), ])
p_avoid_highspeed_fordf <- rep(p_avoid_highspeed, each = nrow(df_long))
df_final <- cbind(df_final, p_avoid_highspeed_fordf, result_df)

##mean number of collisions per boat and summary
df_final$collision_rate1 <- df_final$n_collisions/df_final$boat_number 
df_final$collision_rate2 <- df_final$n_collisions/df_final$area
df_final$collision_rate3 <- df_final$collision_rate1/df_final$area 
df_summary <- df_final %>% group_by(boat_type, p_avoid_highspeed_fordf) %>% 
  summarise(mean_CR = mean(collision_rate, na.rm = T), median_CR = median(collision_rate, na.rm = T), min_CR = min(collision_rate, na.rm = T), max_CR = max(collision_rate, na.rm = T))

##graphs 
##violinplot of collision rate1 ~ p_avoid_highspeed
graph_variable_pavoid_CR1 <- ggplot(df_final, aes(x = as.factor(p_avoid_highspeed_fordf), y = collision_rate1, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  labs(x = "Probability of avoidance with high-speed boats (>=3.5 m/s)", y = "Collision rate (expected no. of collisions per boat)", fill = "Boat type") +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() +
  my_theme 

##violinplot of collision rate2 ~ p_avoid_highspeed
graph_variable_pavoid_CR2 <- ggplot(df_final, aes(x = as.factor(p_avoid_highspeed_fordf), y = collision_rate2, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  labs(x = "Probability of avoidance with high-speed boats (>=3.5 m/s)", y = "Collision rate (expected no. of collisions per m2)", fill = "Boat type") +
  coord_cartesian(ylim = c(0, 0.00005)) +
  theme_bw() +
  my_theme 

##violinplot of collision rate3 ~ p_avoid_highspeed
graph_variable_pavoid_CR3 <- ggplot(df_final, aes(x = as.factor(p_avoid_highspeed_fordf), y = collision_rate3, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  labs(x = "Probability of avoidance with high-speed boats (>=3.5 m/s)", y = "Collision rate (expected no. of collisions per boat per m2)", fill = "Boat type") +
  coord_cartesian(ylim = c(0, 0.00001)) +
  theme_bw() +
  my_theme 

ggsave(path = "Figs",
       filename = "collision_variable_pavoid_CR1.png",
       plot = graph_variable_pavoid_CR1, 
       width = 30, height = 20, units = "cm",dpi = 300)

ggsave(path = "Figs",
       filename = "collision_variable_pavoid_CR2.png",
       plot = graph_variable_pavoid_CR2, 
       width = 30, height = 20, units = "cm",dpi = 300)

ggsave(path = "Figs",
       filename = "collision_variable_pavoid_CR3.png",
       plot = graph_variable_pavoid_CR3, 
       width = 30, height = 20, units = "cm",dpi = 300)

write.csv(df_final, file = paste0("results_variable_pavoid_highspeed_05.03.25.csv"))
write.csv(df_summary, file = paste0("summary_variable_pavoid_highspeed.csv"))

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

##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km from simulation - 2

##estimate number of collisions for all values of p_striking depth and p_avoid_highspeed
idx <- 1  # index tracker
for(x in 1:nrow(new_df)){
  for(i in 1:nrow(df_long)){
    if(df_long$boat_number[i] > 0) {
      #resetting vector length for each i 
      p_avoid <- numeric(length(list_boatspeed[[i]]))
      for(j in 1:length(list_boatspeed[[i]])){
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
      #mean p_avoid
      p_avoid_final[idx] <- mean(p_avoid)
      ##number of collisions follows a binomial distribution where, 
      ##n = number of replications
      ##size = number of trials (i.e number of encounters)
      ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
      n_collisions[idx] <- mean(rbinom(n_sim, 
                                       size = as.integer(df_long$n_encounters[i]), 
                                       p = new_df$p_striking_depth[x] * (1 - p_avoid_final[idx])))
    }
    else {
      p_avoid_final[idx] <- NA
      n_collisions[idx] <- NA
    }
    idx <- idx + 1 
  }
}

result_df <- data.frame(
  p_avoid_final = p_avoid_final,
  n_collisions = n_collisions
)

##prepare dataframe to store results
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = nrow(new_df)), ])
new_df2 <- as.data.frame(new_df[rep(seq_len(nrow(new_df)), each = nrow(df_long)), ])
df_final <- cbind(df_final, new_df2, result_df)

##mean number of collisions per boat and summary
df_final$collision_rate1 <- df_final$n_collisions/df_final$boat_number 
df_final$collision_rate2 <- df_final$n_collisions/df_final$area
df_final$collision_rate3 <- df_final$collision_rate1/df_final$area  
df_summary <- df_final %>% group_by(boat_type, p_striking_depth, p_avoid_highspeed) %>% 
  summarise(mean_CR = mean(collision_rate, na.rm = T), median_CR = median(collision_rate, na.rm = T), 
            min_CR = min(collision_rate, na.rm = T), max_CR = max(collision_rate, na.rm = T))

##violinplot of collision rate1 ~ p_striking depth & p_avoid_highspeed
graph_variable_both_CR1 <- ggplot(df_final, aes(x = factor(p_striking_depth), y = collision_rate1, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  facet_wrap(~ p_avoid_highspeed) +
  scale_y_continuous(position = "left") +
  labs(title = "Collision rate by p-avoidance with high-speed boats (>=3.5 m/s)",
    x = "Probability of being within striking depth", 
    y = "Collision rate (expected no. of collisions per boat)",
    fill = "Boat type") +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() +
  my_theme 

##violinplot of collision rate2 ~ p_striking depth & p_avoid_highspeed
graph_variable_both_CR2 <- ggplot(df_final, aes(x = factor(p_striking_depth), y = collision_rate2, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  facet_wrap(~ p_avoid_highspeed) +
  scale_y_continuous(position = "left") +
  labs(title = "Collision rate by p-avoidance with high-speed boats (>=3.5 m/s)",
       x = "Probability of being within striking depth", 
       y = "Collision rate (expected no. of collisions per m2)",
       fill = "Boat type") +
  coord_cartesian(ylim = c(0, 0.00005)) +
  theme_bw() +
  my_theme 

##violinplot of collision rate3 ~ p_striking depth & p_avoid_highspeed
graph_variable_both_CR3 <- ggplot(df_final, aes(factor(x = p_striking_depth), y = collision_rate3, fill = boat_type)) + 
  geom_violin(position =  position_dodge(0.9)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.9)) +
  facet_wrap(~ p_avoid_highspeed) +
  scale_y_continuous(position = "left") +
  labs(title = "Collision rate by p-avoidance with high-speed boats (>=3.5 m/s)",
       x = "Probability of being within striking depth", 
       y = "Collision rate (expected no. of collisions per boat per m2)",
       fill = "Boat type") +
  coord_cartesian(ylim = c(0, 0.00001)) +
  theme_bw() +
  my_theme 

ggsave(path = "Figs",
       filename = "collision_variable_both_CR1.png",
       plot = graph_variable_both_CR1, 
       width = 40, height = 30, units = "cm",dpi = 300)

ggsave(path = "Figs",
       filename = "collision_variable_both_CR2.png",
       plot = graph_variable_both_CR2, 
       width = 40, height = 30, units = "cm",dpi = 300)

ggsave(path = "Figs",
       filename = "collision_variable_both_CR3.png",
       plot = graph_variable_both_CR3, 
       width = 40, height = 30, units = "cm",dpi = 300)

write.csv(df_final, file = paste0("results_variable_both_05.03.25.csv"))
write.csv(df_summary, file = paste0("summary_variable_both.csv"))




