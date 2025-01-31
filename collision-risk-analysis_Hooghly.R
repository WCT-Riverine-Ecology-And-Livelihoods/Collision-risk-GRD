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

##changing the striking depth parameter alone
p_striking_depth <- seq(0.05, 0.3, by = 0.05)
p_striking_depth <- rep(p_striking_depth, each = nrow(df_long))

##generate speeds for the boats per km using a normal distribution with the following parameters
##small boat speed: mean of 5 kmph and 3 sd, converted into m/s
##medium boat speed: mean of 10 kmph and 5 sd, converted into m/s
##large boat speed: mean of 15 kmph and 5 sd, converted into m/s
set.seed(12)
lam_e <- c()
lam_FD <- c()
p_avoid <- c()
n_sim <- 10000
df_final <- as.data.frame(df_long[rep(seq_len(nrow(df_long)), times = length(p_striking_depth)), ])
for(i in 1:nrow(df_final)){
  if(df_final$boat_number[i] > 0) {
  v.boat <- if(df_final$boat_type[i] == "small_boats") {
    rnorm(df_final$boat_number[i], 1.4, 0.8) 
  }
  else if (df_final$boat_type[i] == "medium_boats"){
    rnorm(df_final$boat_number[i], 2.8, 1.4) 
  }
  else {
    rnorm(df_final$boat_number[i], 4.2, 1.4) 
  }
  ##calculate encounter rate for 1 boat-1 dolphin for all the boats in each km
  ##setting p_avoid to vary according to the speed of boat in each km and taking the mean of p_avoid per km to calculate encounter rate
  for (j in 1:length(v.boat)){
    lam_e[j] <- getLambda(vb = v.boat[j], vm = v.dolph, S = df_final$area[i], rc = df_final$rc_radius[i])  #encounter rate (assumes a fixed animal speed)
    #encounter rate of 1 boat - 1 dolphin when boat is travelling the fixed distance of 1000 m
    lam_FD[j] <- lam_e[j] * (distance/v.boat[j])
    p_avoid[j] <- if(v.boat[j] <= 1.5)
      {
      0.95
      }
    else if(v.boat[j] > 1.5 & v.boat[j] < 3.5){
      0.7
    }
    else {
      0.5
    }
  }
  #encounter rate per km per year
  df_long$lam_final[i] <- sum(lam_FD) * n.days * n.dolph 
  df_long$p_avoid[i] <- mean(p_avoid)
  #number of encounters follows a poisson distribution of final encounter rate
  df_long$n_encounters[i] <- mean(rpois(n_sim, df_long$lam_final[i]))
  ##number of collisions follows a binomial distribution where, 
  ##n = number of replications
  ##size = number of trials (i.e number of encounters)
  ##p = probability of collision, which is the probability of being within striking depth AND not avoiding the boat
  df_long$n_collisions[i] <- mean(rbinom(n_sim, 
                                         size = as.integer(df_long$n_encounters[i]), 
                                         p = p_striking_depth * (1 - df_long$p_avoid[i])))
}
  else {
   df_long$lam_final[i] <- 0
   df_long$p_avoid[i] <- NA
   df_long$n_encounters[i] <- 0
   df_long$n_collisions[i] <- 0
  }
}

write.csv(df_long, file = "collision-results_22.01.25.csv")
