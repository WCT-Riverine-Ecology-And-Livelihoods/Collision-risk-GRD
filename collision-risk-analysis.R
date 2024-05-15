library(dplyr)
library(ggplot2)
options(scipen = 999)
# function for mean encounter rate (lambda_e) for 1 animal-1 boat
# when animal speed is assumed to be fixed
getLambda <- function(vb, vm, rc, S)#vb: boat speed; vm: animal speed; rc= critical distance of encounter 
{
  getJ <- Vectorize(function(vb, vm) {
    alpha <- (2*vm*vb)/(vm^2+vb^2) #alpha 
    f.theta <- function(theta, alpha) sqrt(1 - alpha * cos(theta))/(2 * pi)
    integrate(f.theta, lower = 0, upper = 2 * pi, alpha = alpha)$value})
  
  (2*rc/S)*sqrt(vm^2 + vb^2)*getJ(vb, vm)
}

#parameters required for above function
area <- 300000000 # Hooghly study area in m^2 (300 km. sq.)
distance <- 550000 #distance travelled by a boat on the 550 km of Hooghly
n.days <- 30 #number of days
n.boats <- 20 #number of boats
n.dolph <- 250 #number of dolphins
v.dolph <-  0.9722222 #fixed speed (m/s) for ganges river dolphin based on  Sasaki-Yamamoto et al. (2012) - 3.5 km/hr
l_boat <- c(3, 10, 20) #length of boat (m)
w_boat <- c(1.5, 4, 8) #width of boat (m)
l_dolph <- rep(2.25, 3) #length of dolphin (m)
w_dolph <- rep(1.25, 3) #width of dolphin (m)

#three ways of getting critical distance - rc
df <- data.frame(l_boat, w_boat, l_dolph, w_dolph)
df$radius_boat <- sqrt((l_boat * w_boat/pi)) #calculating radius of encounter using areas of boat and animal
df$radius_dolph <- sqrt((l_dolph * w_dolph/pi))
df$rc_radius <- df$radius_boat + df$radius_dolph 
df$rc_width <- df$w_boat + df$w_dolph #when dolph and boat are parallel to each other
df$rc_lw <- df$w_boat + df$l_dolph #when dolph is perpendicular to boat

v.boat <- seq(2, 4.2, 0.01) #in m/s (based on vboat of 8-15 km/hr)

#data frame with all combinations of rc_radius and boat speed
df_rc_radius <- expand.grid(df$rc_radius, v.boat)
colnames(df_rc_radius) <- c("rc_radius", "v.boat")
df_rc_radius <- df_rc_radius %>% mutate(boat_type = case_when(
  rc_radius == rc_radius[1] ~ "small",
  rc_radius == rc_radius[2] ~ "medium",
  rc_radius == rc_radius[3] ~ "large",
))
  
df_rc_radius$area <- rep(area, nrow(df_rc_radius))
df_rc_radius$v.dolph <- rep(v.dolph, nrow(df_rc_radius))

#encounter rate
df_rc_radius$lambda_e <- getLambda(vb = df_rc_radius$v.boat, vm = df_rc_radius$v.dolph, S = df_rc_radius$area, rc = df_rc_radius$rc_radius)  #encounter rate (assumes a fixed animal speed)

#encounter rate when a boat is travelling a fixed distance
df_rc_radius$lambda_FD <- (df_rc_radius$lambda_e) * (distance/df_rc_radius$v.boat) 
df_rc_radius$lam_e_final <-  df_rc_radius$lambda_FD * n.boats * n.days * n.dolph #lambda re-scaled to account for number of days, number of boats and number of animals

#number of encounters follows a poisson distribution of final encounter rate
n_sim <- 10000
n_encounters <- c()
for(i in 1:nrow(df_rc_radius)){
  n_encounters[i] <- mean(rpois(n_sim, df_rc_radius$lam_e_final[i]))
}
df_rc_radius <- cbind(df_rc_radius, n_encounters)

#keeping probability of avoidance as 0
p_striking_depth <- 0.1
n_collisions <- c()
for(i in 1:nrow(df_rc_radius)){
 n_collisions[i] <- mean(rbinom(n_sim, size = as.integer(df_rc_radius$n_encounters[i]), p = p_striking_depth))
}
df_rc_radius <- cbind(df_rc_radius, n_collisions)
df_rc_radius$n_collisions <- round(df_rc_radius$n_collisions, 0)

##graph
ggplot(df_rc_radius, aes(x = v.boat, y = n_collisions, color = boat_type)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(breaks = seq(from = 100000, to = 600000, by = 50000)) +
  ylab("Number of collisions") +
  xlab("speed of boat (m/s)")
  
