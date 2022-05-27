require(tidyverse)
d <- read.csv("Greenhouse Randomization.csv")
d2 <- d %>% 
  filter(Number != 30,Number !=26,Number !=14,
         VWC_Sensor =="s")
d3 <- d %>% 
  filter(VWC_Sensor !="s")
set.seed(05052022)
sample(d2$Number, size = length(d2$Number), replace = FALSE)         
length(sample(d$Number, replace = FALSE))



###Create treatments ####
set.seed(05052022)
ww <- data.frame(Number = sample(d2$Number, size = 7, replace = FALSE))
ww$treatment <- "WW"
d2 <- left_join(d2, ww,by = "Number")

droughtburial <- data.frame(Number = sample(d2$Number[which(is.na(d2$treatment))], size = 7, replace = FALSE))
droughtburial$treatment <- "droughtburial"
d2 <- left_join(d2, droughtburial, by = "Number")

d2$treatment.x[which(!is.na(d2$treatment.y))]<- d2$treatment.y[which(!is.na(d2$treatment.y))]

drought <- data.frame(Number = sample(d2$Number[which(is.na(d2$treatment.x))], size = 7, replace = FALSE))
drought$treatment <- "drought"
d2 <- left_join(d2, drought, by = "Number")

d2$treatment.x[which(!is.na(d2$treatment))]<- d2$treatment[which(!is.na(d2$treatment))]

d2$treatment.x[which(is.na(d2$treatment.x))]<- "burial"
d2 <- d2[,-c(3,6,7)]
names(d2)[4]<- "Treatment"
d3$Treatment <- sample(c(unique(d2$Treatment),"drought","burial"),size = length(d3$Treatment), replace = FALSE)
dFINAL <- rbind(d2,d3)
#write_csv(dFINAL, path = "DATA/Metadata.csv")
####rand with burial ####
rnd <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/Metadata.csv") %>% 
  filter(Treatment != "burial"& Treatment!= "droughtburial")
set.seed(05252022) #### put the date of measurements to get the randomization
sample(rnd$Number, size = length(rnd$Number), replace = FALSE)    
drought <- rnd %>% filter(Treatment == "drought")
ww <- rnd %>% filter(Treatment == "WW")


sample(ww$Number, size = nrow(ww), replace = FALSE)
sample(drought$Number, size = nrow(drought), replace = FALSE)


all <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/Metadata.csv") 
VWC_Sensors <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/Metadata.csv") %>% 
  filter(VWC_Sensor == "s")
T1 <-  VWC_Sensors %>% filter(Treatment == "WW")
T2 <-  VWC_Sensors %>% filter(Treatment == "drought")
T3 <- VWC_Sensors %>% filter(Treatment =="droughtburial")
T4 <- VWC_Sensors %>% filter(Treatment == "burial")

set.seed(05262022)

sample(T1$Number, size = 7, replace = FALSE)
sample(T2$Number, size = 7, replace = FALSE)
sample(T3$Number, size = 8, replace = FALSE)
sample(T4$Number, size = 8, replace = FALSE)
sample(all$Number, size = nrow(all), replace = FALSE)
