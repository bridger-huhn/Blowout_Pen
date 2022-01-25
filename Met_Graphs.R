require(tidyverse)
d <- read.delim("DATA/Met_Data/Sanddune-Pen-1_PAR_RAD_PRECIP (1).dat", skip = 1,sep = ",")
d <- d[-c(1,2),]
d2 <- d %>% 
  mutate_if(is.character,as.numeric)
d2$TIMESTAMP <- d$TIMESTAMP

plot(d2$Rain_mm_Tot)
sum(d2$Rain_mm_Tot)

soil2 <- read.delim("DATA/Met_Data/Sanddune-Pen-1_TEROS11_Table.dat", skip = 1,sep = ",")
soil2 <- soil2[-c(1,2),]  
soil <- soil2 %>% 
  mutate_if(is.character,as.numeric)
soil$TIMESTAMP <- soil2$TIMESTAMP
rm(soil2)
plot(soil$VWC_N1_Avg,
     ylim = c(0.05,.18))
points(soil$VWC_S2_Avg, col = "red")
points(soil$VWC_S3_Avg, col = "blue")
points(soil$VWC_S4_Avg, col = "green")
