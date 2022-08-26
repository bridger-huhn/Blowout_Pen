require(tidyverse)
require(lubridate)
require(plotly)
source("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/LICOR_6400_Cleaner.R")

#gets all files with the given file type should be .csv

d <- clean_raw_csvs("/Users/bridgerhuhn/Dropbox/2022_Pen_GH/Licor_6400csv/05302022-05312022") %>% 
  separate(comment, into = c("rep","curve"), sep = "_")

# put treatments in
meta <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/Metadata.csv")[,1:4] %>% 
  mutate(rep = as.character(Number)) %>% 
  select(!Number)
d<- left_join(d, meta)

#filter between ACi and light curves
ACi <- d %>% filter(grepl( "ACi", curve)) 
light <- d %>% filter(grepl("light", curve))




### photosynq data clean ####


phot <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/photosynq/spreadsheet (5).csv") %>% select(ID,time)


phot1 <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/photosynq/blowout-pen-greenhouse (4).csv") %>% mutate(ID = Datum.ID)
photosynq <- left_join(phot1, phot, by = "ID") %>% 
  select(time.y, Plant.ID,names(phot1)[which(str_detect(names(phot1), pattern = "Phi2"))]) %>% 
  pivot_longer(!c(time.y, Plant.ID),names_to = 'PAR', values_to = "phi") %>%  
  separate(col = PAR, into = c("trash", "PAR"), sep = "_", remove = FALSE) %>% 
  mutate(PAR = as.numeric(PAR)) %>% 
  select(!trash) %>% 
  mutate(time = as.POSIXct(time.y, format = "%m/%d/%Y %H:%M %p")) %>% 
  select(!time.y) %>% 
  rename(rep = Plant.ID)

photosynq <- left_join(photosynq, meta, by = "rep") %>% 
  select(!c(VWC_Sensor, Size))

photosynqfo <- left_join(phot1, phot, by = "ID") %>% 
  select( time.y, Plant.ID,names(phot1)[which(grepl("Fo",names(phot1)))]) %>% 
  pivot_longer(!c(time.y, Plant.ID),names_to = 'PAR', values_to = "Fo") %>%  
  separate(col = PAR, into = c("trash", "PAR"), sep = "_", remove = FALSE) %>% 
  mutate(PAR = as.numeric(PAR)) %>% 
  select(!trash) %>% 
  mutate(time = as.POSIXct(time.y, format = "%m/%d/%Y %H:%M %p")) %>% 
  select(!time.y) %>% 
  rename(rep = Plant.ID)

photosynqfo <- left_join(photosynqfo, meta, by = "rep") %>% 
  select(!c(VWC_Sensor, Size))
droughtfo <- photosynqfo %>% filter(Treatment == "drought")

photosynqNPQ <- left_join(phot1, phot, by = "ID") %>% 
  select( time.y, Plant.ID,names(phot1)[which(grepl("NPQ",names(phot1)))]) %>% 
  select(!names(phot1)[which(grepl("NPQt",names(phot1)))]) %>% 
  pivot_longer(!c(time.y, Plant.ID),names_to = 'PAR', values_to = "NPQ") %>%  
  separate(col = PAR, into = c("trash", "PAR"), sep = "_", remove = FALSE) %>% 
  mutate(PAR = as.numeric(PAR)) %>% 
  select(!trash) %>% 
  mutate(time = as.POSIXct(time.y, format = "%m/%d/%Y %H:%M %p")) %>% 
  select(!time.y) %>% 
  rename(rep = Plant.ID)

photosynqNPQ <- left_join(photosynqNPQ, meta, by = "rep") %>% 
  select(!c(VWC_Sensor, Size))
#### getting meteorological data in ####
meteor <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/meteorological/Pestemon_GH_Table1.csv", skip = 3)
names(meteor) <- names(read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/meteorological/Pestemon_GH_Table1.csv", skip = 1))

## look up table for what ports in the logger correspond to which reps
## .1-.15 = mux B, 2.1-2.15 = Mux A
lu_meteor <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/lu_meteorological_ports.csv")  %>% mutate("port_Mux"= paste(port,port.1, sep = "_")) %>% 
  select(-port,-port.1)
names(lu_meteor)[2]<- "port"

### this line grabs the port of each measurement
port <- gsub(".*[.]([^.]+)[.].*", "\\1", names(meteor))

### adds what mux each port is from
port[which(port == "1")[1]:which(port == "15")[2]] <- paste(port[which(port == "1")[1]:which(port == "15")[2]], "B", sep = "_")

port[which(port == "1")[1]:which(port == "15")[2]] <- paste(port[which(port == "1")[1]:which(port == "15")[2]], "A", sep = "_")


meas <- gsub("(WP_kPa).*","\\1", names(meteor))
meas <- gsub("(kohms).*","\\1",meas)

namesList <- left_join(data.frame(port = port, meas = meas), lu_meteor)
namesListTemp <- c(rep("",nrow(namesList)))
for (i in 1:nrow(namesList)) {
  if(is.na(namesList$rep[i])){
    namesListTemp[i] <- namesList$port[i]
  }else{namesListTemp[i]<- paste(namesList$meas[i], namesList$rep[i], sep = "_")}
}

names(meteor) <- namesListTemp

rm(meas, port,namesList,namesListTemp, lu_meteor)

meteor$TIMESTAMP <- as.POSIXct(meteor$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S")
#### graphing ####

for (i in 1:length(which(grepl("Pa",names(meteor))))) {
  lines(log10(meteor[,which(grepl("Pa",names(meteor)))[i]]))
}

metlong <- meteor %>% 
  select(TIMESTAMP,which(grepl("Pa",names(meteor)))) %>% 
  pivot_longer(cols = -TIMESTAMP, names_to ="ID", values_to = "kPa") %>% 
  filter(kPa >0 & kPa < 750) %>% 
  separate(ID, into = c("WP","KPa", "rep"), sep = "_") %>% 
  left_join(meta) %>% 
  select(-Size, -VWC_Sensor,-WP, -KPa)

###water content through time
ggplotly(ggplot(metlong, aes(x = TIMESTAMP, y = kPa,shape = rep)) +
  geom_line(aes(colour = Treatment )))

### NPQ late in experiment
ggplotly(ggplot(photosynqNPQ %>% filter(month(time)==5, day(time)==  30|31), aes(x = PAR, y = NPQ, color = Treatment)) + 
           geom_jitter())

#### PHI late in exp
ggplotly(ggplot(photosynq %>% filter(month(time)==5, day(time)==  30), aes(x = PAR, y = phi, color = Treatment)) + 
  geom_jitter())

#### Phi late in exp
ggplotly(ggplot(photosynq %>% filter(month(time)==5, day(time)==  30|31), aes(x = PAR, y = phi, color = Treatment)) + 
  geom_jitter())

### plots Light Curves
ggplotly(ggplot(light, aes(x = PARi, y = Photo, color = Treatment))+
           geom_point())

### plots ACi Curves
ggplotly(ggplot(ACi, aes(x = Ci, y = Photo, color = Treatment))+
           geom_point()+
           xlim(0,2000))

