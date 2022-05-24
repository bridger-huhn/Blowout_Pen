require(tidyverse)
require(plotly)
require(lubridate)
d <- read.delim("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/Met_Data/Sanddune-Pen-1_PAR_RAD_PRECIP (1).dat", skip = 1,sep = ",")
d <- d[-c(1,2),]
d2 <- d %>% 
  mutate_if(is.character,as.numeric)
d2$TIMESTAMP <- d$TIMESTAMP

plot(d2$Rain_mm_Tot)
sum(d2$Rain_mm_Tot)

#####lower site####
soil2 <- read.delim("https://github.com/bridger-huhn/Blowout_Pen/raw/main/DATA/Met_Data/Sanddune-Pen-1_TEROS11_Table.dat", skip = 1,sep = ",")
soil2 <- soil2[-c(1,2),]  
soil <- soil2 %>% 
  mutate_if(is.character,as.numeric)
soil$TIMESTAMP <- soil2$TIMESTAMP
longSoilVWC <- soil %>% 
  mutate(D = as.POSIXct(TIMESTAMP)) %>% 
  mutate(DayOfYear = yday(D)) %>% 
  pivot_longer(cols = c(VWC_S1_Avg:T_N4_Std)) %>% 
  select(c(-TIMESTAMP)) %>% 
  rename(TIMESTAMP = D) %>% 
  mutate(depth = case_when(grepl("1", name) ~ "1",
                         grepl("2", name) ~ "2",
                         grepl("3", name) ~ "3",
                         grepl("4", name) ~ "4")) %>% 
  mutate(pit = case_when(grepl("N", name) ~ "N",
                         grepl("M", name) ~ "M",
                         grepl("S[0-9]", name) ~ "S")) %>% 
  filter(grepl("VWC", name)) %>% 
  filter(grepl("Avg", name)) %>% 
  rename(VWC = value) %>% 
  mutate(site = "lower dune")
  
rm(soil2, soil)

#####top of dune#####
soil2 <- read.delim("https://github.com/bridger-huhn/Blowout_Pen/raw/main/DATA/Met_Data/Sanddune-Pen-2_TEROS11_Table.dat", skip = 1,sep = ",")
soil2 <- soil2[-c(1,2),]  
soil <- soil2 %>% 
  mutate_if(is.character,as.numeric)
soil$TIMESTAMP <- soil2$TIMESTAMP
longSoilVWC2 <- soil %>% 
  mutate(D = as.POSIXct(TIMESTAMP)) %>% 
  mutate(DayOfYear = yday(D)) %>% 
  pivot_longer(cols = c(VWC_S1_Avg:T_N4_Std)) %>% 
  select(c(-TIMESTAMP)) %>% 
  rename(TIMESTAMP = D) %>% 
  mutate(depth = case_when(grepl("1", name) ~ "1",
                           grepl("2", name) ~ "2",
                           grepl("3", name) ~ "3",
                           grepl("4", name) ~ "4")) %>% 
  mutate(pit = case_when(grepl("N", name) ~ "N",
                         grepl("M", name) ~ "M",
                         grepl("S[0-9]", name) ~ "S")) %>% 
  filter(grepl("VWC", name)) %>% 
  filter(grepl("Avg", name)) %>% 
  rename(VWC = value) %>% 
  mutate(site = "upper dune")

rm(soil2, soil)
longSoilVWC<- rbind(longSoilVWC, longSoilVWC2)

### base R plots
# plot(soil2$VWC_M1_Avg,
#      ylim = c(0.05,.18))
# points(soil2$VWC_M2_Avg, col = "red")
# points(soil2$VWC_M3_Avg, col = "blue")
# points(soil2$VWC_M4_Avg, col = "green")

# animate soil moisture for a given pit
# longSoilVWC has columns "DayOfYear", "VWC", "depth", and "pit"
hProfileAnim <- longSoilVWC %>%
  filter(hour(TIMESTAMP) == 0, minute(TIMESTAMP) == 0) %>% 
  mutate(frames =  DayOfYear,
  ) %>% 
  plot_ly(
    x = ~ VWC, 
    y = ~ depth, 
    color = ~ pit,
    frame = ~ DayOfYear, 
    type = 'scatter',
    mode = 'lines',
    linetype = ~ site
  ) %>% 
  layout(xaxis = list(),
         yaxis = list(autorange="reversed")) %>% 
  animation_opts(
    frame = 100,
    transition = 20
  ) %>% 
  layout(
    xaxis = list(title = 'Daily Mean VWC'), 
    yaxis = list(title = 'Depth (cm)')
  ) %>% 
  animation_slider(
    currentvalue = list(prefix = "Timestamp: ")
  )
hProfileAnim


##### examine data #####
PRP<- read.delim("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/Met_Data/Sanddune-Pen-1_ATMOS14_Table.dat", skip = 1,sep = ",")
PRP <- PRP[-c(1:3),] 
PAR <- read.delim()
AT<- read.delim("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/Met_Data/Sanddune-Pen-1_ATMOS22_Table_correctdates.dat", skip = 1,sep = ",")
AT <- AT[-c(1:3),]

d <- merge(AT,PRP,by= c("TIMESTAMP","RECORD"))
d <- merge(d, soil2, by = c("TIMESTAMP","RECORD")) %>% 
  filter(atmosphericPress_Avg >78)%>% 
  filter(relativeHumidity_Avg> 0) %>% 
  select(TIMESTAMP,RECORD, AirTemp_Avg,windSpeed_Avg, vaporPress_Avg, VWC_S1_Avg:T_N4_Avg)

d$Avg_VWC <-  (as.numeric(d$VWC_S2_Avg)+
                          as.numeric(d$VWC_N2_Avg)+
                          as.numeric(d$VWC_M2_Avg))/3
d$wetting <- NA
for (i in 2:nrow(d)) {
  d$wetting[i]<- d$Avg_VWC[i]-d$Avg_VWC[i-1]
}
d$TIMESTAMP <- as.POSIXct(d$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S")
d$day <- yday(d$TIMESTAMP)
write_csv(data.frame(d$day), path = "/Users/bridgerhuhn/Documents/School/objective analysis of field data/HW#4/day.csv")
#write_csv(d, path = "/Users/bridgerhuhn/Documents/School/objective analysis of field data/HW#4/dat.csv")
ggplot(data = d, aes(x =as.numeric(T_S1_Avg), y = as.numeric(T_S2_Avg), col = as.numeric(RECORD)))+
  geom_point()
ggplotly((ggplot(data = d, aes(color= 1:nrow(d) ,y =as.numeric(VWC_M3_Avg), x = as.numeric(PAR))))+
           geom_point()+
           scale_color_continuous(low = "red", high = "green")+
           xlim(.03,.11)+
           ylim(.03,.11))