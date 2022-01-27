require(tidyverse)
require(plotly)
require(lubridate)
d <- read.delim("DATA/Met_Data/Sanddune-Pen-1_PAR_RAD_PRECIP (1).dat", skip = 1,sep = ",")
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
