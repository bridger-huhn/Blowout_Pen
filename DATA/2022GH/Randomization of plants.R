require(tidyverse)
dfilt <- read.csv("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/Greenhouse Randomization.csv")
set.seed(04272022)
irga <- dfilt %>% 
  filter(Size == "L"| Size =="M")
photosynq <- dfilt$Number
sample(irga, length(irga),replace = FALSE)
sample(photosynq, length(photosynq), replace = FALSE)
