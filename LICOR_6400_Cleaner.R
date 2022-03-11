library(tidyverse)

#this sets the wd to whichever file you want to extract files from
setwd("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/DATA/2022GH/LICOR")
#gets all files with the given file type
allFiles<-list.files(pattern = "*.csv")

outDF<- data.frame() ### creates a dataframe to store

for (i in 1:length(allFiles)){
  #for the current file read it in
  dat<- read.csv(allFiles[i], row.names = NULL)
  
  #puts meta data in a column
  meta<- dat[1,1]
  #stores the dat in which this file was created
  dat$meta <- meta 
  
  ## renames columns using by finding "obs"
  names(dat)<- as.character(unlist(dat[min(which(dat[,1] == "Obs")),])) 
  
  ### some irgas have Mch columns, this removes those
  dat <-dat[,-(which(grepl("Mch",names(dat))))]
  dat <- dat[,-82]
  #binds data frames together
  outDF<- rbind(outDF, dat)
}
dat <- outDF
## this function gets wrid of rows that aren't necessary
LC<-function(dat){
  #creates a comments column
  dat$comment<-"kitten"
  
  
  #This is to put comments in comment column####
  for (i in 1:nrow(dat)) {
    if (i>1){
      if (dat[i,1]=="Remark=") {
        dat$comment[i] <- dat[i,2]
      } else {
        dat$comment[i]<- dat$comment[i-1]
      }
    }
  }
  
  #deletes rows where Remarks don't take measuremetns (Remarks that we didn't type in on the machine)
  todelete<- c()
  for (i in 1:nrow(dat)) {
    if (i>1){
      if (dat[i,1]=="Remark=" & grepl(pattern = "\"", dat[i,2])) {## any row with remark, and a " in it are put into a list
        todelete <- append(todelete,i)
      }
    }
  }
  for (i in 1:nrow(dat)) {
    if (i>1){
      if (grepl(pattern = "light|ACi|aci",dat[i,which(names(dat) == "comment")])) {
        dat$comment[i] <- dat[i,2]
      } else {
        dat$comment[i]<- dat$comment[i-1]
      }
    }
  }
  dat <- dat[-todelete,]
  
  ## deletes leading rows that where named "kitten" earlier ^^^^
  dat <- dat[-which(dat$comment=="kitten"),] 
  dat <- dat[which(grepl(pattern = ":", dat[,2])),]
  #makes all numbers numeric and non numbers NAs
  dat[,1]<-as.numeric(as.character(unlist(dat[,1])))
  dat$comment <- sub("^.{1}", "", dat$comment)
  dat$comment <- sub(".{2}$","", dat$comment)
  dat$comment <- sub(".+ ","",dat$comment)
  return(dat)
}

d<- LC(outDF)
rm(outDF)

d$comment[which(d$comment == "n")]<- "7_ACi"
d<- d[-which(d$comment == "23_light"),]
d$comment[which(d$comment =="23_light_redo/actual")]<- "23_light"

require(tidyverse)
ggplot(d %>% filter(grepl("ACi",comment)),aes(x = as.numeric(CO2R),y = as.numeric(PhiPS2), color = comment))+
  geom_point()
ggplot(d %>% filter(grepl("light",comment)),aes(x = as.numeric(PARi),y = as.numeric(PhiPS2), color = comment))+
  geom_point()
plot(d$PARi)
ggplot(d %>% filter(as.numeric(d$PhiCO2)>0), aes(x = as.numeric(PhiCO2),y = as.numeric(PhiPS2)))+
  geom_point()+
  xlim(c(0,.03))
