#Average nearest neighbor calculation
#again, based from the DLC style

#first, we must calculate a new matrix of the distances between each fish to each other
#for an N=5, this will result in 10 distances. 

#distance function:

distance_function <- function(x, y, n, p) {
  sqrt((n-x)^2+(p-y)^2)
}

F1_F2 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X2, trialdata$Y2)
F1_F3 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X3, trialdata$Y3)
F1_F4 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X4, trialdata$Y4)
F1_F5 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X5, trialdata$Y5)
F2_F3 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X3, trialdata$Y3)
F2_F4 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X4, trialdata$Y4)
F2_F5 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X5, trialdata$Y5)
F3_F4 <- distance_function(trialdata$X3, trialdata$Y3, trialdata$X4, trialdata$Y4)
F3_F5 <- distance_function(trialdata$X3, trialdata$Y3, trialdata$X5, trialdata$Y5)
F4_F5 <- distance_function(trialdata$X4, trialdata$Y4, trialdata$X5, trialdata$Y5)

trialdata_coordinatedistance <- bind_cols(F1_F2, F1_F3, F1_F4, F1_F5, F2_F3, F2_F4, F2_F5, F3_F4, F3_F5, F4_F5)
names(trialdata_coordinatedistance) <- c("F1_F2", "F1_F3", "F1_F4", "F1_F5", "F2_F3", "F2_F4", "F2_F5", "F3_F4", "F3_F5", "F4_F5")
rm(F1_F2, F1_F3, F1_F4, F1_F5, F2_F3, F2_F4, F2_F5, F3_F4, F3_F5, F4_F5)

#now we have our distance matrix

#let's pull out the minimum distance, the maximum distance and the total distance from this df:

trialdata_coordinatedistance$distmin <- apply(trialdata_coordinatedistance[1:10], MARGIN =  1, FUN = min, na.rm = T)
trialdata_coordinatedistance$distmin[is.infinite(trialdata_coordinatedistance$distmin)] <- NA  
trialdata_coordinatedistance$distmax <- apply(trialdata_coordinatedistance[1:10], MARGIN =  1, FUN = max, na.rm = T)
trialdata_coordinatedistance$disttotal <- apply(trialdata_coordinatedistance[1:10], MARGIN =  1, FUN = sum, na.rm = T)
trialdata_coordinatedistance$distavg <- apply(trialdata_coordinatedistance[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#calculating avg. distance for randomly shuffled points:

#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[-c(11:14)]

###SEEDS

#seed for test
#set.seed(858)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 adults 1yr 30 min bonsai 040822
#set.seed(3892)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 adults 6mo 30 min bonsai 040822
#set.seed(8321)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 15dpf 30 min bonsai 040822
#set.seed(9809)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 23dpf 30 min bonsai 1 042522
#set.seed(8282)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 23dpf 30 min bonsai 2 042522
#set.seed(7768793)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 17dpf 30 min bonsai in dish 041922
#set.seed(6786)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 19dpf 5 min bonsai in dish 1 041222
#set.seed(38393)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#seed for 5 larva 19dpf 5 min bonsai in dish 2 041222
#set.seed(8383)
#trialdata_coordinatedistance_tmp <- trialdata_coordinatedistance_tmp[sample(1:nrow(trialdata_coordinatedistance_tmp)), sample(1:ncol(trialdata_coordinatedistance_tmp))] 
#trialdata_coordinatedistance_tmp$distavg <- apply(trialdata_coordinatedistance_tmp[1:10], MARGIN =  1, FUN = mean, na.rm = T)

#####and then we bind

#trialdata_coordinatedistance <- cbind(trialdata_coordinatedistance, trialdata_coordinatedistance_tmp$distavg)
#names(trialdata_coordinatedistance)[15]<- "randomdistavg"

#average distance index:

t#rialdata_coordinatedistance$distindex <- trialdata_coordinatedistance$distavg/trialdata_coordinatedistance$randomdistavg

#my calculation:

#trialdata_coordinatedistance$distindex2 <- trialdata_coordinatedistance$distavg/trialdata_coordinatedistance$distmax

#quicknote: if it's exactly equal to 1 that is because there is exactly 1 distance :)

#trialdata_coordinatedistance["distindex2"][trialdata_coordinatedistance["distindex2"] == 1] <- NA

#bind in a frame reference for the chart below:

trialdata_coordinatedistance <- trialdata_coordinatedistance %>%
  mutate(Frame = row_number())

trialdata_coordinatedistance <- trialdata_coordinatedistance %>%
  mutate(Second = Frame/25)

#######plots:

#ggplot(trialdata_coordinatedistance, aes (x=Frame, y=distindex))+
#  geom_line(size=0.3)+
#  ylab("Distance Index (Nearest Neighbor/Nearest Neighbor Shuffled)")+
#  xlab("Time (Frames)")+
#  theme(text = element_text(size = 20))+   
#  theme_classic()

#ggplot(trialdata_coordinatedistance, aes (x=Frame, y=distindex2))+
#  geom_line(size=0.3)+
#  ylab("Distance Index (Nearest Neighbor/Nearest Neighbor Shuffled)")+
#  xlab("Time (Frames)")+
#  theme(text = element_text(size = 20))+   
#  coord_cartesian(ylim=c(0, 1))+
#  theme_classic()

#after suggestion for just average distance over time: 
#note: maximum dist here would be like 540 mm

#ggplot(trialdata_coordinatedistance, aes (x=Frame, y=distavg))+
#  geom_line(size=0.3)+
#  ylab("Average Distance Between Animals (mm)")+
#  xlab("Time (Frames)")+
#  theme(text = element_text(size = 20))+   
#  theme_classic()


