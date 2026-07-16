precast_list_k1to25 = list()
for (k in 1:25){
  precast_list_k1to25[[paste0("k",as.character(k))]] <- PRECAST(PRECASTObj_log30, K = k)
}

##write loop to get all the ICs
###load final PRECAST clusters
aicList<-list()
for(k in 1:25){
  PRECASTObj <- SelectModel(precast_list_k1to25[[k]],criteria='AIC')
  aicList[k]<-PRECASTObj@resList$icMat[2]
  print(aicList[k])
  rm(PRECASTObj)
  gc()
}

saveRDS(precast_list_k1to25, file.path(data_path,"precast_log30_list_aic_k1to25.rds"))

aic<-as.numeric(unlist(aicList))

# Adjust the indices to start from 2
aic_df <- data.frame(Index = 1:25, AIC = aic)

png(file.path(fig_path,"kxaic.png"),width = 600, height = 400)
# Plot using ggplot2 with y-axis labels in scientific notation
ggplot(aic_df, aes(x = Index, y = AIC)) +
  geom_line() + # Add line
  geom_point()+
  scale_y_continuous(labels = scales::scientific)
dev.off()