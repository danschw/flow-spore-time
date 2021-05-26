## Coevolution with a seed bank - analysis of sporulation in evolved lines
library(renv)
# renv::restore()
library(here)
#analysis of flow-cytometry population data using Karava's analysis.R as reference.
source(here("src/functions.R"))
source(here("src/DS_functions.R"))
set.seed(1)
#choose data

   folders <- 
   list.dirs(here("data/FCM/fcs/"),recursive = F)
   

   
   #set this manually to process data of a single day (folder)
   day = "sample_data"
   folders <- folders[grep(paste0(day), folders)]
   folders <- list.dirs(folders)
   # remove the folder of folders
   folders <- folders[-1]

   # Make directories to store the data
   if (! dir.exists(here("fig/gate_plots", day))){
      dir.create(here("fig/gate_plots", day))
   }

   if (! dir.exists(here("data/output", day))){
      dir.create(here("data/output", day))
   }
   
# dilution of culture analyzed in FCM:
dilution <- "x100"

#fields to parse from the .fcs file name (separated by underscore)
sample.var <- c("species","strain","colony","t.culture","medium","well","rep","xt","num") 

for (folder in folders[]){
#### Load data, sample set ####
fcsset <- flowCreateFlowSet(filepath = folder, sample_variables = sample.var,
                            transformation = FALSE,separators = "[-_]")
#transform with arcsine, recpmendded by Karava et al.
fcsset <- Transform.Novocyte(fcsset)

#data frame to collect stats
df.stats <- fcsset%>%
   flowFcsToDf(.)%>%
   select(.,sample.var)%>%
   distinct()
df.stats$dilution <- dilution
df.stats$total.events <- NA
df.stats$volume.nL <- NA
df.stats$limit.events <- NA
df.stats$limit.volume.uL <- NA
df.stats$singlets <- NA
df.stats$neg.rmv <- NA
df.stats$noise.cutoff <- NA


#### Gating for singlets with flowStats ####

# The gate function needs to be applied to each sample seperatly
# get number of samples
n.sample <- nrow(fcsset@phenoData@data)



#initialise list to store singlet plots
plot.list <- vector('list', n.sample)

for (i in 1:n.sample){
   # collect metadata for stats table
   df.stats$volume.nL[i] <- as.numeric(fcsset[[i]]@description$`$VOL`)
   df.stats$total.events[i] <- as.numeric(fcsset[[i]]@description$`$TOT`)
   df.stats$limit.volume.uL[i] <- as.numeric(fcsset[[i]]@description$`#NCVolumeLimits`)
   df.stats$limit.events[i] <- as.numeric(fcsset[[i]]@description$`#NCEventsLimits`)
   

      singlet_gate <- gate_singlet(fcsset[[i]], area = "asinh.FSC.A", height = "asinh.FSC.H", filterId = "Singlets",wider_gate = TRUE )
      df.stats$singlets[i] <- summary(filter(fcsset[[i]],singlet_gate))@true

      #plot gate
      id <- fcsset[[i]]@description$GUID
      plot.list[[i]] <-
         as.ggplot(
            ggcyto(fcsset[[i]], aes(x = `asinh.FSC.A`, y =  `asinh.FSC.H`))+
               geom_hex(bins = 500)+
               geom_gate(singlet_gate, size=0.01)+
               geom_stats(adjust=c(0.2,0.75))+
               theme_cowplot(font_size = 8)+
               scale_y_continuous(labels =  function(x) format(x, scientific = TRUE))+
               scale_x_continuous(labels =  function(x) format(x, scientific = TRUE))+
               facet_null()+
               ggtitle(df.stats$well[i])

            )

      #apply gate
      fcsset[[i]] <- fcsset[[i]] %>%
            Subset(singlet_gate)
      
      
      # filter negatives
      neg.gate <- rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.FSC.A" = c(0, 15),"asinh.SSC.A" = c(0, 15))
      df.stats$neg.rmv[i] <-  df.stats$singlets[i]- summary(filter(fcsset[[i]],neg.gate))@true
      #apply gate
      fcsset[[i]] <- fcsset[[i]] %>%
         Subset(neg.gate) # remove negative
      
      #find cutoff to remove noise by FSC
      noise<- rangeGate(x = fcsset[[i]], "asinh.FSC.A", alpha=0.5,sd=3, absolute = F)
      df.stats$noise.cutoff[i] <- noise@min[[1]]

}

#save plot
ggsave2(filename = paste0(fcsset[[i]]@description$`$SRC`,".pdf" ), 
        plot = plot_grid(plotlist = plot.list),
        path = here("fig/gate_plots/", day))


#### Gating for singlets with norm2Filter ####
# df.set <- fcsset %>%
#    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.A", scale.factor = 10)) %>%
#    Subset(rectangleGate("asinh.BL1.A" = c(0, 15))) %>%
#    flowFcsToDf(.)

# plotting noise gate

noise.plot <-   
ggcyto(fcsset,aes(asinh.FSC.A))+
   geom_density(fill="grey80")+
   geom_vline(data = df.stats, aes(xintercept=noise.cutoff), color="red")+
   facet_wrap(~well)+theme_cowplot()

ggsave2(filename = paste0("noise_",fcsset[[i]]@description$`$SRC`,".pdf" ), 
        plot = noise.plot,
        path = here("fig/gate_plots/", day))

#### transform to dataframe ####
df.set <- fcsset%>%
      flowFcsToDf(.)
# no. of events before noise filter
df.stats <- 
   df.set %>%
      group_by(well) %>%
      summarise(noisy.events=n())%>%
      full_join(df.stats,.)

# # # plot
# df.set %>%
#       ggplot(aes(asinh.SSC.A, asinh.BL1.A))+
#       geom_hex(bins=500)+
#       facet_wrap(~well)

clean.df <- df.set[0,]

# remove noise
wells <- levels(as.factor(df.stats$well))
for (wl in wells){
   clean.df <- rbind(clean.df,
   df.set %>%
      dplyr::filter(well==wl)%>%
      dplyr::filter(asinh.FSC.A>df.stats$noise.cutoff[df.stats$well==wl]))
}


df.stats <- 
   clean.df %>%
      group_by(well) %>%
      summarise(clean.events=n())%>%
      full_join(df.stats,.)

complot <- 
   ggplot(df.set, aes(asinh.SSC.A, asinh.BL1.A))+
   geom_point(size=0.1, color="grey")+
   geom_point(data=clean.df,size=0.1, color="blue")+
   geom_density2d(col = "white",  size = 0.1, alpha = 0.5) +
   scale_alpha_continuous(guide = FALSE) +
   facet_wrap(~well)+
   theme_cowplot()+
   panel_border()

# ggsave2(paste0("fig/gate_plots/scatterNoise_",fcsset[[i]]@description$`$SRC`,".png"), complot)
ggsave2(filename = paste0("scatterNoise_",fcsset[[i]]@description$`$SRC`,".png" ), 
        plot = complot,
        path = here("fig/gate_plots/", day))
#### predict centers of sub-populations ####

# Is host strain WT or mutant?
mutant.host <- ( grepl("SN",df.stats$strain[1]) |
               grepl("dSpo",df.stats$strain[1]))
# Load reference containing all subpopulations
# I pre-compiled  models
# if (mutant.host){
# base::load(here("data/cluster_models","SNCt_cluster_model.Rdata"))
# } else {
# base::load(here("data/cluster_models","WLCt_cluster_model.Rdata"))
# }
# # use this to generate prediction model from data frame
# #WT model based on day 0.5 WLCt_L1
# df.mix <- clean.df%>%
#       # dplyr::filter(phage=="noPHI")%>%
#       dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
#       as.matrix() %>%
#       mclust::Mclust(data = ., G = 2) #I  have only 2 clusters
# save(df.mix, file=here("data/cluster_models","WLCt_cluster_model.Rdata")
# #dSpoIIE model based on day 0.5 SNCt_L1
# df.mix <- clean.df%>%
#    # dplyr::filter(phage=="noPHI")%>%
#    dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
#    as.matrix() %>%
#    mclust::Mclust(data = ., G = 1) #I  have only 1 cluster (veg cells)
# save(df.mix, file=here("data/cluster_models","SNCt_cluster_model.Rdata")
#WT evolved model based on WSCt-L1_D14_T24H_C8 on new novocyte
# df.mix <- clean.df%>%
#       # dplyr::filter(phage=="noPHI")%>%
#       dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
#       as.matrix() %>%
#       mclust::Mclust(data = ., G = 2) #I  have only 2 clusters
# save(df.mix, file=here("data/cluster_models","WSCt-L1_D14_T24H_C8_cluster_model.Rdata"))

# model from data of new novocyte 
base::load(file=here("data/cluster_models","sample_data_cluster_model.Rdata"))

# getting centers for visualization and export

centers.list.df <- t(df.mix$parameters$mean)

if (mutant.host){
   centers.list.df <-
   rbind(centers.list.df, matrix(c(10,5), nrow = 1))
}
# center.locs <- factor(df.mix$classification, levels = c(1, 2))

# assigning cluster to population based on SYBR fluoresence (BL1)
# higher SYBR => veg cell
pop.tbl <- data.frame(cluster=c(1,2), pop=NA)
pop.tbl$pop[which.max(centers.list.df[,"asinh.BL1.A"])] <- "veg"
pop.tbl$pop[which.min(centers.list.df[,"asinh.BL1.A"])] <- "spore"

### Prediction of clusters for all samples ####
cluster.predict <- clean.df %>%
      dplyr::select( asinh.FSC.A, asinh.BL1.A) %>%
      predict.Mclust(df.mix,.)

clean.df$pop <- sapply(cluster.predict$classification,  function(x) {ifelse(x==pop.tbl$cluster[1], pop.tbl$pop[1],pop.tbl$pop[2])})

p <-       
   # data.frame(clean.df, cluster = cluster.predict$classification) %>%
   ggplot(clean.df, aes(asinh.FSC.A, asinh.BL1.A)) +
   geom_hex(aes(fill = pop), bins = 300) + # ,alpha=..ncount.. #order= ?
   geom_density2d(color="black",  size = 0.1) +
   # xlim(c(5, 15)) + ylim(c(2.5, 15)) +
   # scale_fill_viridis(
   #    discrete = TRUE, end = 0.8, name = "", direction = -1,
   #    guide = FALSE, option = 'C'
   # ) +
   # scale_alpha_continuous(guide = FALSE) +
   theme_bw()+ 
   geom_point(aes(centers.list.df[1, 1], centers.list.df[1, 2]), col = "blue", size = 1) +
   geom_point(aes(centers.list.df[2, 1], centers.list.df[2, 2]), col = "blue", size = 1) +
   facet_wrap(~well)
# ggsave(paste0("fig/gate_plots/cluster_",fcsset[[i]]@description$`$SRC`,".png"), plot = p)
ggsave2(filename = paste0("cluster_",fcsset[[i]]@description$`$SRC`,".png" ), 
        plot = p,
        path = here("fig/gate_plots/", day))

#### Get the quantities ####


clust.count <- data.frame(clean.df, cluster = cluster.predict$classification) %>%
   group_by(well, cluster) %>%
   summarize(count = n()) 
for(i in pop.tbl$cluster){
   clust.count$cluster[clust.count$cluster %in%  pop.tbl$cluster[i]] <-  pop.tbl$pop[i]
}

df.stats <-    
   spread(clust.count ,cluster, count)%>%
   full_join(df.stats,.)

# assign a spore column for mutant to maintain compatability of count tables
if (mutant.host | (!"spore" %in% colnames(df.stats))){
   df.stats$spore <- 0
}

   
# calculate concentrations based on volume and dilution
df.stats$spore.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
                     df.stats$spore/(df.stats$volume.nL/1e6)
df.stats$veg.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
   df.stats$veg/(df.stats$volume.nL/1e6)  

# write results to file
write_csv(df.stats,
         file = here("data/output/",day,paste0(fcsset[[i]]@description$`$SRC`,".csv")))
  
print(paste("done",folder))
} #folder loop

