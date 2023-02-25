## Coevolution with a seed bank - analysis of sporulation in evolved lines
library(renv)
# renv::init()
# renv::restore()
library(here)
#analysis of flow-cytometry population data using Karava's analysis.R as reference.
source(here("src/functions.R"))
source(here("src/DS_functions.R"))
library(multimode)
set.seed(1)
#choose data
   folders <- 
   list.dirs(here("data/FCM/fcs/"),recursive = F)

   if (! dir.exists(here("data/output"))){
     dir.create(here("data/output"), recursive = T)
   }
   
#fields to parse from the .fcs file name (separated by underscore and/or hyphen)
sample.var <- c("sample.date","exp.plate","t.exp","bio.rep","dilution","well","rep","xt","num") 

for (folder in folders){
  
# get data from sample name and set up folders
  # "day" in this experiment is a single sample with 3 technical reps.
  day = str_remove(folder, ".*/")
  
  # dilution of culture analyzed in FCM, extracted from sample name:
  dilution <- str_extract(day, "x[\\d]*")
  
  # Make directory to store the plots
  if (! dir.exists(here("fig/gate_plots", day))){
    dir.create(here("fig/gate_plots", day), recursive = T)
  }
  

  
  
#### Load data, sample set ####
fcsset <- flowCreateFlowSet(filepath = folder, sample_variables = sample.var,
                            transformation = FALSE,separators = "[-_\\.]")
#transform with arcsine, recommended by Karava et al.
fcsset <- Transform.Novocyte(fcsset)

#data frame to collect stats
df.stats <- fcsset%>%
   flowFcsToDf(.)%>%
   select(., all_of(sample.var))%>%
   distinct()
df.stats$dilution <- dilution
df.stats$total.events <- NA
df.stats$volume.nL <- NA
df.stats$limit.events <- NA
df.stats$limit.volume.uL <- NA
df.stats$singlets <- NA
df.stats$neg.rmv <- NA
df.stats$noise.cutoff <- NA
df.stats$sybr.modes<- NA
df.stats$sybr.antimode<- NA


#### Gating for singlets with flowStats ####

# The gate function needs to be applied to each sample seperatly
# get number of samples
n.sample <- nrow(fcsset@phenoData@data)



#initialize list to store singlet plots
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
      
      #find cutoff to remove noise by SSC
      noise<- rangeGate(x = fcsset[[i]], "asinh.SSC.A", alpha=0.1,sd=3, absolute = F)
      df.stats$noise.cutoff[i] <- noise@min[[1]]

}

#save plot
ggsave2(filename = paste0("singlets_",fcsset[[i]]@description$`$SRC`,".pdf" ), 
        plot = plot_grid(plotlist = plot.list),
        path = here("fig/gate_plots/", day))


#### Gating for singlets with norm2Filter ####
# df.set <- fcsset %>%
#    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.A", scale.factor = 10)) %>%
#    Subset(rectangleGate("asinh.BL1.A" = c(0, 15))) %>%
#    flowFcsToDf(.)

# plotting noise gate

noise.plot <-   
# ggcyto(fcsset,aes(asinh.FSC.A))+
  ggcyto(fcsset,aes(asinh.SSC.A))+
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

complotFSC <- 
   ggplot(df.set, aes(asinh.FSC.A, asinh.BL1.A))+
   geom_point(size=0.1, color="grey")+
   geom_point(data=clean.df,size=0.1, color="blue")+
   geom_density2d(col = "white",  size = 0.1, alpha = 0.5) +
   scale_alpha_continuous(guide = FALSE) +
   facet_wrap(~well)+
   theme_cowplot()+
   panel_border()

complotSSC <- 
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
        plot = plot_grid(complotFSC, complotSSC, ncol = 1),
        path = here("fig/gate_plots/", day), height = 4, width = 4)




#### predict centers of sub-populations ------------------------------------


# number of sub-populations ##

# test for number of modes on SYBR
# If there is only one mode assume all veg

#list to collect plots
p.list <- list()
for (wl in wells){
 # test if uni- or multi- modal
   unimodal <- 
    clean.df %>% 
    dplyr::filter(well==wl) %>% 
    pull(asinh.BL1.A) %>% 
    diptest::dip.test(.)
  unimodal <- unimodal$p.value >0.05
  df.stats$sybr.modes[df.stats$well==wl] <- if_else(unimodal, 1,2)
  
  modes <- 
    clean.df %>% 
    dplyr::filter(well==wl) %>% 
    pull(asinh.BL1.A) %>%
  locmodes(., display = F, mod0 = if_else(unimodal, 1,2))
  
  # SYBR population split
  df.stats$sybr.antimode[df.stats$well==wl] <- 
    if_else(unimodal, NA_real_,modes$locations[2])
  
  vline.col = "red"
  if(!unimodal) vline.col =  c("red","blue","red")
  
  p.list[[wl]] <- 
    clean.df %>% 
    dplyr::filter(well==wl) %>% 
    ggplot(aes(asinh.BL1.A))+
    geom_histogram()+
    geom_vline(xintercept = modes$locations, color =vline.col)+
    theme_bw()+
    ggtitle(wl) 
  
}

ggsave(here("fig/gate_plots/", day,paste0("sybrModes_",day,".png" )), 
       plot = plot_grid(plotlist=p.list, ncol = 1),
       width = 5, height = 6)

rm(p.list)


# categorize spore-veg ----------------------------------------------------

# separate subpopulations by SYBR -----------------------------------------
sybr.df <- clean.df[0,]

# assign to sub population
wells <- levels(as.factor(df.stats$well))
for (wl in wells){
  sybr.df <- 
    rbind(sybr.df,
          clean.df %>%
            dplyr::filter(well==wl)%>%
            dplyr::mutate(sub.pop = 
              case_when(
                # in case of one SYBR population assign all as veg 
                # (all spores is rare, unless there is no stain)
                (df.stats$sybr.modes[df.stats$well==wl] == 1) ~ "veg",
                # is event SYBR higher than threshold than assign veg
                asinh.BL1.A>df.stats$sybr.antimode[df.stats$well==wl] ~ "veg",
                # whatever is left is a spore.
                TRUE ~ "spore"
              )
            )
          )
}



complot <- 
  ggplot(sybr.df, aes(asinh.SSC.A, asinh.BL1.A))+
  geom_point(aes(color = sub.pop), size=0.1)+
  geom_density2d(col = "white",  size = 0.1, alpha = 0.5) +
  scale_alpha_continuous(guide = FALSE) +
  facet_wrap(~well)+
  theme_cowplot()+
  panel_border()+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3) ) )

# ggsave2(paste0("fig/gate_plots/scatterNoise_",fcsset[[i]]@description$`$SRC`,".png"), complot)
ggsave2(filename = paste0("scatterSYBR_",fcsset[[i]]@description$`$SRC`,".png" ), 
        plot = complot,
        path = here("fig/gate_plots/", day))

#### Get the quantities ####

df.stats <- 
  sybr.df %>%
  count(well, sub.pop) %>% 
  pivot_wider(names_from = sub.pop, values_from = n) %>% 
  full_join(df.stats,.)

# add spores = 0 to veg only populations
if(is.null(df.stats$spore)){
  df.stats$spore <- 0
}
   
# calculate concentrations based on volume and dilution
df.stats$spore.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
                     df.stats$spore/(df.stats$volume.nL/1e6)
df.stats$veg.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
   df.stats$veg/(df.stats$volume.nL/1e6)  

# write results to file
write_csv(df.stats,
         path = here("data/output/",paste0(fcsset[[i]]@description$`$SRC`,".csv")))
  
print(paste("done",folder))
} #folder loop

