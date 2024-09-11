######################################################################################
########## Hierarchical spatial clustering of climate awareness survey 2021 ##########
##########        Zammarchi, G. & Maranzano, P. (2024+)                     ##########
######################################################################################

##### Libraries
library(tidyverse)
library(ggplot2)
library(ClustGeo)
library(sf)
library(ggpubr)
library(NbClust)
library(viridis)
library(moments)
library(openxlsx)
'%notin%' <- Negate('%in%')

##### Setup
cols <- c("Morelli et al. (2024)" = "#FF9933",
          "Chavent et al. (2018)" = "#0000FF")

##### Auxiliary functions
source("H:/Il mio Drive/SpatialClustering/ClusterGeoTS - Raffaele Mattera/AuxFuns - hclustgeo majority rule.R", encoding = 'UTF-8')

##### Folder
setwd("H:/Il mio Drive/Awareness_ZammarchiMaranzano/Analysis")

##### Load data
load("../Data/Dataset_finale.RData")
data <- dataset_red

##### Check on data structure
summary(data)
data %>%
  group_by(Code,Entity) %>%
  summarise(n = n()) %>%
  View(title = "Obs. per country")
data %>%
  group_by(Variable) %>%
  summarise(n = n()) %>%
  View(title = "Obs. per variable")


##### Extract polygons for World's countries
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_red <- world %>%
  filter(type %in% c("Sovereign country","Country","Sovereignty")) %>%
  select(Code = iso_a3_eh)

##### Add geometries to data
data_sf <- left_join(x = world_red, y = data, by = c("Code"))
data_sf <- st_as_sf(data_sf)




####################################################
########## Plot of climate awareness rate ##########
####################################################
Period <- "Avg_2014_2022"
p <- data_sf %>%
  select(Code,Entity,Variable,Period = any_of(Period)) %>%
  pivot_wider(names_from = Variable, values_from = Period) %>%
  ggplot() +
  geom_sf(mapping = aes(fill = R12)) +
  theme_minimal() +
  labs(title = "Share of people with medium-high or high climate awareness (2021)",
       subtitle = "'I know a moderate amount about it' = Medium-high\n'I know a lot about it' = High",
       fill = "Share % of the overall sample") +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 50) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 14),axis.text = element_text(size = 10))
ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = "Awareness_High.png")

p <- data_sf %>%
  select(Code,Entity,Variable,Period = any_of(Period)) %>%
  pivot_wider(names_from = Variable, values_from = Period) %>%
  ggplot() +
  geom_sf(mapping = aes(fill = R12)) +
  theme_minimal() +
  labs(title = "Share of people with medium-low or low climate awareness (2021)",
       subtitle = "'I have never heard of it' = Low\n'I know a little about it' = Medium-Low",
       fill = "Share % of the overall sample") +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 50) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 14),axis.text = element_text(size = 10))
ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = "Awareness_Low.png")










###########################################################################
#################### Main analysis: low awareness only ####################
###########################################################################
AwarenessOnly <- TRUE
Dissim_std <- FALSE

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40))

if (AwarenessOnly == TRUE) {
  data_wide <- data_wide %>%
    select(Code,Entity,geometry,R = R12)
}

##### Compute distance matrix of coordinates
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide)
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

########## Spatial clustering ##########
# Algorithm finding alpha and K
a_ch<-NULL #alpha proposed by Chavent 2018
a_max<-NULL #alpha from the maximum
k_max<-15
d_a<-0.01
wss<-matrix(data=NA,nrow=k_max,ncol=3)
q<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_ch<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_ch[1,]<-0
exp_in_max<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_max[1,]<-0
sil <- du <- ci <- mc <- ch <- matrix(data=NA,nrow=k_max-1,ncol=7)
# total inertia d0 d1
T1 <- inertdiss(as.dist(D0_norm))
T2 <- inertdiss(as.dist(D1_norm))
for (k in 2:k_max){
  # k <- 8
  
  cr <- choicealpha(as.dist(D0_norm),as.dist(D1_norm), range.alpha = seq(0, 1, d_a), K=k,
                    graph=T, scale=FALSE)
  a <- as.data.frame(cr$Qnorm)
  b <- as.data.frame(cr$Q)
  a$tot <- a$Q0norm + a$Q1norm
  b$tot <- (b$Q0*T1 + b$Q1*T2)/(T1 + T2)
  
  ##### Chavent et al. (2018): select alpha such that min{|Q1norm-Q0norm|}
  a_ch[k] <- cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  
  ##### Morelli et al. (2024): select alpha such that max{|Q1norm+Q0norm|}
  a_max[k] <- cr$range.alpha[which.max(b$tot)]
  
  exp_in_ch[k,1]<-b[a_ch[k]==cr$range.alpha,1]
  exp_in_max[k,1]<-b[a_max[k]==cr$range.alpha,1]
  exp_in_ch[k,2]<-b[a_ch[k]==cr$range.alpha,2]
  exp_in_max[k,2]<-b[a_max[k]==cr$range.alpha,2]
  print(k)
}
exp_in_ch[,3]<-(exp_in_ch[,1]*T1+exp_in_ch[,2]*T2)/(T1+T2)
exp_in_max[,3]<-(exp_in_max[,1]*T1+exp_in_max[,2]*T2)/(T1+T2)

##### Optimal alphas given initial set of K
val_a <- sort(unique(c(a_ch[2:k_max], a_max[2:k_max])))
MajRule <- vector(mode = "list", length = length(val_a))

##### Computating the set of indices for each AlphaStar_K
for(i in 1:length(val_a)){
  # i <- 1
  a <- val_a[i]
  
  # For each value of alpha we compute the combined matrix D, the tree from hier clu
  dissim <- (1-a)*(as.dist(D0_norm)) + a*(as.dist(D1_norm))
  if (Dissim_std == TRUE) {
    dissim <- dissim/max(dissim)
  }
  
  # hclustgeo given alpha_star = a for every k = 1,...,k_max (redundant)
  ind <- NbClust(data=cbind(data_wide_mat,D1_centr),diss=dissim,distance=NULL,
                 min.nc = 2, max.nc = k_max, method="ward.D",index="all")
  
  # Extract indicators only for the k associated with alpha_star = a
  # N.B. alpha_star changes for Chavent (ch) and Maximum (max)
  dum_ch <- trimws(a_ch[2:k_max]) == a
  dum_max <- trimws(a_max[2:k_max]) == a
  
  # Silhouette: Rousseeuw (1987)
  sil[dum_ch,1] <- "Sil"
  sil[dum_ch,2] <- a
  sil[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  sil[dum_ch,4] <- ind$All.index[dum_ch,"Silhouette"]
  sil[dum_max,5] <- a
  sil[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  sil[dum_max,7] <- ind$All.index[dum_max,"Silhouette"]
  # Dunn (1974)
  du[dum_ch,1] <- "Dunn"
  du[dum_ch,2] <- a
  du[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  du[dum_ch,4] <- ind$All.index[dum_ch,"Dunn"]
  du[dum_max,5] <- a
  du[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  du[dum_max,7] <- ind$All.index[dum_max,"Dunn"]
  # Hubert and Levin (1976)
  ci[dum_ch,1] <- "Cindex"
  ci[dum_ch,2] <- a
  ci[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  ci[dum_ch,4] <- ind$All.index[dum_ch,"Cindex"]
  ci[dum_max,5] <- a
  ci[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  ci[dum_max,7] <- ind$All.index[dum_max,"Cindex"]
  # Calinski and Harabasz (1974)
  ch[dum_ch,1] <- "CH"
  ch[dum_ch,2] <- a
  ch[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  ch[dum_ch,4] <- ind$All.index[dum_ch,"CH"]
  ch[dum_max,5] <- a
  ch[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  ch[dum_max,7] <- ind$All.index[dum_max,"CH"]
  # McClain and Rao (1975)
  mc[dum_ch,1] <- "McClain"
  mc[dum_ch,2] <- a
  mc[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  mc[dum_ch,4] <- ind$All.index[dum_ch,"McClain"]
  mc[dum_max,5] <- a
  mc[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  mc[dum_max,7] <- ind$All.index[dum_max,"McClain"]
}
colnames(sil) <- colnames(du) <- colnames(mc) <- colnames(ci) <- colnames(ch) <- c(
  "Index",
  "alpha_star_ch","K_ch","Index_ch",
  "alpha_star_max","K_max","Index_max"
)
sil <- data.frame(sil)
mc <- data.frame(mc)
du <- data.frame(du)
ci <- data.frame(ci)
ch <- data.frame(ch)

########## Analysis of the indices ##########
Indices <- bind_rows(sil,du,mc,ci,ch)

##### Plot of the selected indices
Indices_plot <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2024)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))
p2 <- Indices_plot %>%
  ggplot(mapping = aes(x = K, y = Index_val,col = Criterion,group = Criterion)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  geom_point(data = Indices_plot %>% filter(K == 4 & Index %in% c("C-index (Hubert and Levin, 1976) -- Min")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(4) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Calinski and Harabasz (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(5) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Calinski and Harabasz (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(14) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Dunn (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(8) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Dunn (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(2) & Index %in% c("McClain and Rao (1975) -- Min")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(3) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Silhouette (Rousseeuw, 1987) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(5) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Silhouette (Rousseeuw, 1987) -- Max")),
             size = 7, shape = 1, col = "red") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  scale_x_continuous(breaks = seq(from=2,to=15,by=1),limits = c(2,15)) +
  labs(title = "",
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("Index values at $\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
p <- ggarrange(p2,ncol = 1,align = "hv",common.legend = TRUE,legend = "bottom")
p <- annotate_figure(p,top = text_grob("Indices for selecting the optimal number of clusters",
                                       face = "bold",size = 25,vjust = -0.50))
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 1: indices are evaluated for a growing number of clusters at each optimal $\\alpha = \\alpha^*_K$"),
                                       vjust = +2.5,hjust = +0.6,size = 12))
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 2: $\\alpha = \\alpha^*_K$ are either determined according to Chavent et al. (2018) or Morelli et al. (2024)"),
                                       vjust = +5.0,hjust = +0.545,size = 12))
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = "R12/Indices_OptimK.png")

p1 <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Alpha = case_when(grepl("max",Alpha) ~ "Morelli et al. (2024)",
                           grepl("ch",Alpha) ~ "Chavent et al. (2018)",
                           TRUE ~ Alpha),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index)) %>%
  ggplot(mapping = aes(x = K, y = alpha_val,col = Alpha, group = Alpha)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = seq(from=2,to=15,by=1),limits = c(2,15)) +
  labs(title = latex2exp::TeX("Optimal $\\alpha$ values (i.e., $\\alpha^*_K$) conditioning on $K$ clusters"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("$\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
ggpubr::ggexport(p1,width = 3400, height = 2200,res = 250,filename = "R12/Optim_Alpha_K.png")


########## Comparison between indices at alpha_*_K and at alpha_K = 0
Indices <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2024)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))

ind0 <- matrix(data=NA,nrow=k_max-1,ncol=7)
a <- 0
dissim <- (1-a)*(as.dist(D0_norm)) + a*(as.dist(D1_norm))
if (Dissim_std == TRUE) {
  dissim <- dissim/max(dissim)
}
ind <- NbClust(data=cbind(data_wide_mat,D1_centr),diss=dissim,distance=NULL,
               min.nc = 2, max.nc = k_max, method="ward.D",index="all")
ind0[,1] <- a
ind0[,2] <- as.numeric(names(ind$All.index[,"Silhouette"]))
ind0[,3] <- ind$All.index[,"Silhouette"]
ind0[,4] <- ind$All.index[,"Dunn"]
ind0[,5] <- ind$All.index[,"Cindex"]
ind0[,6] <- ind$All.index[,"CH"]
ind0[,7] <- ind$All.index[,"McClain"]
colnames(ind0) <- c("Alpha","K","Sil","Dunn","Cindex","CH","McClain")

ind0 <- ind0 %>%
  as_tibble() %>%
  pivot_longer(cols = c("Sil","Dunn","Cindex","CH","McClain"),names_to = "Index",values_to = "Index_val0") %>%
  mutate(Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index)) %>%
  select(-c(Alpha))

Indices_augm <- full_join(x = Indices, y = ind0, by = c("K","Index"))
Indices_augm <- Indices_augm %>%
  mutate(Diff_Index0 = (Index_val - Index_val0)/Index_val0*100) %>%
  select(Index,K,Criterion,Diff_Index0) %>%
  filter(K %in% c(2,3,4,5,8,14))

up <- Indices_augm[Indices_augm$Diff_Index0 >= 0,]
down <- Indices_augm[Indices_augm$Diff_Index0 < 0,]


########## Plot of percentage gain or loss w.r.t. alpha=0
p <- ggplot() + 
  geom_bar(data = up, 
           mapping = aes(x = as.factor(K), y = Diff_Index0,
                         fill = Criterion,
                         group = Criterion),width = 0.5,
           position = "dodge", stat = "identity") + 
  geom_bar(data = down, 
           mapping = aes(x = as.factor(K), y = Diff_Index0,
                         fill = Criterion,
                         group = Criterion),width = 0.5,
           position = "dodge", stat = "identity") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  labs(title = latex2exp::TeX("Gain/loss of the indices evaluated at $\\alpha^*_K$ w.r.t. the indices evaluated at $\\alpha_K = 0$"),
       subtitle = latex2exp::TeX("Gain/loss is computed as $GL = \\frac{Index_{\\alpha^*_K} - Index_{\\alpha_K = 0}$}{Index_{\\alpha_K = 0}} \\times 100"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("% increase or decrease")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) + 
  scale_fill_manual("Criterion",values = cols)
ggpubr::ggexport(p,width = 3000, height = 2000,res = 210,filename = "R12/GainLoss_wrt_alpha0.png")

##### Optimal combinations
# Comb_optim <- Indices %>%
#   select(K,alpha_val) %>%
#   mutate(K = as.numeric(K),
#          alpha_val = as.numeric(alpha_val)) %>%
#   filter(K %in% c(2,3,4,5,8,14)) %>%
#   distinct()

Comb_optim <- rbind(
  # K = 4 and alpha = alpha_max
  c(4,0.60),
  # K = 5 and alpha = alpha_max
  c(5,0.62),
  # K = 4 and alpha = 0
  c(4,0),
  # K = 5 and alpha = 0
  c(5,0)
)
Comb_optim <- as.data.frame(Comb_optim)

clustering <- clustering0 <- list()
ClusterNames <- character(length = dim(Comb_optim)[1])
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]

  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt=NULL, scale=FALSE)
  clustering[[i]] <- cutree(hc, K_star)
  ClusterNames[i] <- paste0("Cluster_K=",K_star,"Alpha=",alpha_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt=NULL, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  # names(clustering0)[i] <- paste0("Cluster_K=",K_star,"Alpha=0")
  
  data_wide <- bind_cols(data_wide,clustering[[i]])

  Period <- "Avg_2014_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("R12/Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}
colnames(data_wide)[5:dim(data_wide)[2]] <- ClusterNames

##### Coherence among clustering methods: adjusted Rank index
Rand_mat <- matrix(data = NA, nrow = dim(Comb_optim)[1], ncol = dim(Comb_optim)[1])
CellName <- character(length = dim(Comb_optim)[1])
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  CellName[i] <- paste0("K=",K_star," & Alpha=",alpha_star)
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  for (j in 1:dim(Comb_optim)[1]) {
    stats <- fpc::cluster.stats(DistMat, clustering = clustering[[i]], alt.clustering = clustering[[j]])
    Rand_mat[i,j] <- stats$corrected.rand
  }
}
colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
Rand_mat

##### Groups characterization
stat1 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=4Alpha=0.6`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat2 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=4Alpha=0`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat3 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=5Alpha=0.62`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat4 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=5Alpha=0`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))


### Export GreenPeace dataset
##### Creo file Excel
wb <- createWorkbook("ClusterAwareness")
##### Save Excel sheet
addWorksheet(wb,ClusterNames[1])
writeData(wb, sheet = ClusterNames[1], stat1, colNames = T)
addWorksheet(wb,ClusterNames[2])
writeData(wb, sheet = ClusterNames[2], stat1, colNames = T)
addWorksheet(wb,ClusterNames[3])
writeData(wb, sheet = ClusterNames[3], stat1, colNames = T)
addWorksheet(wb,ClusterNames[4])
writeData(wb, sheet = ClusterNames[4], stat1, colNames = T)
##### Salvatggio file Excel
saveWorkbook(wb,"R12/ClusterStats_Awareness.xlsx",overwrite = T)










#############################################################################
#################### Robustness analysis: high awareness ####################
#############################################################################
AwarenessOnly <- TRUE
Dissim_std <- FALSE

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40))

if (AwarenessOnly == TRUE) {
  data_wide <- data_wide %>%
    select(Code,Entity,geometry,R = R34)
}

##### Compute distance matrix of coordinates
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide)
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

########## Spatial clustering ##########
# Algorithm finding alpha and K
a_ch<-NULL #alpha proposed by Chavent 2018
a_max<-NULL #alpha from the maximum
k_max<-15
d_a<-0.01
wss<-matrix(data=NA,nrow=k_max,ncol=3)
q<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_ch<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_ch[1,]<-0
exp_in_max<-matrix(data=NA,nrow=k_max,ncol=3)
exp_in_max[1,]<-0
sil <- du <- ci <- mc <- ch <- matrix(data=NA,nrow=k_max-1,ncol=7)
# total inertia d0 d1
T1 <- inertdiss(as.dist(D0_norm))
T2 <- inertdiss(as.dist(D1_norm))
for (k in 2:k_max){
  # k <- 8
  
  cr <- choicealpha(as.dist(D0_norm),as.dist(D1_norm), range.alpha = seq(0, 1, d_a), K=k,
                    graph=T, scale=FALSE)
  a <- as.data.frame(cr$Qnorm)
  b <- as.data.frame(cr$Q)
  a$tot <- a$Q0norm + a$Q1norm
  b$tot <- (b$Q0*T1 + b$Q1*T2)/(T1 + T2)
  
  ##### Chavent et al. (2018): select alpha such that min{|Q1norm-Q0norm|}
  a_ch[k] <- cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  
  ##### Morelli et al. (2024): select alpha such that max{|Q1norm+Q0norm|}
  a_max[k] <- cr$range.alpha[which.max(b$tot)]
  
  exp_in_ch[k,1]<-b[a_ch[k]==cr$range.alpha,1]
  exp_in_max[k,1]<-b[a_max[k]==cr$range.alpha,1]
  exp_in_ch[k,2]<-b[a_ch[k]==cr$range.alpha,2]
  exp_in_max[k,2]<-b[a_max[k]==cr$range.alpha,2]
  print(k)
}
exp_in_ch[,3]<-(exp_in_ch[,1]*T1+exp_in_ch[,2]*T2)/(T1+T2)
exp_in_max[,3]<-(exp_in_max[,1]*T1+exp_in_max[,2]*T2)/(T1+T2)

##### Optimal alphas given initial set of K
val_a <- sort(unique(c(a_ch[2:k_max], a_max[2:k_max])))
MajRule <- vector(mode = "list", length = length(val_a))

##### Computating the set of indices for each AlphaStar_K
for(i in 1:length(val_a)){
  # i <- 1
  a <- val_a[i]
  
  # For each value of alpha we compute the combined matrix D, the tree from hier clu
  dissim <- (1-a)*(as.dist(D0_norm)) + a*(as.dist(D1_norm))
  if (Dissim_std == TRUE) {
    dissim <- dissim/max(dissim)
  }
  
  # hclustgeo given alpha_star = a for every k = 1,...,k_max (redundant)
  ind <- NbClust(data=cbind(data_wide_mat,D1_centr),diss=dissim,distance=NULL,
                 min.nc = 2, max.nc = k_max, method="ward.D",index="all")
  
  # Extract indicators only for the k associated with alpha_star = a
  # N.B. alpha_star changes for Chavent (ch) and Maximum (max)
  dum_ch <- trimws(a_ch[2:k_max]) == a
  dum_max <- trimws(a_max[2:k_max]) == a
  
  # Silhouette: Rousseeuw (1987)
  sil[dum_ch,1] <- "Sil"
  sil[dum_ch,2] <- a
  sil[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  sil[dum_ch,4] <- ind$All.index[dum_ch,"Silhouette"]
  sil[dum_max,5] <- a
  sil[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  sil[dum_max,7] <- ind$All.index[dum_max,"Silhouette"]
  # Dunn (1974)
  du[dum_ch,1] <- "Dunn"
  du[dum_ch,2] <- a
  du[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  du[dum_ch,4] <- ind$All.index[dum_ch,"Dunn"]
  du[dum_max,5] <- a
  du[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  du[dum_max,7] <- ind$All.index[dum_max,"Dunn"]
  # Hubert and Levin (1976)
  ci[dum_ch,1] <- "Cindex"
  ci[dum_ch,2] <- a
  ci[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  ci[dum_ch,4] <- ind$All.index[dum_ch,"Cindex"]
  ci[dum_max,5] <- a
  ci[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  ci[dum_max,7] <- ind$All.index[dum_max,"Cindex"]
  # Calinski and Harabasz (1974)
  ch[dum_ch,1] <- "CH"
  ch[dum_ch,2] <- a
  ch[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  ch[dum_ch,4] <- ind$All.index[dum_ch,"CH"]
  ch[dum_max,5] <- a
  ch[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  ch[dum_max,7] <- ind$All.index[dum_max,"CH"]
  # McClain and Rao (1975)
  mc[dum_ch,1] <- "McClain"
  mc[dum_ch,2] <- a
  mc[dum_ch,3] <- which(trimws(a_ch[2:k_max]) == a) + 1
  mc[dum_ch,4] <- ind$All.index[dum_ch,"McClain"]
  mc[dum_max,5] <- a
  mc[dum_max,6] <- which(trimws(a_max[2:k_max]) == a) + 1
  mc[dum_max,7] <- ind$All.index[dum_max,"McClain"]
}
colnames(sil) <- colnames(du) <- colnames(mc) <- colnames(ci) <- colnames(ch) <- c(
  "Index",
  "alpha_star_ch","K_ch","Index_ch",
  "alpha_star_max","K_max","Index_max"
)
sil <- data.frame(sil)
mc <- data.frame(mc)
du <- data.frame(du)
ci <- data.frame(ci)
ch <- data.frame(ch)

########## Analysis of the indices ##########
Indices <- bind_rows(sil,du,mc,ci,ch)

##### Plot of the selected indices
Indices_plot <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2024)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))
p2 <- Indices_plot %>%
  ggplot(mapping = aes(x = K, y = Index_val,col = Criterion,group = Criterion)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  geom_point(data = Indices_plot %>% filter(K %in% c(5) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("C-index (Hubert and Levin, 1976) -- Min")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(15) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("C-index (Hubert and Levin, 1976) -- Min")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(3) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Calinski and Harabasz (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(4) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Calinski and Harabasz (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(3) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Dunn (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(4) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Dunn (1974) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(2) & Index %in% c("McClain and Rao (1975) -- Min")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(3) & Criterion %in% c("Chavent et al. (2018)") & Index %in% c("Silhouette (Rousseeuw, 1987) -- Max")),
             size = 7, shape = 1, col = "red") + 
  geom_point(data = Indices_plot %>% filter(K %in% c(4) & Criterion %in% c("Morelli et al. (2024)") & Index %in% c("Silhouette (Rousseeuw, 1987) -- Max")),
             size = 7, shape = 1, col = "red") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  scale_x_continuous(breaks = seq(from=2,to=15,by=1),limits = c(2,15)) +
  labs(title = "",
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("Index values at $\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
p <- ggarrange(p2,ncol = 1,align = "hv",common.legend = TRUE,legend = "bottom")
p <- annotate_figure(p,top = text_grob("Indices for selecting the optimal number of clusters",
                                       face = "bold",size = 25,vjust = -0.50))
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 1: indices are evaluated for a growing number of clusters at each optimal $\\alpha = \\alpha^*_K$"),
                                       vjust = +2.5,hjust = +0.6,size = 12))
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 2: $\\alpha = \\alpha^*_K$ are either determined according to Chavent et al. (2018) or Morelli et al. (2024)"),
                                       vjust = +5.0,hjust = +0.545,size = 12))
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = "R34/Indices_OptimK.png")

p1 <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Alpha = case_when(grepl("max",Alpha) ~ "Morelli et al. (2024)",
                           grepl("ch",Alpha) ~ "Chavent et al. (2018)",
                           TRUE ~ Alpha),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index)) %>%
  ggplot(mapping = aes(x = K, y = alpha_val,col = Alpha, group = Alpha)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  scale_x_continuous(breaks = seq(from=2,to=15,by=1),limits = c(2,15)) +
  labs(title = latex2exp::TeX("Optimal $\\alpha$ values (i.e., $\\alpha^*_K$) conditioning on $K$ clusters"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("$\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
ggpubr::ggexport(p1,width = 3400, height = 2200,res = 250,filename = "R34/Optim_Alpha_K.png")


########## Comparison between indices at alpha_*_K and at alpha_K = 0
Indices <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2024)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))

ind0 <- matrix(data=NA,nrow=k_max-1,ncol=7)
a <- 0
dissim <- (1-a)*(as.dist(D0_norm)) + a*(as.dist(D1_norm))
if (Dissim_std == TRUE) {
  dissim <- dissim/max(dissim)
}
ind <- NbClust(data=cbind(data_wide_mat,D1_centr),diss=dissim,distance=NULL,
               min.nc = 2, max.nc = k_max, method="ward.D",index="all")
ind0[,1] <- a
ind0[,2] <- as.numeric(names(ind$All.index[,"Silhouette"]))
ind0[,3] <- ind$All.index[,"Silhouette"]
ind0[,4] <- ind$All.index[,"Dunn"]
ind0[,5] <- ind$All.index[,"Cindex"]
ind0[,6] <- ind$All.index[,"CH"]
ind0[,7] <- ind$All.index[,"McClain"]
colnames(ind0) <- c("Alpha","K","Sil","Dunn","Cindex","CH","McClain")

ind0 <- ind0 %>%
  as_tibble() %>%
  pivot_longer(cols = c("Sil","Dunn","Cindex","CH","McClain"),names_to = "Index",values_to = "Index_val0") %>%
  mutate(Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index)) %>%
  select(-c(Alpha))

Indices_augm <- full_join(x = Indices, y = ind0, by = c("K","Index"))
Indices_augm <- Indices_augm %>%
  mutate(Diff_Index0 = (Index_val - Index_val0)/Index_val0*100) %>%
  select(Index,K,Criterion,Diff_Index0) %>%
  filter(K %in% c(2,3,4,5,15))

up <- Indices_augm[Indices_augm$Diff_Index0 >= 0,]
down <- Indices_augm[Indices_augm$Diff_Index0 < 0,]


########## Plot of percentage gain or loss w.r.t. alpha=0
p <- ggplot() + 
  geom_bar(data = up, 
           mapping = aes(x = as.factor(K), y = Diff_Index0,
                         fill = Criterion,
                         group = Criterion),width = 0.5,
           position = "dodge", stat = "identity") + 
  geom_bar(data = down, 
           mapping = aes(x = as.factor(K), y = Diff_Index0,
                         fill = Criterion,
                         group = Criterion),width = 0.5,
           position = "dodge", stat = "identity") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  labs(title = latex2exp::TeX("Gain/loss of the indices evaluated at $\\alpha^*_K$ w.r.t. the indices evaluated at $\\alpha_K = 0$"),
       subtitle = latex2exp::TeX("Gain/loss is computed as $GL = \\frac{Index_{\\alpha^*_K} - Index_{\\alpha_K = 0}$}{Index_{\\alpha_K = 0}} \\times 100"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("% increase or decrease")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) + 
  scale_fill_manual("Criterion",values = cols)
ggpubr::ggexport(p,width = 3000, height = 2000,res = 210,filename = "R34/GainLoss_wrt_alpha0.png")

##### Optimal combinations
# Comb_optim <- Indices %>%
#   select(K,alpha_val) %>%
#   mutate(K = as.numeric(K),
#          alpha_val = as.numeric(alpha_val)) %>%
#   filter(K %in% c(2,3,4,5,8,14)) %>%
#   distinct()

Comb_optim <- rbind(
  # K = 4 and alpha = alpha_max
  c(4,0.60),
  # K = 5 and alpha = alpha_max
  c(5,0.62),
  # K = 4 and alpha = 0
  c(4,0),
  # K = 5 and alpha = 0
  c(5,0)
)
Comb_optim <- as.data.frame(Comb_optim)

clustering <- clustering0 <- list()
ClusterNames <- character(length = dim(Comb_optim)[1])
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt=NULL, scale=FALSE)
  clustering[[i]] <- cutree(hc, K_star)
  ClusterNames[i] <- paste0("Cluster_K=",K_star,"Alpha=",alpha_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt=NULL, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  # names(clustering0)[i] <- paste0("Cluster_K=",K_star,"Alpha=0")
  
  data_wide <- bind_cols(data_wide,clustering[[i]])
  
  Period <- "Avg_2014_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("R34/Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}
colnames(data_wide)[5:dim(data_wide)[2]] <- ClusterNames

##### Coherence among clustering methods: adjusted Rank index
Rand_mat <- matrix(data = NA, nrow = dim(Comb_optim)[1], ncol = dim(Comb_optim)[1])
CellName <- character(length = dim(Comb_optim)[1])
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  CellName[i] <- paste0("K=",K_star," & Alpha=",alpha_star)
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  for (j in 1:dim(Comb_optim)[1]) {
    stats <- fpc::cluster.stats(DistMat, clustering = clustering[[i]], alt.clustering = clustering[[j]])
    Rand_mat[i,j] <- stats$corrected.rand
  }
}
colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
Rand_mat

##### Groups characterization
stat1 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=4Alpha=0.6`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat2 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=4Alpha=0`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat3 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=5Alpha=0.62`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))
stat4 <- data_wide %>%
  as_tibble() %>%
  group_by(`Cluster_K=5Alpha=0`) %>%
  summarise(Mean = mean(R,na.rm=T),
            SD = sd(R,na.rm=T),
            Skew = skewness(R,na.rm=T),
            Kurt = kurtosis(R,na.rm=T))


### Export GreenPeace dataset
##### Creo file Excel
wb <- createWorkbook("ClusterAwareness")
##### Save Excel sheet
addWorksheet(wb,ClusterNames[1])
writeData(wb, sheet = ClusterNames[1], stat1, colNames = T)
addWorksheet(wb,ClusterNames[2])
writeData(wb, sheet = ClusterNames[2], stat2, colNames = T)
addWorksheet(wb,ClusterNames[3])
writeData(wb, sheet = ClusterNames[3], stat3, colNames = T)
addWorksheet(wb,ClusterNames[4])
writeData(wb, sheet = ClusterNames[4], stat4, colNames = T)
##### Salvatggio file Excel
saveWorkbook(wb,"R34/ClusterStats_Awareness.xlsx",overwrite = T)










##########################################################################
#################### Robustness analysis: correlation ####################
##########################################################################
# https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-hierarchical-clustering

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40))

for (rc in c(TRUE,FALSE)) {
  RedCorr <- rc
  ##### Compute Pearson's linear correlation
  corr <- data_wide %>%
    sf::st_drop_geometry() %>%
    select(-c(Code,Entity)) %>%
    rename('Medium-high and high climate awareness' = R34,
           'Medium-low and low climate awareness' = R12) %>%
    cor() %>%
    reshape2::melt(na.rm = TRUE) %>%
    rename(Reg1 = Var1, Reg2 = Var2, Corr2021 = value)
  corrmat <- corr %>%
    pivot_longer(cols = 3:last_col()) %>%
    # Standardize for fill/color scaling in tile
    group_by(name) %>%
    mutate(value_std = scale(x = value, center = T, scale = T)) %>%
    ungroup()
  if (RedCorr == TRUE) {
    corrmat <- corrmat %>%
      filter(Reg1 %in% c("Medium-high and high climate awareness",
                         "Medium-low and low climate awareness"))
  }
  p_corr <- corrmat %>%
    ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
    geom_tile(color = "white") +
    scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) + 
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() + 
    geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "",
      plot.title = element_text(size = 20, face = "bold",hjust = +0.75)
    ) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) + 
    labs(title = "Linear correlation among awareness, socio-economic and climate-related variables (average 2014-2021)")
  if (RedCorr == TRUE) {
    ggexport(p_corr,width = 2400, height = 2000, res = 150, filename = "CorrMat_red.png")
  } else {
    ggexport(p_corr,width = 2400, height = 2000, res = 150, filename = "CorrMat.png")
  }
}








######################################################################################################################
#################### Robustness analysis: awareness, socio-economic and climate-related variables ####################
######################################################################################################################

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40))

##### Compute distance matrix of coordinates
data_wide <- data_wide %>%
  select(-c(R34),
         -c(yr_sch,rnna,HCI,TerritorialEmiss_tCO2PC),
         -c(cwd,r50mm,rx5day,hi41,hi35,tr29),
         -contains("Share")) 
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide)
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

Comb_optim <- rbind(
  # K = 4 and alpha = alpha_max
  c(4,0.60),
  # K = 5 and alpha = alpha_max
  c(5,0.62),
  # K = 4 and alpha = 0
  c(4,0),
  # K = 5 and alpha = 0
  c(5,0)
)
Comb_optim <- as.data.frame(Comb_optim)

clustering_AwareClimateSocio <- clustering0 <- list()
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt=NULL, scale=FALSE)
  clustering_AwareClimateSocio[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt=NULL, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareClimateSocio[[i]])
  Period <- "Avg_2014_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareClimateSocio[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Robustness_AwareClimateSocio_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}

##### Coherence among clustering methods: adjusted Rank index
Rand_mat <- matrix(data = NA, nrow = dim(Comb_optim)[1], ncol = dim(Comb_optim)[1])
CellName <- character(length = dim(Comb_optim)[1])
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  CellName[i] <- paste0("K=",K_star," & Alpha=",alpha_star)
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  for (j in 1:dim(Comb_optim)[1]) {
    stats <- fpc::cluster.stats(DistMat, clustering = clustering[[i]], alt.clustering = clustering[[j]])
    Rand_mat[i,j] <- stats$corrected.rand
  }
}
colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
Rand_mat





######################################################################################################
#################### Robustness analysis: awareness and climate-related variables ####################
######################################################################################################

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40))

##### Compute distance matrix of coordinates
data_wide <- data_wide %>%
  select(-c(R34),
         -c(yr_sch,rnna,HCI,TerritorialEmiss_tCO2PC,CarbonIntens_Electr,EnerIntens_PrimEnergy,
            HDI,csh_g,rgdpna,EmpRate,TerritorialEmiss_IntensGDP_KgThs,TradeOpenness),
         -c(cwd,r50mm,rx5day,hi41,hi35,tr29),
         -contains("Share")) 
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide)
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

Comb_optim <- rbind(
  # K = 4 and alpha = alpha_max
  c(4,0.60),
  # K = 5 and alpha = alpha_max
  c(5,0.62),
  # K = 4 and alpha = 0
  c(4,0),
  # K = 5 and alpha = 0
  c(5,0)
)
Comb_optim <- as.data.frame(Comb_optim)

clustering_AwareClimate <- clustering0 <- list()
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt=NULL, scale=FALSE)
  clustering_AwareClimate[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt=NULL, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareClimate[[i]])
  Period <- "Avg_2014_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareClimate[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Robustness_AwareClimateSocio_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}





###########################################################################################
#################### Robustness analysis: awareness and socio-economic ####################
###########################################################################################

##### Dataset in wide form: 2014-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2014_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2014_2022) %>%
  select(-c(hd40)) %>%
  mutate(TerritorialEmiss_tCO2PC = TerritorialEmiss_tCO2PC*1000)

##### Compute distance matrix of coordinates
data_wide <- data_wide %>%
  select(-c(R34),
         -c(yr_sch,rnna,HCI,TerritorialEmiss_tCO2PC),
         -c(cwd,r50mm,rx5day,hi41,hi35,tr29,cdd,hd30,pr,tas,tx84rr,wsdi),
         -contains("Share")) 
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide)
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

Comb_optim <- rbind(
  # K = 4 and alpha = alpha_max
  c(4,0.60),
  # K = 5 and alpha = alpha_max
  c(5,0.62),
  # K = 4 and alpha = 0
  c(4,0),
  # K = 5 and alpha = 0
  c(5,0)
)
Comb_optim <- as.data.frame(Comb_optim)

clustering_AwareSocio <- clustering0 <- list()
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt=NULL, scale=FALSE)
  clustering_AwareSocio[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt=NULL, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareSocio[[i]])
  Period <- "Avg_2014_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareSocio[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Robustness_AwareSocio_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}





###########################################################################################################################
#################### Robustness analysis: Coherence among clustering methods using adjusted Rank index ####################
###########################################################################################################################

Rand_mat_Rob1 <- Rand_mat_Rob2 <- Rand_mat_Rob3 <- Rand_mat_Main <- matrix(data = NA, nrow = length(clustering[[1]]), ncol = 6)
for (i in 1:length(clustering_AwareSocio)) {
  Rand_mat_Rob1[,i+2] <- clustering_AwareClimateSocio[[i]]
}
Rand_mat_Rob1[,1] <- "AwareClimateSocio"
Rand_mat_Rob1[,2] <- 1:length(clustering[[1]])
Rand_mat_Rob1 <- data.frame(Rand_mat_Rob1)
for (i in 1:length(clustering_AwareSocio)) {
  Rand_mat_Rob2[,i+2] <- clustering_AwareClimate[[i]]
}
Rand_mat_Rob2[,1] <- "AwareClimate"
Rand_mat_Rob2[,2] <- 1:length(clustering[[1]])
Rand_mat_Rob2 <- data.frame(Rand_mat_Rob2)
for (i in 1:length(clustering_AwareClimate)) {
  Rand_mat_Rob3[,i+2] <- clustering_AwareSocio[[i]]
}
Rand_mat_Rob3[,1] <- "AwareSocio"
Rand_mat_Rob3[,2] <- 1:length(clustering[[1]])
Rand_mat_Rob3 <- data.frame(Rand_mat_Rob3)
for (i in 1:length(clustering_AwareClimate)) {
  Rand_mat_Main[,i+2] <- clustering[[i]]
}
Rand_mat_Main[,1] <- "Awareness only (main)"
Rand_mat_Main[,2] <- 1:length(clustering[[1]])
Rand_mat_Main <- data.frame(Rand_mat_Main)

for (i in 1:dim(Comb_optim)[1]) {
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  colnames(Rand_mat_Rob1)[i+2] <- colnames(Rand_mat_Rob2)[i+2] <- colnames(Rand_mat_Rob3)[i+2] <- colnames(Rand_mat_Main)[i+2] <- paste0("K=",K_star," & Alpha=",alpha_star)
}
colnames(Rand_mat_Rob1)[1] <- colnames(Rand_mat_Rob2)[1] <- colnames(Rand_mat_Rob3)[1] <- colnames(Rand_mat_Main)[1] <- "Scenario"
colnames(Rand_mat_Rob1)[2] <- colnames(Rand_mat_Rob2)[2] <- colnames(Rand_mat_Rob3)[2] <- colnames(Rand_mat_Main)[2] <- "Idx_cnt"

RobustnessClustering <- bind_rows(Rand_mat_Rob1,Rand_mat_Rob2,Rand_mat_Rob3,Rand_mat_Main)
RobustnessClustering <- RobustnessClustering %>%
  pivot_longer(cols = 3:last_col(), names_to = "Combination", values_to = "Cluster") %>%
  mutate(Key = paste0(Scenario," - ",Combination)) %>%
  dplyr::select(Key,Idx_cnt,Cluster) %>%
  pivot_wider(names_from = c(Key),values_from = Cluster) %>%
  select(-c(Idx_cnt))

CellName <- colnames(RobustnessClustering)
Rand_mat <- matrix(data = NA, nrow = length(CellName), ncol = length(CellName))
alpha_star <- 0
DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
for (i in 1:length(CellName)) {
  for (j in 1:length(CellName)) {
    stats <- fpc::cluster.stats(DistMat,
                                clustering = as.numeric(pull(RobustnessClustering[,i])),
                                alt.clustering = as.numeric(pull(RobustnessClustering[,j])))
    Rand_mat[i,j] <- stats$corrected.rand
  }
}
colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
Rand_mat

### Export adjusted rand indices
##### Creo file Excel
wb <- createWorkbook("RobustnessCoherence")
##### Save Excel sheet
addWorksheet(wb,"AdjRandIndex")
writeData(wb, sheet = "AdjRandIndex", Rand_mat, colNames = T)
##### Salvatggio file Excel
saveWorkbook(wb,"Robustness/RobustnessCoherence.xlsx",overwrite = T)


### Plot adjusted rand indices (corr matrix)
corr <- Rand_mat %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, AdjRanIndex = value)
corrmat <- corr %>%
  pivot_longer(cols = 3:last_col()) %>%
  # Standardize for fill/color scaling in tile
  group_by(name) %>%
  mutate(value_std = scale(x = value, center = T, scale = T)) %>%
  ungroup()
p_corr <- corrmat %>%
  ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "",
    plot.title = element_text(size = 20, face = "bold",hjust = +0.85)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + 
  labs(title = "Adjusted Rand Index for the main clusterings and the robustness clustering")
ggexport(p_corr,width = 2400, height = 2000, res = 200, filename = "Robustness/AdjrandIndex.png")
