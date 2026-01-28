######################################################################################
########## Hierarchical spatial clustering of climate awareness survey 2022 ##########
##########        Zammarchi, G. & Maranzano, P. (2025+)                     ##########
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
library(gt)
'%notin%' <- Negate('%in%')

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
  filter(type %in% c("Sovereign country","Country","Sovereignty") | adm0_a3 == "ISR") %>%
  select(Code = iso_a3_eh, Code_short = iso_a2_eh)

##### Add geometries to data
data_sf <- left_join(x = world_red, y = data, by = c("Code"))
data_sf <- st_as_sf(data_sf)

##### Setup
Period <- "Avg_2019_2022"
cols <- c("Morelli et al. (2025)" = "#FF9933",
          "Chavent et al. (2018)" = "#0000FF",
          "G1" = "#CECECE", "G2" = "#888888", "G3" = "#4F4F4F", "G4" = "#000000")
Weigths <- FALSE
if (Weigths == FALSE) {
  FolderW <- "Unweighted/"
  wt <- NULL
} else {
  FolderW <- "Weighted/"
  wt <- data_sf %>% filter(Variable == "pop") %>% select(Avg_2019_2022) %>% st_drop_geometry() %>% pull()
}




####################################################
########## Plot of climate awareness rate ##########
####################################################
p <- data_sf %>%
  select(Code,Entity,Variable,Period = any_of(Period)) %>%
  pivot_wider(names_from = Variable, values_from = Period) %>%
  ggplot() +
  geom_sf(mapping = aes(fill = R34)) +
  theme_minimal() +
  labs(title = "Share of people with medium-high or high climate awareness (2022)",
       subtitle = "'I know a moderate amount about it' = Medium-high\n'I know a lot about it' = High",
       fill = "Share % of the overall sample") +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 50) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 14),axis.text = element_text(size = 10))
ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = "Revision/Awareness_High.png")

p <- data_sf %>%
  select(Code,Entity,Variable,Period = any_of(Period)) %>%
  pivot_wider(names_from = Variable, values_from = Period) %>%
  ggplot() +
  geom_sf(mapping = aes(fill = R12)) +
  theme_minimal() +
  labs(title = "Share of people with medium-low or low climate awareness (2022)",
       subtitle = "'I have never heard of it' = Low\n'I know a little about it' = Medium-Low",
       fill = "Share % of the overall sample") +
  scale_fill_gradient2(low = "darkgreen", mid = "yellow", high = "red", midpoint = 50) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 14),axis.text = element_text(size = 10))
ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = "Revision/Awareness_Low.png")

data_p <- data_sf %>%
  select(Entity,Variable,Period = any_of(Period)) %>%
  filter(Variable %in% c("FaceSurvey22Resp")) %>%
  mutate(Variable = case_when(Variable == "FaceSurvey22Resp" ~ "Survey 2022: number of respondents (unweighted)")) %>%
  arrange(Entity) %>%
  mutate(quantile_class = cut_number(Period, n = 5))
q_breaks_round <- round(quantile(data_p$Period,probs = seq(0, 1, length.out = 6),na.rm = TRUE),0)
q_labels <- paste0(q_breaks_round[-length(q_breaks_round)]," – ",q_breaks_round[-1])
p <- data_p %>%
  ggplot(mapping = aes(fill = quantile_class)) + 
  geom_sf() + 
  scale_fill_manual(
    values = c("#c6dbef","#6baed6","#2171b5","#08309e","#08306b"),
    name = "",
    guide = guide_legend(nrow = 1),
    labels = q_labels,
  ) +
  labs(title = "Survey 2022: number of respondents (unweighted)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 18),
        axis.title = element_text(size = 14),axis.text = element_text(size = 10)) +
  theme(legend.position = "bottom")
ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = "Revision/RespondentsUnweighted.png")

data_sf %>%
  select(Entity,Variable,Period = any_of(Period)) %>%
  filter(Variable %in% c("FaceSurvey22Resp","R12","R34")) %>%
  mutate(Period = round(Period)) %>%
  pivot_wider(names_from = Variable,values_from = Period) %>%
  arrange(Entity) %>%
  st_drop_geometry() %>%
  xtable::xtable()

library(readr)
read_csv("H:/Il mio Drive/Awareness_ZammarchiMaranzano/Data/Facebook/facebook-users-by-country-2025.csv") %>%
  mutate(FacebookUsers_2025 = FacebookUsers_2025/1000000,
         FacebookUsers_2024 = FacebookUsers_2024/1000000,
         FacebookUsers_2023 = FacebookUsers_2023/1000000) %>%
  select(Code_short = flagCode,contains("FacebookUsers")) %>%
  left_join(x = data_sf %>%
              select(Code_short,Entity,Variable,Period = any_of(Period)) %>%
              filter(Variable %in% c("FaceSurvey22Resp","R12","R34","pop")) %>%
              pivot_wider(names_from = Variable, values_from = Period) %>%
              select(Code_short,Entity,FaceSurvey22Resp,R12,R34,pop), by = c("Code_short")) %>%
  arrange(Entity) %>%
  st_drop_geometry() %>%
  xtable::xtable(digits = c(0,0,0,0,0,0,2,2,2,2,0,0,0,0))










###########################################################################
#################### Main analysis: low awareness only ####################
###########################################################################
AwarenessOnly <- TRUE
Dissim_std <- TRUE

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,any_of(Period)) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Period)

if (AwarenessOnly == TRUE) {
  data_wide <- data_wide %>%
    select(Code,Entity,geometry,R = R12)
}

##### Compute distance matrix of coordinates
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide, which = "Euclidean")
# Range normalization
D1_norm <- D1/max(D1)
# Check coherence among Great circle (Geodetic) and Euclidean distances
D1_geo <- sf::st_distance(data_wide)
Check <- data.frame(Code = data_wide$Code,Entity = data_wide$Entity,Eucl = rowMeans(D1),Geod = rowMeans(D1_geo))
Check <- Check %>%
  mutate(rEucl = rank(Eucl),
         rGeod = rank(Geod),
         rdiff = rEucl - rGeod)
p1 <- Check %>%
  ggplot(mapping = aes(x = rEucl,y = rGeod)) +
  geom_point(size = 3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "orange") +
  ggpmisc::stat_correlation(size = 10, col ="black", family = "bold") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18)) + 
  labs(x = "Rank on average Euclidean distance",
       y = "Rank on average Geodetic distance",
       title = "Rank computed on the average distance between a country and the rest of the world")
p2 <- Check %>%
  ggplot(mapping = aes(x = Eucl,y = Geod)) +
  geom_point(size = 3) + 
  geom_smooth() + 
  geom_smooth(method = "lm", colour = "orange") +
  ggpmisc::stat_correlation(size = 10, col ="black", family = "bold") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18)) + 
  labs(x = "Average Euclidean distance",
       y = "Average Geodetic distance",
       title = "Average distance between a country and the rest of the world")
p <- ggarrange(p2,p1,ncol = 2,align = "hv",common.legend = TRUE,legend = "bottom")
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = "Revision/Distances.png")
  


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
k_max<-10
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
  
  cr <- choicealpha(as.dist(D0_norm),as.dist(D1_norm), range.alpha = seq(0, 1, d_a), K=k, graph=T, scale=FALSE,
                    wt = wt)
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
                 min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
  
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
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))

Indices_plot_opt <- Indices_plot %>%
  group_by(Index,Criterion) %>%
  summarise(Max_K = K[which.max(Index_val)], Max = max(Index_val),
            Min_K = K[which.min(Index_val)], Min = min(Index_val)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    OptK = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                            "Dunn (1974) -- Max",
                            "Calinski and Harabasz (1974) -- Max")) {
      Max_K
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min_K
    },
    OptVal = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                              "Dunn (1974) -- Max",
                              "Calinski and Harabasz (1974) -- Max")) {
      Max
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min
    })

p2 <- Indices_plot %>%
  ggplot(mapping = aes(x = K, y = Index_val,col = Criterion,group = Criterion)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  geom_point(data = Indices_plot_opt, mapping = aes(x = OptK, y = OptVal), size = 7, shape = 1, col = "red") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
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
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 2: $\\alpha = \\alpha^*_K$ are either determined according to Chavent et al. (2018) or Morelli et al. (2025)"),
                                       vjust = +5.0,hjust = +0.545,size = 12))
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"R12/Indices_OptimK.png"))

p1 <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Alpha = case_when(grepl("max",Alpha) ~ "Morelli et al. (2025)",
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
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
  labs(title = latex2exp::TeX("Optimal $\\alpha$ values (i.e., $\\alpha^*_K$) conditioning on $K$ clusters"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("$\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
ggpubr::ggexport(p1,width = 3400, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"R12/Optim_Alpha_K.png"))


########## Comparison between indices at alpha_*_K and at alpha_K = 0
Indices <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
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
               min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
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
  mutate(K = factor(K,levels = c("2","3","4","5","6","7","8","9","10"),
                    labels = c("2","3","4","5","6","7","8","9","10"),ordered = T))

up <- Indices_augm[Indices_augm$Diff_Index0 >= 0,]
down <- Indices_augm[Indices_augm$Diff_Index0 < 0,]


########## Plot of percentage gain or loss w.r.t. alpha=0
p <- Indices_augm %>%
  ggplot() + 
  geom_bar(mapping = aes(x = K, y = Diff_Index0,
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
ggpubr::ggexport(p,width = 3000, height = 2000,res = 210,filename = paste0("Revision/",FolderW,"R12/GainLoss_wrt_alpha0.png"))

##### Optimal combinations
if (Weigths == TRUE) {
  Comb_optim <- rbind(
    c(4,0.21),
    c(3,0.43),
    c(4,0),
    c(3,0)
  )
} else {
  Comb_optim <- rbind(
    c(4,0.33),
    c(3,0.17),
    c(4,0),
    c(3,0)
  )
}
Comb_optim <- as.data.frame(Comb_optim)

clusteringL <- clustering0 <- list()
ClusterNames <- character(length = dim(Comb_optim)[1])
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]

  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt = wt, scale=FALSE)
  clusteringL[[i]] <- cutree(hc, K_star)
  ClusterNames[i] <- paste0("Cluster_K=",K_star,"Alpha=",alpha_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt = wt, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  # names(clustering0)[i] <- paste0("Cluster_K=",K_star,"Alpha=0")
  
  data_wide <- bind_cols(data_wide,clusteringL[[i]])

  Period <- "Avg_2019_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clusteringL[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Revision/",FolderW,"R12/Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}
colnames(data_wide)[5:dim(data_wide)[2]] <- ClusterNames

# ##### Coherence among clustering methods: adjusted Rank index
# Rand_mat <- matrix(data = NA, nrow = dim(Comb_optim)[1], ncol = dim(Comb_optim)[1])
# CellName <- character(length = dim(Comb_optim)[1])
# for (i in 1:dim(Comb_optim)[1]) {
#   # i = 1
#   alpha_star <- Comb_optim$Alpha[i]
#   K_star <- Comb_optim$K[i]
#   CellName[i] <- paste0("K=",K_star," & Alpha=",alpha_star)
#   DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
#   for (j in 1:dim(Comb_optim)[1]) {
#     stats <- fpc::cluster.stats(DistMat, clustering = clustering[[i]], alt.clustering = clustering[[j]])
#     Rand_mat[i,j] <- stats$corrected.rand
#   }
# }
# colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
# Rand_mat

##### Groups characterization
StatPlot <- data_wide %>%
  as_tibble() %>%
  pivot_longer(cols = contains("Cluster_"), names_to = "Setting",values_to = "Cluster") %>%
  group_by(Setting,Cluster) %>%
  summarise(NumObs = n(),
            Mean = mean(R,na.rm=T),
            'Std.Dev.' = sd(R,na.rm=T),
            Skewness = skewness(R,na.rm=T),
            Kurtosis = kurtosis(R,na.rm=T)) %>%
  pivot_longer(cols = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"), names_to = "Stat", values_to = "Value") %>%
  ungroup() %>%
  mutate(
    Setting = paste0(
      "K = ",
      stringi::stri_extract(Setting, regex = "(?<=Cluster_K=)[0-9]+"),
      "  &  Alpha = ",
      stringi::stri_extract(Setting, regex = "(?<=Alpha=)[0-9]+\\.?[0-9]*")
    ),
    Stat = factor(Stat,labels = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"),
                  levels = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"),
                  ordered = TRUE)
  ) %>%
  mutate(Cluster = paste0("G",Cluster),
         Cluster = factor(Cluster,levels = c("G1","G2","G3","G4"),labels = c("G1","G2","G3","G4"),ordered = T))

StatTable <- StatPlot %>%
  pivot_wider(names_from = Stat, values_from = Value)

p <- StatPlot %>%
  ggplot(mapping = aes(x = Cluster, y = Value, group = Cluster, fill = Cluster)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ggh4x::facet_grid2(rows = vars(Stat), cols = vars(Setting), scales = "free") +
  theme_bw() + 
  labs(y = "",
       title = "Descriptive statistics for the share of people with low or medium-low climate change awareness",
       ) + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20, face = "bold",hjust = +0.5),
        plot.subtitle = element_text(size = 17),
        strip.text = element_text(size = 18),axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 14)) + 
  scale_fill_manual(values = cols)
ggexport(p,width = 3000, height = 2400, res = 200, filename = paste0("Revision/",FolderW,"R12/DescrByGroup.png"))

### Export statistics by group
##### Creo file Excel
wb <- createWorkbook("ClusterAwareness")
##### Save Excel sheet
StatTableList <- StatTable %>%
  group_split(Setting, .keep = TRUE)
for (i in 1:length(StatTableList)) {
  addWorksheet(wb,ClusterNames[i])
  writeData(wb, sheet = ClusterNames[i], StatTableList[[i]], colNames = T)
}
##### Salvatggio file Excel
saveWorkbook(wb,paste0("Revision/",FolderW,"R12/ClusterStats_Awareness.xlsx"),overwrite = T)










#############################################################################
#################### Robustness analysis: high awareness ####################
#############################################################################
AwarenessOnly <- TRUE
Dissim_std <- TRUE

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,any_of(Period)) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Period)

if (AwarenessOnly == TRUE) {
  data_wide <- data_wide %>%
    select(Code,Entity,geometry,R = R34)
}

##### Compute distance matrix of coordinates
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide, which = "Euclidean")
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
k_max<-10
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
                    graph=T, scale=FALSE, wt = wt)
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
                 min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
  
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
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))

Indices_plot_opt <- Indices_plot %>%
  group_by(Index,Criterion) %>%
  summarise(Max_K = K[which.max(Index_val)], Max = max(Index_val),
            Min_K = K[which.min(Index_val)], Min = min(Index_val)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    OptK = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                            "Dunn (1974) -- Max",
                            "Calinski and Harabasz (1974) -- Max")) {
      Max_K
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min_K
    },
    OptVal = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                              "Dunn (1974) -- Max",
                              "Calinski and Harabasz (1974) -- Max")) {
      Max
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min
    })

p2 <- Indices_plot %>%
  ggplot(mapping = aes(x = K, y = Index_val,col = Criterion,group = Criterion)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  geom_point(data = Indices_plot_opt, mapping = aes(x = OptK, y = OptVal), size = 7, shape = 1, col = "red") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
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
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 2: $\\alpha = \\alpha^*_K$ are either determined according to Chavent et al. (2018) or Morelli et al. (2025)"),
                                       vjust = +5.0,hjust = +0.545,size = 12))
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"R34/Indices_OptimK.png"))

p1 <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Alpha = case_when(grepl("max",Alpha) ~ "Morelli et al. (2025)",
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
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
  labs(title = latex2exp::TeX("Optimal $\\alpha$ values (i.e., $\\alpha^*_K$) conditioning on $K$ clusters"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("$\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
ggpubr::ggexport(p1,width = 3400, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"R34/Optim_Alpha_K.png"))


########## Comparison between indices at alpha_*_K and at alpha_K = 0
Indices <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
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
               min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
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
  mutate(K = factor(K,levels = c("2","3","4","5","6","7","8","9","10"),
                    labels = c("2","3","4","5","6","7","8","9","10"),ordered = T))

########## Plot of percentage gain or loss w.r.t. alpha=0
p <- Indices_augm %>%
  ggplot() + 
  geom_bar(mapping = aes(x = K, y = Diff_Index0,
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
ggpubr::ggexport(p,width = 3000, height = 2000,res = 210,filename = paste0("Revision/",FolderW,"R34/GainLoss_wrt_alpha0.png"))

##### Optimal combinations
if (Weigths == TRUE) {
  Comb_optim <- rbind(
    c(4,0.21),
    c(3,0.43),
    c(4,0),
    c(3,0)
  )
} else {
  Comb_optim <- rbind(
    c(4,0.33),
    c(3,0.17),
    c(4,0),
    c(3,0)
  )
}
Comb_optim <- as.data.frame(Comb_optim)

clusteringH <- clustering0 <- list()
ClusterNames <- character(length = dim(Comb_optim)[1])
colnames(Comb_optim) <- c("K","Alpha")

##### Map of the clusters evaluated at the optimal parameters
set.seed(12345)
for (i in 1:dim(Comb_optim)[1]) {
  # i = 1
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  
  DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt = wt, scale=FALSE)
  clusteringH[[i]] <- cutree(hc, K_star)
  ClusterNames[i] <- paste0("Cluster_K=",K_star,"Alpha=",alpha_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt = wt, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  # names(clustering0)[i] <- paste0("Cluster_K=",K_star,"Alpha=0")
  
  data_wide <- bind_cols(data_wide,clusteringH[[i]])
  
  Period <- "Avg_2019_2022"
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clusteringH[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Revision/",FolderW,"R34/Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}
colnames(data_wide)[5:dim(data_wide)[2]] <- ClusterNames

# ##### Coherence among clustering methods: adjusted Rank index
# Rand_mat <- matrix(data = NA, nrow = dim(Comb_optim)[1], ncol = dim(Comb_optim)[1])
# CellName <- character(length = dim(Comb_optim)[1])
# for (i in 1:dim(Comb_optim)[1]) {
#   # i = 1
#   alpha_star <- Comb_optim$Alpha[i]
#   K_star <- Comb_optim$K[i]
#   CellName[i] <- paste0("K=",K_star," & Alpha=",alpha_star)
#   DistMat <- (1-alpha_star)*(as.dist(D0_norm)) + alpha_star*(as.dist(D1_norm))
#   for (j in 1:dim(Comb_optim)[1]) {
#     stats <- fpc::cluster.stats(DistMat, clustering = clustering[[i]], alt.clustering = clustering[[j]])
#     Rand_mat[i,j] <- stats$corrected.rand
#   }
# }
# colnames(Rand_mat) <- rownames(Rand_mat) <- CellName
# Rand_mat

##### Groups characterization
StatPlot <- data_wide %>%
  as_tibble() %>%
  pivot_longer(cols = contains("Cluster_"), names_to = "Setting",values_to = "Cluster") %>%
  group_by(Setting,Cluster) %>%
  summarise(NumObs = n(),
            Mean = mean(R,na.rm=T),
            'Std.Dev.' = sd(R,na.rm=T),
            Skewness = skewness(R,na.rm=T),
            Kurtosis = kurtosis(R,na.rm=T)) %>%
  pivot_longer(cols = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"), names_to = "Stat", values_to = "Value") %>%
  ungroup() %>%
  mutate(
    Setting = paste0(
      "K = ",
      stringi::stri_extract(Setting, regex = "(?<=Cluster_K=)[0-9]+"),
      "  &  Alpha = ",
      stringi::stri_extract(Setting, regex = "(?<=Alpha=)[0-9]+\\.?[0-9]*")
    ),
    Stat = factor(Stat,labels = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"),
                  levels = c("NumObs","Mean","Std.Dev.","Skewness","Kurtosis"),
                  ordered = TRUE)
  ) %>%
  mutate(Cluster = paste0("G",Cluster),
         Cluster = factor(Cluster,levels = c("G1","G2","G3","G4"),labels = c("G1","G2","G3","G4"),ordered = T))

StatTable <- StatPlot %>%
  pivot_wider(names_from = Stat, values_from = Value)

p <- StatPlot %>%
  ggplot(mapping = aes(x = Cluster, y = Value, group = Cluster, fill = Cluster)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ggh4x::facet_grid2(rows = vars(Stat), cols = vars(Setting), scales = "free") +
  theme_bw() + 
  labs(y = "",
       title = "Descriptive statistics for the share of people with high or medium-high climate change awareness",
       ) + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20, face = "bold",hjust = +0.5),
        plot.subtitle = element_text(size = 17),
        strip.text = element_text(size = 18),axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 14)) + 
  scale_fill_manual(values = cols)
ggexport(p,width = 3000, height = 2400, res = 200, filename = paste0("Revision/",FolderW,"R34/DescrByGroup.png"))


### Export GreenPeace dataset
##### Creo file Excel
wb <- createWorkbook("ClusterAwareness")
##### Save Excel sheet
StatTableList <- StatTable %>%
  group_split(Setting, .keep = TRUE)
for (i in 1:length(StatTableList)) {
  addWorksheet(wb,ClusterNames[i])
  writeData(wb, sheet = ClusterNames[i], StatTableList[[i]], colNames = T)
}
##### Salvatggio file Excel
saveWorkbook(wb,paste0("Revision/",FolderW,"R34/ClusterStats_Awareness.xlsx"),overwrite = T)





################################################################################################
#################### Robustness analysis: correlation and descriptive stats ####################
################################################################################################
# https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-hierarchical-clustering

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,any_of(Period)) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Period)

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
    labs(title = "Linear correlation among awareness, socio-economic and climate-related variables (average 2019-2022)")
  if (RedCorr == TRUE) {
    ggexport(p_corr,width = 2400, height = 2000, res = 150, filename = paste0("Revision/",FolderW,"CorrMat_red.png"))
  } else {
    ggexport(p_corr,width = 2400, height = 2000, res = 150, filename = paste0("Revision/",FolderW,"CorrMat.png"))
  }
}

##### Descriptive statistics
data_wide %>%
  st_drop_geometry() %>%
  mutate(TerritorialEmiss_IntensGDP_KgThs = TerritorialEmiss_IntensGDP_KgThs*1000000,
         csh_g = csh_g*100,
         HDI = HDI*100,
         TradeOpenness = TradeOpenness*100) %>%
  pivot_longer(cols = 3:last_col(), values_to = "Value", names_to = "Variable") %>%
  group_by(Variable) %>%
  summarise(
    Min = min(Value),
    Average = mean(Value),
    Median = median(Value),
    Max = max(Value),
    SD = sd(Value),
    Skew = moments::skewness(Value),
    Kurt = moments::kurtosis(Value)
    ) %>%
  mutate(across(is.numeric,function(x) round(x,digits = 2))) %>%
  xtable::xtable()








######################################################################################################################
#################### Robustness analysis: awareness, socio-economic and climate-related variables ####################
######################################################################################################################
AwarenessOnly <- TRUE
Dissim_std <- TRUE

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2019_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2019_2022) %>%
  mutate(TerritorialEmiss_IntensGDP_KgThs = TerritorialEmiss_IntensGDP_KgThs*1000000,
         csh_g = csh_g*100,
         HDI = HDI*100,
         TradeOpenness = TradeOpenness*100)
data_wide <- data_wide %>%
  select(-c(R34))

##### Compute distance matrix of coordinates
D1_centr <- sf::st_coordinates(sf::st_centroid(data_wide, of_largest_polygon = TRUE))
D1 <- sf::st_distance(data_wide,which = "Euclidean")
# Range normalization
D1_norm <- D1/max(D1)

##### Compute Euclidean distance matrix of features
data_wide_mat <- data_wide %>%
  sf::st_drop_geometry() %>%
  select(-c("Code","Entity")) %>%
  # Robust normalization
  mutate(across(is.numeric,
                function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))
  # # Robust normalization
  # mutate(across(c(is.numeric, -contains("Lon"),-contains("Lat"),-contains("Intens")),
  #               function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))),
  #        across(contains("Intens"),
  #               function(x) (x - mean(x,na.rm = T))/sd(x,na.rm = T)))
D0 <- data_wide_mat %>%
  as.matrix() %>%
  # Range normalization
  dist(method = "euclidean")
D0_norm <- D0/max(D0)

########## Spatial clustering ##########
# Algorithm finding alpha and K
a_ch<-NULL #alpha proposed by Chavent 2018
a_max<-NULL #alpha from the maximum
k_max<-10
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
                    graph=T, scale=FALSE, wt = wt)
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
  ind <- NbClust(data=cbind(data_wide_mat,D1_centr),diss=dissim,distance=NULL, alphaBeale = 1,
                 min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
  
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
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
                               grepl("ch",Criterion) ~ "Chavent et al. (2018)",
                               TRUE ~ Criterion),
         Index = case_when(Index == "Sil" ~ "Silhouette (Rousseeuw, 1987) -- Max",
                           Index == "Dunn" ~ "Dunn (1974) -- Max",
                           Index == "McClain" ~ "McClain and Rao (1975) -- Min",
                           Index == "Cindex" ~ "C-index (Hubert and Levin, 1976) -- Min",
                           Index == "CH" ~ "Calinski and Harabasz (1974) -- Max",
                           TRUE ~ Index))

Indices_plot_opt <- Indices_plot %>%
  group_by(Index,Criterion) %>%
  summarise(Max_K = K[which.max(Index_val)], Max = max(Index_val),
            Min_K = K[which.min(Index_val)], Min = min(Index_val)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    OptK = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                            "Dunn (1974) -- Max",
                            "Calinski and Harabasz (1974) -- Max")) {
      Max_K
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min_K
    },
    OptVal = if (Index %in% c("Silhouette (Rousseeuw, 1987) -- Max",
                              "Dunn (1974) -- Max",
                              "Calinski and Harabasz (1974) -- Max")) {
      Max
    } else if (Index %in% c("McClain and Rao (1975) -- Min",
                            "C-index (Hubert and Levin, 1976) -- Min")) {
      Min
    })

p2 <- Indices_plot %>%
  ggplot(mapping = aes(x = K, y = Index_val,col = Criterion,group = Criterion)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 3) + 
  geom_point(data = Indices_plot_opt, mapping = aes(x = OptK, y = OptVal), size = 7, shape = 1, col = "red") + 
  facet_wrap(~ Index, scales = "free", ncol = 5) + 
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
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
p <- annotate_figure(p,top = text_grob(latex2exp::TeX("Note 2: $\\alpha = \\alpha^*_K$ are either determined according to Chavent et al. (2018) or Morelli et al. (2025)"),
                                       vjust = +5.0,hjust = +0.545,size = 12))
ggpubr::ggexport(p,width = 4000, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"Combined/Indices_OptimK.png"))

p1 <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Alpha = case_when(grepl("max",Alpha) ~ "Morelli et al. (2025)",
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
  scale_x_continuous(breaks = seq(from=2,to=k_max,by=1),limits = c(2,k_max)) +
  labs(title = latex2exp::TeX("Optimal $\\alpha$ values (i.e., $\\alpha^*_K$) conditioning on $K$ clusters"),
       x = latex2exp::TeX("Number of clusters $K$"), y = latex2exp::TeX("$\\alpha^*_K$")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 25,face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16)) + 
  scale_color_manual("Criterion",values = cols)
ggpubr::ggexport(p1,width = 3400, height = 2200,res = 250,filename = paste0("Revision/",FolderW,"Combined/Optim_Alpha_K.png"))


########## Comparison between indices at alpha_*_K and at alpha_K = 0
Indices <- Indices %>%
  select(Index,K = K_ch,Index_ch,Index_max,alpha_star_ch,alpha_star_max) %>%
  pivot_longer(cols = c("Index_ch","Index_max"),names_to = "Criterion",values_to = "Index_val") %>%
  pivot_longer(cols = c("alpha_star_ch","alpha_star_max"),names_to = "Alpha",values_to = "alpha_val") %>%
  mutate(K = as.numeric(K),
         Index_val = as.numeric(Index_val),
         Criterion = case_when(grepl("max",Criterion) ~ "Morelli et al. (2025)",
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
               min.nc = 2, max.nc = k_max, method="ward.D2",index="all")
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
  mutate(K = factor(K,levels = c("2","3","4","5","6","7","8","9","10"),
                    labels = c("2","3","4","5","6","7","8","9","10"),ordered = T))

########## Plot of percentage gain or loss w.r.t. alpha=0
p <- Indices_augm %>%
  ggplot() + 
  geom_bar(mapping = aes(x = K, y = Diff_Index0,
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
ggpubr::ggexport(p,width = 3000, height = 2000,res = 210,filename = paste0("Revision/",FolderW,"Combined/GainLoss_wrt_alpha0.png"))

##### Optimal combinations
if (Weigths == TRUE) {
  Comb_optim <- rbind(
    c(4,0.21),
    c(3,0.43),
    c(4,0),
    c(3,0)
  )
} else {
  Comb_optim <- rbind(
    c(4,0.33),
    c(3,0.17),
    c(4,0),
    c(3,0)
  )
}
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
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt = wt, scale=FALSE)
  clustering_AwareClimateSocio[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt = wt, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareClimateSocio[[i]])
  Period <- "Avg_2019_2022"
  colnames(data_wide)[which(colnames(data_wide) == "cluster")] <- paste0("Cluster_K",K_star,"_Alpha",alpha_star)
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareClimateSocio[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Revision/",FolderW,"Robustness/Robustness_AwareClimateSocio_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}
colnames(data_wide)[(dim(data_wide)[2]-length(ClusterNames)+1):dim(data_wide)[2]] <- ClusterNames

##### Groups characterization
StatPlot <- data_wide %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  pivot_longer(cols = contains("Cluster_"), names_to = "Setting",values_to = "Cluster") %>%
  select(-c(Code,Entity)) %>%
  pivot_longer(cols = c(everything(),-c("Setting","Cluster")), names_to = "Var",values_to = "Value") %>%
  group_by(Setting,Cluster,Var) %>%
  summarise(NumObs = n(),
            q25 = quantile(Value,probs = 0.25, na.rm=T),
            Mean = mean(Value,na.rm=T),
            Median = median(Value,na.rm=T),
            'Std.Dev.' = sd(Value,na.rm=T),
            q75 = quantile(Value,probs = 0.75, na.rm=T)) %>%
  pivot_longer(cols = c("NumObs","q25","Mean","Median","Std.Dev.","q75"), names_to = "Stat", values_to = "Value") %>%
  ungroup() %>%
  mutate(
    Setting = paste0(
      "K = ",
      stringi::stri_extract(Setting, regex = "(?<=Cluster_K=)[0-9]+"),
      "  &  Alpha = ",
      stringi::stri_extract(Setting, regex = "(?<=Alpha=)[0-9]+\\.?[0-9]*")
    ),
    Stat = factor(Stat,labels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  levels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  ordered = TRUE),
    Cluster = paste0("G",Cluster),
    Cluster = factor(Cluster,levels = c("G1","G2","G3","G4"),labels = c("G1","G2","G3","G4"),ordered = T)
  )

gt_tbl1 <- StatPlot %>% 
  filter(grepl(pattern = "K = 3",x = Setting),
         Var != "FaceSurvey22Resp") %>%
  mutate(Setting = case_when(
    grepl(pattern = "K = 3",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 3  &  Alpha = p",
    grepl(pattern = "K = 4",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 4  &  Alpha = p",
    TRUE ~ Setting
  )) %>%
  pivot_wider(names_from = c(Setting,Var),values_from = Value, names_sep = "__") %>% 
  mutate(across(is.numeric,function(x) round(x,digits = 2))) %>%
  gt::gt(groupname_col = "Cluster",
         rowname_col = c("Stat"),
         caption = "Awareness + climate + socioeconomic variables (robustness): descriptive statistics by group and setting"
         ) %>%
  tab_spanner(
    label = "K = 3  &  Alpha = 0",
    columns = starts_with("K = 3  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 3  &  Alpha = ",Comb_optim[Comb_optim$K == 3 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 3  &  Alpha = p__")
  ) %>%
  cols_label(
    # 
    `K = 3  &  Alpha = 0__R12` = "R12",
    `K = 3  &  Alpha = 0__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 3  &  Alpha = 0__csh_g` = "csh_g",
    `K = 3  &  Alpha = 0__EmpRate` = "EmpRate",
    `K = 3  &  Alpha = 0__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 3  &  Alpha = 0__HDI` = "HDI",
    `K = 3  &  Alpha = 0__pop` = "pop",
    `K = 3  &  Alpha = 0__rgdpna_pc` = "rgdpna_pc",
    `K = 3  &  Alpha = 0__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 3  &  Alpha = 0__TradeOpenness` = "TrdOpen",
    # 
    `K = 3  &  Alpha = p__R12` = "R12",
    `K = 3  &  Alpha = p__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 3  &  Alpha = p__csh_g` = "csh_g",
    `K = 3  &  Alpha = p__EmpRate` = "EmpRate",
    `K = 3  &  Alpha = p__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 3  &  Alpha = p__HDI` = "HDI",
    `K = 3  &  Alpha = p__pop` = "pop",
    `K = 3  &  Alpha = p__rgdpna_pc` = "rgdpna_pc",
    `K = 3  &  Alpha = p__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 3  &  Alpha = p__TradeOpenness` = "TrdOpen",
    # 
    `K = 3  &  Alpha = 0__cdd` = "cdd",
    `K = 3  &  Alpha = 0__hd30` = "hd30",
    `K = 3  &  Alpha = 0__pr` = "pr",
    `K = 3  &  Alpha = 0__tas` = "tas",
    `K = 3  &  Alpha = 0__tx84rr` = "tx84rr",
    `K = 3  &  Alpha = 0__wsdi` = "wsdi",
    # 
    `K = 3  &  Alpha = p__cdd` = "cdd",
    `K = 3  &  Alpha = p__hd30` = "hd30",
    `K = 3  &  Alpha = p__pr` = "pr",
    `K = 3  &  Alpha = p__tas` = "tas",
    `K = 3  &  Alpha = p__tx84rr` = "tx84rr",
    `K = 3  &  Alpha = p__wsdi` = "wsdi"
  )
gt::gtsave(gt_tbl1, paste0("Revision/",FolderW,"Robustness/DescrByGroup_AwareSocioClimate_1.png"),
           vwidth  = 3500, vheight = 2000)

gt_tbl2 <- StatPlot %>% 
  filter(grepl(pattern = "K = 4",x = Setting),
         Var != "FaceSurvey22Resp") %>%
  mutate(Setting = case_when(
    grepl(pattern = "K = 3",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 3  &  Alpha = p",
    grepl(pattern = "K = 4",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 4  &  Alpha = p",
    TRUE ~ Setting
  )) %>%
  pivot_wider(names_from = c(Setting,Var),values_from = Value, names_sep = "__") %>% 
  mutate(across(is.numeric,function(x) round(x,digits = 2))) %>%
  gt::gt(groupname_col = "Cluster",
         rowname_col = c("Stat"),
         caption = "Awareness + climate + socioeconomic variables (robustness): descriptive statistics by group and setting"
         ) %>%
  tab_spanner(
    label = "K = 4  &  Alpha = 0",
    columns = starts_with("K = 4  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 4  &  Alpha = ",Comb_optim[Comb_optim$K == 4 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 4  &  Alpha = p__")
  ) %>%
  cols_label(
    # 
    `K = 4  &  Alpha = 0__R12` = "R12",
    `K = 4  &  Alpha = 0__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 4  &  Alpha = 0__csh_g` = "csh_g",
    `K = 4  &  Alpha = 0__EmpRate` = "EmpRate",
    `K = 4  &  Alpha = 0__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 4  &  Alpha = 0__HDI` = "HDI",
    `K = 4  &  Alpha = 0__pop` = "pop",
    `K = 4  &  Alpha = 0__rgdpna_pc` = "rgdpna_pc",
    `K = 4  &  Alpha = 0__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 4  &  Alpha = 0__TradeOpenness` = "TrdOpen",
    # 
    `K = 4  &  Alpha = p__R12` = "R12",
    `K = 4  &  Alpha = p__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 4  &  Alpha = p__csh_g` = "csh_g",
    `K = 4  &  Alpha = p__EmpRate` = "EmpRate",
    `K = 4  &  Alpha = p__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 4  &  Alpha = p__HDI` = "HDI",
    `K = 4  &  Alpha = p__pop` = "pop",
    `K = 4  &  Alpha = p__rgdpna_pc` = "rgdpna_pc",
    `K = 4  &  Alpha = p__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 4  &  Alpha = p__TradeOpenness` = "TrdOpen",
    # 
    `K = 4  &  Alpha = 0__cdd` = "cdd",
    `K = 4  &  Alpha = 0__hd30` = "hd30",
    `K = 4  &  Alpha = 0__pr` = "pr",
    `K = 4  &  Alpha = 0__tas` = "tas",
    `K = 4  &  Alpha = 0__tx84rr` = "tx84rr",
    `K = 4  &  Alpha = 0__wsdi` = "wsdi",
    # 
    `K = 4  &  Alpha = p__cdd` = "cdd",
    `K = 4  &  Alpha = p__hd30` = "hd30",
    `K = 4  &  Alpha = p__pr` = "pr",
    `K = 4  &  Alpha = p__tas` = "tas",
    `K = 4  &  Alpha = p__tx84rr` = "tx84rr",
    `K = 4  &  Alpha = p__wsdi` = "wsdi"
  )
gt::gtsave(gt_tbl2, paste0("Revision/",FolderW,"Robustness/DescrByGroup_AwareSocioClimate_2.png"),
           vwidth  = 3500, vheight = 2000)





######################################################################################################
#################### Robustness analysis: awareness and climate-related variables ####################
######################################################################################################

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2019_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2019_2022) %>%
  mutate(TerritorialEmiss_IntensGDP_KgThs = TerritorialEmiss_IntensGDP_KgThs*1000000,
         csh_g = csh_g*100,
         HDI = HDI*100,
         TradeOpenness = TradeOpenness*100)

##### Compute distance matrix of coordinates
data_wide <- data_wide %>%
  select(Code,Entity,cdd,hd30,pr,tas,tx84rr,wsdi,R12) 
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

if (Weigths == TRUE) {
  Comb_optim <- rbind(
    c(4,0.21),
    c(3,0.43),
    c(4,0),
    c(3,0)
  )
} else {
  Comb_optim <- rbind(
    c(4,0.33),
    c(3,0.17),
    c(4,0),
    c(3,0)
  )
}
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
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt = wt, scale=FALSE)
  clustering_AwareClimate[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt = wt, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareClimate[[i]])
  Period <- "Avg_2019_2022"
  colnames(data_wide)[which(colnames(data_wide) == "cluster")] <- paste0("Cluster_K",K_star,"_Alpha",alpha_star)
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareClimate[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Revision/",FolderW,"Robustness/Robustness_AwareClimate_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}

##### Groups characterization
StatPlot <- data_wide %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  pivot_longer(cols = contains("Cluster_"), names_to = "Setting",values_to = "Cluster") %>%
  select(-c(Code,Entity)) %>%
  pivot_longer(cols = c(everything(),-c("Setting","Cluster")), names_to = "Var",values_to = "Value") %>%
  mutate(Value = case_when(Var == "rgdpna_pc" ~ Value/1000,
                           Var == "pr" ~ Value/100,
                           TRUE ~ Value)) %>%
  group_by(Setting,Cluster,Var) %>%
  summarise(NumObs = n(),
            q25 = quantile(Value,probs = 0.25, na.rm=T),
            Mean = mean(Value,na.rm=T),
            Median = median(Value,na.rm=T),
            'Std.Dev.' = sd(Value,na.rm=T),
            q75 = quantile(Value,probs = 0.75, na.rm=T)) %>%
  pivot_longer(cols = c("NumObs","q25","Mean","Median","Std.Dev.","q75"), names_to = "Stat", values_to = "Value") %>%
  ungroup() %>%
  mutate(
    Setting = paste0(
      "K = ",
      stringi::stri_extract(Setting, regex = "(?<=Cluster_K)[0-9]+"),
      "  &  Alpha = ",
      stringi::stri_extract(Setting, regex = "(?<=Alpha)[0-9]+\\.?[0-9]*")
    ),
    Stat = factor(Stat,labels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  levels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  ordered = TRUE),
    Cluster = paste0("G",Cluster),
    Cluster = factor(Cluster,levels = c("G1","G2","G3","G4"),labels = c("G1","G2","G3","G4"),ordered = T)
    )

gt_tbl <- StatPlot %>% 
  mutate(Setting = case_when(
    grepl(pattern = "K = 3",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 3  &  Alpha = p",
    grepl(pattern = "K = 4",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 4  &  Alpha = p",
    TRUE ~ Setting
  )) %>%
  pivot_wider(names_from = c(Setting,Var),values_from = Value, names_sep = "__") %>% 
  mutate(across(is.numeric,function(x) round(x,digits = 2))) %>%
  gt::gt(groupname_col = "Cluster",
         rowname_col = c("Stat"),
         caption = "Awareness + climate variables (robustness): descriptive statistics by group and setting") %>%
  tab_spanner(
    label = "K = 3  &  Alpha = 0",
    columns = starts_with("K = 3  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 3  &  Alpha = ",Comb_optim[Comb_optim$K == 3 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 3  &  Alpha = p__")
  ) %>%
  tab_spanner(
    label = "K = 4  &  Alpha = 0",
    columns = starts_with("K = 4  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 4  &  Alpha = ",Comb_optim[Comb_optim$K == 4 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 4  &  Alpha = p__")
  ) %>%
  cols_label(
    # 
    `K = 3  &  Alpha = 0__R12` = "R12",
    `K = 3  &  Alpha = 0__cdd` = "cdd",
    `K = 3  &  Alpha = 0__hd30` = "hd30",
    `K = 3  &  Alpha = 0__pr` = "pr",
    `K = 3  &  Alpha = 0__tas` = "tas",
    `K = 3  &  Alpha = 0__tx84rr` = "tx84rr",
    `K = 3  &  Alpha = 0__wsdi` = "wsdi",
    # 
    `K = 3  &  Alpha = p__R12` = "R12",
    `K = 3  &  Alpha = p__cdd` = "cdd",
    `K = 3  &  Alpha = p__hd30` = "hd30",
    `K = 3  &  Alpha = p__pr` = "pr",
    `K = 3  &  Alpha = p__tas` = "tas",
    `K = 3  &  Alpha = p__tx84rr` = "tx84rr",
    `K = 3  &  Alpha = p__wsdi` = "wsdi",
    # 
    `K = 4  &  Alpha = 0__R12` = "R12",
    `K = 4  &  Alpha = 0__cdd` = "cdd",
    `K = 4  &  Alpha = 0__hd30` = "hd30",
    `K = 4  &  Alpha = 0__pr` = "pr",
    `K = 4  &  Alpha = 0__tas` = "tas",
    `K = 4  &  Alpha = 0__tx84rr` = "tx84rr",
    `K = 4  &  Alpha = 0__wsdi` = "wsdi",
    # 
    `K = 4  &  Alpha = p__R12` = "R12",
    `K = 4  &  Alpha = p__cdd` = "cdd",
    `K = 4  &  Alpha = p__hd30` = "hd30",
    `K = 4  &  Alpha = p__pr` = "pr",
    `K = 4  &  Alpha = p__tas` = "tas",
    `K = 4  &  Alpha = p__tx84rr` = "tx84rr",
    `K = 4  &  Alpha = p__wsdi` = "wsdi"
  )
# cols_width(
#   Cluster  ~ px(90),
#   Stat      ~ px(100),
#   everything() ~ px(155)
# )
gt::gtsave(gt_tbl, paste0("Revision/",FolderW,"Robustness/DescrByGroup_AwareClimate.png"),
           vwidth  = 3000, vheight = 2000)





###########################################################################################
#################### Robustness analysis: awareness and socio-economic ####################
###########################################################################################

##### Dataset in wide form: 2019-2022
data_wide <- data_sf %>%
  select(Code,Entity,Variable,Avg_2019_2022) %>%
  filter(!is.na(Entity)) %>%
  pivot_wider(names_from = Variable, values_from = Avg_2019_2022) %>%
  mutate(TerritorialEmiss_IntensGDP_KgThs = TerritorialEmiss_IntensGDP_KgThs*1000000,
         csh_g = csh_g*100,
         HDI = HDI*100,
         TradeOpenness = TradeOpenness*100)

##### Compute distance matrix of coordinates
data_wide <- data_wide %>%
  select(Code,Entity,CarbonIntens_Electr,EnerIntens_PrimEnergy,HDI,csh_g,pop,rgdpna_pc,
         EmpRate,TerritorialEmiss_IntensGDP_KgThs,TradeOpenness,R12) 
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

if (Weigths == TRUE) {
  Comb_optim <- rbind(
    c(4,0.21),
    c(3,0.43),
    c(4,0),
    c(3,0)
  )
} else {
  Comb_optim <- rbind(
    c(4,0.33),
    c(3,0.17),
    c(4,0),
    c(3,0)
  )
}
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
  hc <- hclustgeo(as.dist(D0_norm,upper = T), as.dist(D1_norm), alpha=alpha_star, wt = wt, scale=FALSE)
  clustering_AwareSocio[[i]] <- cutree(hc, K_star)
  
  DistMat0 <- as.dist(D0_norm)
  hc0 <- hclustgeo(as.dist(D0_norm,upper = T), wt = wt, scale=FALSE)
  clustering0[[i]] <- cutree(hc0, K_star)
  data_wide <- bind_cols(data_wide,cluster = clustering_AwareSocio[[i]])
  Period <- "Avg_2019_2022"
  colnames(data_wide)[which(colnames(data_wide) == "cluster")] <- paste0("Cluster_K",K_star,"_Alpha",alpha_star)
  p <- data_wide %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = as.factor(clustering_AwareSocio[[i]]))) +
    theme_minimal() +
    labs(title = latex2exp::TeX(paste0("Optimal clustering results: $\\alpha^*=$",alpha_star,
                                       " & $K^*=$",K_star)),
         fill = "Clusters") +
    theme(plot.title = element_text(face = "bold",size = 20))
  ggpubr::ggexport(p,width = 2400, height = 1200,res = 250,filename = paste0("Revision/",FolderW,"Robustness/Robustness_AwareSocio_Map_Cluster_OptComb_K",K_star,"_Alpha",alpha_star,".png"))
}

##### Groups characterization
StatPlot <- data_wide %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  pivot_longer(cols = contains("Cluster_"), names_to = "Setting",values_to = "Cluster") %>%
  select(-c(Code,Entity)) %>%
  pivot_longer(cols = c(everything(),-c("Setting","Cluster")), names_to = "Var",values_to = "Value") %>%
  mutate(Value = case_when(Var == "rgdpna_pc" ~ Value/1000,
                           TRUE ~ Value)) %>%
  group_by(Setting,Cluster,Var) %>%
  summarise(NumObs = n(),
            q25 = quantile(Value,probs = 0.25, na.rm=T),
            Mean = mean(Value,na.rm=T),
            Median = median(Value,na.rm=T),
            'Std.Dev.' = sd(Value,na.rm=T),
            q75 = quantile(Value,probs = 0.75, na.rm=T)) %>%
  pivot_longer(cols = c("NumObs","q25","Mean","Median","Std.Dev.","q75"), names_to = "Stat", values_to = "Value") %>%
  ungroup() %>%
  mutate(
    Setting = paste0(
      "K = ",
      stringi::stri_extract(Setting, regex = "(?<=Cluster_K)[0-9]+"),
      "  &  Alpha = ",
      stringi::stri_extract(Setting, regex = "(?<=Alpha)[0-9]+\\.?[0-9]*")
    ),
    Stat = factor(Stat,labels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  levels = c("NumObs","q25","Mean","Median","Std.Dev.","q75"),
                  ordered = TRUE),
    Cluster = paste0("G",Cluster),
    Cluster = factor(Cluster,levels = c("G1","G2","G3","G4"),labels = c("G1","G2","G3","G4"),ordered = T)
  )

gt_tbl <- StatPlot %>% 
  mutate(Setting = case_when(
    grepl(pattern = "K = 3",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 3  &  Alpha = p",
    grepl(pattern = "K = 4",x = Setting) & grepl(pattern = "Alpha = 0.",x = Setting) ~ "K = 4  &  Alpha = p",
    TRUE ~ Setting
  )) %>%
  pivot_wider(names_from = c(Setting,Var),values_from = Value, names_sep = "__") %>% 
  mutate(across(is.numeric,function(x) round(x,digits = 2))) %>%
  gt::gt(groupname_col = "Cluster",
         rowname_col = c("Stat"),
         caption = "Awareness + socioeconomic variables (robustness): descriptive statistics by group and setting") %>%
  tab_spanner(
    label = "K = 3  &  Alpha = 0",
    columns = starts_with("K = 3  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 3  &  Alpha = ",Comb_optim[Comb_optim$K == 3 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 3  &  Alpha = p__")
  ) %>%
  tab_spanner(
    label = "K = 4  &  Alpha = 0",
    columns = starts_with("K = 4  &  Alpha = 0__")
  ) %>%
  tab_spanner(
    label = paste0("K = 4  &  Alpha = ",Comb_optim[Comb_optim$K == 4 & Comb_optim$Alpha > 0,]$Alpha),
    columns = starts_with("K = 4  &  Alpha = p__")
  ) %>%
  cols_label(
    # 
    `K = 3  &  Alpha = 0__R12` = "R12",
    `K = 3  &  Alpha = 0__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 3  &  Alpha = 0__csh_g` = "csh_g",
    `K = 3  &  Alpha = 0__EmpRate` = "EmpRate",
    `K = 3  &  Alpha = 0__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 3  &  Alpha = 0__HDI` = "HDI",
    `K = 3  &  Alpha = 0__pop` = "pop",
    `K = 3  &  Alpha = 0__rgdpna_pc` = "rgdpna_pc",
    `K = 3  &  Alpha = 0__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 3  &  Alpha = 0__TradeOpenness` = "TrdOpen",
    # 
    `K = 3  &  Alpha = p__R12` = "R12",
    `K = 3  &  Alpha = p__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 3  &  Alpha = p__csh_g` = "csh_g",
    `K = 3  &  Alpha = p__EmpRate` = "EmpRate",
    `K = 3  &  Alpha = p__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 3  &  Alpha = p__HDI` = "HDI",
    `K = 3  &  Alpha = p__pop` = "pop",
    `K = 3  &  Alpha = p__rgdpna_pc` = "rgdpna_pc",
    `K = 3  &  Alpha = p__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 3  &  Alpha = p__TradeOpenness` = "TrdOpen",
    # 
    `K = 4  &  Alpha = 0__R12` = "R12",
    `K = 4  &  Alpha = 0__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 4  &  Alpha = 0__csh_g` = "csh_g",
    `K = 4  &  Alpha = 0__EmpRate` = "EmpRate",
    `K = 4  &  Alpha = 0__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 4  &  Alpha = 0__HDI` = "HDI",
    `K = 4  &  Alpha = 0__pop` = "pop",
    `K = 4  &  Alpha = 0__rgdpna_pc` = "rgdpna_pc",
    `K = 4  &  Alpha = 0__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 4  &  Alpha = 0__TradeOpenness` = "TrdOpen",
    # 
    `K = 4  &  Alpha = p__R12` = "R12",
    `K = 4  &  Alpha = p__CarbonIntens_Electr` = "CarbInt_Ele",
    `K = 4  &  Alpha = p__csh_g` = "csh_g",
    `K = 4  &  Alpha = p__EmpRate` = "EmpRate",
    `K = 4  &  Alpha = p__EnerIntens_PrimEnergy` = "EnInt_PrimEn",
    `K = 4  &  Alpha = p__HDI` = "HDI",
    `K = 4  &  Alpha = p__pop` = "pop",
    `K = 4  &  Alpha = p__rgdpna_pc` = "rgdpna_pc",
    `K = 4  &  Alpha = p__TerritorialEmiss_IntensGDP_KgThs` = "TerEm_IntGDP",
    `K = 4  &  Alpha = p__TradeOpenness` = "TrdOpen",
  )
  # cols_width(
  #   Cluster  ~ px(90),
  #   Stat      ~ px(100),
  #   everything() ~ px(155)
  # )
gt::gtsave(gt_tbl, paste0("Revision/",FolderW,"Robustness/DescrByGroup_AwareSocio.png"),
           vwidth  = 3500, vheight = 2000)





###########################################################################################################################
#################### Robustness analysis: Coherence among clustering methods using adjusted Rank index ####################
###########################################################################################################################

Rand_mat_Rob1 <- Rand_mat_Rob2 <- Rand_mat_Rob3 <- Rand_mat_Main <- Rand_mat_Main2 <- matrix(data = NA, nrow = length(clustering_AwareClimateSocio[[1]]), ncol = 6)
# 
for (i in 1:length(clustering_AwareClimateSocio)) {
  Rand_mat_Rob1[,i+2] <- clustering_AwareClimateSocio[[i]]
}
Rand_mat_Rob1[,1] <- "Robustness: Awareness + Climate + Socio-Economic"
Rand_mat_Rob1[,2] <- 1:length(clusteringH[[1]])
Rand_mat_Rob1 <- data.frame(Rand_mat_Rob1)
# 
for (i in 1:length(clustering_AwareClimate)) {
  Rand_mat_Rob2[,i+2] <- clustering_AwareClimate[[i]]
}
Rand_mat_Rob2[,1] <- "Robustness: Awareness + Climate"
Rand_mat_Rob2[,2] <- 1:length(clusteringH[[1]])
Rand_mat_Rob2 <- data.frame(Rand_mat_Rob2)
# 
for (i in 1:length(clustering_AwareSocio)) {
  Rand_mat_Rob3[,i+2] <- clustering_AwareSocio[[i]]
}
Rand_mat_Rob3[,1] <- "Robustness: Awareness + Socio"
Rand_mat_Rob3[,2] <- 1:length(clusteringH[[1]])
Rand_mat_Rob3 <- data.frame(Rand_mat_Rob3)
# 
for (i in 1:length(clusteringL)) {
  Rand_mat_Main[,i+2] <- clusteringL[[i]]
}
Rand_mat_Main[,1] <- "Main: Low awareness only"
Rand_mat_Main[,2] <- 1:length(clusteringH[[1]])
Rand_mat_Main <- data.frame(Rand_mat_Main)
# 
for (i in 1:length(clustering_AwareClimateSocio)) {
  Rand_mat_Main2[,i+2] <- clusteringH[[i]]
}
Rand_mat_Main2[,1] <- "Robustness: High awareness only"
Rand_mat_Main2[,2] <- 1:length(clusteringH[[1]])
Rand_mat_Main2 <- data.frame(Rand_mat_Main2)

for (i in 1:dim(Comb_optim)[1]) {
  alpha_star <- Comb_optim$Alpha[i]
  K_star <- Comb_optim$K[i]
  colnames(Rand_mat_Rob1)[i+2] <- colnames(Rand_mat_Rob2)[i+2] <- colnames(Rand_mat_Rob3)[i+2] <- colnames(Rand_mat_Main)[i+2] <- colnames(Rand_mat_Main2)[i+2] <- paste0("K=",K_star," & Alpha=",alpha_star)
}
colnames(Rand_mat_Rob1)[1] <- colnames(Rand_mat_Rob2)[1] <- colnames(Rand_mat_Rob3)[1] <- colnames(Rand_mat_Main)[1] <- colnames(Rand_mat_Main2)[1] <- "Scenario"
colnames(Rand_mat_Rob1)[2] <- colnames(Rand_mat_Rob2)[2] <- colnames(Rand_mat_Rob3)[2] <- colnames(Rand_mat_Main)[2] <- colnames(Rand_mat_Main2)[2] <- "Idx_cnt"

RobustnessClustering <- bind_rows(Rand_mat_Rob1,Rand_mat_Rob2,Rand_mat_Rob3,Rand_mat_Main2,Rand_mat_Main)
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
saveWorkbook(wb,paste0("Revision/",FolderW,"Robustness/RobustnessCoherence.xlsx"),overwrite = T)


### Plot adjusted rand indices (corr matrix)
corr <- Rand_mat %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, AdjRanIndex = value)
corrmat <- corr %>%
  pivot_longer(cols = 3:last_col()) %>%
  # Standardize for fill/color scaling in tile
  group_by(name) %>%
  mutate(value_std = scale(x = value, center = T, scale = T)) %>%
  ungroup() %>%
  arrange(desc(Reg1),desc(Reg2))
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
ggexport(p_corr,width = 2400, height = 2000, res = 200, filename = paste0("Revision/",FolderW,"AdjrandIndex.png"))




