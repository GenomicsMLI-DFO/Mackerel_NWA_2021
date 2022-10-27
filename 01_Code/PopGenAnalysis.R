# Info --------------------------------------------------------------------

# Population genetics analysis on mackerel data
# Mackerel ddrad dataset - only adults
# Audrey Bourret
# 2021-07-06
#

# Library -----------------------------------------------------------------

library(here)

library(tidyverse)
library(readxl)

#library(reshape)
#library(vcfR)
library(adegenet)
#library(parallel)
library(hierfstat)


# Change the behaviours of remotes
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

# Load
library(remotes)
remotes::install_github("biodray/QuickPop")
library(QuickPop)


# PCA
#library("adegraphics")
#library("lattice")
#library("gplots")
#library("ggmap")

# RDA
#library(codep)
#library(adespatial)
#library(adegraphics)
#library(vegan)
#library(ape)
#library(car)

# Internal functions
#for(i in 1:length( list.files("./02_Functions") )){
#  source(file.path("./02_Functions",  list.files("./02_Functions")[i]))  
#}

# Define initial working directory
#current.wd <- getwd()

# Functions ---------------------------------------------------------------

`%nin%` = Negate(`%in%`)

# Data --------------------------------------------------------------------

#load(file.path(get.value("filter.data"),"SNPs.data"))

pop.data <- read_excel(here("00_Data/ProjectInfo.xlsx"), sheet ="Data")

pop.data  <- pop.data  %>% mutate(Sample = paste(Sequencage, ID_IBIS, sep ="_"),
                                  Plaque_ID = paste(Sequencage, "_P",Plaque, sep =""),
                        new.NAFO = ifelse(is.na(NAFO), Region,
                                          ifelse(
                                          NAFO == "NA", Region,
                                          ifelse(NAFO %in% c("4W", "4X"), "4WX", 
                                          ifelse(NAFO %in% c("4R", "4S"), "4RS",
                                                 ifelse(NAFO %in% c("3Ps", "4Vn"), "3Ps4Vn", 
                                                        ifelse(NAFO %in% c("5Y", "5Ze"), "5YZe",
                                                               ifelse(NAFO %in% c("6A", "6B"), "6AB",
                                                                       NAFO)))))))
)


#pop.data %>% View()

pop.data %>% pull(SNP_panel) %>% unique()

# Basic Sample Stats ------------------------------------------------------

# Table 1 summary stats

tab1 <- pop.data %>% filter(Cat_Sample == "Sample",
                    Dev_Stage == "Adult") %>%
  mutate(Country = factor(Country, levels = c("Canada", "US", "Europe")),
         NAFO = factor(NAFO, levels = c("6A","6B","5Zw","5Ze", "5Y",  "4W", "4X", "3Ps", "4Vn", "4T", "4S", "4R", "3K")),
         Month = as.numeric(as.character(Month)),
         Age = as.numeric(as.character(Age))) %>% 
  
  group_by(Country, NAFO,TypeSurvey,Year) %>% summarise(Month = paste(min(Month, na.rm = T), 
                                                                      max(Month, na.rm = T),
                                                                      sep = " - "),
                                                        MeanAge = mean(Age, na.rm = T) %>% round(digits = 1),
                                                        minAge = min(Age, na.rm = T),
                                                        maxAge = max(Age, na.rm = T),
                                                        N.female = length(Sample[Sex == "F"]),
                                                        N.male = length(Sample[Sex == "M"]),
                                                        N.NA = n() - N.female - N.male,
                                                        SexRatioFMNA = paste(N.female, N.male, N.NA,
                                                                          sep = ":"),
                                                        N.all = n(),
                                                        N.panel = length(Sample[SNP_panel == "Yes"])) %>% 
  select(-c(N.female, N.male, N.NA))  

tab1 %>% View()

#write.table(as.data.frame(tab1), "clipboard", col.names = T, row.names = F)

# How many Adults were send to GQ

pop.data %>% group_by(Cat_Sample, Dev_Stage, Sequencage) %>% summarise(N = n())

# How many samples by panel

pop.data %>% filter(Cat_Sample == "Sample",
                    Dev_Stage == "Adult",
                    SNP_panel == "Yes") %>%  
            group_by(Country) %>% summarise(N = n())

# How many ref samples

pop.data %>% pull(REF_assign) %>% unique()

pop.data %>% filter(Cat_Sample == "Sample",
                    Dev_Stage == "Adult",
                    SNP_panel == "Yes",
                    REF_assign %in% c("CAN", "US")) %>%  
  mutate(Age = as.numeric(as.character(Age))) %>% 
  group_by(Country, Month, NAFO) %>% summarise(N = n(),
                                               MeanAge = mean(Age, na.rm = T) %>% round(digits = 1),
                                               minAge = min(Age, na.rm = T),
                                               maxAge = max(Age, na.rm = T))

pop.data %>% filter(Cat_Sample == "Sample",
                    Dev_Stage == "Adult",
                    SNP_panel == "Yes",
                    REF_assign %in% c("CAN", "US")) %>%  
  mutate(Age = as.numeric(as.character(Age))) %>% 
  group_by(Country) %>% summarise(N = n(),
                                               MeanAge = mean(Age, na.rm = T) %>% round(digits = 1),
                                               minAge = min(Age, na.rm = T),
                                               maxAge = max(Age, na.rm = T))
# Basic raw reads stats ---------------------------------------------------

NSrun <- read_excel(here("00_Data/Raw_GQ_Stats/WESTPOPMAC_NovaSeqReadSet_2021-07-09.xlsx"))

HIrun <- read_excel(here("00_Data/Raw_GQ_Stats/MAQUEREAUX_HiSeqReadSet_2021-07-09.xlsx"))

names(HIrun) == names(NSrun)

GQrun <- bind_rows(NSrun, HIrun)

GQrun %>% names()

GQrun %>% mutate(Nreads = as.character(`Nombre de lectures`),
                 Nreads = map_chr(Nreads, function(d){
                   str_extract_all(d, "[:digit:]") %>% unlist() %>% paste(collapse = "")
                 })
             ,
                 Nreads = as.numeric(as.character(Nreads))) %>% 
          group_by(`Type de séquençage`) %>% 
          summarise(Nlibrary = n(),
                    Nsample = Nlibrary*96,
                    Nreads = sum(Nreads),
                    Mean.reads = Nreads/Nsample)

# Loading VCF data --------------------------------------------------------

list.files(here("00_Data"))

load(here("00_Data/gen.america.878samples.10703snps.data")) 
load(here("00_Data/gen.global.878samples.10832snps.data")) 

# Subset new data sets ----------------------------------------------------

pop.data %>% head()
pop.data %>% filter(Cat_Sample == "Sample",
                Dev_Stage == "Adult",
                SNP_panel == "Yes") %>% 
  group_by(Country, new.NAFO) %>% summarise(N = n(), 
                                        
                                        
                                        Nuni = length(unique(ID_GQ)))
# Gobal data
ID.samples.glob <- pop.data %>% filter(Cat_Sample == "Sample",
                                       Dev_Stage == "Adult",
                                       SNP_panel == "Yes") %>% 
                    pull(ID_GQ)

gl.global.samples <- gl.global.10[indNames(gl.global.10) %in% ID.samples.glob]

pop(gl.global.samples) <- data.frame(ID_GQ = indNames(gl.global.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = ifelse(Country %in% c("Canada", "US"), "NWA", "NEA")
  ) %>% pull(new.pop) 

pop(gl.global.samples) %>% table()

# Idem but for gi
gi.global.samples <- gi.global.10[indNames(gi.global.10) %in% ID.samples.glob]

pop(gi.global.samples) <- data.frame(ID_GQ = indNames(gi.global.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = ifelse(Country %in% c("Canada", "US"), "NWA", "NEA")
  ) %>% pull(new.pop) 

pop(gi.global.samples) %>% table()

# American data
ID.samples.ame <- pop.data %>% filter(Cat_Sample == "Sample",
                                  Dev_Stage == "Adult",
                                  SNP_panel == "Yes",
                                  Country %in% c("Canada", "US")) %>% 
  pull(ID_GQ)

gl.america.samples <- gl.america.10[indNames(gl.america.10) %in% ID.samples.ame]


pop(gl.america.samples) <- data.frame(ID_GQ = indNames(gl.america.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO
  ) %>% pull(new.pop) 

pop(gl.america.samples) %>% table()


gi.america.samples <- gi.america.10[indNames(gi.america.10) %in% ID.samples.ame]
pop(gi.america.samples) <-  data.frame(ID_GQ = indNames(gi.america.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO
  ) %>% pull(new.pop) 

pop(gi.america.samples) %>% table()


# Ref data

ID.samples.ref <- pop.data %>% filter(Cat_Sample == "Sample",
                                  Dev_Stage == "Adult",
                                  SNP_panel == "Yes",
                                  Country %in% c("Canada", "US"),
                                  REF_assign %in% c("CAN", "US")) %>% 
  pull(ID_GQ)

gl.ref <- gl.america.10[indNames(gl.america.10) %in% ID.samples.ref]


pop(gl.ref) <- data.frame(ID_GQ = indNames(gl.ref)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = REF_assign) %>% pull(new.pop) 

pop(gl.ref) %>% table()

# Need a gi format for latter transfo
gi.ref <- gi.america.10[indNames(gi.america.10) %in% ID.samples.ref]


pop(gi.ref) <- data.frame(ID_GQ = indNames(gi.ref)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = REF_assign) %>% pull(new.pop) 

pop(gi.ref) %>% table()

# No ref data (for assignment)

ID.samples.no.ref <- pop.data %>% filter(Cat_Sample == "Sample",
                                     Dev_Stage == "Adult",
                                     SNP_panel == "Yes",
                                     Country %in% c("Canada", "US"),
                                     REF_assign == "NA") %>% 
  pull(ID_GQ)


gl.no.ref <- gl.america.10[indNames(gl.america.10) %in% ID.samples.no.ref]


pop(gl.no.ref) <- data.frame(ID_GQ = indNames(gl.no.ref)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO) %>% pull(new.pop) 

pop(gl.no.ref) %>% table()


# Check duplicate plate 

ID.dup <-  data.frame(ID = sapply(str_split(indNames(gl.global.10), "_"), `[`, 2)) %>% 
              group_by(ID) %>% summarise(N = n()) %>% filter(N == 2) %>% pull(ID)

ID.dup.full <- data.frame(ID_GQ = indNames(gl.global.10) %>% str_subset(paste(ID.dup, collapse = "|"))) %>% 
  left_join(pop.data) %>%
  filter(Plaque == 2) %>% pull(ID_GQ) 
ID.dup.full

# Create a structure file only for duplicated individuals

gl.global.dup <-   gl.global.10[indNames(gl.global.10) %in% ID.dup.full]
pop(gl.global.dup) <-  sapply(str_split(indNames(gl.global.dup), "_"), `[`, 1)
pop(gl.global.dup) %>% table()


gi.global.dup <-   gi.global.10[indNames(gi.global.10) %in% ID.dup.full]
pop(gi.global.dup) <-  sapply(str_split(indNames(gi.global.dup), "_"), `[`, 1)
pop(gi.global.dup) %>% table()

#hf.global.dup <-  genind2hierfstat(gi.global.dup)

# Create ref + NEA for diversity analysis

ID.div <- c(pop.data %>% filter(Cat_Sample == "Sample",
                              Dev_Stage == "Adult",
                              SNP_panel == "Yes",
                              Country %in% c("Canada", "US"),
                              REF_assign %in% c("CAN", "US")) %>% 
                              pull(ID_GQ),
            pop.data %>% filter(Cat_Sample == "Sample",
                                Dev_Stage == "Adult",
                                SNP_panel == "Yes",
                                Country %in% c("Europe")) %>% 
              pull(ID_GQ)
      )

gi.global.div <-   gi.global.10[indNames(gi.global.10) %in% ID.div]
pop(gi.global.div) <-   data.frame(ID_GQ = indNames(gi.global.div)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = ifelse(is.na(REF_assign) | REF_assign == "NA", Country, REF_assign)) %>% pull(new.pop) 

pop(gi.global.div) %>% table()

gi.global.div4 <-   gi.global.10[indNames(gi.global.10) %in% ID.div]
pop(gi.global.div4) <-   data.frame(ID_GQ = indNames(gi.global.div4)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = ifelse(is.na(REF_assign) | REF_assign == "NA", Region, REF_assign)) %>% pull(new.pop) 

pop(gi.global.div4) %>% table()

# Check duplicate individuals ---------------------------------------------

# With a DAPC
daPop.dup.prelim <- dapc(gl.global.dup, pop = pop(gl.global.dup), n.da=100, n.pca=80)
temp <- optim.a.score(daPop.dup.prelim)

scatter.dapc(daPop.dup.prelim, xax=1, scree.pca=T, posi.pca = "topleft")
#scatter.dapc(daPop.ame.opti, xax=1, yax=2, scree.pca=T, posi.pca="bottomleft")

# With an AMOVA
# 
 sequencer.vec <- pop(gi.global.dup)
 ind.dup.vec   <- sapply(str_split(indNames(gi.global.dup), "_"), `[`, 2)
# loci <- hf.global.dup [, -1] 
# 
# sequencer.vec %>% table()
# ind.dup.vec %>% table()
# 
# amova.seq <- test.g(loci, level = sequencer.vec) 
# amova.seq2 <- varcomp.glob(levels = data.frame(sequencer.vec), loci, diploid = TRUE) 
# amova.seq3 <- varcomp.glob(levels = data.frame(ind.dup.vec), loci, diploid = TRUE) 
# 
# amova.seq4 <- test.between(loci, test.lev = sequencer.vec, rand.unit = ind.dup.vec, nperm = 100) 
# amova.seq5 <- test.between(loci, test.lev = ind.dup.vec, rand.unit = sequencer.vec, nperm = 100) 
# 
# amova.seq
# amova.seq2
# amova.seq3
# amova.seq4
# amova.seq5

library(poppr)

strata(gi.global.dup) <- data.frame(SEQ = sequencer.vec, IND = ind.dup.vec)

# amova.result <- poppr.amova(gi.global.dup, ~IND/SEQ, method = "pegas")  
# amova.result
# 
# 
# amova.result.ind <- poppr.amova(gi.global.dup, ~IND+SEQ, within = T)  
# amova.result.ind
# 
# amova.test <- randtest(amova.result) # Test for significance
# plot(amova.test)
# 
# amova.test


amova.result.seq <- poppr.amova(gi.global.dup, ~SEQ)  
amova.result.seq

amova.test.seq <- randtest(amova.result.seq) # Test for significance
plot(amova.test)

amova.test.seq
  
# Basic stats -------------------------------------------------------------

na.gl.ind<- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}

# Missing value distribution among NAFO division

gl.global.samples
gl.america.samples

graph.NAa <- data.frame(ID_GQ = indNames(gl.global.samples)) %>% mutate(nNA = na.gl.ind(gl.global.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "Greenland",
                              ifelse(Country == "Canada", "Canada",
                                     ifelse(Country == "US", "US", NA)))),
         Site = factor(Site, levels = c("Canada", "US", "Greenland", "Bay of Biscay"))) %>% 
  ggplot(aes(x = Site, y = nNA * 100)) +
  geom_violin()+
  geom_jitter(height = 0, aes(col = Sequencage), alpha = 0.5) +
  labs(x = "", y = "Individual missing value (%)") + 
  guides(colour = guide_legend(title = "Sequencing strategy"))+
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

graph.NAa

graph.NAb <- data.frame(ID_GQ = indNames(gl.america.samples)) %>% mutate(nNA = na.gl.ind(gl.america.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))),
         new.NAFO = factor(new.NAFO, levels = c("6AB","5Zw","5YZe",  "4WX", "3Ps4Vn", "4T", "4RS", "3K"))) %>% 
  ggplot(aes(x = new.NAFO, y = nNA * 100)) +
  geom_violin()+
  geom_jitter(height = 0, aes(col = Sequencage), alpha = 0.5) +
  #facet_grid(.~Country, scales = "free", space = "free") +
  labs(x = "", y = "Individual missing value (%)") + 
  guides(colour = guide_legend(title = "Sequencing strategy"))+
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

graph.NAb


figNA <- ggpubr::ggarrange(graph.NAa + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                   graph.NAb + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  labels = c("A", "B"),
                  common.legend = T, legend = "top",
                  ncol = 2, widths = c(2,3), align = "hv")

figNA

ggsave(filename = here("02_Results/figSupp_NAdistribution.png"), plot = figNA, 
       width = 7, height = 4 , units = "in",  bg = "white",
       dpi = 300)


# Stats on NA  

## NWA-NEA
data.frame(ID_GQ = indNames(gl.global.samples)) %>% mutate(nNA = na.gl.ind(gl.global.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(US = ifelse(Country == "US", "US", "OTHERS")) %>% 
  group_by(US) %>% summarise(MeanNA = mean(nNA),
                             SdNA = sd(nNA))
## NWA
data.frame(ID_GQ = indNames(gl.america.samples)) %>% mutate(nNA = na.gl.ind(gl.america.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(US = ifelse(Country == "US", "US", "OTHERS")) %>% 
  group_by(US) %>% summarise(MeanNA = mean(nNA),
                             SdNA = sd(nNA))


library(poppr)

div.res <- poppr(gi.global.div, total = F, sample = 9, plot = T)
div.res


div.rar <- diversity_ci(gi.global.div, ci = 95, total = F)  

div.rar

div4.rar <- diversity_ci(gi.global.div4, ci = 95, total = F)  

div4.rar


# Should be revised

hf.america.samples <- genind2hierfstat(gi.america.samples)

basicstat.america.samples <- basic.stats(hf.america.samples, diploid = TRUE, digits = 4)


str(basicstat.america.samples)

Ho <- basicstat.america.samples$Ho %>% as.data.frame() %>% 
  mutate(Key = row.names(.)) %>%   
  pivot_longer(cols = names(.) %>% str_subset("Key", negate = T), names_to = "Pop", values_to = "Ho")

str(Hs)

Ho %>% ggplot(aes(x = Pop, y = Ho)) + 
  geom_boxplot() +
  theme_bw()




# PCA ---------------------------------------------------------------------

pca.america <- glPca(gl.america.samples, center = TRUE, scale = FALSE, 
                     parallel = F, n.core = 1, nf = 1000)


pca.america.scale <- glPca(gl.america.samples, center = TRUE, scale = TRUE, 
                            parallel = T, n.core = 40, nf = 5)


# Need to be rerun

pca.global <- glPca(gl.global.samples, center = TRUE, scale = FALSE, 
                    parallel = F, n.core = 1, nf = 1000)

#save(list = c("pca.global","pca.america"),
#     file = here("02_Results/01_PCA/PCA_june2021.data"))

load(here("02_Results/01_PCA/PCA_june2021.data"))

pca.global$eig[1]/sum(pca.global$eig) # proportion of variation explained by 1st axis
pca.global$eig[2]/sum(pca.global$eig)# proportion of variation explained by 2nd axis 
pca.global$eig[3]/sum(pca.global$eig) # proportion of variation explained by 3rd axis


#res.PCA.global.table %>% nrow()

res.PCA.global.table <- data.frame(Sample = row.names(pca.global$scores),
                                   score = pca.global$scores[,1:10]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                ifelse(Country == "Canada", "NWA - Canada",
                ifelse(Country == "US", "NWA - US", NA)))),
         nNA = na.gl.ind(gl.global.samples))



mean.axis <- res.PCA.global.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                   Mean.score.Axis.2 = mean(score.PC2))

res.PCA.global.table <- res.PCA.global.table %>% left_join(mean.axis )

fig2a <- res.PCA.global.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Site)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  
  #geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
   #            size = 0.2, alpha = 0.5)+
  geom_point(aes(), alpha = 0.5, size = 2) + 
  
  
  scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PC1 - 0.7 %", y = "PC2 - 0.2 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

fig2a

fig2a.NA <- res.PCA.global.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
  #             size = 0.2, alpha = 0.5)+
  scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  
  
  geom_point(aes(shape = Site, fill = nNA), alpha = 0.75, size = 2) + 
 # scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  scale_fill_distiller(palette = "Spectral", limits = c(0, 0.4)) +
  scale_shape_manual(values = c(22, 23, 24, 25)) +
  labs(x = "PC1 - 0.7 %", y = "PC2 - 0.2 %") + 
  guides(fill = guide_colorbar(title = "% missing values", title.position = "top"),
        shape = guide_legend(title = element_blank(), nrow = 2))+
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom")
#axis.title.x = element_blank())

fig2a.NA


pca.america$eig[1]/sum(pca.america$eig) # proportion of variation explained by 1st axis
pca.america$eig[2]/sum(pca.america$eig)# proportion of variation explained by 2nd axis 



res.PCA.america.table <- data.frame(Sample = row.names(pca.america$scores),
                                    score = pca.america$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))),
         nNA = na.gl.ind(gl.america.samples))

mean.axis <- res.PCA.america.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                   Mean.score.Axis.2 = mean(score.PC2))

res.PCA.america.table <- res.PCA.america.table %>% left_join(mean.axis )


fig2b <- res.PCA.america.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
  #             size = 0.2, alpha = 0.5)+
    geom_point(aes(col = Site), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("firebrick2", "dodgerblue1")) +
    labs(x = "PC1 - 0.3 %", y = "PC2 - 0.3 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

fig2b


fig2b.NA <- res.PCA.america.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(shape = Site, fill = nNA), alpha = 0.75, size = 2) + 
  # scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  scale_fill_distiller(palette = "Spectral", limits = c(0,.4)) +
  scale_shape_manual(values = c(24, 25)) + 
  labs(x = "PC1 - 0.3 %", y = "PC2 - 0.3 %") + 
  guides(fill = guide_colorbar(title = "% missing values", title.position = "top"),
         shape = guide_legend(title = element_blank(), nrow = 2))+
  
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom")
#axis.title.x = element_blank())

fig2b.NA


ggsave(filename = here("02_Results/fig2a_PCA_global.png"), plot = fig2a, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/fig2b_PCA_NWA.png"), plot = fig2b, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

# Joining fig4 a and b
library(ggpubr)

fig2 <- ggpubr::ggarrange(fig2a + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  fig2b + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                                    labels = c("A", "B"),
                 common.legend = T, legend = "bottom",
                  ncol = 2)

fig2

ggsave(filename = here("02_Results/fig2_PCA_wLines.png"), plot = fig2, 
       width = 7, height = 4 , units = "in",
       dpi = 300)


fig2.NA <- ggarrange(fig2a.NA + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  fig2b.NA + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  labels = c("A", "B"),
                  common.legend = T, legend = "bottom",
                  ncol = 2)

fig2.NA

ggsave(filename = here("02_Results/fig2_PCA_NAs.png"), plot = fig2.NA, 
       width = 7, height = 4 , units = "in", bg = "white",
       dpi = 300)


# PCA - ref samples -------------------------------------------------------

pca.ref <- glPca(gl.ref, center = TRUE, scale = FALSE, 
                 parallel = F, n.core = 40, nf = 1000)

pca_var(pca.ref, naxe = 10)

pca.ref$eig[1]/sum(pca.ref$eig) 
pca.ref$eig[2]/sum(pca.ref$eig) 


#save(list = c("pca.ref","pca.ref.outliers", "pca.ref.neutral"),
#     file = here("02_Results/01_PCA/PCA_ref.data"))

load(here("02_Results/01_PCA/PCA_ref.data"))

pca.ref.neutral$eig[1]/sum(pca.ref.neutral$eig) 
pca.ref.neutral$eig[2]/sum(pca.ref.neutral$eig) 

pca.ref.outliers$eig[1]/sum(pca.ref.outliers$eig) 
pca.ref.outliers$eig[2]/sum(pca.ref.outliers$eig) 

res.PCA.ref.table <- data.frame(Sample = row.names(pca.ref$scores),
                                 score = pca.ref$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.ref.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                 Mean.score.Axis.2 = mean(score.PC2))

res.PCA.ref.table <- res.PCA.ref.table %>% left_join(mean.axis ) %>% mutate(nNA = na.gl.ind(gl.ref))



figS4.PCA <- res.PCA.ref.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
    stat_ellipse(aes(col =REF_assign),  level = 0.70) +
  #ggforce::geom_mark_ellipse(aes(label = REF_assign, col =REF_assign, filter = !is.na(REF_assign)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2, col = REF_assign),
  #             size = 0.2, alpha = 0.5)+
  
  geom_point(aes(col = REF_assign), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"), labels = c("Northern", "Southern")) +
  labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS4.PCA

figS4.PCA.NA <- res.PCA.ref.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
#  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
 #              size = 0.2, alpha = 0.5)+
  geom_point(aes(shape = REF_assign, fill = nNA), alpha = 0.75, size = 2) + 
  # scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  scale_fill_distiller(palette = "Spectral", limits = c(0,.4)) +
  scale_shape_manual(values = c(24, 25),  labels = c("Northern", "Southern")) + 
  labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  guides(fill = guide_colorbar(title = "% missing values", title.position = "top"),
         shape = guide_legend(title = element_blank(), nrow = 2))+
  
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom")
#axis.title.x = element_blank())
figS4.PCA.NA



ggsave(filename = here("02_Results/figS4_PCAref.png"), plot = figS4.PCA, 
       width = 4, height = 4 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figS4_PCAref_NA.png"), plot = figS4.PCA.NA, 
       width = 4, height = 4 , units = "in",
       dpi = 300)


res.PCA.ref.neutral.table <- data.frame(Sample = row.names(pca.ref.neutral$scores),
                                score = pca.ref.neutral$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.ref.neutral.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                Mean.score.Axis.2 = mean(score.PC2))

res.PCA.ref.neutral.table <- res.PCA.ref.neutral.table %>% left_join(mean.axis ) %>% mutate(nNA = na.gl.ind(gl.ref.neutral))



figS4.PCA.neutral <- res.PCA.ref.neutral.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2, col = REF_assign),
  #             size = 0.2, alpha = 0.5)+
  
  geom_point(aes(col = REF_assign), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"), labels = c("Northern", "Southern")) +
  labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS4.PCA.neutral



res.PCA.ref.outliers.table <- data.frame(Sample = row.names(pca.ref.outliers$scores),
                                        score = pca.ref.outliers$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.ref.outliers.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                        Mean.score.Axis.2 = mean(score.PC2))

res.PCA.ref.outliers.table <- res.PCA.ref.outliers.table %>% left_join(mean.axis ) %>% mutate(nNA = na.gl.ind(gl.ref.outliers))



figS4.PCA.outliers <- res.PCA.ref.outliers.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
#  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2, col = REF_assign),
 #              size = 0.2, alpha = 0.5)+
  
  geom_point(aes(col = REF_assign), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"), labels = c("Northern", "Southern")) +
  labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS4.PCA.outliers

ggsave(filename = here("02_Results/figS4_PCAref_neutral.png"), plot = figS4.PCA.neutral, 
       width = 4, height = 4 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figS4_PCAref_outliers.png"), plot = figS4.PCA.outliers, 
       width = 4, height = 4 , units = "in",
       dpi = 300)


figS.ref <- ggpubr::ggarrange(figS4.PCA + theme(axis.title = element_text(size = 10 ),
                                      
                                      plot.margin = margin(10, 10, 10, 20, "pt")),
                      figS4.PCA.neutral + theme(axis.title = element_text(size = 10 ),
                                      
                                      plot.margin = margin(10, 10, 10, 20, "pt")),
                     labels = LETTERS,
                     common.legend = T, legend = "bottom",
                     ncol = 2)

figS.ref

ggsave(filename = here("02_Results/figS_REF_AllvsNeutral.png"), plot = fig2.NA, 
       width = 7, height = 4 , units = "in", bg = "white",
       dpi = 300)


fig2a$layers[[3]]$aes_params$size <- 1
fig2b$layers[[3]]$aes_params$size <- 1
figS4.PCA$layers[[3]]$aes_params$size <- 1

fig2v2 <- ggpubr::ggarrange(fig2a + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")) +
                            guides(colour = guide_legend(override.aes = list(size=3)))
                              
                          ,
                  fig2b + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  figS4.PCA + theme(axis.title = element_text(size = 10 ),
                                    
                                    plot.margin = margin(10, 10, 10, 20, "pt")),
                  labels = LETTERS,
                  common.legend = T, legend = "bottom",
                  ncol = 3)

fig2v2

ggsave(filename = here("02_Results/fig2_PCA_withREF.png"), plot = fig2v2, 
       width = 8, height = 3 , units = "in", bg = "white",
       dpi = 300)



figNAv2 <- ggpubr::ggarrange(fig2a.NA + theme(axis.title = element_text(size = 10 ),
                                  
                                  plot.margin = margin(10, 10, 10, 20, "pt"))
                    ,
                    fig2b.NA + theme(axis.title = element_text(size = 10 ),
                                  
                                  plot.margin = margin(10, 10, 10, 20, "pt")),
                    figS4.PCA.NA + theme(axis.title = element_text(size = 10 ),
                                      
                                      plot.margin = margin(10, 10, 10, 20, "pt")),
                    labels = LETTERS,
                    common.legend = T, legend = "bottom",
                    ncol = 3)

figNAv2

ggsave(filename = here("02_Results/fig2_PCA_withREF_NA.png"), plot = figNAv2, 
       width = 8, height = 3 , units = "in", bg = "white",
       dpi = 300)



# PCoA --------------------------------------------------------------------
library(dartR)

# Function to count NA




df.america.samples        <- tab(gl.america.samples,  NA.method = c("asis"))
df.america.samples.NAmax        <- tab(gl.america.samples,  NA.method = c("mean"))
df.america.samples.NAmax  <- apply(df.america.samples , 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))


euc.america.samples.NAmax <- dist(df.america.samples.NAmax, method = "euclidean", diag = T, upper = T)

pcoa.america <-  gl.pcoa(euc.america.samples.NAmax, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

pcoa.america$eig[1]/sum(pcoa.america$eig)
pcoa.america$eig[2]/sum(pcoa.america$eig)
pcoa.america$eig[3]/sum(pcoa.america$eig)

df.global.samples        <- tab(gl.global.samples,  NA.method = c("asis"))
df.global.samples.NAmax        <- tab(gl.global.samples,  NA.method = c("mean"))
df.global.samples.NAmax  <- apply(df.global.samples , 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
euc.global.samples.NAmax <- dist(df.global.samples.NAmax, method = "euclidean", diag = T, upper = T)

pcoa.global <-  gl.pcoa(euc.global.samples.NAmax, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

gl.pcoa.plot(pcoa.global, gl.global.samples, xaxis = 1 , yaxis = 3)


pcoa.global$eig[1]/sum(pcoa.global$eig)
pcoa.global$eig[2]/sum(pcoa.global$eig)
pcoa.global$eig[3]/sum(pcoa.global$eig)

str(pcoa.global)

res.PCoA.global.table <- data.frame(Sample = row.names(pcoa.global$scores),
                                   score = pcoa.global$scores[,1:05]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCoA.global.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.Axis.1),
                                                                     Mean.score.Axis.2 = mean(score.Axis.2))

res.PCoA.global.table <- res.PCoA.global.table %>% left_join(mean.axis ) %>% mutate(nNA = na.gl.ind(gl.global.samples))



fig2a.pcoa <- res.PCoA.global.table  %>% 
  ggplot(aes(x = score.Axis.1, y = score.Axis.2)) +
  #stat_ellipse(aes(col = Site), level = 0.95, type = "norm") +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.Axis.1, yend= score.Axis.2, col = Site),
               size = 0.2, alpha = 0.5)+
  
  geom_point(aes(col = Site), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PCo 1 (0.7 %)", y = "PCo 2 (0.4 %)") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
     legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

fig2a.pcoa


res.PCoA.global.table  %>% 
  ggplot(aes(x = score.Axis.1, y = score.Axis.2, col = nNA)) +
  #stat_ellipse(aes(col = Site), level = 0.95, type = "norm") +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.Axis.1, yend= score.Axis.2, col = Site),
#               size = 0.2, alpha = 0.5)+
  scale_colour_distiller(palette = "Spectral") +
  geom_point(alpha = 0.5, size = 2) + 
 # scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PCo 1 (0.7 %)", y = "PCo 2 (0.4 %)") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
 facet_wrap(~Country) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank())



res.PCoA.america.table <- data.frame(Sample = row.names(pcoa.america$scores),
                                    score = pcoa.america$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCoA.america.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.Axis.1),
                                                        Mean.score.Axis.2 = mean(score.Axis.2))

res.PCoA.america.table <- res.PCoA.america.table %>% left_join(mean.axis )

fig2b.pcoa <- res.PCoA.america.table  %>% 
  ggplot(aes(x = score.Axis.1, y = score.Axis.2)) +
 # stat_ellipse(aes(col = Site), level = 0.95, type = "norm") +
  #stat_conf_ellipse(aes(color = Site), level = 0.95) +
#  ggforce::geom_mark_ellipse(aes(label = Site, col =Site))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +

  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.Axis.1, yend= score.Axis.2, col = Site),
               size = 0.2, alpha = 0.5)+
  geom_point(aes(col = Site), alpha = 0.5, size = 2) + 
  scale_colour_manual(values = c("firebrick2", "dodgerblue1")) +
  labs(x = "PCo 1 (0.4 %)", y = "PCo 2 (0.3 %)") + 
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank())

  fig2b.pcoa


ggsave(filename = here("02_Results/fig2a_PCoA_global.png"), plot = fig2a.pcoa, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/fig2b_PCoA_NWA.png"), plot = fig2b.pcoa, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

# Joining fig4 a and b
library(ggpubr)

fig2.pcoa <- ggarrange(fig2a.pcoa + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  fig2b.pcoa + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")),
                  labels = c("A", "B"),
                  common.legend = T, legend = "bottom", 
                  ncol = 2)

fig2.pcoa

ggsave(filename = here("02_Results/fig2_PCoA.png"), plot = fig2.pcoa, 
       width = 7, height = 4 , units = "in", bg = "white",
       dpi = 300)


df.ref        <- tab(gl.ref,  NA.method = c("mean"))
df.ref.NAmax  <- apply(df.ref , 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
euc.ref.NAmax <- dist(df.ref.NAmax, method = "euclidean", diag = T, upper = T)

pcoa.ref <-  gl.pcoa(euc.ref.NAmax, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

gl.pcoa.plot(pcoa.ref, gl.ref, xaxis = 1 , yaxis = 2)


pcoa.ref$eig[1]/sum(pcoa.ref$eig)
pcoa.ref$eig[2]/sum(pcoa.ref$eig)
pcoa.ref$eig[3]/sum(pcoa.ref$eig)



res.PCoA.ref.table <- data.frame(Sample = row.names(pcoa.ref$scores),
                                     score = pcoa.ref$scores[,1:5]) %>% 
  left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCoA.ref.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.Axis.1),
                                                                 Mean.score.Axis.2 = mean(score.Axis.2))

res.PCoA.ref.table <- res.PCoA.ref.table %>% left_join(mean.axis ) %>% mutate(nNA = na.gl.ind(gl.ref))

hist(na.gl.ind(gl.ref))

figS4.PCoA <- res.PCoA.ref.table  %>% 
  ggplot(aes(x = score.Axis.1, y = score.Axis.2)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.Axis.1, yend= score.Axis.2, col = REF_assign),
               size = 0.2, alpha = 0.5)+
  
  geom_point(aes(col = REF_assign), alpha = 0.5, size = 2) + 
 scale_colour_manual(values = c("firebrick2", "dodgerblue1"), labels = c("Northern", "Southern")) +
  labs(x = "PCo 1 (1.2 %)", y = "PCo 2 (1.1 %)") + 
  #labs(x = "PC1 - 1.1 %", y = "PC2 - 1.1 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS4.PCoA



res.PCoA.ref.table  %>% 
  ggplot(aes(x = score.Axis.1, y = score.Axis.2, col = nNA)) +
  #stat_ellipse(aes(col = Site), level = 0.95, type = "norm") +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.Axis.1, yend= score.Axis.2, col = Site),
  #               size = 0.2, alpha = 0.5)+
  scale_colour_gradient(low = "blue", high = "red") +
  geom_point(alpha = 0.5, size = 2) + 
  # scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PCo 1 (0.7 %)", y = "PCo 2 (0.4 %)") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  facet_wrap(~Country) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank())


  ggsave(filename = here("02_Results/figS4_PCoAref.png"), plot = figS4.PCoA, 
         width = 4, height = 4 , units = "in",
         dpi = 300)

# Structure analysis ------------------------------------------------------

# The Monte Carlo Markov Chain was run for 100,000 steps, following a burn-in period of 10,000 steps
# Parallel_structure() relies on the mclapply() function from the Parallel package

#install.packages("ParallelStructure", repos="http://R-Forge.R-project.org")

library(ParallelStructure)
library(dartR)

# Convertion to structure

pop.df <- data.frame(pop = c("6AB","5Zw","5YZe",  "4WX", "3Ps4Vn", "4T", "4RS", "3K", "BOB", "Groenland")) %>% 
             mutate(num = 1:length(unique(pop)))
pop.df

# Function to save in the right format
gl2str.noMarker <- function(gl.int, file.out, path.out){
  ploidy(gl.int) <- 2 #  Doesn't work if it is not specify
  
  ind.pop <-  data.frame(pop = pop(gl.int)) %>% 
    left_join(pop.df) %>% pull(num)
  
  gl2structure(x = gl.int,   addcolumns = ind.pop,  outfile = file.out, outpath = path.out)
  
  #remove marker
  convert.int <- readLines(file.path(path.out, file.out))
  
  write(convert.int[-1], file=file.path(path.out, file.out))
  
  cat("\nThere's", nInd(gl.int), "individuals and", nLoc(gl.int), "loci in this dataset",
      "\nPopulations are : ", paste(sort(unique(ind.pop)), collapse=", "))
  
}

# Create a subset with less missing data

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")]

  names(res) <- names(res) %>% str_remove("[.]0")

  return(res)

}


# Function to create a list of loci, from a genind object

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
  # Create vectors for each loci
  MAF.res <- adegenet::minorAllele(gi)
  NA.res  <- na.gi.count(gi) # Internal function

  # Filter by threshold
  MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
  cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")

  NA.loc <- names(NA.res[NA.res <= NA.trs])
  cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")

  # LOCI with both conditions
  LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
  LOCI.res %>% length()

  cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")

  return(LOCI.res)
}

# Create vectors of good loci

loc.global.MAF10.NA.05 <- filter.MAF.NA(gi.global.samples, MAF.trs = 0.10, NA.trs = 0.05)
loc.america.MAF10.NA.05 <- filter.MAF.NA(gi.america.samples, MAF.trs = 0.10, NA.trs = 0.05)

# Change de pop name for global dataset
pop(gl.global.samples) <- data.frame(ID_GQ = indNames(gl.global.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO
  ) %>% pull(new.pop) 

gl2str.noMarker(gl.global.samples, file.out = "gl.global.samples.str", path.out = "./02_Results/02_Structure/"  )
gl2str.noMarker(gl.america.samples, file.out = "gl.america.samples.str", path.out = "./02_Results/02_Structure/"  )

# Reduce subset
gl2str.noMarker(gl.global.samples[, loc.global.MAF10.NA.05], file.out = "gl.global.samples.MAF10.NA05.str", path.out = "./02_Results/02_Structure/"  )
gl2str.noMarker(gl.america.samples[, loc.america.MAF10.NA.05], file.out = "gl.america.samples.MAF10.NA05.str", path.out = "./02_Results/02_Structure/"  )

create.jobs <- function(nk = 1:4, nrep = 1, popnum = 1, burn = 1000, iteration = 10000){
  
  # Prevent scientific notation
  burn <- format(burn, scientific=F)
  iteration <- format(iteration, scientific=F)
  
  cat("\nCreating", length(nk) * nrep, "jobs:",
      "\nK =", min(nk), "-", max(nk),
      "\nN rep =", nrep,
      "\nBurning =", burn,
      "\nMCMC =", iteration,
      "\n\n"
  )
  
  # Job list
  jobs_df <- data.frame(jobID = character(),
                        pop =  character(),
                        k = numeric(),
                        burnin = character(),
                        MCMC = character(),
                        stringsAsFactors = F)
  
  
  for(x in 1:nrep){
    jobs_df.int <- data.frame(jobID = paste(paste0("k",nk), x, sep="_r"),
                              pop = paste(popnum, collapse = ","),
                              k = nk,
                              burnin = burn,
                              MCMC = iteration) # doit être en texte
    
    jobs_df <- rbind(jobs_df, jobs_df.int)
  }
  
  return(jobs_df)
}

# Global job
jobs_df_1 <- create.jobs(nk = 1:10, nrep = 10, popnum = 1:10, burn = 100000, iteration = 200000)
# NWA jobs
jobs_df_2 <- create.jobs(nk = 1:8, nrep = 10, popnum = 1:8, burn = 100000, iteration = 200000)

write(t(jobs_df_1), ncol=length(jobs_df_1[1,]), file='02_Results/02_Structure/01_Global/joblist_final.txt')
write(t(jobs_df_2), ncol=length(jobs_df_2[1,]), file='02_Results/02_Structure/02_America/joblist_final.txt')

# RUN

STR_path='~/Documents/Programs/Structure/'

# Run structure

parallel_structure(structure_path=STR_path, joblist='02_Results/02_Structure/01_Global/joblist_final.txt', n_cpu=8, infile='02_Results/02_Structure/gl.global.samples.MAF10.NA05.str', outpath='02_Results/02_Structure/01_Global/', numinds=629, numloci=4140, printqhat=1, plot_output=0)

file.copy(from = here("results_summary.csv"),
          to = here("02_Results/02_Structure/01_Global/results_summary_2021-08-20.csv"))

parallel_structure(structure_path=STR_path, joblist='02_Results/02_Structure/02_America/joblist_final.txt', n_cpu=8, infile='02_Results/02_Structure/gl.america.samples.MAF10.NA05.str', outpath='02_Results/02_Structure/02_America/', numinds=557, numloci=3930, printqhat=1, plot_output=0)
#parallel_structure(structure_path=STR_path, joblist='02_Results/02_Structure/02_America/joblist_2.txt', n_cpu=8, infile='02_Results/02_Structure/gl.america.samples.MAF10.NA05.str', outpath='02_Results/02_Structure/02_America/', numinds=557, numloci=3930, printqhat=1, plot_output=1)

file.copy(from = here("results_summary.csv"),
          to = here("02_Results/02_Structure/02_America/results_summary_2021-08-20.csv"))




# Try with pophelper - but to hard to install
# http://www.royfrancis.com/pophelper/index.html


# Compute Evanno by hand

summary.global.str <- QuickPop::strc_readsummary(here("02_Results/02_Structure/01_Global/results_summary_2021-08-20.csv"))
summary.global.str

graph.global.Evanno <- QuickPop::strc_evanno(summary.global.str, plot = "resume", plot.max = TRUE)
graph.global.Evanno

ggsave(filename = file.path(here::here(), "02_Results", "02_Structure", "Evanno.global.png"), 
       plot = graph.global.Evanno$plot.full,
       height = 5, width = 5, units = "in")          


summary.america.str <- QuickPop::strc_readsummary(here("02_Results/02_Structure/02_America/results_summary_2021-08-20.csv"))
summary.america.str

graph.america.Evanno <- QuickPop::strc_evanno(summary.america.str, plot = "resume", plot.max = TRUE)
graph.america.Evanno

ggsave(filename = file.path(here::here(), "02_Results", "02_Structure", "Evanno.america.png"), 
       plot = graph.america.Evanno$plot.full,
       height = 5, width = 5, units = "in")  

# Create figures

# RUN CLUMPP

# Create param file

clumpp.param <- function(folder.path = ".", k){
  # Read files
  struct.files <- list.files(folder.path, full.names = T, pattern = paste0("results_job_k", k, "_")) %>% str_subset("_q")
  
  n.ind <- length(readLines(struct.files[1]))
  n.run <- length(struct.files)
  
  # New file names
  pop.file   <- file.path(folder.path, paste0("clumpp.popfile.k", 2))
  perm.file  <- file.path(folder.path, paste0("clumpp.perm.k", 2)) 
  misc.file  <- file.path(folder.path, paste0("clumpp.misc.k", 2))
  out.file   <- file.path(folder.path, paste0("clumpp.out.k", 2))
  param.file <- file.path(folder.path, paste0("clumpp.param.k", 2))
  
  
  cat("", file = pop.file, append = FALSE, sep = "") 
  
  for(x in seq_along(struct.files)){
    
    
    struct.int <- utils::read.table(struct.files[x], header = F)
    names(struct.int) <- c("ID", "PopNum", paste0("Q", 1:(base::ncol(struct.int) - 
                                                            2)))
    
    struct.int <- struct.int %>% mutate(ID = 1:n.ind,
                                        ID = paste0(ID, ":")) %>% 
      relocate(PopNum, .after = last_col())
    
    res <- vector()
    
    for(y in 1:nrow(struct.int)){
    
      res[y] <- paste(struct.int[y,] , collapse= " ")
    }
    
    
    cat(res, "\n", file = pop.file, append = TRUE, sep = "\n") 
  
    
    
    
      
  }
  
  cat("DATATYPE 1",
      paste("INDFILE ", pop.file),
      paste("POPFILE", pop.file),
      paste("OUTFILE", out.file),
      paste("MISCFILE", misc.file),
      paste("K", k),
      paste("C", n.ind),
      paste("R", n.run),
      "M 1",
      "W 1",
      "S 2",
      #Additional options for greedy
      "GREEDY_OPTION 2",
      "REPEATS 1000",
      "PERMUTATIONFILE",
      #Optional outputs
      "PRINT_PERMUTED_DATA 1",
      paste("PERMUTED_DATAFILE", perm.file),
      "PRINT_EVERY_PERM 0",
      "EVERY_PERMFILE",
      "PRINT_RANDOM_INPUTORDER 0",
      "RANDOM_INPUTORDERFILE",
      #Advanced options
      "OVERRIDE_WARNINGS 0",
      "ORDER_BY_RUN 1",
      
      sep = "\n",
      file = param.file, append = FALSE)
  
  
}
  


clumpp <- function(folder.path = ".", k){
  
  param.file   <- file.path(folder.path, paste0("clumpp.param.k", 2))
  out.file     <- file.path(folder.path, paste0("clumpp.out.k", 2))
  final.file   <- file.path(folder.path, paste0("clumpp.result.k", 2))
  
  # Extract.ID
  struct.files <- list.files(folder.path, full.names = T, pattern = paste0("results_job_k", k, "_")) %>% str_subset("_q")
    struct.int <- utils::read.table(struct.files[1], header = F)
  
  ID <- struct.int[,1]
  
  cmd <- paste(param.file)

  A <- system2(command = "CLUMPP", args = cmd, stdout = TRUE)
  
  print(A)
  
  res  <-  utils::read.table(out.file, header = F)

  names(res) <-  c("ID", paste0("Q", 1:(base::ncol(struct.int) - 
                                                    2)), "PopNum")

  res$ID <- ID

  write_csv(res, file = final.file)    
}  
  
clumpp.readq <- function (file.path){
  struct.int <- readr::read_csv(file.path)
  struct.int$K <- (ncol(struct.int) - 2)
  res <- tidyr::pivot_longer(struct.int, stringr::str_subset(names(struct.int), 
                                                             "ID|PopNum|K", negate = T), names_to = "Qpop", values_to = "Qvalue")
  return(res)
} 

      


clumpp.param(folder.path = "./02_Results/02_Structure/02_America", k = 2)

clumpp(folder.path = "./02_Results/02_Structure/02_America", k = 2)

K2.america <- clumpp.readq("./02_Results/02_Structure/02_America/clumpp.result.k2")

clumpp.param(folder.path = "./02_Results/02_Structure/01_Global", k = 2)

clumpp(folder.path = "./02_Results/02_Structure/01_Global", k = 2)

K2.global <- clumpp.readq("./02_Results/02_Structure/01_Global/clumpp.result.k2")


str_graph.1 <- bind_rows(K2.global) %>% left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "Greenland", new.NAFO)),
         Site = factor(Site, levels = c("6AB","5Zw","5YZe", "4WX", "3Ps4Vn", "4T", "4RS", "3K", "Greenland", "Bay of Biscay"))) %>% 
  ggplot(aes(x =  reorder(ID, Qvalue, FUN = function(x) min(x)), y = Qvalue, fill = Qpop)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("dodgerblue1", "chartreuse4")) +
  
  #scale_fill_manual(values=col1)+
  #scale_fill_manual(breaks = c("European", "Canadian", "Undefined 1", "Undefined 2"), values = c("royalblue2", "firebrick2", "goldenrod2", "purple2"), aesthetics = "fill") +
  #scale_y_continuous(limits = c(0,1)) +
  facet_grid(K ~ Site, space = "free", scale = "free", labeller = labeller(K = label_both)) + 
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )

str_graph.1 

ggsave(filename = file.path(here::here(), "02_Results", "02_Structure", "Structure.global.K2.png"), 
       plot = str_graph.1,
       height = 5, width = 12, units = "in")   


id.order <- K2.america %>% filter(Qpop == "Q2") %>% arrange(Qvalue) %>% pull(ID)

str_graph.2 <- bind_rows(K2.america) %>% left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "Greenland", new.NAFO)),
         Site = factor(Site, levels = c("6AB","5Zw","5YZe", "4WX", "3Ps4Vn", "4T", "4RS", "3K", "Greenland", "Bay of Biscay")),
         ID_GQ = fct_relevel(ID, id.order)) %>% 
  ggplot(aes(x =  ID_GQ, y = Qvalue, fill = Qpop)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("dodgerblue1", "firebrick2")) +
  
  #scale_fill_manual(values=col1)+
  #scale_fill_manual(breaks = c("European", "Canadian", "Undefined 1", "Undefined 2"), values = c("royalblue2", "firebrick2", "goldenrod2", "purple2"), aesthetics = "fill") +
  #scale_y_continuous(limits = c(0,1)) +
  facet_grid(K ~ Site, space = "free", scale = "free", labeller = labeller(K = label_both)) + 
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )

str_graph.2 


# Joining fig4 a and b
library(ggpubr)

figS3 <- ggarrange(ggarrange(graph.global.Evanno$plot.resume + theme(plot.margin = margin(55, 10, 10, 0, "pt")),
                    str_graph.1 + theme(plot.margin = margin(0, 10, 40, 10, "pt")),
                  labels = c("A", ""), ncol = 2, nrow = 1, widths = c(2,3)),
          ggarrange(graph.america.Evanno$plot.resume  + theme(plot.margin = margin(35, 10, 10, 0, "pt")),
          str_graph.2 + theme(plot.margin = margin(0, 10, 40, 10, "pt")),
          labels = c("B", ""), ncol = 2, nrow = 1, widths = c(2,3)),
                
                  common.legend = T, legend = "none", 
                  ncol = 1, nrow = 2)

figS3


ggsave(filename = here("02_Results/figS3_Structure.png"), plot = figS3, 
       width = 10, height = 6 , units = "in",
       dpi = 300)

# DAPC ---------------------------------------------------------------------

pop(gl.global.samples) %>% table

pop.global <- data.frame(ID_GQ = indNames(gl.global.samples)) %>% left_join(pop.data) %>% 
mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                     ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                            ifelse(Country == "Canada", "NWA - Canada",
                                   ifelse(Country == "US", "NWA - US", NA)))))
pop.global %>% group_by(Site) %>% summarise(N = n())



daPop.global.prelim <- dapc(gl.global.samples, pop = factor(pop.global$Site), n.da=100, n.pca=200, glPca= pca.global)
scatter.dapc(daPop.global.prelim, xax=1, yax=2, scree.pca=T)

# Chose the best number of PC
temp <- optim.a.score(daPop.global.prelim, n.sim=20) # 115 - 152

daPop.global.opti <- dapc(gl.global.samples, pop = factor(pop.global$Site), n.da=3, n.pca=50, glPca= pca.global)
scatter.dapc(daPop.global.opti, xax=1, yax=1, scree.pca=T)

dapc.global.res <- data.frame(ID_GQ = row.names( daPop.global.opti$ind.coord),
                   LD1 =  daPop.global.opti$ind.coord[,1]) %>% left_join(pop.data) %>% 
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))


library(ggridges)

ggdapc1.global <- dapc.global.res %>% ggplot(aes(x = LD1, y = Site, fill = Site)) + 
  geom_density_ridges(alpha = 0.75, scale = 1.5) + 
  #scale_fill_manual(values =  c("firebrick2", "dodgerblue1")) +
  scale_fill_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  scale_y_discrete(labels = c("Bay of Biscay", "Greenland", "Canada", "US")) +
  #geom_vline(xintercept = c( -1.6479146, 0.7751274), lty = "dashed", col =  c("firebrick2", "dodgerblue1")) +
  
  labs(x = "DA1", y = "") +
  theme_bw() + theme( legend.title = element_blank())
ggdapc1.global

ggdapc2.global <-data.frame(eig = daPop.global.opti$eig) %>% mutate(order = 1:nrow(.)) %>% 
  ggplot(aes(x = order, y = eig)) + 
  geom_bar(stat = "identity", fill = c("darkgray", rep("white", length(daPop.global.opti$eig) - 1)), col = "black") +
  
  labs(y = "DA eigenvalue") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
ggdapc2.global            

ggdapc3.global <- data.frame(eig = cumsum(daPop.global.opti$pca.eig) / sum(daPop.global.opti$pca.eig))  %>% 
  mutate(order = 1:nrow(.))  %>% ggplot(aes(x = order, y = eig)) + 
  
  
  geom_bar(stat = "identity", fill = c(rep("darkgray", daPop.global.opti$n.pca ), rep("white", length(daPop.global.opti$pca.eig) - daPop.global.opti$n.pca)), col = c(rep("black",daPop.global.opti$n.pca ), rep("gray",length(daPop.global.opti$pca.eig) - daPop.global.opti$n.pca))) +
  labs(y = "PCA eigenvalue") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
ggdapc3.global


# NWA

pop(gl.america.samples) %>% table

#vec.pop <- data.frame(ID_GQ = indNames(gl.america.samples)) %>% left_join(pop.data) %>% 
#  mutate(new.pop = ifelse(REF_assign %in% c("CAN", "US"), paste0("Ref-", REF_assign), new.NAFO)) %>% 
#  pull(new.pop)

daPop.ame.prelim <- dapc(gl.america.samples, pop = pop(gl.america.samples), n.da=10, n.pca=200, glPca= pca.america)
temp <- optim.a.score(daPop.ame.prelim, n.sim=20) # 115 - 152

daPop.ame.opti <- dapc(gl.america.samples, pop = pop(gl.america.samples), 
                       n.da=10, n.pca=120,  glPca= pca.america)

adegenet::scatter.dapc(daPop.ame.opti, xax=1, yax=2, scree.pca=T, posi.pca="bottomleft", legend = T)


library(ggridges)

dapc.ame.res <- data.frame(ID_GQ = row.names( daPop.ame.opti$ind.coord),
                  LD1 =  daPop.ame.opti$ind.coord[,1]) %>% left_join(pop.data) %>% 
 mutate(#new.pop = ifelse(REF_assign %in% c("CAN", "US"), paste0("Ref-", REF_assign), new.NAFO),
    
    new.NAFO = factor(new.NAFO, levels = c("Ref-US", "6AB","5Zw","5YZe",  "4WX", "3Ps4Vn", "4T", "4RS", "3K", "Ref-CAN")))

ggdapc1.ame <- dapc.ame.res %>% ggplot(aes(x = LD1, y = new.NAFO, fill = Country)) + 
  geom_density_ridges(alpha = 0.75, scale = 1.5) + 
  scale_fill_manual(values =  c("firebrick2", "dodgerblue1")) +
  #geom_vline(xintercept = c( -1.6479146, 0.7751274), lty = "dashed", col =  c("firebrick2", "dodgerblue1")) +
  labs(x = "DA1", y = "") +
  theme_bw()
ggdapc1.ame

ggdapc2.ame <-data.frame(eig = daPop.ame.opti$eig) %>% mutate(order = 1:nrow(.)) %>% 
              ggplot(aes(x = order, y = eig)) + 
             geom_bar(stat = "identity", fill = c("darkgray", rep("white", length(daPop.ame.opti$eig) - 1)), col = "black") +
            
  labs(y = "DA eigenvalue") +
  theme_bw() +
             theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())
ggdapc2.ame            

ggdapc3.ame <- data.frame(eig = cumsum(daPop.ame.opti$pca.eig) / sum(daPop.ame.opti$pca.eig))  %>% 
  mutate(order = 1:nrow(.))  %>% ggplot(aes(x = order, y = eig)) + 

  
    geom_bar(stat = "identity", fill = c(rep("darkgray", daPop.ame.opti$n.pca ), rep("white", length(daPop.ame.opti$pca.eig) - daPop.ame.opti$n.pca)), col = c(rep("black",daPop.ame.opti$n.pca ), rep("gray",length(daPop.ame.opti$pca.eig) - daPop.ame.opti$n.pca))) +
  labs(y = "PCA eigenvalue") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
ggdapc3.ame


ggdapc.ame <-  ggpubr::ggarrange(ggdapc1.ame , ggpubr::ggarrange(ggdapc2.ame, ggdapc3.ame, nrow = 2, ncol = 1,  align = "hv"),
                  nrow = 1, ncol = 2, common.legend = T, align = "h", widths = c(3,2), legend = "none")

ggdapc.ame

ggdapc.global <-  ggpubr::ggarrange(ggdapc1.global , ggpubr::ggarrange(ggdapc2.global, ggdapc3.global, nrow = 2, ncol = 1,  align = "hv"),
                                 nrow = 1, ncol = 2, common.legend = T, align = "h", widths = c(3,2))

ggdapc.global


ggdapc.all <- ggpubr::ggarrange(ggdapc.global, ggdapc.ame, labels = c("A", "B"), nrow = 2, ncol = 1,
                  common.legend = T, align = "hv")


ggdapc.all


ggsave(filename = here("02_Results/figS_DAPC_ALL.png"), plot = ggdapc.all, 
       width = 6, height = 8 , units = "in", bg = "white",
       dpi = 300)



# Fst ---------------------------------------------------------------------

library(dartR)

# Fst between NWA and NEA
table(pop(gl.global.samples))


global.FST <- gl.fst.pop(gl.global.samples, nboots = 999, percent = 95, nclusters = 1)
global.FST$Fsts
global.FST$Bootstraps$`Lower bound CI limit`
global.FST$Bootstraps$`Upper bound CI limit`

# Fst among NWA 
table(pop(gl.america.samples))

america.FST <- gl.fst.pop(gl.america.samples, nboots = 999, percent = 95, nclusters = 1)
america.FST

#america.FST <- gl.fst.pop(gl.america.samples[indNames(gl.america.samples) %nin% bad.samples], nboots = 100, percent = 95, nclusters = 1)
#america.FST
#gl.america.samples[indNames(gl.america.samples) %nin% bad.samples]

# Fst value between REF
table(pop(gl.ref))

ref.FST <- gl.fst.pop(gl.ref, nboots = 999, percent = 95, nclusters = 1)
ref.FST$Fsts["US","CAN"]
ref.FST$Bootstraps$`Lower bound CI limit`
ref.FST$Bootstraps$`Upper bound CI limit`
ref.FST$Bootstraps$`p-value`

# Save Fst results

save(list = c("global.FST","america.FST", "ref.FST"),
     file = file.path(here("02_Results/04_Fst/Fst.data")))

load(file.path(here("02_Results/04_Fst/Fst.data")))

# Figure

res.america.FST <- america.FST$Bootstraps %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")



library(ggplot2)

res.america.FST <- res.america.FST %>% mutate(Comparison = paste(Population2, Population1, sep = "-"),
                                              Comparison = ifelse(Comparison == "3Ps4Vn-4WX", "4WX-3Ps4Vn",
                                                                  #ifelse(Comparison == "4T-4WX", "4WX-4T",
                                                                  #        ifelse(Comparison == "3K-4WX", "4WX-3K",
                                                                  ifelse(Comparison == "5Zw-6AB", "6AB-5Zw",
                                                                         ifelse(Comparison == "5YZe-6AB", "6AB-5YZe",
                                                                                ifelse(Comparison == "5YZe-6AB", "6AB-5YZe",       
                                                                                       
                                                                                       
                                                                                       
                                                                                       Comparison)))),
                                              
                                              Comparison = factor(Comparison, levels = c("6AB-5Zw", "6AB-5YZe", "5Zw-5YZe", 
                                                                                         "4WX-3Ps4Vn", "4WX-4T", "4WX-4RS", "4WX-3K", "3Ps4Vn-4T", "3Ps4Vn-4RS", "3Ps4Vn-3K", "4T-4RS", "4T-3K", "4RS-3K",
                                                                                         "6AB-4WX","6AB-3Ps4Vn", "6AB-4T", "6AB-4RS", "6AB-3K",
                                                                                         "5Zw-4WX","5Zw-3Ps4Vn", "5Zw-4T", "5Zw-4RS", "5Zw-3K",
                                                                                         "5YZe-4WX","5YZe-3Ps4Vn", "5YZe-4T", "5YZe-4RS", "5YZe-3K")),
                                              Levels = ifelse(Comparison %in% c("3K-3Ps4Vn","3K-4RS", "3K-4T", "3K-4WX", "3Ps4Vn-3K", 
                                                                                "3Ps4Vn-4RS", "3Ps4Vn-4T", "3Ps4Vn-4WX",
                                                                                "4RS-4T", "4RS-4WX", "4RS-3K", "4T-4WX",
                                                                                "4WX-3Ps4Vn", "4T-3Ps4Vn", "4RS-3Ps4Vn", "4T-3K", "4T-4RS", "4WX-4RS", "4WX-4T", "4WX-3K"), "Intra-Canada", 
                                                              ifelse(Comparison %in% c("5YZe-5Zw", "5YZe-6AB", "5Zw-6AB", "6AB-5YZe", "6AB-5Zw", "5Zw-5YZe"), "Intra-US",    
                                                                     "Inter-Country")),
                                              Levels = factor(Levels, levels = c("Inter-Country", "Intra-Canada", "Intra-US"))
)

res.america.FST

res.america.FST %>% group_by(Levels) %>% summarise(MinFST = min(Fst),
                                                   MaxFST = max(Fst))

res.america.FST %>% filter(Levels == "Intra-Canada") %>% arrange(desc(`Lower bound CI limit`))

write_csv(res.america.FST, file.path(here::here(), "02_Results", "04_Fst", "FST_NWA_results.csv"))

ggFst <- res.america.FST %>%
  ggplot(aes(x = Comparison, y = Fst)) + 
  geom_hline(yintercept = 0, col = "gray") +
  facet_grid(. ~ Levels, space = "free", scale = "free") +
  labs(x = "") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))


ggFst 

fig3 <- ggFst +
  geom_hline(yintercept = ref.FST$Fsts["US","CAN"], lty = "dashed", col = "red", cex = 1) +  
  geom_point(size = 4) +
  geom_pointrange(aes(ymin = `Lower bound CI limit`, ymax = `Upper bound CI limit`)) 


fig3.v2 <- ggFst + geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin = ref.FST$Bootstraps$`Lower bound CI limit`, ymax = ref.FST$Bootstraps$`Upper bound CI limit`),
                  fill= "lemonchiffon", col = "lemonchiffon2", cex = 0.2)+
   geom_hline(yintercept = ref.FST$Fsts["US","CAN"], lty = "dashed", col = "gold2", cex = 0.5) +  
  geom_point(size = 4) +
  geom_pointrange(aes(ymin = `Lower bound CI limit`, ymax = `Upper bound CI limit`)) 
fig3.v2

ggsave(filename = here("02_Results/fig3_Fst_NWA.png"), plot = fig3, 
       width = 7, height = 4 , units = "in",
       dpi = 300)

#ggsave(filename = here("02_Results/fig3.v2_Fst_NWA.png"), plot = fig3.v2, 
#       width = 7, height = 4 , units = "in",
#       dpi = 300)

ggsave(filename = here("02_Results/fig3.v3_Fst_NWA.png"), plot = fig3.v2, 
       width = 7, height = 4 , units = "in",
       dpi = 300)
# Assignation to references : MC validation -----------------------------------------------

# Save external in structure format

library(dartR)

# RDA

gl2structure(x = gl.ref,   addcolumns =  pop(gl.ref),  outfile = "gl.ref_june2021.str",   outpath = file.path("./02_Results/03_assignPOP/"))

gl2structure(x = gl.no.ref,   addcolumns =  pop(gl.no.ref),  outfile = "gl.no.ref_june2021.str",   outpath = file.path("./02_Results/03_assignPOP/"))


#bad.ID

gl2structure(x = gl.ref[indNames(gl.ref) %nin% bad.ID],   addcolumns =  pop(gl.ref[indNames(gl.ref) %nin% bad.ID]),  outfile = "gl.ref_june2021_15.str",   outpath = file.path("./02_Results/03_assignPOP/"))

gl2structure(x = gl.no.ref[indNames(gl.no.ref) %nin% bad.ID],   addcolumns =  pop(gl.no.ref[indNames(gl.no.ref) %nin% bad.ID]),  outfile = "gl.no.ref_june2021_15.str",   outpath = file.path("./02_Results/03_assignPOP/"))

gl.ref[indNames(gl.ref) %nin% bad.ID] %>% pop() %>% table()

# Import structure

model.assign <- c("lda", "svm", "naiveBayes", "tree", "randomForest") 
#model.assign <- c("svm") 

ref.group <- c("gl.ref_june2021.str")

# MONTE-CARLO CROSS VALIDATION

library(assignPOP)


# Loop over ref.group
for(x in ref.group){
  
  print(x)
  
  data.str <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/", x))
  
  # Loop over the models
  for(m in model.assign){
    
    print(m)
    
    assign.MC(data.str, train.inds=c(0.7), train.loci=c(0.1, 0.25, 0.5, 1),
              loci.sample="fst", iterations=30, model=m, dir=paste0("./02_Results/03_assignPOP/", x %>% str_remove("str"),m,"/"))
    # Unload the memory  
    gc()  
  }
  
}


# Create a big ref

accuMC.all <- data.frame()

for(x in ref.group){
  
  print(x)
  
  # Loop over the models
  for(m in model.assign){
    
    print(m)
    
    accuMC <- accuracy.MC(paste0("./02_Results/03_assignPOP/", x %>% str_remove("str"),m,"/"))
    accuMC <- accuMC %>% mutate(ref = x, model = m)
    
    accuMC.all <- bind_rows(accuMC.all, accuMC)
    
  }
  
}

accuMC.all <- accuMC.all %>% pivot_longer(names(accuMC.all) %>% str_subset("assign.rate"),
                                          names_to = "group", values_to = "assign.rate") %>% 
  mutate(group = group %>% str_remove("assign.rate."),
         assign.rate = as.numeric(as.character(assign.rate)),
         ref = ref %>% str_remove("gl.america.10.") %>% str_remove(".str") %>% str_remove("gl."),
         model = as.factor(model)
  ) %>% 
  filter(!is.na(assign.rate),
         group != "all"#,
         #str_detect(ref, "2REF")
  ) 


# Change levels of model
levels(accuMC.all$model)
levels(accuMC.all$model) <- c("LDA", "Naive Bayes", "Random Forest", "SVM", "Decision tree")         

# Stats

accuMC.all %>% group_by(model, train.loci, group) %>% 
               summarise(mean = mean(assign.rate),
                         sd = sd(assign.rate))


fig4a <- accuMC.all %>% mutate(P.train.loci = as.numeric(as.character(train.loci)) * 100) %>% 
  ggplot(aes(x = group, y = assign.rate, fill = factor(P.train.loci))) +
  geom_hline(yintercept = c(0.5), lty = "dashed") +
  geom_boxplot() +
  facet_grid(. ~ model, space = "free_x", scale = "free_x") +
  labs(x = "Reference contingent", y = "Assigment accuracy") +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  scale_fill_manual(values = c("gray100", "gray85", "gray70", "gray40"))+
  guides(fill=guide_legend(title="% of SNPs")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top")

fig4a

ggsave(filename = here("02_Results/fig4a_MCvalidation.png"), plot = fig4a, 
       width = 7, height = 3 , units = "in",
       dpi = 300)

fig4a.v2 <- accuMC.all %>% mutate(P.train.loci = as.numeric(as.character(train.loci)) * 100) %>% 
  filter(model == "SVM",
         P.train.loci == 100) %>% 
  ggplot(aes(x = group, y = assign.rate)) +
  geom_hline(yintercept = c(0.5), lty = "dashed") +
  geom_boxplot() +
  #facet_grid(. ~ model, space = "free_x", scale = "free_x") +
  labs(x = "Reference contingent", y = "Assigment accuracy") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  scale_fill_manual(values = c("firebrick2", "dodgerblue1"))+
    theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

fig4a.v2

ggsave(filename = here("02_Results/fig4a.v2_MCvalidation.png"), plot = fig4a.v2, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

# Assignation to references : Simulated Hybrids -----------------------------------------------

devtools::install_github("bwringe/parallelnewhybrid") #required
devtools::install_github("rystanley/genepopedit") #required
devtools::install_github("bwringe/hybriddetective") #This package



library("parallelnewhybrid")


ref.repro.5.str    <- assignPOP::read.Structure(file.path("./03_Results/07_Analysis/07d_Assignations/", "gl.america.10.repro5.2REF.str"))

ref.repro.5.str %>% str()

ref.repro.5.str$SampleID

gl.repro.5.REF <-   gl.america.10[indNames(gl.america.10) %in% ref.repro.5.str$SampleID]

gl.repro.5.REF[,]

# Conversion in a genpop format to create hybrids

library(dartR)

df.ref <- genind2df(gi.ref)

genpop.transfo <- function(x) {
  ifelse(is.na(x), "000000",
         ifelse(x == "00", "100100",
                ifelse(x == "01", "100200",
                       ifelse(x == "11", "200200", "000000"
                       ))))
} 


#unflatten.ref <- df.ref %>% mutate(pop = paste0(pop, row.names(.) %>% str_remove("HI|NS")))


unflatten.ref <- bind_cols(df.ref %>% mutate(pop = paste0(pop, row.names(.) %>% str_remove("HI|NS"))) %>% select(pop) ,
                                    df.ref %>% select(-pop) %>% mutate_all(genpop.transfo))

unflatten.ref[1:10, 1:10]

# Then save as genpop

genepopedit::genepop_unflatten(df = unflatten.ref, path = file.path("./02_Results/03_assignPOP/NewHybrids/", "gi.ref.txt"))

# Check that the new genpop file is OK

genepopedit::genepop_detective("./02_Results/03_assignPOP/NewHybrids/gi.ref.txt",
                               variable = "PopNum")

genepopedit::genepop_detective("./02_Results/03_assignPOP/NewHybrids/gi.ref.txt",
                               variable = "Allele")
# Create Hybrids

# FEEL LIKE THE FOLLOWING LINES WERE JUST TESTS
#hybriddetective::freqbasedsim_GTFreq(GenePopData = "./02_Results/03_assignPOP/NewHybrids/gi.ref.txt")

#write.table(x = sim_PurePops, file = "./03_Results/07_Analysis/07e_Hybrids/test.txt",
#            row.names = F, col.names = F, quote = F)

#hybriddetective::freqbasedsim_GTFreq(GenePopData = "./03_Results/07_Analysis/07e_Hybrids/test.txt")

hybriddetective::freqbasedsim_AlleleSample(GPD = "./02_Results/03_assignPOP/NewHybrids/gi.ref.txt",
                                           pop.groups = c("CAN", "US"),
                                           NumSims = 3, NumReps = 3, prop.sample = 0.9)

# Converting back to genpop

files.to.convert <- list.files("./02_Results/03_assignPOP/NewHybrids/", pattern = "NH", full.names = T)
files.to.convert

for(x in files.to.convert){
  
  print(x)
  
  temp <- read.table(x, skip =  4, header = T, colClasses = "character")
  
  names(temp) <- names(temp) %>% str_remove("X")
  # then save as genpop
  
  genepopedit::genepop_unflatten(df = temp, path = x %>% str_remove("_NH"))
  
}


# Do the assignments ...
ref.gp <- assignPOP::read.Genepop(file.path("./02_Results/03_assignPOP/NewHybrids/", "gi.ref.txt"))

# Perform the assignation

for(s in 1:3){
  
  for(r in 1:3){
    
    print(paste0("S",s,"R", r))
    
    hyb.gp  <- assignPOP::read.Genepop(paste0("./02_Results/03_assignPOP/NewHybrids/gi.ref_S",s,"R", r,".txt"))
    
    assignPOP::assign.X( x1=ref.gp , x2=hyb.gp , dir=paste0("./02_Results/03_assignPOP/NewHybrids/Assignation_S",s,"R", r,"/"), model="svm", mplot = F)
    
  }
  
}

# Collect the results

# Sample name - this file is created 
ID.key <- read.table("./02_Results/03_assignPOP/NewHybrids/gi.ref_individuals.txt",
                     header = F, col.names = "ID") 

ID.key %>% head()

# Empty data.frame
RES <- data.frame(Ind.ID = character(), 
                  pred.pop = character(), 
                  pop.1 = numeric(), 
                  pop.2 = numeric(),
                  Analyse = character())

for(s in 1:3){
  
  for(r in 1:3){
    
    temp <- read.table(file.path(paste0("./02_Results/03_assignPOP/NewHybrids/Assignation_S",s,"R", r,"/"), "AssignmentResult.txt"), header = T) 
    
    temp$Ind.ID <- ID.key$ID
    temp$Analyse <- paste0("S",s,"R", r)
    
    RES <- bind_rows(RES, temp)
    
  }
  
}

RES %>% head()

RES.final <- RES %>% mutate(Num = sapply(str_split(Ind.ID, "_"), `[`, 3),
                            Num = ifelse(is.na(Num),sapply(str_split(Ind.ID, "_"), `[`, 2), Num ),
                            CAT = Ind.ID %>% str_remove(paste0("_", Num)),
                            Serie = str_sub(Analyse, 1,2),
                            Rep = str_sub(Analyse,3,4)
)

RES.final %>% head()

# Stats
RES.final %>% group_by(Ind.ID, CAT, Serie) %>% 
  summarise(freq.pop.1 = mean(pop.1)) %>% 
  group_by(CAT) %>% summarise(min = min(freq.pop.1),
                              max = max(freq.pop.1),
                              mean = mean(freq.pop.1),
                              sd = sd(freq.pop.1),
                              diff.min = 1 - min,
                              diff.max = 1 - max,
                              N = n(),
                              Nnorth = length(Ind.ID[freq.pop.1 >0.5]))


fig4b <- RES.final %>% group_by(Ind.ID, CAT, Serie) %>% 
  summarise(freq.pop.1 = mean(pop.1)) %>% #View()
  mutate(CAT = factor(CAT, levels = c("Pure_US","BC_US","F1.out", "F2.out", "BC_CAN","Pure_CAN"))) %>% 
  ggplot(aes(x=freq.pop.1, y = CAT)) +
  geom_vline(xintercept = 0.5, lty = "dashed") +
  #geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0.5,xmax=Inf),
  #          fill=scales::alpha("firebrick2", 1/150))+
  #geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=-Inf,xmax=0.5),
  #          fill= scales::alpha("dodgerblue1", 1/150))+
  geom_boxplot() +
  scale_y_discrete(labels = c("Pure southern", "Backcross southern", "F1 hybrids", "F2 hybrids", "Backcross northern", "Pure northern")) +
  labs(y = "Simulated genotype", x = "Membership probability") +
  #geom_jitter(aes(col = Serie)) +
  theme_bw() 

fig4b

ggsave(filename = here("02_Results/fig4b_Hybridvalidation.png"), plot = fig4b, 
       width = 7, height = 3 , units = "in",
       dpi = 300)

# Joining fig4 a and b
library(ggpubr)

fig4 <- ggarrange(fig4a + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")), 
                  ggarrange(fig4b + theme(axis.title = element_text(size = 10 ),
                                          plot.margin = margin(10, 10, 10, 20, "pt")),ncol = 2),
          labels = c("A", "B"),
          heights = c(5,4),
          nrow = 2)

fig4

ggsave(filename = here("02_Results/fig4_AssignmentValidation.png"), plot = fig4, 
       width = 7, height = 6 , units = "in",
       dpi = 300)

fig4.v2 <- ggarrange(fig4a.v2 + theme(axis.title = element_text(size = 10 ),
                                
                                plot.margin = margin(10, 10, 10, 20, "pt")), 
                  fig4b + theme(axis.title = element_text(size = 10 ),
                                          plot.margin = margin(10, 10, 10, 20, "pt")),
                  labels = c("A", "B"),
                  widths =c(3,5),
                  ncol = 2, align = "h")

fig4.v2

ggsave(filename = here("02_Results/fig4.v2_AssignmentValidation.png"), plot = fig4.v2, 
       width = 6, height = 3 , units = "in",
       dpi = 300)

# Assignation to references : random samples  -----------------------------------------------
library(assignPOP)

for(i in 21:100){

print(i)

# Randomized pop  
  
gl.rand <-  gl.ref   
pop(gl.rand) <- sample(pop(gl.ref), size = nInd(gl.rand)) 

gl2structure(x = gl.rand,   addcolumns =  pop(gl.rand),  outfile = paste0("gl.rand_", i, ".str"),   outpath = file.path("./02_Results/03_assignPOP/Random_ref"))

# Perform MCMC

rand.str <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/Random_ref",paste0("gl.rand_", i, ".str")))

assign.MC(rand.str, train.inds=c(0.7), train.loci=c(1),
            loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/03_assignPOP/Random_ref/",  paste0("gl.rand_", i),"/"),
          processors = 40)
# unload the memory
gc() 
}

# Create a big ref

accuMC.rand <- data.frame()

for(i in 1:100){
  
    print(i)
    
    accuMC <- accuracy.MC(paste0("./02_Results/03_assignPOP/Random_ref/", paste0("gl.rand_", i),"/"))
    accuMC <- accuMC %>% mutate(repetition = i)
    
    accuMC.rand <- bind_rows(accuMC.rand, accuMC)
    
  }
  

accuMC.rand <- accuMC.rand %>% pivot_longer(names(accuMC.rand) %>% str_subset("assign.rate"),
                                          names_to = "group", values_to = "assign.rate") %>% 
  mutate(group = group %>% str_remove("assign.rate."),
         assign.rate = as.numeric(as.character(assign.rate))
         
  ) %>% 
  filter(!is.na(assign.rate),
         group != "all"#,
         #str_detect(ref, "2REF")
  ) 


accuMC.rand %>% group_by(repetition) %>% 
  summarise(mean = mean(assign.rate)) %>% 
  ggplot(aes(y = mean, x = 1)) +
  geom_hline(yintercept = c(0.5), lty = "dashed") +
  geom_boxplot() +
  geom_jitter(height = 0) +
  #facet_grid(. ~ repetition, space = "free_x", scale = "free_x") +
  labs(x = "", y = "Assigment accuracy") +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  scale_fill_manual(values = c("gray100", "gray85", "gray70", "gray40"))+
  scale_y_continuous(limits = c(0,1)) +
  guides(fill=guide_legend(title="% of SNPs")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top")


accuMC.rand %>% group_by(repetition) %>% 
  summarise(mean.assign = mean(assign.rate)) %>% 
  ungroup() %>% summarise(mean = mean(mean.assign),
                          sd = sd(mean.assign),
                          min = min(mean.assign),
                          max = max(mean.assign),
                          median = median(mean.assign))

# Assignation to references : Adults assigments -----------------------------------------------

# Load the datasets
others.str    <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/", "gl.no.ref_june2021.str"))
ref.str       <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/", "gl.ref_june2021.str"))


others.str    <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/", "gl.no.ref_june2021_15.str"))
ref.str       <- assignPOP::read.Structure(file.path("./02_Results/03_assignPOP/", "gl.ref_june2021_15.str"))


# Perform the assigments
assign.X(x1=ref.str , x2=others.str, dir="./02_Results/03_assignPOP/Assignations_SVM/", model="svm")

#assign.X(x1=ref.str , x2=others.str, dir="./02_Results/03_assignPOP/Assignations_SVM_15/", model="svm")

#assignPOP::assign.X(x1=ref.str , x2=others.str, dir="./02_Results/03_assignPOP/Assignations_LDA/", model="lda")

# Load the result
assign.others <- read.table(file.path("./02_Results/03_assignPOP/Assignations_SVM/", "AssignmentResult.txt"), header = T) 
#assign.others <- read.table(file.path("./02_Results/03_assignPOP/Assignations_SVM_15/", "AssignmentResult.txt"), header = T) 

#assign.others <- read.table(file.path("./02_Results/03_assignPOP/Assignations_LDA/", "AssignmentResult.txt"), header = T) 


# graph

# Rapid check - no more 5Y in 5YZe
assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  group_by(new.NAFO, NAFO) %>% summarise(N = n())

# Stats

df.na <- data.frame(Sample = indNames(gl.america.samples),
                    NNA = na.gl.ind(gl.america.samples))

write_csv(df.na, "na.test.2021-11-17.csv")

bad.ID <-df.na %>% filter(NNA>=0.10) %>% pull(Sample)

df.na %>% View()

test20 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA, y = CAN, col = Sequencage) )+ geom_point() + 
  facet_wrap(~new.NAFO, nrow = 2)

test20

test30 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA, y = CAN, col = Sequencage) )+ geom_point() + 
  facet_wrap(~new.NAFO, nrow = 2)

test30

test15 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA, y = CAN, col = Sequencage) )+ geom_point() + 
  facet_wrap(~new.NAFO, nrow = 2)

test15



assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  group_by(new.NAFO) %>% summarise(NCAN = length(Ind.ID[pred.pop == "CAN"]),
                                   NUS = length(Ind.ID[pred.pop == "US"]),
                                   PropCAN = NCAN / (NCAN + NUS),
                                   PropUS = NUS / (NCAN + NUS),
                                   NCAN70 = length(Ind.ID[CAN >= .70]),
                                   NUS70 = length(Ind.ID[US >= .70]),
                                   PropCAN70 = NCAN70 / (NCAN70 + NUS70),
                                   PropUS70 = NUS70 / (NCAN70 + NUS70))

assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  group_by(Country) %>% summarise(NCAN = length(Ind.ID[pred.pop == "CAN"]),
                                   NUS = length(Ind.ID[pred.pop == "US"]),
                                   PropCAN = NCAN / (NCAN + NUS),
                                   PropUS = NUS / (NCAN + NUS),
                                  NCAN70 = length(Ind.ID[CAN >= .70]),
                                  NUS70 = length(Ind.ID[US >= .70]),
                                  PropCAN70 = NCAN70 / (NCAN70 + NUS70),
                                  PropUS70 = NUS70 / (NCAN70 + NUS70))

assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  group_by(new.NAFO) %>% summarise(NCAN = length(Ind.ID[pred.pop == "CAN"]),
                                   NUS = length(Ind.ID[pred.pop == "US"]),
                                   NCAN70 = length(Ind.ID[CAN >= .70]),
                                   NUS70 = length(Ind.ID[US >= .70]),
                                   PropCAN70 = NCAN70 / (NCAN70+NUS70),
                                   PropUS70 = NUS70 / (NCAN70+NUS70),
                                   Prop70 = ((NCAN70+NUS70)/(NCAN+NUS)))

assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  summarise(NCAN = length(Ind.ID[pred.pop == "CAN"]),
                                   NUS = length(Ind.ID[pred.pop == "US"]),
                                   NCAN75 = length(Ind.ID[CAN >= .75]),
                                   NUS75 = length(Ind.ID[US >= .75]),
                                   PropCAN75 = NCAN75 / (NCAN75+NUS75),
                                   PropUS75 = NUS75 / (NCAN75+NUS75),
                                   Prop75 = ((NCAN75+NUS75)/(NCAN+NUS)))


assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  group_by(Country) %>% 
  summarise(NCAN = length(Ind.ID[pred.pop == "CAN"]),
            NUS = length(Ind.ID[pred.pop == "US"]),
            NCAN75 = length(Ind.ID[CAN >= .75]),
            NUS75 = length(Ind.ID[US >= .75]),
            PropCAN75 = NCAN75 / (NCAN),
            PropUS75 = NUS75 / (NUS),
            Prop75 = ((NCAN75+NUS75)/(NCAN+NUS)))


fig5 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
    mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
           new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN", "Northern", 
                               ifelse(pred.pop == "US", "Southern", "ERROR"))) %>% 
  # Faire les stats pour les graphiques
  group_by(new.NAFO, new.pred.pop) %>% 
  summarise(N = n()) %>% 
  mutate(SUM = sum(N),
         freq = N / sum(N)) %>% 
  #Graphique
  ggplot(aes(x = new.NAFO, y = freq, fill = new.pred.pop)) +
  geom_bar(stat = "identity", position = "stack") +
  #coord_polar("y", start = 0) +
  # Ajouter les N
  scale_fill_manual(values = c("firebrick2", "dodgerblue1")) +
  geom_text(aes(y = -0.01, label = paste0("n=",SUM)), vjust = 1, col = "black", cex = 3) +
  labs(x = "NAFO zone", y = "Frequency") +
  scale_y_continuous(limit = c(-0.02,1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

fig5

ggsave(filename = here("02_Results/fig5_AdultAssignments.png"), plot = fig5, 
       width = 4, height = 3.5 , units = "in",
       dpi = 300)


 fig5a <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
   mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
          new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
          new.pred.pop = ifelse(pred.pop == "CAN", "Northern", 
                                ifelse(pred.pop == "US", "Southern", "ERROR"))) %>% 
   # Faire les stats pour les graphiques
   group_by(new.NAFO, new.pred.pop) %>% 
   summarise(N = n()) %>% 
   mutate(SUM = sum(N),
          freq = N / sum(N)) %>% 
   ggplot(aes(x = 1, y = freq, fill = new.pred.pop)) +
   geom_bar(width = , stat = "identity", color = "gray10", cex = 0.2) +
   coord_polar("y", start = 0) + 
   scale_fill_manual(values = c("firebrick2", "dodgerblue1")) +
   #scale_y_continuous(limits = c(-1,1)) +
   geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 4, col = "black", cex = 3) +
   facet_grid(.~new.NAFO) + theme_void() +
   theme(legend.title = element_blank(),
         legend.position = "bottom",
         plot.background = element_rect(fill = "white", colour = NA),
         plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))

fig5a

fig5a.v1 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
         new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN" & CAN >=.7, "Northern 70%", 
                        ifelse(pred.pop == "CAN" & CAN >=.5, "Northern 50%",
                        ifelse(pred.pop == "US" & US >=.7, "Southern 70%", 
                        ifelse(pred.pop == "US" & US >=.5, "Southern 50%",
                                "ERROR")))),
         new.pred.pop = factor(new.pred.pop, levels = c("Northern 50%", "Northern 70%", "Southern 70%", "Southern 50%"))) %>% 
  # Faire les stats pour les graphiques
  group_by(new.NAFO, new.pred.pop) %>% 
  summarise(N = n()) %>% 
  mutate(SUM = sum(N),
         freq = N / sum(N)) %>% 
  ggplot(aes(x = 1, y = freq, fill = new.pred.pop)) +
  geom_bar(width = , stat = "identity", color = "gray10", cex = 0.2) +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = c("salmon1","firebrick2",  "dodgerblue1", "skyblue2")) +
  #scale_y_continuous(limits = c(-1,1)) +
  geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 4, col = "black", cex = 3) +
  facet_grid(.~new.NAFO) + theme_void() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))

fig5a.v1

fig5a.v2 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
         new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN" & CAN >=.7, "Northern 70%", 
                               ifelse(pred.pop == "CAN" & CAN >=.5, "Unassigned",
                                      ifelse(pred.pop == "US" & US >=.7, "Southern 70%", 
                                             ifelse(pred.pop == "US" & US >=.5, "Unassigned",
                                                    "ERROR")))),
         new.pred.pop = factor(new.pred.pop, levels = c("Northern 70%", "Southern 70%", "Unassigned"))) %>% 
  # Faire les stats pour les graphiques
  group_by(new.NAFO, new.pred.pop) %>% 
  summarise(N = n()) %>% 
  mutate(SUM = sum(N),
         freq = N / sum(N)) %>% 
  ggplot(aes(x = 1, y = freq, fill = new.pred.pop)) +
  geom_bar(width = , stat = "identity", color = "gray10", cex = 0.2) +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = c("firebrick2",  "dodgerblue1", "gray")) +
  #scale_y_continuous(limits = c(-1,1)) +
  geom_text(aes(y = 0, label = paste0("n=",SUM)), vjust = 4, col = "black", cex = 3) +
  facet_grid(.~new.NAFO) + theme_void() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))

fig5a.v2

id.order <- assign.others  %>% arrange(CAN) %>% pull(Ind.ID)

fig5b <- assign.others %>% pivot_longer(c("CAN", "US"), names_to = "cluster", values_to = "q") %>% 
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
         new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN", "Northern", 
                               ifelse(pred.pop == "US", "Southern", "ERROR")),
         cluster = factor(cluster, levels = c("US", "CAN"))) %>% 


  #Graphique
  ggplot(aes(x = fct_relevel(Ind.ID, id.order), y = q, fill = cluster)) +
  #ggplot(aes(x =  reorder(ID_GQ, Qvalue, FUN = function(x) min(x)), y = Qvalue, fill = Qpop)) +
  geom_bar(stat = "identity", width = 1) +
  geom_hline(yintercept = 0.5, col = "gray10", lty = "dashed") +
  scale_fill_manual(values = c( "dodgerblue1", "firebrick2")) +

  #scale_fill_manual(values=col1)+
  #scale_fill_manual(breaks = c("European", "Canadian", "Undefined 1", "Undefined 2"), values = c("royalblue2", "firebrick2", "goldenrod2", "purple2"), aesthetics = "fill") +
  #scale_y_continuous(limits = c(0,1)) +
  facet_grid(. ~ new.NAFO, scale = "free", labeller = labeller(K = label_both)) + 
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        #strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        #strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(5, "pt"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10, unit = "pt") )

fig5b

library(ggpubr)

fig5.v2 <- ggarrange(fig5a + theme( plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
                                    panel.spacing = unit(5, "pt")), 
                     fig5b,
                  labels = c("A", "B"),
                  heights = c(1,2),
                  align = "v",
                  nrow = 2,
                  common.legend = TRUE,
                  legend = "none")


fig5.v2

ggsave(filename = here("02_Results/fig5.v2_AdultAssignments.png"), plot = fig5.v2, 
       width = 8, height = 4 , units = "in",
       dpi = 300)

#

fig5c <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
         new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN", "Northern", 
                               ifelse(pred.pop == "US", "Southern", "ERROR"))) %>% 
  
  
  #Graphique
  ggplot(aes(x = CAN, fill = pred.pop)) +
     geom_vline(xintercept = 0.5, col = "gray10", lty = "dashed") +
  geom_histogram(binwidth = 0.1, closed = "right") +
 # geom_density(alpha = 0.5)+
   geom_hline(yintercept = 0, col = "gray10") +
  scale_fill_manual(values = c( "firebrick2", "dodgerblue1")) +
  
  #scale_fill_manual(values=col1)+
  #scale_fill_manual(breaks = c("European", "Canadian", "Undefined 1", "Undefined 2"), values = c("royalblue2", "firebrick2", "goldenrod2", "purple2"), aesthetics = "fill") +
  #scale_y_continuous(limits = c(0,1)) +
  facet_grid(. ~ new.NAFO, labeller = labeller(K = label_both)) + 
  labs(x="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.ticks = element_line(colour = "black"),
        #strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        #strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(5, "pt"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.title = element_blank(),
        #axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10, unit = "pt") )

fig5c


fig5c.v1 <- assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  mutate(new.NAFO = ifelse(new.NAFO == "5YZe", "5Ze", new.NAFO) ,
         new.NAFO = factor(new.NAFO, levels = c("6AB", "5Zw", "5Ze", "4WX", "3Ps4Vn", "4T", "4RS", "3K")),
         new.pred.pop = ifelse(pred.pop == "CAN" & CAN >=.7, "Northern 70%", 
                               ifelse(pred.pop == "CAN" & CAN >=.5, "Northern 50%",
                                      ifelse(pred.pop == "US" & US >=.7, "Southern 70%", 
                                             ifelse(pred.pop == "US" & US >=.5, "Southern 50%",
                                                    "ERROR")))),
         new.pred.pop = factor(new.pred.pop, levels = c("Northern 50%", "Northern 70%", "Southern 70%", "Southern 50%"))) %>% 
  
  
  #Graphique
  ggplot(aes(x = CAN, fill = new.pred.pop)) +
  geom_vline(xintercept =c(0.5), col = c("gray10"), lty = "dashed") +
  geom_vline(xintercept =c(0.3, 0.7), col = c("gray70"), lty = "dashed") +
    geom_histogram(binwidth = 0.1, closed = "right") +
  # geom_density(alpha = 0.5)+
  geom_hline(yintercept = 0, col = "gray10") +
  scale_fill_manual(values = c("salmon1","firebrick2",  "dodgerblue1", "skyblue2")) +
  
  #scale_fill_manual(values=col1)+
  #scale_fill_manual(breaks = c("European", "Canadian", "Undefined 1", "Undefined 2"), values = c("royalblue2", "firebrick2", "goldenrod2", "purple2"), aesthetics = "fill") +
  #scale_y_continuous(limits = c(0,1)) +
  facet_grid(. ~ new.NAFO, labeller = labeller(K = label_both)) + 
  labs(x="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.ticks = element_line(colour = "black"),
        #strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        #strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(5, "pt"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.title = element_blank(),
        #axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10, unit = "pt") )

fig5c.v1

fig5.v3 <- ggpubr::ggarrange(fig5a + theme( #plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
                                  panel.spacing = unit(20, "pt")), 
                     fig5c,
                     labels = c("A", "B"),
                     heights = c(1,2),
                     align = "v",
                     nrow = 2,
                     common.legend = TRUE,
                     legend = "none")


fig5.v3


fig5.v4 <- ggpubr::ggarrange(fig5a.v1 + theme( #plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
  panel.spacing = unit(20, "pt")), 
  fig5c.v1,
  labels = c("A", "B"),
  heights = c(1,2),
  align = "v",
  nrow = 2,
  common.legend = TRUE,
  legend = "none")


fig5.v4
ggsave(filename = here("02_Results/fig5.v4_AdultAssignments.png"), plot = fig5.v4, 
       width = 8, height = 4 , units = "in",
       dpi = 300)


# Is there a correlation between structure and assigment? NO! 

test <- assign.others %>% left_join(K2.america %>% filter(Qpop == "Q1"),
                            by = c("Ind.ID" = "ID_GQ"))


cor.test(test$US, test$Qvalue)

# PCAdapt -----------------------------------------------------------------

library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")

# Convertion to plink .bed format

ID.global <- indNames(gl.global.samples) 
loc.global <- data.frame(ID = locNames(gl.global.samples))

write.csv(loc.global, file.path(here::here(),  "00_Data", "Loc.global.csv"), 
          row.names = F, quote = F)

vcf.path.global <- file.path(here::here(), "../../../Extra_Storage/Projets/Maquereaux_RAD/00_Data/06a_Filtering.denovo/global.r0.9.MAF0.1/06e_AllIndividuals/populations.snps.vcf")

cmd1a <- paste("--vcf", vcf.path.global, 
              #"--recode",
              "--plink-tped",
              paste("--indv", ID.global, collapse = " "),
             # "--snps", file.path(here::here(),  "00_Data", "Loc.global.csv"),
              "--out", file.path(here::here(), "00_Data/",  paste0("populations.global.", nrow(loc.global), "snps.", length(ID.global), "ind.vcf"))
)

cmd1a

A1 <- system2("vcftools", cmd1a, stdout=T, stderr=T)
A1 %>% tail()


ID.america <- indNames(gl.america.samples) 
loc.america <- data.frame(ID = locNames(gl.america.samples))

write.csv(loc.america, file.path(here::here(),  "00_Data", "Loc.america.csv"), 
          row.names = F, quote = F)

vcf.path.america <- file.path(here::here(), "../../../Extra_Storage/Projets/Maquereaux_RAD/00_Data/06a_Filtering.denovo/america.r0.9.MAF0.1/06e_AllIndividuals/populations.snps.vcf")

cmd1b <- paste("--vcf", vcf.path.america, 
              #"--recode",
              "--plink-tped",
              paste("--indv", ID.america, collapse = " "),
              # "--snps", file.path(here::here(),  "00_Data", "Loc.america.csv"),
              "--out", file.path(here::here(), "00_Data/",  paste0("populations.america.", nrow(loc.america), "snps.", length(ID.america), "ind.vcf"))
)

cmd1

A1 <- system2("vcftools", cmd1b, stdout=T, stderr=T)
A1 %>% tail()


# Correct the map file to remove chromosomes (change to 0)
global.tped.plink <- read.delim("00_Data/populations.global.10832snps.629ind.vcf.tped",
                         header = F)
global.tped.plink[1:6, 1:6]
global.tped.plink$V1 <- 0

write.table(global.tped.plink, 
            file = "00_Data/populations.global.10832snps.629ind.vcf.tped",
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


america.tped.plink <- read.delim("00_Data/populations.america.10703snps.557ind.vcf.tped",
                         header = F)
head(america.tped.plink)
america.ped.plink$V1 <- 0

write.table(america.tped.plink, 
            file = "00_Data/populations.america.10703snps.557ind.vcf.tped",
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Make BED

cmd2a <- paste("--tfam", "./00_Data/populations.america.10703snps.557ind.vcf.tfam", 
               "--tped", "./00_Data/populations.america.10703snps.557ind.vcf.tped", 
               "--make-bed", 
               "--out", "./00_Data/populations.america.10703snps.557ind"

)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a


cmd2b <- paste("--tfam", "./00_Data/populations.global.10832snps.629ind.vcf.tfam", 
               "--tped", "./00_Data/populations.global.10832snps.629ind.vcf.tped", 
               "--make-bed", 
               "--out", "./00_Data/populations.global.10832snps.629ind"
               
)

A2b <- system2("plink", cmd2b, stdout=T, stderr=T)
A2b

# Read .bed in PCAadapt
pca.america_genotype  <- read.pcadapt("00_Data/populations.america.10703snps.557ind.bed",
                                          type = "bed")

pca.global_genotype   <- read.pcadapt("00_Data/populations.global.10832snps.629ind.bed",
                                         type = "bed")

# file to get the real snps name

pcadapt.america.snp <- read.delim("00_Data/populations.america.10703snps.557ind.bim",
                                header = F) %>% pull(V2)

pcadapt.global.snp <- read.delim("00_Data/populations.global.10832snps.629ind.bim",
                                 header = F) %>% pull(V2)


# Run pcadapt

K.init <- 30

pcadapt.america <- pcadapt(pca.america_genotype , K =K.init)
pcadapt.global  <- pcadapt(pca.global_genotype , K = K.init)

# Check screeplot

plot(pcadapt.america, option = "screeplot") 
plot(pcadapt.global, option = "screeplot")

# Check structure

america.pop <-  data_frame(ID_GQ = indNames(gl.america.samples)) %>% 
  left_join(pop.data) %>% 
  pull(Country)
  
plot(pcadapt.america, option = "scores", pop = factor(america.pop)) 


global.pop <-  data_frame(ID_GQ = indNames(gl.global.samples)) %>% 
  left_join(pop.data) %>% 
  pull(Country)

plot(pcadapt.global, option = "scores", pop = factor(global.pop)) 

# K = 2 pour america et K = 1 pour global

pcadapt.america.k2  <- pcadapt(pca.america_genotype , K = 2)
pcadapt.global.k1   <- pcadapt(pca.global_genotype , K = 1)

plot(pcadapt.america.k2, option = "manhattan")
plot(pcadapt.global.k1, option = "manhattan")

hist(pcadapt.america.k2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
hist(pcadapt.global.k1$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(pcadapt.america.k2, option = "qqplot")
plot(pcadapt.global.k1, option = "qqplot")

# Statistics
#x$pvalues 
alpha <- 0.05

qval.america.k2 <- qvalue::qvalue(pcadapt.america.k2$pvalues)$qvalues
outliers.america.k2 <- which(qval.america.k2 < alpha)
length(outliers.america.k2)

qval.global.k1 <- qvalue::qvalue(pcadapt.global.k1$pvalues)$qvalues
outliers.global.k1 <- which(qval.global.k1 < alpha)
length(outliers.global.k1)

out.america.PCA <- pcadapt.america.snp[outliers.america.k2] %>% str_replace(":", "_")
out.global.PCA <- pcadapt.global.snp[outliers.global.k1]%>% str_replace(":", "_")

# Bayescan ----------------------------------------------------------------

library(dartR)
?gl2bayescan

# Data conversion

gl.america.samples %>% pop() %>% table()
gl.america.samples@other$loc.metrics.flags$monomorphs <- FALSE

gl2bayescan(gl.america.samples, outfile = "00_Data/bayescan_america.txt", outpath = ".")
gl2bayescan(gl.america.samples, outfile = "00_Data/bayescan_america_NAFO.txt", outpath = ".")


gl.global.samples %>% pop() %>% table()
gl.global.samples@other$loc.metrics.flags$monomorphs <- FALSE


gl2bayescan(gl.global.samples, outfile = "00_Data/bayescan_global.txt", outpath = ".")

# Run bayescan

cmd <- paste("./00_Data/bayescan_america.txt",
             "-snp",
             "-od", "./02_Results/05_bayescan", #Output directory file
             "-threads", 20
             )

cmd

A <- system2("bayescan", cmd, stdout = T, stderr = T)


cmd <- paste("./00_Data/bayescan_america_NAFO.txt",
             "-snp",
             "-od", "./02_Results/05_bayescan", #Output directory file
             "-threads", 20
)

cmd

A <- system2("bayescan", cmd, stdout = T, stderr = T)


cmd <- paste("./00_Data/bayescan_global.txt",
             "-snp",
             "-od", "./02_Results/05_bayescan", #Output directory file
             "-threads", 20
)

cmd

system2("bayescan", cmd, stdout = T, stderr = T)


# Extracting results


source("~/Documents/Programs/BayeScan/R_functions/plot_R.r")

plot_bayescan("02_Results/05_bayescan/bayescan_america_fst.txt", FDR = 0.05)
plot_bayescan("02_Results/05_bayescan/bayescan_america_NAFO_fst.txt", FDR = 0.05)
plot_bayescan("02_Results/05_bayescan/bayescan_global_fst.txt", FDR = 0.05)


plot_bayescan %>% View()

bayes.global <- read.table("02_Results/05_bayescan/bayescan_global_fst.txt")
head(bayes.global )
hist(bayes.global$fst)

outliers.global.fst <- bayes.global %>% filter(qval < 0.05) %>% row.names() %>% as.numeric()
out.global.fst <- locNames(gl.global.samples)[outliers.global.fst]

intersect(out.global.PCA, out.global.fst)

bayes.america <- read.table("02_Results/05_bayescan/bayescan_america_NAFO_fst.txt")
head(bayes.america)
hist(bayes.america$fst)

outliers.america.fst <- bayes.america %>% filter(qval < 0.05) %>% row.names() %>% as.numeric()

out.america.fst <- locNames(gl.america.samples)[outliers.america.fst]

intersect(out.america.PCA, out.america.fst)

# Neutral analysis --------------------------------------------------------

#out.global <-  c(outt.global.fst) %>% unique()
out.global <-  c(out.global.fst, out.global.PCA) %>% unique()
length(out.global)

gl.global.samples.neutral  <- gl.global.samples[,locNames(gl.global.samples) %nin% out.global]
gl.global.samples.outliers <- gl.global.samples[,locNames(gl.global.samples) %in% out.global]


out.america <-  c(out.america.fst, out.america.PCA) %>% unique()
length(out.america)

gl.america.samples.neutral  <- gl.america.samples[,locNames(gl.america.samples) %nin% out.america]
gl.america.samples.outliers  <- gl.america.samples[,locNames(gl.america.samples) %in% out.america]


#PCA

#save(list = c("pca.america.neutral","pca.america.outliers", 
#              "pca.global.neutral","pca.global.outliers",
#              "pca.ref.neutral", "pca.ref.outliers",
#                "out.america", "out.global"),
#     file = here("02_Results/01_PCA/PCA_neutralvsoutliers.data"))

load(here("02_Results/01_PCA/PCA_neutralvsoutliers.data"))

# On neutral also

gl.ref.neutral  <- gl.ref[,locNames(gl.ref) %nin% out.america]
gl.ref.outliers  <- gl.ref[,locNames(gl.ref) %in% out.america]


pca.ref.neutral <- glPca(gl.ref.neutral, center = TRUE, scale = FALSE, 
                         parallel = T, n.core = 40, nf = 1000)


pca.ref.outliers <- glPca(gl.ref.outliers, center = TRUE, scale = FALSE, 
                          parallel = T, n.core = 40, nf = 1000)





pca.global.outliers <- glPca(gl.global.samples.outliers, center = TRUE, scale = FALSE, 
                              parallel = T, n.core = 40, nf = 1000)

pca.global.outliers$eig[1]/sum(pca.global.outliers$eig) # proportion of variation explained by 1st axis
pca.global.outliers$eig[2]/sum(pca.global.outliers$eig)# proportion of variation explained by 2nd axis 

res.PCA.global.outliers.table <- data.frame(Sample = row.names(pca.global.outliers$scores),
                                   score = pca.global.outliers$scores[,1:10]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.global.outliers.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                   Mean.score.Axis.2 = mean(score.PC2))

res.PCA.global.outliers.table <- res.PCA.global.outliers.table %>% left_join(mean.axis )

figS.global.outliers.pca <- res.PCA.global.outliers.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Site)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  
 # geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
  #             size = 0.2, alpha = 0.5)+
  geom_point(aes(), alpha = 0.5, size = 2) + 
  
  
  scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PC1 - 4.7 %", y = "PC2 - 1.6 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS.global.outliers.pca


pca.global.neutral <- glPca(gl.global.samples.neutral, center = TRUE, scale = FALSE, 
                             parallel = T, n.core = 40, nf = 1000)

pca.global.neutral$eig[1]/sum(pca.global.neutral$eig) # proportion of variation explained by 1st axis
pca.global.neutral$eig[2]/sum(pca.global.neutral$eig)# proportion of variation explained by 2nd axis 

res.PCA.global.neutral.table <- data.frame(Sample = row.names(pca.global.neutral$scores),
                                            score = pca.global.neutral$scores[,1:10]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.global.neutral.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                            Mean.score.Axis.2 = mean(score.PC2))

res.PCA.global.neutral.table <- res.PCA.global.neutral.table %>% left_join(mean.axis )

figS.global.neutral.pca <- res.PCA.global.neutral.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Site)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  
#  geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
 #              size = 0.2, alpha = 0.5)+
  geom_point(aes(), alpha = 0.5, size = 2) + 
  
  
  scale_colour_manual(values = c("chocolate4", "chartreuse4", "firebrick2", "dodgerblue1")) +
  labs(x = "PC1 - 0.6 %", y = "PC2 - 0.3 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS.global.neutral.pca



pca.america.outliers <- glPca(gl.america.samples.outliers, center = TRUE, scale = FALSE, 
                             parallel = T, n.core = 40, nf = 1000)

pca.america.outliers$eig[1]/sum(pca.america.outliers$eig) # proportion of variation explained by 1st axis
pca.america.outliers$eig[2]/sum(pca.america.outliers$eig)# proportion of variation explained by 2nd axis 

res.PCA.america.outliers.table <- data.frame(Sample = row.names(pca.america.outliers$scores),
                                            score = pca.america.outliers$scores[,1:10]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.america.outliers.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                            Mean.score.Axis.2 = mean(score.PC2))

res.PCA.america.outliers.table <- res.PCA.america.outliers.table %>% left_join(mean.axis )

figS.america.outliers.pca <- res.PCA.america.outliers.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Site)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  
  #geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
  #             size = 0.2, alpha = 0.5)+
  geom_point(aes(), alpha = 0.5, size = 2) + 
  
  
  scale_colour_manual(values = c("firebrick2", "dodgerblue1")) +
  labs(x = "PC1 - 2.9 %", y = "PC2 - 2.7 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS.america.outliers.pca

pca.america.neutral <- glPca(gl.america.samples.neutral, center = TRUE, scale = FALSE, 
                            parallel = T, n.core = 40, nf = 1000)

pca.america.neutral$eig[1]/sum(pca.america.neutral$eig) # proportion of variation explained by 1st axis
pca.america.neutral$eig[2]/sum(pca.america.neutral$eig) # proportion of variation explained by 2nd axis 

res.PCA.america.neutral.table <- data.frame(Sample = row.names(pca.america.neutral$scores),
                                           score = pca.america.neutral$scores[,1:10]) %>% 
  left_join(pop.data) %>% #View()
  mutate(Site = ifelse(new.NAFO == "BOB" , "NEA - Bay of Biscay", 
                       ifelse(new.NAFO == "Groenland" , "NEA - Greenland",
                              ifelse(Country == "Canada", "NWA - Canada",
                                     ifelse(Country == "US", "NWA - US", NA)))))

mean.axis <- res.PCA.america.neutral.table %>% group_by(Site) %>% summarise(Mean.score.Axis.1 = mean(score.PC1),
                                                                           Mean.score.Axis.2 = mean(score.PC2))

res.PCA.america.neutral.table <- res.PCA.america.neutral.table %>% left_join(mean.axis )

figS.america.neutral.pca <- res.PCA.america.neutral.table  %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Site)) +
  #  stat_ellipse(aes(col = REF)) +
  #ggforce::geom_mark_ellipse(aes(label = Pop.RefAdults, col = Pop.RefAdults, filter = !is.na(Ref.adults)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  
  #geom_segment(aes(x = Mean.score.Axis.1, y = Mean.score.Axis.2, xend = score.PC1, yend= score.PC2),
  #             size = 0.2, alpha = 0.5)+
  geom_point(aes(), alpha = 0.5, size = 2) + 
  
  
  scale_colour_manual(values = c("firebrick2", "dodgerblue1")) +
  labs(x = "PC1 - 0.3 %", y = "PC2 - 0.3 %") + 
  #scale_x_continuous(limits = c(-2, 10))
  #directlabels::geom_dl(aes(label = pop.final)), method = "smart.grid" ) +
  #facet_grid(dataset ~ nloci.MEM, scale = "free") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #strip.text = element_text(angle = 90),
    #panel.grid = element_blank(),
    #panel.spacing = unit(0, "cm"),
    #panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
#axis.title.x = element_blank())

figS.america.neutral.pca





ggsave(filename = here("02_Results/figS_PCA_global_neutral.png"), plot = figS.global.neutral.pca, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figS_PCA_global_outliers.png"), plot = figS.global.outliers.pca, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figS_PCA_america_neutral.png"), plot = figS.america.neutral.pca, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figS_PCA_america_outliers.png"), plot = figS.america.neutral.pca, 
       width = 3, height = 3 , units = "in",
       dpi = 300)


# Joining fig4 a and b

figS.pca <- ggpubr::ggarrange(figS.global.neutral.pca+ theme(axis.title = element_text(size = 10 ),
                                          
                                          plot.margin = margin(10, 10, 10, 20, "pt")) + labs (title = "Neutral SNPs"),
                                            
                      figS.america.neutral.pca + theme(axis.title = element_text(size = 10 ),
                                                        
                                                        plot.margin = margin(10, 10, 10, 20, "pt"))  + labs (title = "Neutral SNPs"),
                       
                      figS.global.outliers.pca + theme(axis.title = element_text(size = 10 ),
                                                       
                                                       plot.margin = margin(10, 10, 10, 20, "pt"))  + labs (title = "Outlier SNPs"),
                      
                      
                      figS.america.outliers.pca + theme(axis.title = element_text(size = 10 ),
                                                                           
                                                                           plot.margin = margin(10, 10, 10, 20, "pt"))  + labs (title = "Outliers SNPs"),
                                          
                       labels = LETTERS ,common.legend =T, legend = "bottom" , ncol = 2,nrow = 2
                       )

figS.pca

figS2.pca <- ggpubr::ggarrange(figS.global.neutral.pca + theme(axis.title = element_text(size = 10 ),
                             
                             plot.margin = margin(10, 10, 10, 20, "pt")),
                        figS.america.neutral.pca + theme(axis.title = element_text(size = 10 ),
                             
                             plot.margin = margin(10, 10, 10, 20, "pt")),
          labels = c("A", "B"),
          common.legend = T, legend = "bottom", 
          ncol = 2)
figS2.pca

figS.global.neutral.pca$layers[[3]]$aes_params$size <- 1
figS.america.neutral.pca$layers[[3]]$aes_params$size <- 1
figS4.PCA.neutral$layers[[3]]$aes_params$size <- 1


figS3.pca <- ggpubr::ggarrange(figS.global.neutral.pca + theme(axis.title = element_text(size = 10 ),
                                                       
                                                       plot.margin = margin(10, 10, 10, 20, "pt")) +
                                 guides(colour = guide_legend(override.aes = list(size=3))),
                       figS.america.neutral.pca + theme(axis.title = element_text(size = 10 ),
                                                        
                                                        plot.margin = margin(10, 10, 10, 20, "pt")),
                       figS4.PCA.neutral + theme(axis.title = element_text(size = 10 ),
                                                        
                                                        plot.margin = margin(10, 10, 10, 20, "pt")),
                       labels = LETTERS,
                      # labels = c("A", "B"),
                       common.legend = T, legend = "bottom", 
                       ncol = 3)
figS3.pca


ggsave(filename = here("02_Results/figSoutliers_PCA_v3.png"), plot = figS3.pca, 
       width = 8, height = 3 , units = "in", bg = "white",
       dpi = 300)



# Others ------------------------------------------------------------------



df.ref[,outliers_pcadapt]

gi.ref[,outliers_pcadapt] %>% 

nLoc(gi.transcripto)

indNames(gi.)

Mydata2 <- genind2hierfstat(gi.transcripto) 

basicstat <- basic.stats(Mydata2, diploid = TRUE, digits = 2) 
names(basicstat)  

x <- indpca(Mydata2) 
plot(x, cex = 0.7)


#

## Library dartR

pop(gl.america.samples) <- data.frame(ID_GQ = indNames(gl.america.samples)) %>% 
  left_join(pop.data) %>%  pull(Country) 


library(dartR)

pca.dartR <- gl.pcoa(gl.america.samples, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

gl.pcoa.plot(pca.dartR,gl.america.samples, xaxis = 1 , yaxis = 2)


gl.america.ind <- gl.america.samples
pop(gl.america.ind) <- indNames(gl.america.ind)


res1 <- gl.dist.pop(gl.america.samples, method = "euclidean")

str(res1)
str(res2)

res2 <- stats::dist(tab(gl.america.samples,  NA.method = c("asis")), method = "euclidean", diag = F, upper = F)

cor(res1, res2)

#dist.dart.1 <- gl.dist.pop(gl.america.ind[1:10], method="euclidean")

pcoa.dart.5 <-  gl.pcoa(gl.america.samples, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5, correction = "lingoes")

gl.pcoa.plot(pcoa.dart.5,gl.america.samples, xaxis = 1 , yaxis = 2)

#


dist1 <- dist(tab(gl.america.samples,  NA.method = c("asis")), method = "euclidean", diag = T, upper = T)

dist2 <- dist(tab(gl.america.samples,  NA.method = c("mean")), method = "euclidean", diag = T, upper = T)

dist3 <- dist(tab(gl.america.samples,  NA.method = c("zero")), method = "euclidean", diag = T, upper = T)

mat4 <- tab(gl.global.samples,  NA.method = c("asis"))

mat4 

mat4 <- apply(mat4, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

vecNA <- apply(mat4, 1, function(x) length( is.na(x)[is.na(x) == TRUE]))

data.frame(ID_GQ = row.names(mat4 ), n.NA = vecNA) %>% left_join(pop.data) %>% 
  ggplot(aes(x= new.NAFO, y = n.NA)) + geom_violin()


sum(is.na(dist4)) # No NAs

length( is.na(mat4[1,])[is.na(mat4[1,]) == TRUE])

table(is.na(mat4[1,]))

dist4

dist4 <- dist(mat4, method = "euclidean", diag = T, upper = T)


pcoa1.dartR <- gl.pcoa(dist1, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

pcoa1.dartR %>% str()

gl.pcoa.plot(pcoa1.dartR ,gl.america.samples, xaxis = 1 , yaxis = 2)


pcoa2.dartR <- gl.pcoa(dist2, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

gl.pcoa.plot(pcoa2.dartR ,gl.america.samples, xaxis = 1 , yaxis = 2)


pcoa3.dartR <- gl.pcoa(dist3, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)


gl.pcoa.plot(pcoa3.dartR ,gl.america.samples, xaxis = 1 , yaxis = 2)


pcoa4.dartR <- gl.pcoa(dist4, nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

cov4 <- cov(mat4)

str(cov4)

eig4 <- eigen(cov4)


gl.pcoa.plot(pcoa4.dartR ,gl.global.samples, xaxis = 1 , yaxis = 2)


cor.test(dist2, dist3)
cor.test(dist2, dist4)
cor.test(dist3, dist4)

cov

mme.pca <- eigen(adegenet::tab(gl.america.samples)) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments



gi.america.samples

library(QuickPop)

loc.NA05  <- filter.MAF.NA(gi.america.samples, MAF.trs = 0.05, NA.trs = 0.05)


pca.dartR.NA05 <- gl.pcoa(gl.america.samples[,loc.NA05], nfactors = 5, parallel = TRUE, n.cores = 16, verbose = 5)

gl.pcoa.plot(pca.dartR,gl.america.samples, xaxis = 1 , yaxis = 2)


hist(adegenet::glNA(gl.america.samples, alleleAsUnit = F))

pop(glind7) <- indNames(glind7)# redefine the population information
gl.dist.pop(glind7[1:7,], method="euclidean")


?glNA

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

library(LEA)


data("tutorial")
str(tutorial.R)

tutorial.R[1:10, 1:10]

df.america.samples[1:10, 1:10]
?LEA::pca


plot(pc$projections)

data.frame(ID_GQ =  indNames(gl.america.samples),
           pc1 = pc$projections[,1],
           pc2 = pc$projections[,2]
           ) %>% left_join(pop.data) %>% 
  ggplot(aes(x = pc1, y = pc2, col = Country)) + geom_point()
  geom_bar(stat = "identity") +
  facet_grid(. ~ new.NAFO, space = "free", scale = "free")
str(pc$projections[,1])

write.lfmm(df.america.samples, "genotypes.lfmm")

write.geno(df.america.samples, "genotypes.geno")

project = NULL
project = snmf("genotypes.geno",
               K = 1:8,
               entropy = TRUE,
               repetitions = 10,
               project = "new",
               CPU = 20)

plot(project, col = "blue", pch = 19, cex = 1.2)


best = which.min(cross.entropy(project, K = 2))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)
res <- LEA::Q(project, K = 2, run = best)

str(res)

data.frame(ID_GQ =  indNames(gl.america.samples),
           Qvalue = res[,1]) %>% left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Qvalue, fill = Country)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ new.NAFO, space = "free", scale = "free")
  
str(project)



# RDA on sex

pop.data %>% pull(Sex) %>% unique()

ID.withSex <- pop.data %>% filter(Sex %in% c("M", "F")) %>% pull(ID_GQ)

gl.sex <- gl.america.samples[indNames(gl.america.samples) %in% ID.withSex,]
pop(gl.sex) <- data.frame(ID_GQ = indNames(gl.sex)) %>% left_join(pop.data) %>% pull(Sex)

pop(gl.sex) %>% table()

tab(gl.sex)

## Null model
RDA0 <- rda(freq.MAF.final.Gen_ZONE_FG ~ 1,  Variables) 

library(vegan)

## Full model
RDAfull <- vegan::rda(tab(gl.sex) ~ pop(gl.sex))
RDAfull

anova(RDAfull)


res <- as.data.frame(scores(RDAfull, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDAfull, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

arrow.df <- data.frame(name = scores(RDAfull, display="bp", choices=1, scaling=1) %>% row.names(),
                       x0 = 0,
                       y0 = 0,
                       xmax = scores(RDAfull, display="bp", choices=1, scaling=1) %>% as.vector(),
                       ymax = scores(RDAfull, display="bp", choices=2, scaling=1) %>% as.vector()) %>% 
  mutate(name = ifelse(name == "`env$T_Bott_Winter`", "T_Bott_Winter", name))

arrow.factor = 2

library(ggrepel)

gRDA1 <- res %>% ggplot(aes(x = RDA1)) +
  #facet_wrap(~RegionAssesment) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  #geom_point(data = res2, colour = "grey70", cex = 1) +
  geom_density()+
  #scale_color_manual(values = c("blue", "orange", "red"))+
  #geom_segment(data = arrow.df, aes(x = x0, y = y0, xend = xmax * arrow.factor, yend= ymax * arrow.factor),
  #             arrow = arrow(), size = 0.5)+
  #ggrepel::geom_label_repel(data = arrow.df , aes(label = name, x = xmax *arrow.factor, y = ymax * arrow.factor), size = 2, max.overlaps = 20
  #) +
  #geom_text(aes(label = Gen_ZONE_FG %>% str_remove("SFA-|NAFO-"), col = RegionAssesment, cex = 1)) +
  #geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-|NAFO-"),  col = RegionAssesment), size = 3, max.overlaps = 20
  #)# +
  
  labs(x=c("RDA 1"), y=c("RDA 2")) + 
  annotate("text", label = paste("R²-adjusted =", RsquareAdj(RDAfull)$adj.r.squared %>% as.numeric() %>%  round(digits = 6)),
           x = 0.4, y =0.4, vjust = "inward", hjust = "inward") +
  
  theme_bw()  
gRDA1

hist(RDAfull$CCA$v)
dimnames(RDAfull$CCA$v)[[1]][RDAfull$CCA$v > 0.05]

# 75281_80


