# Info --------------------------------------------------------------------

# Test on impacts of HI vs NS 
# Mackerel ddrad dataset - only adults
# Audrey Bourret
# 2022-05-31
#

# Library -----------------------------------------------------------------

# library
library(here)
library(tidyverse)
library(readxl)
library(adegenet)
library(dartR)
library(assignPOP)

`%nin%` = Negate(`%in%`)

#install.packages("BiocManager")
#BiocManager::install(c("SNPRelate", "qvalue")) 

# MetaData --------------------------------------------------------------------

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

# Gen data

# ddRad data --------------------------------------------------------------

load(here("00_Data/gen.america.878samples.10703snps.data")) 

pop.data %>% pull(Sequencage) %>% unique()

ID.samples.NS.ame <- pop.data %>% filter(Sequencage == "NS",
                                      Dev_Stage == "Adult",
                                      SNP_panel == "Yes",
                                      Country %in% c("Canada", "US")) %>% 
  pull(ID_GQ)

ID.samples.NS.ame %>% str_subset("NS_SS0259")

gl.america.NS.samples <- gl.america.10[indNames(gl.america.10) %in% ID.samples.NS.ame]

pop(gl.america.NS.samples) <- data.frame(ID_GQ = indNames(gl.america.NS.samples)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO
  ) %>% pull(new.pop) 

pop(gl.america.NS.samples) %>% table()

gl.america.NS.samples

# Ref data

ID.samples.ref <- pop.data %>% filter(Cat_Sample == "Sample",
                                      Dev_Stage == "Adult",
                                      SNP_panel == "Yes",
                                      Country %in% c("Canada", "US"),
                                      REF_assign %in% c("CAN", "US")) %>% 
  pull(ID_IBIS)

paste0("NS_",ID.samples.ref)

gl.NS.ref <- gl.america.10[indNames(gl.america.10) %in% paste0("NS_",ID.samples.ref)]


pop(gl.NS.ref) <- data.frame(ID_GQ = indNames(gl.NS.ref)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = Country) %>% pull(new.pop) %>% str_replace("Canada", "CAN")

pop(gl.NS.ref) %>% table()

# No ref data (for assignment)


gl.NS.no.ref <- gl.america.NS.samples[indNames(gl.america.NS.samples) %nin% paste0("NS_",ID.samples.ref)]


pop(gl.NS.no.ref) <- data.frame(ID_GQ = indNames(gl.NS.no.ref)) %>% 
  left_join(pop.data) %>% 
  mutate(new.pop = new.NAFO) %>% pull(new.pop) 

pop(gl.NS.no.ref) %>% table()


# ASSIGN ------------------------------------------------------------------

  # Save as structure
  
  gl2structure(x = gl.NS.ref,   addcolumns =  pop(gl.NS.ref),  outfile = paste0("gl.ref.NS.str"),   outpath = file.path("./02_Results/07_NS_tests/"))
  
  gl2structure(x = gl.NS.no.ref,   addcolumns =  pop(gl.NS.no.ref),  outfile =  paste0("gl.no.ref.NS.str"),   outpath = file.path("./02_Results/07_NS_tests/"))
 
  # Load files to perform assignPOP
  
  others.str    <- assignPOP::read.Structure(file.path("./02_Results/07_NS_tests/", paste0("gl.no.ref.NS.str")))

  ref.str  <- assignPOP::read.Structure(file.path("./02_Results/07_NS_tests/",  paste0("gl.ref.NS.str")))
 
  # Perform MCMC

  assign.MC(ref.str, train.inds=c(0.7), train.loci=c(1),
            loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/07_NS_tests/",  "gl.ref.NS_MCMC/"),
            processors = 30)
  # unload the memory
  gc() 
  
  # Perform the assigments
  
  assign.X(x1=ref.str , x2=others.str, dir= paste0("./02_Results/07_NS_tests/", "gl.NS_assign/"), model="svm", mplot = T)


  
  
# Load MCMC results  ---------------------------------------------------------

accuMC.all <- accuracy.MC(paste0("./02_Results/07_NS_tests/",  "gl.ref.NS_MCMC/"))
  

accuMC.all <- accuMC.all %>% pivot_longer(names(accuMC.all) %>% str_subset("assign.rate"),
                                            names_to = "group", values_to = "assign.rate") %>% 
    mutate(group = group %>% str_remove("assign.rate."),
           assign.rate = as.numeric(as.character(assign.rate))#,
           #ref = ref %>% str_remove("gl.america.10.") %>% str_remove(".str") %>% str_remove("gl."),
           #model = as.factor(model)
    ) %>% 
    filter(!is.na(assign.rate),
           group != "all"#,
           #str_detect(ref, "2REF")
    ) 
  
  
#  # Change levels of model
#  levels(accuMC.all$model)
#  levels(accuMC.all$model) <- c("LDA", "Naive Bayes", "Random Forest", "SVM", "Decision tree")         
  
  # Stats
  
  accuMC.all %>%# group_by(group) %>% 
    summarise(mean = mean(assign.rate),
              sd = sd(assign.rate))
  
  
  fig.MC <- accuMC.all %>% #mutate(P.train.loci = as.numeric(as.character(train.loci)) * 100) %>% 
    ggplot(aes(x = group)) +
    geom_hline(yintercept = c(0.5), lty = "dashed") +
    geom_hline(yintercept = c(0.858), col = "firebrick2", lty = "dotted", cex = 1) +
    geom_hline(yintercept = c(0.882), col = "dodgerblue1", lty = "dotted",cex = 1) +
    geom_boxplot(aes(y = assign.rate, fill = group)) +
   # facet_grid(. ~ model, space = "free_x", scale = "free_x") +
    scale_y_continuous(limit = c(0,1)) +
    labs(x = "Reference contingent", y = "Assigment accuracy") +
    scale_x_discrete(labels = c("Northern", "Southern")) +
    scale_fill_manual(values = c("firebrick2",  "dodgerblue1"))+
    #guides(fill=guide_legend(title="% of SNPs")) +
    geom_text(data = data.frame(), aes(y = 0.1, x = c("CAN", "US"), label = c("N ref = 36", "N ref = 63")), vjust = 1, col = "black", cex = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none")
  
  fig.MC
  

  
# Assignment --------------------------------------------------------------

  assign.others <- read.table(file.path("./02_Results/07_NS_tests/", "gl.NS_assign/", "AssignmentResult.txt"), header = T) 

  # graph
  
  # Rapid check - no more 5Y in 5YZe
  assign.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
    group_by(new.NAFO, NAFO) %>% summarise(N = n())
  

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
  
  library(ggpubr)
  

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
  ggsave(filename = here("02_Results/figSNS.AdultAssignments.png"), plot = fig5.v4, 
         width = 8, height = 4 , units = "in",
         dpi = 300)
 
  
  fig5.v5 <-ggpubr::ggarrange(fig.MC,
                              ggpubr::ggarrange(fig5a.v1 + theme( #plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.spacing = unit(20, "pt")), 
    fig5c.v1,
    labels = c("B", "C"),
    heights = c(1,2),
    align = "v",
    nrow = 2,
    common.legend = TRUE,
    legend = "none"),
    labels = c("A", ""), widths = c(1,3)
    )
  
  fig5.v5
  
  ggsave(filename = here("02_Results/figSNS.All.png"), plot = fig5.v5, 
         width = 10, height = 4 , units = "in", bg = "white",
         dpi = 300)
   
  