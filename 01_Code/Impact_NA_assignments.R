# Info --------------------------------------------------------------------

# Test on impact of missing values
# Mackerel ddrad dataset - only adults
# Audrey Bourret
# 2021-11-27
#

# Library -----------------------------------------------------------------


# library
library(here)
library(tidyverse)
library(readxl)
library(dartR)
library(assignPOP)

`%nin%` = Negate(`%in%`)

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

gl.america.samples

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


# Functions to subset datasets --------------------------------------------

# to compute the number of missing value
count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}

add.ind.na.gl <- function(gl, threshold){
  df.gl <- adegenet::tab(gl, NA.method = "asis") 
  df.gl.NA <- apply(df.gl, 1, function(l){
    p.na <- threshold
    n.na <-  length(l[is.na(l) == T])
    n.snps <- length(l)
    n.sampled <- round(n.snps*p.na) - n.na
    if(n.sampled > 0){
      idx.vec <- sample(which(!is.na(l)), size = n.sampled, replace = FALSE)
    } else(idx.vec <- NULL)
    replace(l, idx.vec, NA)
  } 
  )
  gl.NA <- as.genlight(t(df.gl.NA))
  pop(gl.NA) <- pop(gl)
  ploidy(gl.NA) <- 2 #
  return(gl.NA)
}

subset.gl.na <- function(gl, threshold){
  # remove individual with too much missing value
  na.info <- data.frame(Sample = indNames(gl),
                      NNA = count.ind.na.gl(gl))
  
  bad.ID <-na.info %>% filter(NNA > threshold) %>% pull(Sample)
  
  new.gl <- add.ind.na.gl(gl = gl[indNames(gl) %nin% bad.ID], threshold)
  
  return(new.gl)
}

# Only remove individual above a threshold (no NA replacement)
subset.gl.na.soft <- function(gl, threshold){
  # remove individual with too much missing value
  na.info <- data.frame(Sample = indNames(gl),
                        NNA = count.ind.na.gl(gl))
  
  bad.ID <-na.info %>% filter(NNA > threshold) %>% pull(Sample)
  
  new.gl <- gl[indNames(gl) %nin% bad.ID]
  
  return(new.gl)
}


# Create 100 datasets ------------------------------------------------------------

for(j in c(0.10, 0.20, 0.30)){

  print(j)
  
for(i in 1:100){
  
  print(i)
  
  # Subset gl
  
  gl.ref.sub <- subset.gl.na(gl.ref, threshold = j)
  gl.no.ref.sub <- subset.gl.na(gl.no.ref, threshold = j) 
  
  # Save as structure
  
  gl2structure(x = gl.ref.sub,   addcolumns =  pop(gl.ref.sub),  outfile = paste0("gl.ref.NA",j, "_", i, ".str"),   outpath = file.path("./02_Results/06_NA_tests/Datasets_100/"))
  
  gl2structure(x = gl.no.ref.sub,   addcolumns =  pop(gl.no.ref.sub),  outfile =  paste0("gl.no.ref.NA",j, "_", i, ".str"),   outpath = file.path("./02_Results/06_NA_tests/Datasets_100/"))
 
  # Load files to perform assignPOP
  
  others.str    <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Datasets_100/", paste0("gl.no.ref.NA",j, "_", i, ".str")))

  ref.str  <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Datasets_100/",  paste0("gl.ref.NA",j, "_", i, ".str")))
 
  # Perform MCMC

  assign.MC(ref.str, train.inds=c(0.7), train.loci=c(1),
            loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/06_NA_tests/MCMC_100/",  paste0("gl.ref.NA",j, "_", i),"/"),
            processors = 30)
  # unload the memory
  gc() 
  
  # Perform the assigments
  
  assign.X(x1=ref.str , x2=others.str, dir= paste0("./02_Results/06_NA_tests/Assignments_100/", paste0("gl.NA",j, "_", i),"/"), model="svm", mplot = T)

  }

}



# Load MCMC results for 100 datasets ---------------------------------------------------------

# CAN JUST READ RESULT BELOW
accuMC.test <- data.frame()

for(j in c(0.10, 0.20, 0.30)){
  
  print(j)
  
  for(i in 1:100){
    
    print(i)
  
  accuMC <- accuracy.MC(paste0("./02_Results/06_NA_tests/MCMC_100/",  paste0("gl.ref.NA",j, "_", i),"/"))
  accuMC <- accuMC %>% mutate(threshold = j, repetition = i)
  
  accuMC.test <- bind_rows(accuMC.test, accuMC)
  
}

}

#write_csv(accuMC.test, file.path(here::here(), "02_Results/06_NA_tests", "MCMC.test.csv"))

accuMC.test<- read_csv( file.path(here::here(), "02_Results/06_NA_tests", "MCMC.test.csv"))

head(accuMC.test)

accuMC.test.long <- accuMC.test %>% pivot_longer(names(accuMC.test) %>% str_subset("assign.rate"),
                                            names_to = "group", values_to = "assign.rate") %>% 
  mutate(group = group %>% str_remove("assign.rate."),
         assign.rate = as.numeric(as.character(assign.rate))
         
  ) 


accuMC.test.long %>% group_by(threshold, repetition, group) %>% 
  summarise(mean.assign.rate = mean(assign.rate)) %>% 
  filter(group != "all") %>% 
  ggplot(aes(x = group, y = mean.assign.rate, fill = group)) +
  geom_hline(yintercept = c(0.5), lty = "dashed") +
  geom_boxplot() + 
  geom_jitter(height = 0, alpha = 0.2) +
  labs(x = "", y = "Mean mssigment accuracy", title = "Mean accuracy of 100 datasets, at 3 NA thresholds") +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  facet_grid(~ factor(threshold), space = "free", scale = "free") + 
  theme_bw()

# Load assignments results ------------------------------------------------


# CAN JUST READ RESULT BELOW
assign.test.others <- data.frame()

for(j in c(0.10, 0.20, 0.30)){
  
  print(j)
  
  for(i in 1:100){
    
    print(i)
    
    assign.other <- read.table(paste0("./02_Results/06_NA_tests/Assignments_100/", paste0("gl.NA",j, "_", i),"/", "AssignmentResult.txt"), header = T)
    assign.other <- assign.other %>% mutate(threshold = j, repetition = i)
    
    assign.test.others <- bind_rows(assign.test.others, assign.other)
    
  }
  
}


#write_csv(assign.test.others, file.path(here::here(), "02_Results/06_NA_tests", "assignment.test.csv"))
assign.test.others <- read_csv(file.path(here::here(), "02_Results/06_NA_tests", "assignment.test.csv"))

assign.test.others 

assign.others <- read.table(file.path("./02_Results/03_assignPOP/Assignations_SVM/", "AssignmentResult.txt"), header = T)  %>% mutate(Test = "ori")

df.stat.na <- data.frame(Sample = indNames(gl.america.samples),
                         NNA.before = count.ind.na.gl(gl.america.samples)
)


assign.test.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
 # filter(repetition %in% 1:10) %>% 
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= CAN.ori, y = CAN, group = factor(repetition) ))+
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "gam", se = F, col = "gray", lwd = 0.1) + 
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Membership probability full dataset", y =  "Membership probability subset dataset", title = "100 datasets, at 3 NA thresholds") +
  facet_grid(.~threshold) + 
  theme_bw()


assign.test.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
   #filter(repetition %in% 1:10) %>% 
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA.before, y = CAN, group = factor(repetition), col = Country))+
  #geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.01)+
  geom_smooth(data = assign.others %>%  left_join(df.stat.na, by = c("Ind.ID" = "Sample")), 
              aes(group = 1),
              method = "gam", se = F, col = "black", lwd = 0.5) + 
  geom_point(alpha = 0.01)+
  #scale_x_continuous(limits = c(0,0.3)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "Membership probability subset dataset", x =  "% NA full dataset", title = "100 datasets, at 3 NA thresholds") +
  facet_grid(.~threshold) + 
  theme_bw()



# Real results ------------------------------------------------------------

assign.others <- read.table(file.path("./02_Results/03_assignPOP/Assignations_SVM/", "AssignmentResult.txt"), header = T)  %>% mutate(Test = "ori")

accuMC <- accuracy.MC(paste0("./02_Results/03_assignPOP/gl.ref_june2021.svm/"))

accuMC %>% group_by(train.loci) %>% summarise(N = n(),
                                              med.CAN = median(assign.rate.CAN),
                                              med.US = median(assign.rate.US),
                                              med.Tot = median(assign.rate.all))


assign.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA.before, y = CAN, col = Country))+
  #geom_abline(slope = 1, intercept = 0) +
  geom_point()+
  geom_smooth(data = assign.others %>%  left_join(df.stat.na, by = c("Ind.ID" = "Sample")), 
              aes(group = 1),
              method = "gam", se = F, col = "black", lwd = 0.5) + 
  geom_point(alpha = 0.01)+
  #scale_x_continuous(limits = c(0,0.3)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "Membership probability", x =  "% NA", title = "Full datasets") +
 # facet_grid(.~threshold) + 
  theme_bw()


# Truncation --------------------------------------------------------------

gl.ref.NA30 <- subset.gl.na.soft(gl.ref, threshold = 0.3)
gl.ref.NA20 <- subset.gl.na.soft(gl.ref, threshold = 0.2)
gl.ref.NA15 <- subset.gl.na.soft(gl.ref, threshold = 0.15)
gl.ref.NA10 <- subset.gl.na.soft(gl.ref, threshold = 0.1)

gl.no.ref.NA30 <- subset.gl.na.soft(gl.no.ref, threshold = 0.3)
gl.no.ref.NA20 <- subset.gl.na.soft(gl.no.ref, threshold = 0.2)
gl.no.ref.NA15 <- subset.gl.na.soft(gl.no.ref, threshold = 0.15)
gl.no.ref.NA10 <- subset.gl.na.soft(gl.no.ref, threshold = 0.1)

# Do structure files

gl2structure(x = gl.ref.NA30,   addcolumns =  pop(gl.ref.NA30),  outfile = "gl.ref.NA30.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.ref.NA20,   addcolumns =  pop(gl.ref.NA20),  outfile = "gl.ref.NA20.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.ref.NA15,   addcolumns =  pop(gl.ref.NA15),  outfile = "gl.ref.NA15.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.ref.NA10,   addcolumns =  pop(gl.ref.NA10),  outfile = "gl.ref.NA10.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))

gl2structure(x = gl.no.ref.NA30,   addcolumns =  pop(gl.no.ref.NA30),  outfile = "gl.no.ref.NA30.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.no.ref.NA20,   addcolumns =  pop(gl.no.ref.NA20),  outfile = "gl.no.ref.NA20.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.no.ref.NA15,   addcolumns =  pop(gl.no.ref.NA15),  outfile = "gl.no.ref.NA15.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))
gl2structure(x = gl.no.ref.NA10,   addcolumns =  pop(gl.no.ref.NA10),  outfile = "gl.no.ref.NA10.str",   outpath = file.path("./02_Results/06_NA_tests/Subset_threshold_only/"))

# Load files to perform assignPOP

others.NA30.str    <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.no.ref.NA30.str"))
others.NA20.str    <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.no.ref.NA20.str"))
others.NA15.str    <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.no.ref.NA15.str"))
others.NA10.str    <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.no.ref.NA10.str"))


ref.NA30.str       <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.ref.NA30.str"))
ref.NA20.str       <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.ref.NA20.str"))
ref.NA15.str       <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.ref.NA15.str"))
ref.NA10.str       <- assignPOP::read.Structure(file.path("./02_Results/06_NA_tests/Subset_threshold_only/", "gl.ref.NA10.str"))

# Compute MCMC

assign.MC(ref.NA30.str, train.inds=c(0.7), train.loci=c(1),
          loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_30/"),
          processors = 30)

assign.MC(ref.NA20.str, train.inds=c(0.7), train.loci=c(1),
          loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_20/"),
          processors = 30)

assign.MC(ref.NA15.str, train.inds=c(0.7), train.loci=c(1),
          loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_15/"),
          processors = 30)

assign.MC(ref.NA10.str, train.inds=c(0.7), train.loci=c(1),
          loci.sample="fst", iterations=30, model="svm", dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_10/"),
          processors = 30)

# Perform the assignments

assign.X(x1=ref.NA30.str , x2=others.NA30.str, dir="./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA30/", model="svm", mplot = T)
assign.X(x1=ref.NA20.str , x2=others.NA20.str, dir="./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA20/", model="svm", mplot = T)
assign.X(x1=ref.NA15.str , x2=others.NA15.str, dir="./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA15/", model="svm", mplot = T)
assign.X(x1=ref.NA10.str , x2=others.NA10.str, dir="./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA10/", model="svm", mplot = T)

# Load MCMC results - truncations tests -------------------------------------------------------

    accuMC.subset <- bind_rows(accuracy.MC(dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_30/")) %>%  mutate(threshold = 0.3, repetition = "ori"),
                               accuracy.MC(dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_20/")) %>%  mutate(threshold = 0.2, repetition = "ori"),
                               accuracy.MC(dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_15/")) %>%  mutate(threshold = 0.15, repetition = "ori"),
                               accuracy.MC(dir=paste0("./02_Results/06_NA_tests/Subset_threshold_only/MCMC_10/")) %>%  mutate(threshold = 0.1, repetition = "ori")
                              )
    
#write_csv(accuMC.subset, file.path(here::here(), "02_Results/06_NA_tests", "MCMC.subset.csv"))

accuMC.subset <- read_csv(file.path(here::here(), "02_Results/06_NA_tests", "MCMC.subset.csv"))

accuMC.subset.long <- accuMC.subset %>% pivot_longer(names(accuMC.subset) %>% str_subset("assign.rate"),
                                                 names_to = "group", values_to = "assign.rate") %>% 
  mutate(group = group %>% str_remove("assign.rate."),
         assign.rate = as.numeric(as.character(assign.rate))
         
  ) 


accuMC.subset.long %>%  filter(group != "all") %>% 
  ggplot(aes(x = group, y = assign.rate, fill = group)) +
  geom_hline(yintercept = c(0.5), lty = "dashed") +
  geom_boxplot() + 
  #geom_jitter(height = 0, alpha = 0.2) +
  labs(x = "", y = "Assignment accuracy", title = "Accuracy of subset datasets, at 3 NA thresholds") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  facet_grid(~ factor(threshold), space = "free", scale = "free") + 
  theme_bw()

# Load assignment results -------------------------------------------------

assign.subset.others <- bind_rows(read.table(file.path("./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA30/", "AssignmentResult.txt"), header = T) %>% mutate(threshold = 0.3, repetition = "ori"), 
                                  read.table(file.path("./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA20/", "AssignmentResult.txt"), header = T) %>% mutate(threshold = 0.2, repetition = "ori"),
                                  read.table(file.path("./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA15/", "AssignmentResult.txt"), header = T) %>% mutate(threshold = 0.15, repetition = "ori"),
                                  read.table(file.path("./02_Results/06_NA_tests/Subset_threshold_only/Assignations_NA10/", "AssignmentResult.txt"), header = T) %>% mutate(threshold = 0.1, repetition = "ori") 

                                                            ) 
# Stats



#write_csv(assign.subset.others, file.path(here::here(), "02_Results/06_NA_tests", "assignment.subset.csv"))

assign.subset.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
                         left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA.before, y = CAN, col = Country) )+ geom_point(alpha = 0.5) + 
  facet_wrap(~factor(threshold), nrow = 1)


assign.subset.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  #left_join(df.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= CAN.ori, y = CAN, col = Country) )+ geom_point() + 
  facet_wrap(~factor(threshold), nrow = 1)


# Figures -----------------------------------------------------------------

subset.ref.stat <- data.frame(threshold = rep(c(0.1, 0.15, 0.2, 0.3), 2),
                             group = c(rep("CAN", 4), rep("US", 4)),
                             N = c(47, 52, 53, 53, 38, 51, 56, 62)
                             ) %>% 
  mutate(test = paste(round(threshold * 100), "%"))

gl.no.ref.NA15

ID <- assign.subset.others %>% filter(threshold == 0.20) %>% pull(Ind.ID)

cor((assign.subset.others %>% filter(threshold == 0.20) %>% pull(CAN)), 
    (assign.others %>% filter(Ind.ID %in% ID) %>% pull(CAN)), method = "spearman")

subset.no.ref.stat <- data.frame(threshold = c(0.1,0.15, 0.2, 0.3),
                                 N = c(363, 407, 434, 441),
                                 P.NA = c(3.7, 4.6, 5.4, 5.6),
                                 r.NA = c("rho = 0.96", "rho > 0.99", "rho > 0.99", "rho > 0.99"),
                                 r.test.NA = c("rho = 0.75","r = ?",  "rho = 0.86", "rho = 0.88")
) %>% 
  mutate(test = paste(round(threshold * 100), "%"))
                          
assign.test.others %>%  group_by(Ind.ID, threshold) %>% 
  summarise(mean.CAN = mean(CAN),
            sd.CAN = sd(CAN)) %>% 
  left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>% 
  filter(threshold == 0.30) %>% ungroup() %>% select(mean.CAN, CAN.ori) %>% cor(method = "spearman")


subset.ref.stat 


cor.test(assign.others %>%  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% pull(CAN),
    assign.others %>%  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% pull(NNA.before), method = "spearman")

g0 <- assign.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  ggplot(aes(x= NNA.before * 100 , y = CAN, col = Country))+
  #geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.5)+
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"))+
  #geom_smooth(data = assign.others %>%  left_join(df.stat.na, by = c("Ind.ID" = "Sample")), 
  #            aes(group = 1),
  #            method = "gam", se = F, col = "black", lwd = 0.5) + 
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_text(aes(y = 1, x = 30, label = paste0("N individuals = ", 441)), vjust = 1, hjust = 1,  col = "black", cex = 3) +
  geom_text(aes(y = 0.9, x = 30, label = paste0("Overall NA = ", 5.6, "%")), vjust = 1, hjust = 1, col = "black", cex = 3) +
  labs(y = "Membership probability", x =  "% of missing values") +
  # facet_grid(.~threshold) + 
  theme_bw() +   theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")
g0

gl.no.ref


g1.a <- accuMC.subset.long %>%  filter(group != "all") %>%
  filter(threshold %in% c(0.1,0.15,0.2)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x = group, y = assign.rate, fill = group)) +
  geom_hline(yintercept = c(0.858), col = "firebrick2", lty = "dotted", cex = 1) +
  geom_hline(yintercept = c(0.882), col = "dodgerblue1", lty = "dotted",cex = 1) +
  geom_hline(yintercept = c(0.5),  lty = "dashed") +  
  geom_boxplot(alpha = .75) + 
  #geom_jitter(height = 0, alpha = 0.2) +
  labs(x = "", y = "Assignment accuracy") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  facet_grid(~ test, space = "free", scale = "free") + 
  scale_fill_manual(values = c("firebrick2", "dodgerblue1"))+
  geom_text(data = subset.ref.stat %>%   filter(threshold %in% c(0.1,0.15,0.2)), aes(y = 0.1, label = paste0("N ref = ",N)), vjust = 1, col = "black", cex = 3) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")
g1.a

g1.b <- assign.subset.others %>%  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  filter(threshold %in% c(0.1,0.15,0.2)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x= NNA.before * 100, y = CAN, col = Country) )+ 
  geom_point(alpha = 0.5) + 
  facet_grid(~ test) + 
  labs(x = "% of missing values", y = "Membership probability") + 
  geom_text(data = subset.no.ref.stat  %>% filter(threshold %in% c(0.1,0.15,0.2)), aes(y = 1, x = 30, label = paste0("N individuals = ",N)), vjust = 1, hjust = 1,  col = "black", cex = 3) +
  geom_text(data = subset.no.ref.stat  %>%   filter(threshold %in% c(0.1,0.15,0.2)), aes(y = 0.9, x = 30, label = paste0("Overall NA = ",P.NA, "%")), vjust = 1, hjust = 1, col = "black", cex = 3) +
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"))+
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")

g1.b

g1.c <- assign.subset.others %>% left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  filter(threshold %in% c(0.1,0.15,0.2)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x= CAN.ori, y = CAN, col = Country) )+ 
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_point(alpha = 0.5) + 
  facet_grid(~ test, space = "free", scale = "free") + 
  labs(x = "Full membership probability", y = "Filtered membership probability") + 
  geom_text(data = subset.no.ref.stat %>%  filter(threshold %in% c(0.1,0.15,0.2)), aes(y = 1, x = 0, label = paste0(r.NA), vjust = 1, hjust = 0), col = "black" , cex = 3)+
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"))+
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")
g1.c


g1 <- ggpubr::ggarrange(g1.a, g1.b, g1.c,
                  labels = LETTERS,
                  nrow = 3)

g1


g2.a <- accuMC.test.long %>%  group_by(threshold, repetition, group) %>% 
  summarise(mean.assign.rate = mean(assign.rate)) %>% 
  filter(group != "all") %>% 
  filter(threshold %in% c(0.1,0.2, 0.3)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x = group, y = mean.assign.rate, fill = group)) +
  geom_hline(yintercept = c(0.858), col = "firebrick2", lty = "dotted", cex = 1) +
  geom_hline(yintercept = c(0.882), col = "dodgerblue1", lty = "dotted",cex = 1) +
  geom_hline(yintercept = c(0.5),  lty = "dashed") +  
  geom_boxplot(alpha = .75) + 
  labs(x = "", y = "Mean assignment accuracy") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c("Northern", "Southern")) +
  facet_grid(~ test, space = "free", scale = "free") + 
  scale_fill_manual(values = c("firebrick2", "dodgerblue1"))+
  geom_text(data = subset.ref.stat %>%  filter(threshold %in% c(0.1,0.2, 0.3)) , aes(y = 0.1, label = paste0("N ref = ",N)), vjust = 1, col = "black", cex = 3) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")
g2.a


g2.b <- assign.test.others %>%  
  group_by(Ind.ID, threshold) %>% summarise(mean.CAN = mean(CAN),
                                            sd.CAN = sd(CAN)) %>% 
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  left_join(df.stat.na, by = c("Ind.ID" = "Sample")) %>% 
  filter(threshold %in% c(0.1,0.2, 0.3)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x= NNA.before * 100, y = mean.CAN, col = Country) )+ 
 # geom_point(alpha = 0.5) + 
  geom_pointrange(aes(ymin = mean.CAN - sd.CAN, ymax = mean.CAN + sd.CAN), alpha = 0.5, cex = 0.2) +
  facet_grid(~ test) + 
  labs(x = "% of missing values in full dataset", y = "Mean membership probability") + 
  geom_text(data = subset.no.ref.stat %>%    filter(threshold %in% c(0.1,0.2, 0.3)), aes(y = 1, x = 30, label = paste0("N individuals = ",N)), vjust = 1, hjust = 1,  col = "black", cex = 3) +
  geom_text(data = subset.no.ref.stat %>%    filter(threshold %in% c(0.1,0.2, 0.3)), aes(y = 0.9, x = 30, label = paste0("Overall NA = ", test)), vjust = 1, hjust = 1, col = "black", cex = 3) +
  scale_x_continuous(limits = c(0,30)) +
  #scale_y_continuous(limits = c(0,1)) +
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"))+
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")

g2.b


g2.c <- assign.test.others %>%  group_by(Ind.ID, threshold) %>% 
                               summarise(mean.CAN = mean(CAN),
                                           sd.CAN = sd(CAN)) %>% 
    left_join(assign.others %>% select(Ind.ID, CAN.ori = CAN)) %>%
  left_join(pop.data, by = c("Ind.ID" = "Sample")) %>% 
  filter(threshold %in% c(0.1,0.2, 0.3)) %>% 
  mutate(test = paste(round(threshold * 100), "%")) %>% 
  ggplot(aes(x= CAN.ori, y = mean.CAN, col = Country) )+
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_point(alpha = 0.5) + 
  facet_grid(~ test, space = "free", scale = "free") + 
  labs(x = "Full membership probability", y = "Mean filtered membership probability") + 
  geom_text(data = subset.no.ref.stat %>%  filter(threshold %in% c(0.1,0.2, 0.3)), aes(y = 1, x = 0, label = paste0(r.test.NA), vjust = 1, hjust = 0), col = "black" , cex = 3)+
  
  scale_colour_manual(values = c("firebrick2", "dodgerblue1"))+
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none")
g2.c



g1 <- ggpubr::ggarrange(g1.a, g1.b, g1.c,
                  labels = LETTERS,
                  nrow = 3)


g2 <- ggpubr::ggarrange(g2.a, g2.b, g2.c,
                  labels = LETTERS,
                  nrow = 3)

g1
g2

ggpubr::ggarrange(g1, g2)

ggsave(filename = here("02_Results/figSNA_original.png"), plot = g0, 
       width = 3, height = 3 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figSNA_Subset.png"), plot = g1, 
       width = 6, height = 8 , units = "in",
       dpi = 300)

ggsave(filename = here("02_Results/figSNA_Simulations.png"), plot = g2, 
       width = 6, height = 8 , units = "in",
       dpi = 300)

