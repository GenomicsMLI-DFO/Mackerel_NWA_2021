
# Info --------------------------------------------------------------------

# Small script for data conversion
# VCF + metadata limited to the right individuals
# Audrey Bourret
# 2022-03-14
#

# Library -----------------------------------------------------------------

library(here)

library(tidyverse)
library(readxl)

library(adegenet)
library(hierfstat)

`%nin%` = Negate(`%in%`)

# Original MetaData --------------------------------------------------------------------

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


# Original panel ----------------------------------------------------------

load(here("00_Data/gen.america.878samples.10703snps.data")) 
load(here("00_Data/gen.global.878samples.10832snps.data")) 

# Subset new data sets ----------------------------------------------------

# ID Gobal data
ID.samples.glob <- pop.data %>% filter(Cat_Sample == "Sample",
                                       Dev_Stage == "Adult",
                                       SNP_panel == "Yes") %>% 
  pull(ID_GQ)

gl.global.samples <- gl.global.10[indNames(gl.global.10) %in% ID.samples.glob]

Metadata <- data.frame(ID = indNames(gl.global.samples)) %>% 
  left_join(pop.data %>% select(ID = ID_GQ, Country, Lat, Long, Year, Month, Day, NAFO = new.NAFO, Length, Age, Sex))

Metadata %>% View()

write_csv(Metadata, file = "00_Data/Metadata_Mackerel629_NWA2021.csv"
          )


# American data
ID.samples.ame <- pop.data %>% filter(Cat_Sample == "Sample",
                                      Dev_Stage == "Adult",
                                      SNP_panel == "Yes",
                                      Country %in% c("Canada", "US")) %>% 
  pull(ID_GQ)

gl.america.samples <- gl.america.10[indNames(gl.america.10) %in% ID.samples.ame]


SNP.ame <- locNames(gl.america.samples)
SNP.glo <- locNames(gl.global.samples)


vcf.global.path = "/media/genobiwan/Extra_Storage/Projets/Maquereaux_RAD/00_Data/06a_Filtering.denovo/global.r0.9.MAF0.1/06e_AllIndividuals/populations.snps.vcf"
vcf.america.path = "/media/genobiwan/Extra_Storage/Projets/Maquereaux_RAD/00_Data/06a_Filtering.denovo/america.r0.9.MAF0.1/06e_AllIndividuals/populations.snps.vcf"

vcf.test <- vcfR::read.vcfR("00_Data/gen.america.557samples.10703snps.recode.vcf")

cmd.ame <- paste("--vcf", vcf.america.path, 
              "--recode",
              paste("--indv", ID.samples.ame, collapse = " "),
              "--out", "00_Data/gen.america.557samples.10703snps")
cmd.ame

A1 <- system2("vcftools", cmd.ame, stdout=T, stderr=T)
A1

cmd.global <- paste("--vcf", vcf.global.path, 
                 "--recode",
                 paste("--indv", ID.samples.glob , collapse = " "),
                 "--out", "00_Data/gen.global.629samples.10832snps")


cmd.global

A2 <- system2("vcftools", cmd.global, stdout=T, stderr=T)
A2

# Check that the conversion were equivalent

vcf.ame <- vcfR::read.vcfR("00_Data/gen.america.557samples.10703snps.recode.vcf")
gl.ame.test <- vcfR::vcfR2genlight(vcf.ame)

table(locNames(gl.ame.test) == locNames(gl.america.samples))
table(indNames(gl.ame.test) == indNames(gl.america.samples))

vcf.glob <- vcfR::read.vcfR("00_Data/gen.global.629samples.10832snps.recode.vcf")
gl.glob.test <- vcfR::vcfR2genlight(vcf.glob)

table(locNames(gl.glob.test) == locNames(gl.global.samples))
table(indNames(gl.glob.test) == indNames(gl.global.samples))


