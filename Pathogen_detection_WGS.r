###################################################
################## Pathogen-detection pipeline in WGS DLBCL samples
##############

srun --mem-per-cpu=32GB --time=6:00:00 --pty --cpus-per-task=8 bash
module load hpc
module load anaconda3
source activate agshin_condaenv
module load r/4.1.1
R 

library("ggpubr")
library("ComplexHeatmap")
library("reshape2")
library("dplyr", warn.conflicts = FALSE)
require("gridExtra")
library("viridis")
library("ggtext")
library("circlize")
library("readxl")
library("caret")

library(future)
plan("multicore", workers=32)
options(future.globals.maxSize= 32000 * 1024^2)


####### Exploratory data analysis
### Load dataset
parent.directory <- "/home/cluster/abalay/scratch/moved_from_home/WGS_RAW_DLBCL/Kraken2_output"

kraken.file <- read.delim(file.path(parent.directory, "WGS.DLBCL.combined.mpa.txt"), sep="\t", row.names=1)

### change column names
colnames(kraken.file) <- gsub("_dbGaP.28434", "", unlist(strsplit(colnames(kraken.file), ".breport")))

kraken.file <- t(kraken.file)
dim(kraken.file)
# [1]  183 4099

##### create figures dir
dir.figures <- paste(parent.directory, "Figures", sep="/")
if (!dir.exists(dir.figures)){
    dir.create(dir.figures, recursive = TRUE)
}


### plot density distribution of total number of aligned reads
kingdoms <- c("Eukaryota", "Bacteria", "Viruses", "Archaea")
temp.kraken <- kraken.file[, paste("k__", kingdoms, sep="")]

mini <- data.frame(row.names=rownames(temp.kraken), Total_reads=rowSums(temp.kraken))

p <- ggplot(mini, aes(x=log10(Total_reads))) + geom_density() + xlab("Total aligned reads (log10)") + ylab("Density of reads") +
    ggtitle("Total aligned reads distribution across organisms") + 
    geom_vline(data=mini, aes(xintercept = log10(Total_reads), color=rownames(mini)), linetype = "dashed")

png(file.path(dir.figures, "Total_aligned_taxonomic_reads.WGS.DLBCL.png"), res=200, unit="in", height=8, width=11)
p + guides(color=guide_legend(title="Sample", size=5))
dev.off()


### plot stacked plot of viridiplantae, fungi, bacteria, virus, archaea
kingdoms <- c("Eukaryota", "Bacteria", "Viruses", "Archaea")
subkingdoms <- c("Viridiplantae", "Fungi")
temp.kraken <- cbind(kraken.file[, paste("k__", kingdoms, sep="")], kraken.file[, paste("k__Eukaryota|k__", subkingdoms, sep="")])
temp.kraken <- temp.kraken[,2:ncol(temp.kraken)]/sum(temp.kraken[, paste("k__", kingdoms, sep="")])
colnames(temp.kraken) <- c("Bacteria", "Viruses", "Archaea", "Viridiplantae", "Fungi")

temp.kraken <- melt(temp.kraken); colnames(temp.kraken) <- c("Sample", "Kingdom", "Fraction_reads")

p <- ggplot(temp.kraken, aes(x=Sample, y=Fraction_reads, fill=Kingdom)) + geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.title.x=element_text(size=6), axis.text.x = element_text(angle = 90)) +
    scale_color_manual(values=c("red", "blue", "green", "orange", "purple")) +
    xlab("Sample") + ylab("Fraction of reads") + ggtitle("Per sample fraction of reads across organisms") 

png(file.path(dir.figures, "Fraction_of_taxonomic_reads.WGS.DLBCL.png"), res=200, unit="in", height=10, width=15)
p
dev.off()


bacteria.list <- list()
##### bacteria
### subset any fields with species detected
temp.kraken <- subset(kraken.file, select = grepl("k__Bacteria", colnames(kraken.file)))
dim(temp.kraken)
# [1]  183 3289

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1]  183 2094

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

bacteria.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
#     Xanthomonas_euvesicatoria         Staphylococcus_aureus 
#                  4.716681e-05                  2.779526e-05 
#          Pseudomonas_tolaasii    Mycobacterium_tuberculosis 
#                  1.037315e-05                  2.974635e-06 
#        Legionella_pneumophila        Legionella_israelensis 
#                  1.481139e-06                  1.390893e-06 
# Limosilactobacillus_fermentum           Cutibacterium_acnes 
#                  9.931678e-07                  8.965001e-07 
#     Flavobacterium_johnsoniae         Legionella_antarctica 
#                  7.429383e-07                  6.983349e-07

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
#      Xanthomonas_euvesicatoria     Mycobacterium_tuberculosis 
#                   5.910731e-04                   1.279230e-04 
#            Streptococcus_mitis           Pseudomonas_tolaasii 
#                   7.068051e-05                   4.985022e-05 
#          Staphylococcus_aureus       Streptococcus_pneumoniae 
#                   4.790851e-05                   3.131839e-05 
#            Veillonella_atypica            Ralstonia_pickettii 
#                   2.874041e-05                   1.789553e-05 
#  Limosilactobacillus_fermentum         Ralstonia_solanacearum 
#                   1.495122e-05                   1.469552e-05 
#           Streptococcus_oralis          Klebsiella_pneumoniae 
#                   1.371660e-05                   1.320426e-05 
#  Keratinibaculum_paraultunense            Ralstonia_insidiosa 
#                   1.039804e-05                   1.023015e-05 
#            Veillonella_parvula         Xanthomonas_campestris 
#                   1.009485e-05                   7.832734e-06 
#            Cutibacterium_acnes          Staphylococcus_cohnii 
#                   6.355217e-06                   6.345746e-06 
# Streptococcus_pseudopneumoniae      Ralstonia_mannitolilytica 
#                   6.304103e-06                   4.995011e-06 


virus.list <- list()
##### viruses
temp.kraken <- subset(kraken.file, select = grepl("k__Viruses", colnames(kraken.file)))
dim(temp.kraken)
# [1] 183 148

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 183  66

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

virus.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
#   Human_gammaherpesvirus_4      Proteus_virus_Isfahan 
#               6.547868e-06               3.416780e-06 
#    Human_betaherpesvirus_5     Human_mastadenovirus_C 
#               1.367662e-06               7.110777e-08 
# Enterococcus_phage_nattely Papiine_gammaherpesvirus_1 
#               2.275536e-08               2.258685e-08 
#       Escherichia_virus_T4   Torque_teno_mini_virus_5 
#               2.210094e-08               1.863426e-08 
#  Torque_teno_mini_virus_16    Salmonella_phage_SAP012 
#               1.792477e-08               1.589126e-08 

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
#   Human_gammaherpesvirus_4    Human_betaherpesvirus_5 
#               9.072077e-04               2.476609e-04 
#      Proteus_virus_Isfahan     Human_mastadenovirus_C 
#               7.022618e-05               1.055709e-05 
# Papiine_gammaherpesvirus_1  Torque_teno_mini_virus_16 
#               4.035253e-06               3.235720e-06 
#   Torque_teno_mini_virus_5        TTV-like_mini_virus 
#               2.599510e-06               2.471937e-06 
#   Torque_teno_mini_virus_3 Enterococcus_phage_nattely 
#               2.370873e-06               2.124628e-06


fungus.list <- list()
##### fungi
temp.kraken <- subset(kraken.file, select = grepl("k__Eukaryota|k__Fungi", colnames(kraken.file), fixed = TRUE))
dim(temp.kraken)
# [1] 183 160

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 183  66

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

fungus.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
#            Botrytis_cinerea       Naumovozyma_castellii 
#                1.574936e-05                4.684728e-06 
#        Candida_dubliniensis     Kluyveromyces_marxianus 
#                2.908921e-06                1.927210e-06 
#         Pichia_kudriavzevii       Ustilaginoidea_virens 
#                1.898350e-06                1.167466e-06 
#      Aspergillus_luchuensis        Fusarium_graminearum 
#                1.037956e-06                9.579086e-07 
# Colletotrichum_higginsianum          Fusarium_oxysporum 
#                8.539683e-07                6.037838e-07 

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
    #        Botrytis_cinerea       Naumovozyma_castellii 
    #            2.973845e-05                7.097557e-06 
    #    Candida_dubliniensis     Kluyveromyces_marxianus 
    #            4.654923e-06                2.987488e-06 
    #   Ustilaginoidea_virens         Pichia_kudriavzevii 
    #            2.806638e-06                2.745954e-06 
    #  Aspergillus_luchuensis        Fusarium_graminearum 
    #            2.183013e-06                1.880186e-06 
    #        Candida_albicans Colletotrichum_higginsianum 
    #            1.519805e-06                1.504336e-06 


protozoa.list <- list()
##### protozoa
temp.kraken <- subset(kraken.file, select = grepl("k__Eukaryota|p__", colnames(kraken.file), fixed = TRUE))
dim(temp.kraken)
# [1] 183  84

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 183  37

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

protozoa.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Cyanidioschyzon_merolae    Giardia_intestinalis  Cryptosporidium_parvum 
#            3.042289e-07            0.000000e+00            0.000000e+00 
#       Toxoplasma_gondii      Besnoitia_besnoiti        Neospora_caninum 
#            0.000000e+00            0.000000e+00            0.000000e+00 
#      Theileria_annulata         Theileria_parva    Theileria_orientalis 
#            0.000000e+00            0.000000e+00            0.000000e+00 
#          Theileria_equi 
#            0.000000e+00 

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
# Cyanidioschyzon_merolae    Giardia_intestinalis  Cryptosporidium_parvum 
#            8.684557e-07            0.000000e+00            0.000000e+00 
#       Toxoplasma_gondii      Besnoitia_besnoiti        Neospora_caninum 
#            0.000000e+00            0.000000e+00            0.000000e+00 
#      Theileria_annulata         Theileria_parva    Theileria_orientalis 
#            0.000000e+00            0.000000e+00            0.000000e+00 
#          Theileria_equi 
#            0.000000e+00

############################################### Step 1: get fractions of major infectious bacteria, viruses and fungi
##########################################
#################################

##### create figures dir
dir.figures <- "/home/cluster/abalay/scratch/moved_from_home/Kraken2_Figures"
if (!dir.exists(dir.figures)){
    dir.create(dir.figures, recursive = TRUE)
}

##################### Density of proportion of common infectious bacteria species
#### List of selected infectious bacteria species
list.inf.bacteria <- c("Staphylococcus_aureus", "Mycobacterium_tuberculosis", "Legionella_pneumophila", "Legionella_israelensis", "Streptococcus_mitis", 
    "Streptococcus_pneumoniae", "Veillonella_atypica", "Ralstonia_pickettii", "Limosilactobacillus_fermentum", 
    "Streptococcus_oralis", "Klebsiella_pneumoniae", "Ralstonia_insidiosa", 
    "Veillonella_parvula")

#### Combine all samples in one dataframe
names <- names(bacteria.list)
bacteria.df <- do.call(dplyr::bind_rows, bacteria.list)
bacteria.df[is.na(bacteria.df)] <- 0

#### Check if all bacteria in a dataframe
all(list.inf.bacteria %in% colnames(bacteria.df))
# [1] TRUE

write.table(bacteria.df, file=file.path(dir.figures, "Common.bacteria.kraken.dataset.txt"), sep="\t")

#### Density plot for each common infectious bacteria species
ggobj <- list()
for (i in list.inf.bacteria){
    j <- which(list.inf.bacteria == i)
    p <- ggplot(bacteria.df, aes_string(x=i)) + geom_density() + 
    xlab("Fraction of reads") + ylab("Density") + theme(plot.title = element_text(face = "bold")) + 
    ggtitle(i)
    ggobj[[j]] <- p
    }

png(file.path(dir.figures, paste("Distribution_of_aligned_fraction_of_reads_across_bacteria_species.png", sep="")), res=200, unit="in", height=10, width=15)
do.call(grid.arrange, c(ggobj, list(ncol = 4, nrow = 6)))
dev.off()


##########################################
#################################
##################### Density of proportion of common infectious virus species
#### List of selected infectious virus species
list.inf.virus <- c("Human_gammaherpesvirus_4", "Human_betaherpesvirus_5", "Human_mastadenovirus_C")

#### Combine all samples in one dataframe
names <- names(virus.list)
virus.df <- do.call(dplyr::bind_rows, virus.list)
virus.df[is.na(virus.df)] <- 0

#### Check if all bacteria in a dataframe
all(list.inf.virus %in% colnames(virus.df))
# [1] TRUE

write.table(virus.df, file=file.path(dir.figures, "Common.virus.kraken.dataset.txt"), sep="\t")

#### Density plot for each common infectious bacteria species
ggobj <- list()
for (i in list.inf.virus){
    j <- which(list.inf.virus == i)
    p <- ggplot(virus.df, aes_string(x=i)) + geom_density() + 
    xlab("Fraction of reads") + ylab("Density") + theme(plot.title = element_text(face = "bold")) + 
    ggtitle(i)
    ggobj[[j]] <- p
    }

png(file.path(dir.figures, paste("Distribution_of_aligned_fraction_of_reads_across_virus_species.png", sep="")), res=200, unit="in", height=10, width=15)
do.call(grid.arrange, c(ggobj, list(ncol = 3, nrow = 2)))
dev.off()


##########################################
#################################
##################### Density of proportion of common infectious fungi species
#### List of selected infectious fungi species
list.inf.fungi <- c("Candida_dubliniensis", "Kluyveromyces_marxianus", "Pichia_kudriavzevii")

#### Combine all samples in one dataframe
names <- names(fungus.list)
fungi.df <- do.call(dplyr::bind_rows, fungus.list)
fungi.df[is.na(fungi.df)] <- 0

#### Check if all bacteria in a dataframe
all(list.inf.fungi %in% colnames(fungi.df))
# [1] TRUE

write.table(fungi.df, file=file.path(dir.figures, "Common.fungi.kraken.dataset.txt"), sep="\t")

#### Density plot for each common infectious bacteria species
ggobj <- list()
for (i in list.inf.fungi){
    j <- which(list.inf.fungi == i)
    p <- ggplot(fungi.df, aes_string(x=i)) + geom_density() + 
    xlab("Fraction of reads") + ylab("Density") + theme(plot.title = element_text(face = "bold")) + 
    ggtitle(i)
    ggobj[[j]] <- p
    }

png(file.path(dir.figures, paste("Distribution_of_aligned_fraction_of_reads_across_fungi_species.png", sep="")), res=200, unit="in", height=10, width=15)
do.call(grid.arrange, c(ggobj, list(ncol = 3, nrow = 2)))
dev.off()




############################################### Step 2: estimate load of common infectious bacteria, virus and fungi
###################################
##########################
################ Download the genomes

##### create dir for genomes
dir.genomes <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Files_for_scripts/Pathogen_genomes"
if (!dir.exists(dir.genomes)){
    dir.create(dir.genomes, recursive = TRUE)
}

##### create dir for cds sequences of genomes
dir.cds.genomes <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Files_for_scripts/Pathogen_genome_cds"
if (!dir.exists(dir.cds.genomes)){
    dir.create(dir.cds.genomes, recursive = TRUE)
}

##### create dir for annotation files of genomes
dir.anno <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Files_for_scripts/Pathogen_genome.Anno"
if (!dir.exists(dir.anno)){
    dir.create(dir.anno, recursive = TRUE)
}

##### download genomes from ncbi website
full.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/753/085/GCF_001753085.1_ASM175308v1/GCF_001753085.1_ASM175308v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/361/795/GCF_007361795.1_ASM736179v1/GCF_007361795.1_ASM736179v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/960/005/GCF_000960005.1_ASM96000v1/GCF_000960005.1_ASM96000v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/076/835/GCF_002076835.1_ASM207683v1/GCF_002076835.1_ASM207683v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/375/065/GCF_902375065.1_MGYG-HGUT-01444/GCF_902375065.1_MGYG-HGUT-01444_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/819/245/GCF_022819245.1_ASM2281924v1/GCF_022819245.1_ASM2281924v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/025/GCF_900637025.1_46338_H01/GCF_900637025.1_46338_H01_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/186/885/GCF_900186885.1_48903_D01/GCF_900186885.1_48903_D01_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14518/GCF_000845085.1_ViralProj14518_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/945/GCF_000026945.1_ASM2694v1/GCF_000026945.1_ASM2694v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/054/445/GCF_003054445.1_ASM305444v1/GCF_003054445.1_ASM305444v1_genomic.fna.gz")

cds.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/753/085/GCF_001753085.1_ASM175308v1/GCF_001753085.1_ASM175308v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/361/795/GCF_007361795.1_ASM736179v1/GCF_007361795.1_ASM736179v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/960/005/GCF_000960005.1_ASM96000v1/GCF_000960005.1_ASM96000v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/076/835/GCF_002076835.1_ASM207683v1/GCF_002076835.1_ASM207683v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/375/065/GCF_902375065.1_MGYG-HGUT-01444/GCF_902375065.1_MGYG-HGUT-01444_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/819/245/GCF_022819245.1_ASM2281924v1/GCF_022819245.1_ASM2281924v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/025/GCF_900637025.1_46338_H01/GCF_900637025.1_46338_H01_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/186/885/GCF_900186885.1_48903_D01/GCF_900186885.1_48903_D01_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14518/GCF_000845085.1_ViralProj14518_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/945/GCF_000026945.1_ASM2694v1/GCF_000026945.1_ASM2694v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/054/445/GCF_003054445.1_ASM305444v1/GCF_003054445.1_ASM305444v1_cds_from_genomic.fna.gz")

anno.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/753/085/GCF_001753085.1_ASM175308v1/GCF_001753085.1_ASM175308v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/361/795/GCF_007361795.1_ASM736179v1/GCF_007361795.1_ASM736179v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/960/005/GCF_000960005.1_ASM96000v1/GCF_000960005.1_ASM96000v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/076/835/GCF_002076835.1_ASM207683v1/GCF_002076835.1_ASM207683v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/375/065/GCF_902375065.1_MGYG-HGUT-01444/GCF_902375065.1_MGYG-HGUT-01444_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/819/245/GCF_022819245.1_ASM2281924v1/GCF_022819245.1_ASM2281924v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/025/GCF_900637025.1_46338_H01/GCF_900637025.1_46338_H01_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/186/885/GCF_900186885.1_48903_D01/GCF_900186885.1_48903_D01_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14518/GCF_000845085.1_ViralProj14518_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/945/GCF_000026945.1_ASM2694v1/GCF_000026945.1_ASM2694v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/054/445/GCF_003054445.1_ASM305444v1/GCF_003054445.1_ASM305444v1_genomic.gff.gz")

##### check if all species have full genome, cds and annotation files
length(full.genome.urls)
# [1] 19
length(cds.genome.urls)
# [1] 19
length(anno.genome.urls)
# [1] 19


##### Upload files
for (i in full.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    if (file.exists(file.path(dir.genomes, filename), recursive = TRUE) == FALSE){
        system(paste("wget ", i, " -O ", file.path(dir.genomes, filename), sep=""))
    }
}

for (i in cds.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    if (file.exists(file.path(dir.genomes, filename), recursive = TRUE) == FALSE){
        system(paste("wget ", i, " -O ", file.path(dir.cds.genomes, filename), sep=""))
    }
}

for (i in anno.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    if (file.exists(file.path(dir.genomes, filename), recursive = TRUE) == FALSE){
        system(paste("wget ", i, " -O ", file.path(dir.anno, filename), sep=""))
    }
}

##### check if files are uploaded
length(list.files(dir.genomes, recursive = TRUE))
# [1] 40
length(list.files(dir.cds.genomes, recursive = TRUE))
# [1] 40
length(list.files(dir.anno, recursive = TRUE))
# [1] 40

##### classify all genomes to respective organism folders
subfolders <- c("bacteria", "virus", "fungi")

gmap_build <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Packages/gmap-2020-06-01/util/gmap_build"
gmap.dir <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Packages/gmap-2020-06-01"
nthreads <- 4

iit_store <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Packages/gmap-2020-06-01/src/iit_store"
### Note: if GSNAP was not originally installed in the environment run:
# system(paste("cp ", gmap.dir, "/util/fa_coords ", gmap.dir, "/src", sep=""))
# system(paste("cp ", gmap.dir, "/util/md_coords ", gmap.dir, "/src", sep=""))
# system(paste("cp ", gmap.dir, "/util/gmap_process ", gmap.dir, "/src", sep=""))
############ Create genome indices, annotation file and splicesites file for each species
for (j in subfolders){
    genomes <- list.files(file.path(dir.genomes, j))
    if (length(genomes) != 0){
        if (j %in% c("bacteria", "virus")){
            for (i in 1:length(genomes)){
                species.id <- head(strsplit(genomes[i], "[.]")[[1]],1)
                file.genome <- file.path(dir.genomes, j, genomes[i])
                file.anno <- file.path(dir.anno, j, list.files(file.path(dir.anno, j))[i])
                chr.names <- c()
                for (k in system(paste("zcat", file.genome, "| grep '>'", sep=" "), intern = TRUE)){
                    chr.names <- append(chr.names, strsplit(strsplit(k, " ")[[1]][1], ">")[[1]][2])
                }
                chr.names <- paste(unlist(chr.names), collapse=",")
                system(paste(gmap_build, " -d ", species.id, " -g ", file.genome, " -D ", gmap.dir, " -t ", nthreads, " -B ", gmap.dir, "/src", 
                " --circular ", chr.names, sep = ""))
                system(paste("zcat ", file.anno, "| ", iit_store, " -G -l -F -o ", species.id, sep = ""))
            }
        } else {
            for (i in 1:length(genomes)){
                species.id <- head(strsplit(genomes[i], "[.]")[[1]],1)
                file.genome <- file.path(dir.genomes, j, genomes[i])
                file.anno <- file.path(dir.anno, j, list.files(file.path(dir.anno, j))[i])
                system(paste(gmap_build, " -d ", species.id, " -g ", file.genome, " -D ", gmap.dir, " -t ", nthreads, " -B ", gmap.dir, "/src", 
                sep = ""))
                system(paste("zcat ", file.anno, "| ", iit_store, " -G -l -F -o ", species.id, sep = ""))
            }
        }
    }
}

############# Align samples against genomes
##### create dir for pathogen detection
dir.patho.detect <- "/home/cluster/abalay/scratch/moved_from_home/Pathogen_detection_WGS"
if (!dir.exists(dir.patho.detect)){
    dir.create(dir.patho.detect, recursive = TRUE)
}
#### create list with all FASTQ files unmapped to human genome
fastq.vector <- c()
for (j in system(paste("find /home/cluster/abalay/scratch/moved_from_home/WGS_RAW_DLBCL -name '*_unmapped.fastq*'", sep=" "), intern = TRUE)){
    fastq.vector <- append(fastq.vector, j)
}


main.file <- file(paste(dir.patho.detect , "/scheduler_pathogen_gsnap.sh", sep=""))
for (i in 1:length(list.files(dir.genomes, recursive = TRUE))){
    species.id <- basename(head(strsplit(list.files(dir.genomes, recursive = TRUE)[i], "[.]")[[1]],1))
    file.genome <- file.path(dir.genomes, list.files(dir.genomes, recursive = TRUE)[i])
    file.anno <- paste("/home/cluster/abalay/", species.id, ".iit", sep="")
    sub.species.dir <- paste(dir.patho.detect, "/", species.id, sep="")
    if (!dir.exists(sub.species.dir)){
        dir.create(sub.species.dir, recursive = TRUE)
    }
    for (j in fastq.vector){
        if (length(grep('_r1_', j)) > 0){
            sub.sample.dir <- paste(sub.species.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], sep="")
            if (!dir.exists(sub.sample.dir)){
                dir.create(sub.sample.dir, recursive = TRUE)
                dir.create(paste(sub.sample.dir, "/Scripts", sep=""))
            }
            file.r1 <- j
            file.r2 <- paste(strsplit(j, "_r1_")[[1]][1], "_r2_", strsplit(j, "_r1_")[[1]][2],sep="")
            cat("#!/bin/bash\n", "#SBATCH --time=8:00:00\n", "#SBATCH --cpus-per-task=8\n", "#SBATCH --mem-per-cpu=16G\n", "#SBATCH --job-name=", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".fastq")[[1]][1], ".", species.id, "_gsnap.sh", " -o ", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".fastq")[[1]][1],".", species.id, 
                "_gsnap.out -e ", file.path(sub.sample.dir, "/Scripts/"), strsplit(basename(j), ".fastq")[[1]][1], ".", species.id,
                "_gsnap.err\n", 
                "export TMPDIR=", sub.sample.dir, "/Scripts/tmp\n",
                "mkdir -p $TMPDIR\n",
                "cd $TMPDIR\n",
                "export script_dir=/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files\n",
                # "export samtools=$script_dir/Packages/samtools-1.10/samtools\n",
                "export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap\n",
                "export NTHREADS=4\n",
                "$gsnap --gunzip -d ", species.id, " -D $script_dir/Packages/gmap-2020-06-01  -t $NTHREADS  -n 10 -g ", file.anno,
                " -A sam ", file.r1, " ", file.r2, " |samtools view -hb|  samtools sort -@ $NTHREADS > ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam\n", 
                "samtools index ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam.bai\n",
                "samtools stats ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], "_map.stats.txt\n",
                "rm -rf $TMPDIR\n", "chmod +x ", paste(sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh", sep=""), file=paste(sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh", sep=""), sep="", append=TRUE)
            cat("sbatch ", sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh\n", file=paste(dir.patho.detect , "/scheduler_pathogen_gsnap.sh", sep=""), sep="", 
                append=TRUE)
        }
        if (length(grep('_r1_', j)) == 0 & length(grep('_r2_', j)) == 0){
            sub.sample.dir <- paste(sub.species.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], sep="")
            if (!dir.exists(sub.sample.dir)){
            dir.create(sub.sample.dir, recursive = TRUE)
            dir.create(paste(sub.sample.dir, "/Scripts", sep=""))
            }
            file <- j
            cat("#!/bin/bash\n", "#SBATCH --time=8:00:00\n", "#SBATCH --cpus-per-task=8\n", "#SBATCH --mem-per-cpu=16G\n", "#SBATCH --job-name=", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".fastq")[[1]][1], ".", species.id, "_gsnap.sh", " -o ", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".fastq")[[1]][1],".", species.id, 
                "_gsnap.out -e ", file.path(sub.sample.dir, "/Scripts/"), strsplit(basename(j), ".fastq")[[1]][1], ".", species.id,
                "_gsnap.err\n", 
                "export TMPDIR=", sub.sample.dir, "/Scripts/tmp\n",
                "mkdir -p $TMPDIR\n",
                "cd $TMPDIR\n",
                "export script_dir=/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files\n",
                # "export samtools=$script_dir/Packages/samtools-1.10/samtools\n",
                "export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap\n",
                "export NTHREADS=4\n",
                "$gsnap --gunzip -d ", species.id, " -D $script_dir/Packages/gmap-2020-06-01  -t $NTHREADS  -n 10 -g ", file.anno,
                " -A sam ", file, " |samtools view -hb|  samtools sort -@ $NTHREADS > ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam\n", 
                "samtools index ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam.bai\n",
                "samtools stats ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], "_map.stats.txt\n",
                "rm -rf $TMPDIR\n", "chmod +x ", paste(sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh", sep=""), file=paste(sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh", sep=""), sep="", append=TRUE)
            cat("sbatch ", sub.sample.dir, "/Scripts/", strsplit(basename(j), ".fastq")[[1]][1],
                ".", species.id, "_gsnap.sh\n", file=paste(dir.patho.detect , "/scheduler_pathogen_gsnap.sh", sep=""), sep="", 
                append=TRUE)
        }
    }
}


################ get coverage values across whole genome and cds
#### create list with all BAM files mapped to specific pathogen
bam.vector <- c()
for (j in system(paste("find ", dir.patho.detect, "/", list.files(dir.patho.detect)[1], " -name '*.bam'", sep=""), intern = TRUE)){
    bam.vector <- append(bam.vector, j)
}

for (i in 1:length(list.files(dir.genomes, recursive = TRUE))){
    species.id <- basename(head(strsplit(list.files(dir.genomes, recursive = TRUE)[i], "[.]")[[1]],1))
    file.genome <- file.path(dir.genomes, list.files(dir.genomes, recursive = TRUE)[i])
    file.cds.genome <- file.path(dir.cds.genomes, list.files(dir.cds.genomes, recursive = TRUE)[i])
    sub.species.dir <- paste(dir.patho.detect, "/", species.id, sep="")
    for (j in bam.vector){
        sub.sample.dir <- paste(sub.species.dir, "/", strsplit(basename(j), ".bam")[[1]][1], sep="")
        file <- file.path(sub.sample.dir, basename(j))
        cat("#!/bin/bash\n", "#SBATCH --time=8:00:00\n", "#SBATCH --cpus-per-task=8\n", "#SBATCH --mem-per-cpu=16G\n", "#SBATCH --job-name=", file.path(sub.sample.dir, "/Scripts/"),
        strsplit(basename(j), ".bam")[[1]][1], ".", species.id, "_bedtools.sh", " -o ", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".bam")[[1]][1],".", species.id, 
                "_bedtools.out -e ", file.path(sub.sample.dir, "/Scripts/"), strsplit(basename(j), ".bam")[[1]][1], ".", species.id,
                "_bedtools.err\n", 
                "export script_dir=/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Packages\n",
                "export bedtools=$script_dir/bedtools2/bin/bedtools\n",
                "$bedtools genomecov -ibam ", file, " > ",
                sub.sample.dir, "/",
                strsplit(basename(j), ".bam")[[1]][1], ".", species.id, "_bedtools_genomecov.txt\n\n", 
                "$bedtools genomecov -ibam ", file, " -d > ",
                sub.sample.dir, "/",
                strsplit(basename(j), ".bam")[[1]][1], ".", species.id, "_bedtools_cdscov.txt\n",
                file = paste(sub.sample.dir, "/Scripts/",
                strsplit(basename(j), ".bam")[[1]][1], ".", species.id, "_bedtools.sh", sep=""), sep="", append=TRUE)
        cat("sbatch ", sub.sample.dir, "/Scripts/", strsplit(basename(j), ".bam")[[1]][1],
            ".", species.id, "_bedtools.sh\n", file=paste(dir.patho.detect , "/scheduler_pathogen_bedtools.sh", sep=""), sep="", 
            append=TRUE)
    }
}





################ estimate alignment criteria, pathogen load and genome coverage

### to estimate pathogen load we need to retrieve alignment info on Homo Sapiens first
dbgap_wgs.stats.file <- "/home/cluster/abalay/scratch/moved_from_home/WGS_RAW_DLBCL/multiqc_data_1/multiqc_samtools_stats.txt"


### to estimate cds coverage upload features tables
### create dir with features table
dir.feature.table <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Files_for_scripts/Pathogen_feature_tables"
if (!dir.exists(dir.feature.table)){
    dir.create(dir.feature.table, recursive = TRUE)
}

homo.stats.df <- data.frame(read.table(dbgap_wgs.stats.file, sep="\t", header=TRUE))

features.tables.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_feature_table.txt.gz",            
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_feature_table.txt.gz",          
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/753/085/GCF_001753085.1_ASM175308v1/GCF_001753085.1_ASM175308v1_feature_table.txt.gz",        
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/361/795/GCF_007361795.1_ASM736179v1/GCF_007361795.1_ASM736179v1_feature_table.txt.gz",        
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/960/005/GCF_000960005.1_ASM96000v1/GCF_000960005.1_ASM96000v1_feature_table.txt.gz",          
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/076/835/GCF_002076835.1_ASM207683v1/GCF_002076835.1_ASM207683v1_feature_table.txt.gz",        
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/375/065/GCF_902375065.1_MGYG-HGUT-01444/GCF_902375065.1_MGYG-HGUT-01444_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_feature_table.txt.gz",      
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/819/245/GCF_022819245.1_ASM2281924v1/GCF_022819245.1_ASM2281924v1_feature_table.txt.gz",      
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/025/GCF_900637025.1_46338_H01/GCF_900637025.1_46338_H01_feature_table.txt.gz",           
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_feature_table.txt.gz",      
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_feature_table.txt.gz",        
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/186/885/GCF_900186885.1_48903_D01/GCF_900186885.1_48903_D01_feature_table.txt.gz",            
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_feature_table.txt.gz",        
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_feature_table.txt.gz",  
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/085/GCF_000845085.1_ViralProj14518/GCF_000845085.1_ViralProj14518_feature_table.txt.gz",  
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/945/GCF_000026945.1_ASM2694v1/GCF_000026945.1_ASM2694v1_feature_table.txt.gz",            
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/885/GCF_001417885.1_Kmar_1.0/GCF_001417885.1_Kmar_1.0_feature_table.txt.gz",              
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/054/445/GCF_003054445.1_ASM305444v1/GCF_003054445.1_ASM305444v1_feature_table.txt.gz")

# num <- c()
# for (i in select){
#     n <- strsplit(basename(i), "[.]")[[1]][1]
#     num <- append(num, grep(n, list.dirs(dir.patho.detect, recursive = FALSE)))
# }


# select <- c()
# for ( i in anno.genome.urls){
#      n <- strsplit(basename(i), "[.]")[[1]][1]
#      if (length(grep(n, features.tables.urls, value=TRUE)) != 0){
#          select <- append(select, features.tables.urls[grep(n, features.tables.urls)])
#      }
# }

# for (i in setdiff(list.dirs(dir.patho.detect, recursive = FALSE), list.dirs(dir.patho.detect, recursive = FALSE)[num])){
#     system(paste0("rm -r ", i))
# }
##### Upload files
for (i in features.tables.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    if (file.exists(file.path(dir.feature.table, filename), recursive = TRUE) == FALSE){
        system(paste("wget ", i, " -O ", file.path(dir.feature.table, filename), sep=""))
    }
}

### to add kraken2 fractions we need to match NCBI ID to biological pathogen names
pathogen.names <- c("Staphylococcus_aureus", "Candida_dubliniensis", "Mycobacterium_tuberculosis", 
    "Human_mastadenovirus_C", "Human_betaherpesvirus_5", "Streptococcus_mitis", "Kluyveromyces_marxianus",
    "Ralstonia_insidiosa", "Legionella_pneumophila", "Streptococcus_pneumoniae", "Human_gammaherpesvirus_4",
    "Pichia_kudriavzevii", "Legionella_israelensis", "Klebsiella_pneumoniae", "Ralstonia_pickettii",
    "Limosilactobacillus_fermentum", "Veillonella_parvula", "Streptococcus_oralis", "Veillonella_atypica")

##### upload files with kraken2 bacteria and virus fractions
common.bacteria.file <- read.table(file.path(dir.figures, "Common.bacteria.kraken.dataset.txt"), sep="\t", header=TRUE)
common.virus.file <- read.table(file.path(dir.figures, "Common.virus.kraken.dataset.txt"), sep="\t", header=TRUE)
common.fungi.file <- read.table(file.path(dir.figures, "Common.fungi.kraken.dataset.txt"), sep="\t", header=TRUE)
common.pathogen.file <- cbind(common.bacteria.file, common.virus.file, common.fungi.file)

### estimate # of mapped reads, # of accurately mapped bases, error rate = # of mismatches / # of accurately mapped bases,
### average quality of mapped bases to pathogen species, pathogen load, whole genome coverage, cds coverage and kraken2 fraction
big.pathogen.stats.list <- list()

for (i in 1:length(list.dirs(dir.patho.detect, recursive = FALSE))){

    sample.names <- c(); nreads <- c(); nmapped.bases <- c(); err.rate <- c(); ave.qual <- c(); patho.load <- c(); 
    whole.gen.cov <- c()
    cds.gen.cov <- c()
    kraken2.frac <- c()

    full.species.path <- list.dirs(dir.patho.detect, recursive = FALSE)[i]
    species.id <- basename(full.species.path)
    cat("Processing species ", species.id, "\n")

    feature.table.path <- strsplit(file.path(dir.feature.table, grep(species.id, list.files(dir.feature.table), value=TRUE)), ".gz")[[1]][1]
    if (file.exists(feature.table.path) == FALSE){
        system(paste("gzip -d ", file.path(dir.feature.table, grep(species.id, list.files(dir.feature.table), value=TRUE)), sep=""), intern=TRUE)
    }


    feature.table <- read.table(feature.table.path, sep="\t", quote="")
    colnames(feature.table) <- c("# feature", "class", "assembly", "assembly_unit", "seq_type", "chromosome", "genomic_accession", "start", "end", "strand", "product_accession", "non-redundant_refseq", "related_accession", "name", "symbol", "GeneID", "locus_tag", "feature_interval_length", "product_length", "attributes")

    file.names.vector <- system(paste("find ",full.species.path, " -name *_map.stats.txt", sep=""), intern=TRUE)
    for (j in file.names.vector){
        cat("Processing sample ", j, "\n")

        # sample name
        sample.names <- append(sample.names, strsplit(basename(j), ".stats.txt")[[1]][1])

        # nreads
        Nreads <- as.numeric(strsplit(grep("reads mapped:", readLines(j), value=TRUE), "\t")[[1]][3])
        nreads <- append(nreads, Nreads)

        #nmapped.bases
        Nmapped.bases <- as.numeric(strsplit(grep("SN\tbases mapped \\(cigar\\)", readLines(j), value=TRUE), "\t")[[1]][3])
        nmapped.bases <- append(nmapped.bases, Nmapped.bases)

        # error rate
        error.rate <- as.numeric(strsplit(grep("error rate:", readLines(j), value=TRUE), "\t")[[1]][3])
        err.rate <- append(err.rate, error.rate)

        # average quality
        ave.qual <- append(ave.qual, as.numeric(strsplit(grep("average quality:", readLines(j), value=TRUE), "\t")[[1]][3]))

        # pathogen load
        diff.naming <- c("HIV4_2_unmapped_map.stats.txt", "HIV4_1_unmapped_map.stats.txt", "EBV_GM_Ad6_unmapped_map.stats.txt",
        "EBV_CP33_Ad5_unmapped_map.stats.txt", "EBV_CP8_Ad4_unmapped_map.stats.txt", "EBV_AP33_Ad3_unmapped_map.stats.txt",
        "HIV1_1_unmapped_map.stats.txt")

        if (!basename(j) %in% diff.naming){
            nreads.homo <- homo.stats.df[grep(strsplit(basename(j), "_")[[1]][1], homo.stats.df$Sample), "reads_mapped"]
            nmapped.bases.homo <- homo.stats.df[grep(strsplit(basename(j), "_")[[1]][1], homo.stats.df$Sample), "bases_mapped_.cigar."]
        }
        if (basename(j) %in% diff.naming){
            nreads.homo <- homo.stats.df[grep(strsplit(basename(j), "_unmapped")[[1]][1], homo.stats.df$Sample), "reads_mapped"]
            nmapped.bases.homo <- homo.stats.df[grep(strsplit(basename(j), "_unmapped")[[1]][1], homo.stats.df$Sample), "bases_mapped_.cigar."]
        }
        
        ploidy <- 2; purity <- 1; homo.transcriptome.size <- 2897310462

        genome.cov.file <- paste(full.species.path, "/", strsplit(basename(j), "_map.stats.txt")[[1]][1], "/", 
            strsplit(basename(j), "_map.stats.txt")[[1]][1], ".", species.id, 
            "_bedtools_genomecov.txt", sep="")
        genome.cov.file <- read.table(genome.cov.file, sep="\t")

        pathogen.genome.size <- genome.cov.file[which(genome.cov.file$V1 == "genome" & genome.cov.file$V2 == 0), "V4"]
        viral.portion <- (Nreads * Nmapped.bases)/pathogen.genome.size 
        host.portion <- (nreads.homo * nmapped.bases.homo * purity)/(homo.transcriptome.size * ploidy)
        pathogen.load <- viral.portion/host.portion
        patho.load <- append(patho.load, pathogen.load)
        print(pathogen.load)

        # whole genome coverage
        whole.genome.cov <- 1 - genome.cov.file[which(genome.cov.file$V1 == "genome" & genome.cov.file$V2 == 0), "V5"]
        whole.gen.cov <- append(whole.gen.cov, whole.genome.cov)

        # cds coverage
        cds.table <- feature.table %>% filter(`# feature` == "CDS") %>% select(c("# feature", "genomic_accession", "start", "symbol", "feature_interval_length")) %>%
            mutate(end = start + feature_interval_length) %>% as.data.frame
        ### keep duplicate symbol with longer feature interval length
        cds.table <- cds.table %>% arrange(symbol, -feature_interval_length) %>% filter(duplicated(symbol) == FALSE)  
        coordinates <- list()
        for (k in 1:nrow(cds.table)){
            coordinates[[k]] <- seq(cds.table$start[k], cds.table$end[k], by=1)
            len <- length(seq(cds.table$start[k], cds.table$end[k], by=1))
            names(coordinates)[k] <- cds.table$genomic_accession[k]
        }
        cds.cov.file <- paste(full.species.path, "/", strsplit(basename(j), "_map.stats.txt")[[1]][1], "/", 
            strsplit(basename(j), "_map.stats.txt")[[1]][1], ".", species.id, 
            "_bedtools_cdscov.txt", sep="")
        cds.cov.file <- read.table(cds.cov.file, sep="\t")
        cds.cov.file <- cds.cov.file[cds.cov.file$V3 > 0, ]
        covered.bases <- 0; total <- 0
        for (l in unique(cds.table$genomic_accession)){
            coords.subset <- unlist(coordinates[which(names(coordinates) == l)])
            coords.subset <- coords.subset[!duplicated(coords.subset)]
            cds.cov.file.subset <- cds.cov.file[cds.cov.file$V1 == l,]
            covered.bases <- covered.bases + length(intersect(cds.cov.file.subset$V2, coords.subset))
            total <- total + length(coords.subset)
        }
        cds.cov <- covered.bases/total
        cds.gen.cov <- append(cds.gen.cov, cds.cov)

        # kraken2 fraction
        pathogen.name <- pathogen.names[grep(species.id, features.tables.urls)]
        if (!basename(j) %in% diff.naming){
            kraken2.fraction <- common.pathogen.file[grep(strsplit(basename(j), "_")[[1]][1], rownames(common.pathogen.file)), pathogen.name]
        }
        if (basename(j) %in% diff.naming){
            kraken2.fraction <- common.pathogen.file[grep(strsplit(basename(j), "_unmapped")[[1]][1], rownames(common.pathogen.file)), pathogen.name]
        }
        kraken2.frac <- append(kraken2.frac, kraken2.fraction)
        print(kraken2.fraction)
        print(length(kraken2.fraction))
    }
    species.df <- data.frame(Sample.ID = sample.names, Nreads = nreads, Nmapped.bases = nmapped.bases, Error.rate = err.rate, 
        Ave.quality = ave.qual, Pathogen.load = patho.load, Whole.genome.coverage = whole.gen.cov, CDS.coverage = cds.gen.cov,
        Kraken2.load = kraken2.frac)
    big.pathogen.stats.list[[i]] <- species.df
    names(big.pathogen.stats.list)[i] <- species.id
}


saveRDS(big.pathogen.stats.list, file.path(dir.patho.detect, "Pathogen.sample.info.RData"))


big.pathogen.stats.list <- readRDS(file.path(dir.patho.detect, "Pathogen.sample.info.RData"))

# for (i in 1:length(big.pathogen.stats.list)){
#     sample.ids <- unlist(big.pathogen.stats.list[[i]]["Sample.ID"]) %in% a
#     big.pathogen.stats.list[[i]] <- big.pathogen.stats.list[[i]][!sample.ids, ]
# }

actual.pathogen.names <- c("Staphylococcus aureus", "Candida dubliniensis", "Mycobacterium tuberculosis", 
    "Human mastadenovirus C", "Human betaherpesvirus 5", "Streptococcus mitis", "Kluyveromyces marxianus",
    "Ralstonia insidiosa", "Legionella pneumophila", "Streptococcus pneumoniae", "Human gammaherpesvirus 4 (EBV)",
    "Pichia kudriavzevii", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica")

species <- c("Bacteria", "Fungi", "Bacteria", "Virus", "Virus", "Bacteria", "Fungi", "Bacteria",
    "Bacteria", "Bacteria", "Virus", "Fungi", "Bacteria", "Bacteria", "Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria")

### subset error rates and average mapping quality for all pathogens
error.rate <- c(); ave.map.qual <- c(); pathogen.anno <- c(); species.anno <- c()
for (i in 1:length(big.pathogen.stats.list)){
    error.rate <- append(error.rate, big.pathogen.stats.list[[i]]$Error.rate)
    ave.map.qual <- append(ave.map.qual, big.pathogen.stats.list[[i]]$Ave.quality)
    pathogen.anno <- append(pathogen.anno, rep(actual.pathogen.names[i], length(big.pathogen.stats.list[[i]]$Error.rate)))
    species.anno <- append(species.anno, rep(species[i], length(big.pathogen.stats.list[[i]]$Error.rate)))
}

df <- data.frame(Error.rate = error.rate, Ave.Map.Qual = ave.map.qual, Pathogen = pathogen.anno, Species = species.anno)

# axis.x.colors <- c("red", "blue", "springgreen"); names(axis.x.colors) <- unique(df$Species); df$color <- as.vector(axis.x.colors[df$Species])

df$Pathogen <- factor(x = df$Pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica",
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5",
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

df$Species <- factor(x = df$Species, levels = c("Bacteria", "Virus", "Fungi"))

### average mapping quality for all reads across all pathogens
unique(df$Ave.Map.Qual)
# [1] 255

### create directory for the pathogen detection figures
dir.final.figures <- paste(dir.patho.detect, "/Figures", sep="")
if (!dir.exists(dir.final.figures)){
    dir.create(dir.final.figures, recursive = TRUE)
}

### plot error rates per pathogen across samples
p <- ggplot(df, aes(x = Error.rate, y = Pathogen, fill = Pathogen, color=Species)) + 
    geom_violin(trim = FALSE) +
    scale_fill_manual(values=magma(length(unique(df$Pathogen)))) +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.05, color = "brown4", linetype="dashed") +
    xlab("Error rate") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("Error rate across pathogens")

png(paste(dir.final.figures, "/Error.rate.dist.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
p + guides(fill="none")
dev.off()



### subset pathogen load, kraken2 fractions, whole genome coverage and cds coverage
pathogen.load.list <- list()
kraken2.fraction.list <- list()
wgcov.list <- list()
cdscov.list <- list()

for (i in 1:length(big.pathogen.stats.list)){
    pathogen.load.list[[i]] <- big.pathogen.stats.list[[i]]$Pathogen.load
    names(pathogen.load.list)[i] <- actual.pathogen.names[i]
    kraken2.fraction.list[[i]] <- big.pathogen.stats.list[[i]]$Kraken2.load
    names(kraken2.fraction.list)[i] <- actual.pathogen.names[i]
    wgcov.list[[i]] <- big.pathogen.stats.list[[i]]$Whole.genome.coverage
    names(wgcov.list)[i] <- actual.pathogen.names[i]
    cdscov.list[[i]] <- big.pathogen.stats.list[[i]]$CDS.coverage
    names(cdscov.list)[i] <- actual.pathogen.names[i]
}

sample.ids <- big.pathogen.stats.list[[1]]$Sample.ID

pathogen.load.df <- data.frame(do.call(cbind, pathogen.load.list), row.names = sample.ids)
dim(pathogen.load.df)
# [1] 183  19
kraken2.fraction.df <- data.frame(do.call(cbind, kraken2.fraction.list), row.names = sample.ids)
dim(kraken2.fraction.df)
# [1] 183  19
wgcov.df <- data.frame(do.call(cbind, wgcov.list), row.names = sample.ids)
dim(wgcov.df)
# [1] 183  19
cdscov.df <- data.frame(do.call(cbind, cdscov.list), row.names = sample.ids)
dim(cdscov.df)
# [1] 183  19
 
### Heatmap: pathogen load of samples across pathogens
colnames(pathogen.load.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "chocolate4", "Virus" = "blue", "Fungi" = "limegreen"))

h1_top <- HeatmapAnnotation(`Pathogen type` = species, col = species.colors)

viridis.colors <- colorRamp2(seq(0, 0.0004, 0.00004), viridis(11)) ### set color scheme

h1 <- Heatmap(as.matrix(log(1+pathogen.load.df)), name = "log(1+Pathogen load)", col = viridis.colors, 
column_title = "Pathogen",
column_names_rot = 45,
cluster_columns = FALSE,
column_names_side = "top",
column_names_gp = gpar(fontsize = 10),
column_split = species,
row_title = NULL,  
show_row_names=FALSE,
cluster_rows = FALSE,   
top_annotation=h1_top, 
border = TRUE)

png(paste(dir.final.figures, "/Pathogen.load.heatmap.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
h1
dev.off()

### Heatmap: kraken2 fraction of samples across pathogens
colnames(kraken2.fraction.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "chocolate4", "Virus" = "blue", "Fungi" = "limegreen"))

h2_top <- HeatmapAnnotation(`Pathogen type` = species, col = species.colors)

viridis.colors <- colorRamp2(seq(0, 0.05, by=0.001), viridis(51)) ### set color scheme

h2 <- Heatmap(as.matrix(kraken2.fraction.df), name = "Pathogen fraction", col = viridis.colors, 
column_title = "Pathogen",
column_names_rot = 45,
cluster_columns = FALSE,
column_names_side = "top",
column_names_gp = gpar(fontsize = 10),
column_split = species,
row_title = NULL,  
show_row_names=FALSE,
cluster_rows = FALSE,   
top_annotation=h2_top, 
border = TRUE)

png(paste(dir.final.figures, "/Kraken2.fraction.heatmap.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
h2
dev.off()

### Heatmap: whole genome coverage of samples across pathogens
colnames(wgcov.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "chocolate4", "Virus" = "blue", "Fungi" = "limegreen"))

h3_top <- HeatmapAnnotation(`Pathogen type` = species, col = species.colors)

viridis.colors <- colorRamp2(seq(0, 1, by=0.05), viridis(21)) ### set color scheme

h3 <- Heatmap(as.matrix(wgcov.df), name = "WG Cov", col = viridis.colors, 
column_title = "Pathogen",
column_names_rot = 45,
cluster_columns = FALSE,
column_names_side = "top",
column_names_gp = gpar(fontsize = 10),
column_split = species,
row_title = NULL,  
show_row_names=FALSE,
cluster_rows = FALSE,   
top_annotation=h3_top, 
border = TRUE)

png(paste(dir.final.figures, "/WG.Cov.heatmap.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
h3
dev.off()


### Heatmap: cds coverage of samples across pathogens
colnames(cdscov.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "chocolate4", "Virus" = "blue", "Fungi" = "limegreen"))

h4_top <- HeatmapAnnotation(`Pathogen type` = species, col = species.colors)

viridis.colors <- colorRamp2(seq(0, 0.5, by=0.05), viridis(11)) ### set color scheme

h4 <- Heatmap(as.matrix(cdscov.df), name = "CDS Cov", col = viridis.colors, 
column_title = "Pathogen",
column_names_rot = 45,
cluster_columns = FALSE,
column_names_side = "top",
column_names_gp = gpar(fontsize = 10),
column_split = species,
row_title = NULL,  
show_row_names=FALSE,
cluster_rows = FALSE,   
top_annotation=h4_top, 
border = TRUE)

png(paste(dir.final.figures, "/CDS.Cov.heatmap.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
h4
dev.off()




### subset wg cov and cds cov again to plot as point distribution
wgs.cov <- c(); cds.cov <- c(); pathogen.anno <- c(); species.anno <- c()
for (i in 1:length(big.pathogen.stats.list)){
    wgs.cov <- append(wgs.cov, big.pathogen.stats.list[[i]]$Whole.genome.coverage)
    cds.cov <- append(cds.cov, big.pathogen.stats.list[[i]]$CDS.coverage)
    pathogen.anno <- append(pathogen.anno, rep(actual.pathogen.names[i], length(big.pathogen.stats.list[[i]]$Whole.genome.coverage)))
    species.anno <- append(species.anno, rep(species[i], length(big.pathogen.stats.list[[i]]$Whole.genome.coverage)))
}

df <- data.frame(WG.cov = wgs.cov, CDS.cov = cds.cov, Pathogen = pathogen.anno, Species = species.anno)

# axis.x.colors <- c("red", "blue", "springgreen"); names(axis.x.colors) <- unique(df$Species); df$color <- as.vector(axis.x.colors[df$Species])

df$Pathogen <- factor(x = df$Pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

df$Species <- factor(x = df$Species, levels = c("Bacteria", "Virus", "Fungi"))

### plot wg cov distribution for all samples
wg.cov.p <- ggplot(df, aes(x = WG.cov, y = Pathogen, fill = Pathogen, color=Species)) + 
    geom_violin(trim = FALSE) +
    scale_fill_manual(values=magma(length(unique(df$Pathogen)))) +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.1, color = "brown4", linetype="dashed") +
    xlab("WG Coverage") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("WG Coverage across pathogens")

png(paste(dir.final.figures, "/WG.Cov.dist.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
wg.cov.p + guides(fill="none")
dev.off()

### plot cds cov distribution for all samples
cds.cov.p <- ggplot(df, aes(x = CDS.cov, y = Pathogen, fill = Pathogen, color=Species)) + 
    geom_violin(trim = FALSE) +
    scale_fill_manual(values=magma(length(unique(df$Pathogen)))) +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.3, color = "brown4", linetype="dashed") +
    xlab("CDS Coverage") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("CDS Coverage across pathogens")

png(paste(dir.final.figures, "/CDS.Cov.dist.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
cds.cov.p + guides(fill="none")
dev.off()

### Correlation heatmaps
pl.kf.rho <- c(); pl.kf.pval <- c()
pl.wg.rho <- c(); pl.wg.pval <- c()
pl.cds.rho <- c(); pl.cds.pval <- c()
kf.wg.rho <- c(); kf.wg.pval <- c()
kf.cds.rho <- c(); kf.cds.pval <- c()
wg.cds.rho <- c(); wg.cds.pval <- c()


for (i in 1:ncol(pathogen.load.df)){
    temp.df <- data.frame(`Pathogen load` = pathogen.load.df[,i], `Kraken2 fraction` = kraken2.fraction.df[,i],
        `WG Cov` = wgcov.df[,i], `CDS Cov` = cdscov.df[,i])
    spearman.test <- cor.test(temp.df$Pathogen.load, temp.df$Kraken2.fraction, method="spearman")
    pl.kf.rho <- append(pl.kf.rho, spearman.test$estimate)
    pl.kf.pval <- append(pl.kf.pval, spearman.test$p.value)

    spearman.test <- cor.test(temp.df$Pathogen.load, temp.df$WG.Cov, method="spearman")
    pl.wg.rho <- append(pl.wg.rho, spearman.test$estimate)
    pl.wg.pval <- append(pl.wg.pval, spearman.test$p.value)

    spearman.test <- cor.test(temp.df$Pathogen.load, temp.df$CDS.Cov, method="spearman")
    pl.cds.rho <- append(pl.cds.rho, spearman.test$estimate)
    pl.cds.pval <- append(pl.cds.pval, spearman.test$p.value)

    spearman.test <- cor.test(temp.df$Kraken2.fraction, temp.df$WG.Cov, method="spearman")
    kf.wg.rho <- append(kf.wg.rho, spearman.test$estimate)
    kf.wg.pval <- append(kf.wg.pval, spearman.test$p.value)

    spearman.test <- cor.test(temp.df$Kraken2.fraction, temp.df$CDS.Cov, method="spearman")
    kf.cds.rho <- append(kf.cds.rho, spearman.test$estimate)
    kf.cds.pval <- append(kf.cds.pval, spearman.test$p.value)

    spearman.test <- cor.test(temp.df$WG.Cov, temp.df$CDS.Cov, method="spearman")
    wg.cds.rho <- append(wg.cds.rho, spearman.test$estimate)
    wg.cds.pval <- append(wg.cds.pval, spearman.test$p.value)
}

### Line plot of spearman correlation coefficients
#### Pathogen load versus kraken2 fraction
temp.df <- data.frame(R = as.vector(pl.kf.rho), p.val = as.vector(pl.kf.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.6, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(labels = c("Bacteria", "Virus", "Fungi"), values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ Kraken2 fraction")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.kraken2_fraction.corr.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

#### Pathogen load vs whole genome coverage
temp.df <- data.frame(R = as.vector(pl.wg.rho), p.val = as.vector(pl.wg.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.1, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(labels = c("Bacteria", "Virus", "Fungi"), values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ WG Coverage")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.wg_coverage.corr.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Pathogen load vs cds coverage
temp.df <- data.frame(R = as.vector(pl.cds.rho), p.val = as.vector(pl.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.1, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(labels = c("Bacteria", "Virus", "Fungi"), values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.cds_coverage.corr.all.pathogens.WGS.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Kraken2 fraction vs whole genome coverage
temp.df <- data.frame(R = as.vector(kf.wg.rho), p.val = as.vector(kf.wg.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.6, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("Kraken2 fraction ~ WG Coverage")

png(paste(dir.final.figures, "/Line_plot.kraken2_fraction.vs.wg_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

#### Kraken2 fraction vs cds coverage
temp.df <- data.frame(R = as.vector(kf.cds.rho), p.val = as.vector(kf.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.6, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("Kraken2 fraction ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.kraken2_fraction.vs.cds_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Whole genome coverage vs cds coverage
temp.df <- data.frame(R = as.vector(wg.cds.rho), p.val = as.vector(wg.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Staphylococcus aureus", "Mycobacterium tuberculosis", 
    "Streptococcus mitis", "Ralstonia insidiosa", "Legionella pneumophila", 
    "Streptococcus pneumoniae", "Legionella israelensis", "Klebsiella pneumoniae", "Ralstonia pickettii",
    "Limosilactobacillus fermentum", "Veillonella parvula", "Streptococcus oralis", "Veillonella atypica", 
    "Human gammaherpesvirus 4 (EBV)", "Human mastadenovirus C", "Human betaherpesvirus 5", 
    "Pichia kudriavzevii", "Candida dubliniensis", "Kluyveromyces marxianus"))

temp.df$species <- factor(x = temp.df$species, levels = c("Bacteria", "Virus", "Fungi"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue", "limegreen")) + 
    theme_classic() +
    ggtitle("WG Coverage ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.wg_coverage.vs.cds_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

### Kraken2 and Pathogen Load didn't correlate with WG and CDS coverage

### Try classification method: WG coverage > 10%, CDS coverage > 30% and error rate < 5%
for (i in 1:length(big.pathogen.stats.list)){
    temp.vect <- rep(0, nrow(big.pathogen.stats.list[[i]]))
    coordinates <- which(big.pathogen.stats.list[[i]]$Error.rate > 0.05 
    | big.pathogen.stats.list[[i]]$Whole.genome.coverage < 0.1 
    | big.pathogen.stats.list[[i]]$CDS.coverage < 0.3) ### coordinates of no infection

    temp.vect[-coordinates] <- 1
    big.pathogen.stats.list[[i]]$Classification_by_cover_only <- temp.vect
    cat(actual.pathogen.names[i], table(big.pathogen.stats.list[[i]]$Classification_by_cover_only), "\n")
}

# Staphylococcus aureus 183 
# Candida dubliniensis 183 
# Mycobacterium tuberculosis 183 
# Human mastadenovirus C 181 2 
# Human betaherpesvirus 5 181 2 
# Streptococcus mitis 182 1 
# Kluyveromyces marxianus 183 
# Ralstonia insidiosa 182 1 
# Legionella pneumophila 183 
# Streptococcus pneumoniae 182 1 
# Human gammaherpesvirus 4 (EBV) 177 6 
# Pichia kudriavzevii 183 
# Legionella israelensis 183 
# Klebsiella pneumoniae 183 
# Ralstonia pickettii 182 1 
# Limosilactobacillus fermentum 183 
# Veillonella parvula 182 1 
# Streptococcus oralis 182 1 
# Veillonella atypica 182 1

### save the file with pathogen classification
saveRDS(big.pathogen.stats.list, file.path(dir.patho.detect, "Pathogen.sample.info.RData"))












