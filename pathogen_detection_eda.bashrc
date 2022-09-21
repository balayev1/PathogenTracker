###################################################
################## Pathogen-detection pipeline in cell line samples RNA-Seq
##############

module load hpc
module load anaconda3
source activate agshin_condaenv
srun --mem-per-cpu=4GB --time=6:00:00 --pty --cpus-per-task=8 bash
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

####### Exploratory data analysis
### Load dataset
parent.directory <- "/home/cluster/abalay/scratch/Cell_Lines/Kraken2_output"

kraken.file <- read.delim(file.path(parent.directory, "Cell_Lines.patho.combined.mpa.txt"), sep="\t", row.names=1)

### change column names
colnames(kraken.file) <- unlist(strsplit(colnames(kraken.file), ".breport"))

kraken.file <- t(kraken.file)

kraken.file <- kraken.file[grep("ERR", rownames(kraken.file)),]
dim(kraken.file)
# [1]   79 2471

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

png(file.path(dir.figures, "Total_aligned_taxonomic_reads.cell_lines.png"), res=200, unit="in", height=8, width=11)
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

png(file.path(dir.figures, "Fraction_of_taxonomic_reads.cell_lines.png"), res=200, unit="in", height=10, width=15)
p
dev.off()

bacteria.list <- list()
##### bacteria
### subset any fields with species detected
temp.kraken <- subset(kraken.file, select = grepl("k__Bacteria", colnames(kraken.file)))
dim(temp.kraken)
# [1]   79 1976

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1]   79 1253

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

bacteria.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Bacillus_cereus          Ralstonia_insidiosa 
#                 7.645930e-04                 7.422299e-04 
#   Bradyrhizobium_sp._PSBB068 Stenotrophomonas_maltophilia 
#                 3.317304e-04                 2.719918e-04 
#     Bradyrhizobium_sp._BTAi1         Mesorhizobium_terrae 
#                 2.568610e-04                 2.307383e-04 
#          Serratia_marcescens      Pseudomonas_fluorescens 
#                 1.310492e-04                 8.416719e-05 
#          Pseudomonas_sp._TKP        Staphylococcus_aureus 
#                 8.184029e-05                 5.520428e-05 

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
    #      Ralstonia_insidiosa         Mesorhizobium_terrae 
    #             0.0085879335                 0.0030409110 
    # Bradyrhizobium_sp._BTAi1   Bradyrhizobium_sp._PSBB068 
    #             0.0029152663                 0.0026343770 
    #          Bacillus_cereus Stenotrophomonas_maltophilia 
    #             0.0012012454                 0.0011592917 
    #      Serratia_marcescens          Ralstonia_pickettii 
    #             0.0005560365                 0.0005517577 
    #   Ralstonia_solanacearum    Ralstonia_mannitolilytica 
    #             0.0004879332                 0.0003653718 


virus.list <- list()
##### viruses
temp.kraken <- subset(kraken.file, select = grepl("k__Viruses", colnames(kraken.file)))
dim(temp.kraken)
# [1] 79 53


temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 79 11

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

virus.list[[1]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
#    Human_gammaherpesvirus_8    Human_gammaherpesvirus_4 
#                6.905337e-03                2.008352e-03 
#       Proteus_virus_Isfahan      Human_mastadenovirus_C 
#                8.675749e-04                2.300045e-05 
#  Papiine_gammaherpesvirus_1                Aichivirus_B 
#                1.555491e-05                2.109898e-06 
#                 Rotavirus_A  Ralstonia_virus_Raharianne 
#                1.722270e-07                1.352409e-07 
# Staphylococcus_virus_Andhra   Burkholderia_virus_Bcep22 
#                4.655464e-08                6.776647e-09

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
#    Human_gammaherpesvirus_8    Human_gammaherpesvirus_4 
#                4.324526e-02                4.039980e-03 
#       Proteus_virus_Isfahan      Human_mastadenovirus_C 
#                2.341409e-03                8.541751e-05 
#  Papiine_gammaherpesvirus_1                Aichivirus_B 
#                6.049077e-05                2.262823e-05 
#                 Rotavirus_A  Ralstonia_virus_Raharianne 
#                3.787474e-06                3.371310e-06 
# Staphylococcus_virus_Andhra   Burkholderia_virus_Bcep22 
#                7.813794e-07                2.890503e-07 

###################################################
################## Pathogen-detection pipeline in RNA-Seq DLBCL samples
##############
####### Exploratory data analysis
### Load dataset
parent.directory <- "/home/cluster/abalay/scratch/Dbgap_RNA_Seq_DLBCL/Kraken2_output"

kraken.file <- read.delim(file.path(parent.directory, "RNA-Seq.DLBCL.combined.mpa.txt"), sep="\t", row.names=1)

### change column names
colnames(kraken.file) <- unlist(strsplit(colnames(kraken.file), "_dbGaP.28434.breport"))

kraken.file <- t(kraken.file)

dim(kraken.file)
# [1]   71 6726

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

png(file.path(dir.figures, "Total_aligned_taxonomic_reads.RNA-SeqDLBCL.png"), res=200, unit="in", height=8, width=11)
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

png(file.path(dir.figures, "Fraction_of_taxonomic_reads.RNA-SeqDLBCL.png"), res=200, unit="in", height=10, width=15)
p
dev.off()

##### bacteria
### subset any fields with species detected
temp.kraken <- subset(kraken.file, select = grepl("k__Bacteria", colnames(kraken.file)))
dim(temp.kraken)
# [1]   71 6030

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1]   71 4170

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

bacteria.list[[2]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Acinetobacter_lwoffii        Acinetobacter_baumannii 
#                   0.0362780251                   0.0066589663 
#            Salmonella_enterica       Acinetobacter_sp._NEB149 
#                   0.0046891086                   0.0037917901 
#        Acinetobacter_johnsonii          Acinetobacter_sp._WY4 
#                   0.0023772871                   0.0016224480 
# Acinetobacter_sp._FDAARGOS_560               Escherichia_coli 
#                   0.0008923647                   0.0008053546 
#                Yersinia_pestis                Bacillus_cereus 
#                   0.0007783648                   0.0007309256     

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
# Acinetobacter_lwoffii        Acinetobacter_baumannii 
#                    0.198472193                    0.027189165 
#       Acinetobacter_sp._NEB149            Salmonella_enterica 
#                    0.021147821                    0.017678294 
#        Acinetobacter_johnsonii Acinetobacter_sp._FDAARGOS_560 
#                    0.013512632                    0.011130434 
#          Acinetobacter_sp._WY4          Pseudomonas_yamanorum 
#                    0.009013482                    0.008641744 
#                Yersinia_pestis               Escherichia_coli 
#                    0.004091680                    0.002841670

##### viruses
temp.kraken <- subset(kraken.file, select = grepl("k__Viruses", colnames(kraken.file)))
dim(temp.kraken)
# [1] 71 95


temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 71 38

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

virus.list[[2]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Proteus_virus_Isfahan            Human_gammaherpesvirus_4 
#                        1.508351e-04                        3.015329e-05 
# Snyder-Theilen_feline_sarcoma_virus       Drosophila_innubila_nudivirus 
#                        1.029788e-06                        7.530056e-07 
#                Escherichia_virus_T4       Abelson_murine_leukemia_virus 
#                        5.652180e-07                        3.650213e-07 
#           Acinetobacter_virus_Acj61               Yersinia_virus_PYPS2T 
#                        6.463877e-08                        4.190717e-08 
#                   Y73_sarcoma_virus                        Aichivirus_A 
#                        2.440640e-08                        2.355283e-08


colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
# Human_gammaherpesvirus_4               Proteus_virus_Isfahan 
#                        8.998724e-04                        4.234799e-04 
#                Escherichia_virus_T4 Snyder-Theilen_feline_sarcoma_virus 
#                        2.422097e-05                        8.556338e-06 
#       Drosophila_innubila_nudivirus       Abelson_murine_leukemia_virus 
#                        5.096291e-06                        2.978787e-06 
#               Yersinia_virus_PYPS2T                        Aichivirus_A 
#                        2.975409e-06                        1.672251e-06 
#            Escherichia_virus_EcoMF1      Human_immunodeficiency_virus_1 
#                        1.613587e-06                        1.287437e-06

###################################################
################## Pathogen-detection pipeline in RNA-Seq Sammir's DLBCL samples
##############
####### Exploratory data analysis
### Load dataset
parent.directory <- "/home/cluster/abalay/scratch/Sammirs_RNA_Seq/Cutadapt/Kraken2_output"

kraken.file <- read.delim(file.path(parent.directory, "Sammirs_RNA_Seq.patho.combined.mpa.txt"), sep="\t", row.names=1)

### change column names
colnames(kraken.file) <- unlist(strsplit(colnames(kraken.file), ".breport"))

kraken.file <- t(kraken.file)

dim(kraken.file)
# [1]   18 6824

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

png(file.path(dir.figures, "Total_aligned_taxonomic_reads.SammirsRNA-SeqDLBCL.png"), res=200, unit="in", height=8, width=11)
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

png(file.path(dir.figures, "Fraction_of_taxonomic_reads.SammirsRNA-SeqDLBCL.png"), res=200, unit="in", height=10, width=15)
p
dev.off()

##### bacteria
### subset any fields with species detected
temp.kraken <- subset(kraken.file, select = grepl("k__Bacteria", colnames(kraken.file)))
dim(temp.kraken)
# [1]   18 6162

temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1]   18 4225

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

bacteria.list[[3]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Salmonella_enterica          Escherichia_coli       Cutibacterium_acnes 
#              2.403691e-03              2.337013e-03              1.108338e-03 
#     Staphylococcus_aureus           Bacillus_cereus           Yersinia_pestis 
#              8.028474e-04              7.520303e-04              2.371718e-04 
#        Salmonella_bongori     Klebsiella_pneumoniae Burkholderia_pseudomallei 
#              1.438216e-04              1.287629e-04              1.177960e-04 
#       Priestia_megaterium 
#              9.459584e-05    

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
# Salmonella_enterica           Escherichia_coli 
#               0.0222493114               0.0119918778 
#        Cutibacterium_acnes      Staphylococcus_aureus 
#               0.0054590971               0.0032292335 
#            Yersinia_pestis            Bacillus_cereus 
#               0.0027541135               0.0027137077 
#         Salmonella_bongori  Burkholderia_pseudomallei 
#               0.0017537018               0.0011682919 
#  Citrobacter_sp._RHB21-C01 Staphylococcus_epidermidis 
#               0.0008999928               0.0006632246

##### viruses
temp.kraken <- subset(kraken.file, select = grepl("k__Viruses", colnames(kraken.file)))
dim(temp.kraken)
# [1]  18 118


temp.kraken <- temp.kraken[, 2:ncol(temp.kraken)]/mini$Total_reads
temp.kraken[is.na(temp.kraken)] <- 0

temp.kraken <- temp.kraken[, grepl("s__", colnames(temp.kraken))]
dim(temp.kraken)
# [1] 18 43

colnames(temp.kraken) <- gsub(".*s__", "", colnames(temp.kraken))

virus.list[[3]] <- data.frame(temp.kraken)

colmeans <- colMeans(temp.kraken)
colmeans <- colmeans[order(colmeans, decreasing = TRUE)]
head(colmeans, 10)
# Human_gammaherpesvirus_4               Proteus_virus_Isfahan 
#                        1.317625e-02                        3.712793e-03 
#           Mason-Pfizer_monkey_virus              Vibrio_phage_ValB1MD-2 
#                        5.859275e-04                        1.009296e-04 
#          Papiine_gammaherpesvirus_1 Snyder-Theilen_feline_sarcoma_virus 
#                        9.668625e-06                        2.912730e-06 
#      Human_immunodeficiency_virus_1                  Parvovirus_NIH-CQV 
#                        2.784347e-06                        2.374106e-06 
#         Columbid_alphaherpesvirus_1                Escherichia_virus_G4 
#                        2.345646e-06                        1.947459e-06

colmax <- apply(temp.kraken, 2, max)
colmax <- colmax[order(colmax, decreasing = TRUE)]
head(colmax, 10)
#  Human_gammaherpesvirus_4               Proteus_virus_Isfahan 
#                        5.646169e-02                        1.605952e-02 
#           Mason-Pfizer_monkey_virus              Vibrio_phage_ValB1MD-2 
#                        6.005273e-03                        9.715246e-04 
#          Papiine_gammaherpesvirus_1 Snyder-Theilen_feline_sarcoma_virus 
#                        6.542914e-05                        2.572858e-05 
#         Columbid_alphaherpesvirus_1               Simbu_orthobunyavirus 
#                        1.890369e-05                        1.855736e-05 
#      Human_immunodeficiency_virus_1                  Parvovirus_NIH-CQV 
#                        1.315189e-05                        1.285914e-05


############################################### Step 1: get fractions of major infectious bacteria and viruses
##########################################
#################################

##### create figures dir
dir.figures <- "/home/cluster/abalay/scratch/Kraken2_Figures"
if (!dir.exists(dir.figures)){
    dir.create(dir.figures, recursive = TRUE)
}
##################### Density of proportion of common infectious bacteria species
#### List of selected infectious bacteria species
list.inf.bacteria <- c("Bacillus_cereus", "Salmonella_enterica", "Salmonella_bongori", "Escherichia_coli", "Cutibacterium_acnes", 
    "Staphylococcus_aureus", "Staphylococcus_epidermidis", "Yersinia_pestis", "Ralstonia_insidiosa", 
    "Ralstonia_pickettii", "Stenotrophomonas_maltophilia", "Serratia_marcescens", "Acinetobacter_lwoffii", 
    "Acinetobacter_baumannii", "Acinetobacter_sp._NEB149", "Acinetobacter_johnsonii", "Pseudomonas_aeruginosa",
    "Klebsiella_pneumoniae", "Helicobacter_pylori", "Burkholderia_pseudomallei")

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
list.inf.virus <- c("Human_gammaherpesvirus_4", "Human_gammaherpesvirus_8", "Human_immunodeficiency_virus_1")

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


############################################### Step 2: estimate load of common infectious bacteria and virus 
###################################
##########################
################ Download the genomes

##### create dir for genomes
dir.genomes <- "/home/cluster/abalay/scratch/Pipelines_and_Files/Files_for_scripts/Pathogen_genomes"
if (!dir.exists(dir.genomes)){
    dir.create(dir.genomes, recursive = TRUE)
}

##### create dir for cds sequences of genomes
dir.cds.genomes <- "/home/cluster/abalay//scratch/Pipelines_and_Files/Files_for_scripts/Pathogen_genome_cds"
if (!dir.exists(dir.cds.genomes)){
    dir.create(dir.cds.genomes, recursive = TRUE)
}

##### create dir for annotation files of genomes
dir.anno <- "/home/cluster/abalay//scratch/Pipelines_and_Files/Files_for_scripts/Pathogen_genome.Anno"
if (!dir.exists(dir.anno)){
    dir.create(dir.anno, recursive = TRUE)
}

##### download genomes from ncbi website
full.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/220/285/GCF_002220285.1_ASM222028v1/GCF_002220285.1_ASM222028v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/355/GCF_001558355.2_ASM155835v2/GCF_001558355.2_ASM155835v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/439/255/GCF_000439255.1_ASM43925v1/GCF_000439255.1_ASM43925v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/136/215/GCF_004136215.1_ASM413621v1/GCF_004136215.1_ASM413621v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/094/375/GCF_006094375.1_ASM609437v1/GCF_006094375.1_ASM609437v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/798/225/GCF_003798225.1_ASM379822v1/GCF_003798225.1_ASM379822v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/885/GCF_002208885.2_ASM220888v2/GCF_002208885.2_ASM220888v2_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/516/165/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/343/495/GCF_019343495.1_ASM1934349v1/GCF_019343495.1_ASM1934349v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/135/235/GCF_002135235.1_ASM213523v1/GCF_002135235.1_ASM213523v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/337/595/GCF_004337595.1_ASM433759v1/GCF_004337595.1_ASM433759v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/050/665/GCF_003050665.1_ASM305066v1/GCF_003050665.1_ASM305066v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/756/125/GCF_000756125.1_ASM75612v1/GCF_000756125.1_ASM75612v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/GCF_000838265.1_ViralProj14158_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/765/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.fna.gz")

cds.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/220/285/GCF_002220285.1_ASM222028v1/GCF_002220285.1_ASM222028v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/355/GCF_001558355.2_ASM155835v2/GCF_001558355.2_ASM155835v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/439/255/GCF_000439255.1_ASM43925v1/GCF_000439255.1_ASM43925v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/136/215/GCF_004136215.1_ASM413621v1/GCF_004136215.1_ASM413621v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/094/375/GCF_006094375.1_ASM609437v1/GCF_006094375.1_ASM609437v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/798/225/GCF_003798225.1_ASM379822v1/GCF_003798225.1_ASM379822v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/885/GCF_002208885.2_ASM220888v2/GCF_002208885.2_ASM220888v2_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/516/165/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/343/495/GCF_019343495.1_ASM1934349v1/GCF_019343495.1_ASM1934349v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/135/235/GCF_002135235.1_ASM213523v1/GCF_002135235.1_ASM213523v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/337/595/GCF_004337595.1_ASM433759v1/GCF_004337595.1_ASM433759v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/050/665/GCF_003050665.1_ASM305066v1/GCF_003050665.1_ASM305066v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/756/125/GCF_000756125.1_ASM75612v1/GCF_000756125.1_ASM75612v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/GCF_000838265.1_ViralProj14158_cds_from_genomic.fna.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/765/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_cds_from_genomic.fna.gz")

anno.genome.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/220/285/GCF_002220285.1_ASM222028v1/GCF_002220285.1_ASM222028v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/355/GCF_001558355.2_ASM155835v2/GCF_001558355.2_ASM155835v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/439/255/GCF_000439255.1_ASM43925v1/GCF_000439255.1_ASM43925v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/136/215/GCF_004136215.1_ASM413621v1/GCF_004136215.1_ASM413621v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/094/375/GCF_006094375.1_ASM609437v1/GCF_006094375.1_ASM609437v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/798/225/GCF_003798225.1_ASM379822v1/GCF_003798225.1_ASM379822v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/885/GCF_002208885.2_ASM220888v2/GCF_002208885.2_ASM220888v2_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/516/165/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/343/495/GCF_019343495.1_ASM1934349v1/GCF_019343495.1_ASM1934349v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/135/235/GCF_002135235.1_ASM213523v1/GCF_002135235.1_ASM213523v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/337/595/GCF_004337595.1_ASM433759v1/GCF_004337595.1_ASM433759v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/050/665/GCF_003050665.1_ASM305066v1/GCF_003050665.1_ASM305066v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/756/125/GCF_000756125.1_ASM75612v1/GCF_000756125.1_ASM75612v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/GCF_000838265.1_ViralProj14158_genomic.gff.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/765/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.gff.gz")

##### check if all species have full genome, cds and annotation files
length(full.genome.urls)
# [1] 23
length(cds.genome.urls)
# [1] 23
length(anno.genome.urls)
# [1] 23

##### Upload files
for (i in full.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    system(paste("wget ", i, " -O ", file.path(dir.genomes, filename), sep=""))
}

for (i in cds.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    system(paste("wget ", i, " -O ", file.path(dir.cds.genomes, filename), sep=""))
}

for (i in anno.genome.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    system(paste("wget ", i, " -O ", file.path(dir.anno, filename), sep=""))
}

##### check if files are uploaded
length(list.files(dir.genomes))
# [1] 23
length(list.files(dir.cds.genomes))
# [1] 23
length(list.files(dir.anno))
# [1] 23

gmap_build <- "/home/cluster/abalay/scratch/Pipelines_and_Files/Packages/gmap-2020-06-01/util/gmap_build"
gmap.dir <- "/home/cluster/abalay/scratch/Pipelines_and_Files/Packages/gmap-2020-06-01"
nthreads <- 4

iit_store <- "/home/cluster/abalay/scratch/Pipelines_and_Files/Packages/gmap-2020-06-01/src/iit_store"
### Note: if GSNAP was not originally installed in the environment run:
# system(paste("cp ", gmap.dir, "/util/fa_coords ", gmap.dir, "/src", sep=""))
# system(paste("cp ", gmap.dir, "/util/md_coords ", gmap.dir, "/src", sep=""))
# system(paste("cp ", gmap.dir, "/util/gmap_process ", gmap.dir, "/src", sep=""))
############ Create genome indices, annotation file and splicesites file for each species
for (i in 1:length(list.files(dir.genomes))){
    species.id <- head(strsplit(list.files(dir.genomes)[i], "[.]")[[1]],1)
    file.genome <- file.path(dir.genomes, list.files(dir.genomes)[i])
    file.anno <- file.path(dir.anno, list.files(dir.anno)[i])
    chr.names <- c()
    for (j in system(paste("zcat", file.genome, "| grep '>'", sep=" "), intern = TRUE)){
        chr.names <- append(chr.names, strsplit(strsplit(j, " ")[[1]][1], ">")[[1]][2])
    }
    chr.names <- paste(unlist(chr.names), collapse=",")
    system(paste(gmap_build, " -d ", species.id, " -g ", file.genome, " -D ", gmap.dir, " -t ", nthreads, " -B ", gmap.dir, "/src", 
    " --circular ", chr.names, sep = ""))
    system(paste("zcat ", file.anno, "| ", iit_store, " -G -l -F -o ", species.id, sep = ""))
}


############# Align samples against genomes
##### create dir for pathogen detection
dir.patho.detect <- "/home/cluster/abalay/scratch/moved_from_home/Pathogen_detection"
if (!dir.exists(dir.patho.detect)){
    dir.create(dir.patho.detect, recursive = TRUE)
}
#### create list with all FASTQ files unmapped to human genome
fastq.vector <- c()
for (j in system(paste("find /home/cluster/abalay/scratch/ -name '*_unmapped.fastq*'", sep=" "), intern = TRUE)){
    fastq.vector <- append(fastq.vector, j)
}

main.file <- file(paste(dir.patho.detect , "/scheduler_pathogen_gsnap.sh", sep=""))
for (i in 1:length(list.files(dir.genomes))){
    species.id <- head(strsplit(list.files(dir.genomes)[i], "[.]")[[1]],1)
    file.genome <- file.path(dir.genomes, list.files(dir.genomes)[i])
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
                "export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files\n",
                "export samtools=$script_dir/Packages/samtools-1.10/samtools\n",
                "export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap\n",
                "export NTHREADS=4\n",
                "$gsnap --gunzip -d ", species.id, " -D $script_dir/Packages/gmap-2020-06-01  -t $NTHREADS  -n 10 -g ", file.anno,
                " -A sam ", file.r1, " ", file.r2, " |$samtools view -hb|  $samtools sort -@ $NTHREADS > ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam\n", 
                "$samtools index ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam.bai\n",
                "$samtools stats ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
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
                "export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files\n",
                "export samtools=$script_dir/Packages/samtools-1.10/samtools\n",
                "export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap\n",
                "export NTHREADS=4\n",
                "$gsnap --gunzip -d ", species.id, " -D $script_dir/Packages/gmap-2020-06-01  -t $NTHREADS  -n 10 -g ", file.anno,
                " -A sam ", file, " |$samtools view -hb|  $samtools sort -@ $NTHREADS > ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam\n", 
                "$samtools index ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
                sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam.bai\n",
                "$samtools stats ", sub.sample.dir, "/", strsplit(basename(j), ".fastq")[[1]][1], ".bam > ",
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
for (j in system(paste("find ", dir.patho.detect, "/GCF_004136215 -name '*.bam'", sep=""), intern = TRUE)){
    bam.vector <- append(bam.vector, j)
}

for (i in 1:length(list.files(dir.genomes))){
    species.id <- head(strsplit(list.files(dir.genomes)[i], "[.]")[[1]],1)
    file.genome <- file.path(dir.genomes, list.files(dir.genomes)[i])
    file.cds.genome <- file.path(dir.cds.genomes, list.files(dir.cds.genomes)[i])
    sub.species.dir <- paste(dir.patho.detect, "/", species.id, sep="")
    for (j in bam.vector){
        sub.sample.dir <- paste(sub.species.dir, "/", strsplit(basename(j), ".bam")[[1]][1], sep="")
        file <- file.path(sub.sample.dir, basename(j))
        cat("#!/bin/bash\n", "#SBATCH --time=8:00:00\n", "#SBATCH --cpus-per-task=8\n", "#SBATCH --mem-per-cpu=16G\n", "#SBATCH --job-name=", file.path(sub.sample.dir, "/Scripts/"),
        strsplit(basename(j), ".bam")[[1]][1], ".", species.id, "_bedtools.sh", " -o ", file.path(sub.sample.dir, "/Scripts/"), 
                strsplit(basename(j), ".bam")[[1]][1],".", species.id, 
                "_bedtools.out -e ", file.path(sub.sample.dir, "/Scripts/"), strsplit(basename(j), ".bam")[[1]][1], ".", species.id,
                "_bedtools.err\n", 
                "export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files/Packages\n",
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
sammirs.rna_seq.stats.file <- "/home/cluster/abalay/scratch/Sammirs_RNA_Seq/multiqc_reports/cell_lines_fastqc_step3/multiqc_data/multiqc_samtools_stats.txt"
cell_lines.stats.file <- "/home/cluster/abalay/scratch/Cell_Lines/multiqc_reports/cell_lines_fastqc_step3/multiqc_data/multiqc_samtools_stats.txt"
dbgap_rna_seq.stats.file <- "/home/cluster/abalay/scratch/Dbgap_RNA_Seq_DLBCL/multiqc_data/multiqc_samtools_stats.txt"

homo.stats.df <- data.frame(rbind(read.table(sammirs.rna_seq.stats.file, sep="\t", header=TRUE), 
    read.table(cell_lines.stats.file, sep="\t", header=TRUE), 
    read.table(dbgap_rna_seq.stats.file, sep="\t", header=TRUE)))


### to estimate cds coverage upload features tables
### create dir with features table
dir.feature.table <- "/home/cluster/abalay/scratch/Pipelines_and_Files/Files_for_scripts/Pathogen_feature_tables"
if (!dir.exists(dir.feature.table)){
    dir.create(dir.feature.table, recursive = TRUE)
}

features.tables.urls <- c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/220/285/GCF_002220285.1_ASM222028v1/GCF_002220285.1_ASM222028v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/558/355/GCF_001558355.2_ASM155835v2/GCF_001558355.2_ASM155835v2_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/439/255/GCF_000439255.1_ASM43925v1/GCF_000439255.1_ASM43925v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/136/215/GCF_004136215.1_ASM413621v1/GCF_004136215.1_ASM413621v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/094/375/GCF_006094375.1_ASM609437v1/GCF_006094375.1_ASM609437v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/798/225/GCF_003798225.1_ASM379822v1/GCF_003798225.1_ASM379822v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/653/935/GCF_001653935.1_ASM165393v1/GCF_001653935.1_ASM165393v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/466/415/GCF_016466415.2_ASM1646641v2/GCF_016466415.2_ASM1646641v2_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/208/885/GCF_002208885.2_ASM220888v2/GCF_002208885.2_ASM220888v2_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/516/165/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/343/495/GCF_019343495.1_ASM1934349v1/GCF_019343495.1_ASM1934349v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/135/235/GCF_002135235.1_ASM213523v1/GCF_002135235.1_ASM213523v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/337/595/GCF_004337595.1_ASM433759v1/GCF_004337595.1_ASM433759v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/045/595/GCF_011045595.1_ASM1104559v1/GCF_011045595.1_ASM1104559v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/050/665/GCF_003050665.1_ASM305066v1/GCF_003050665.1_ASM305066v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/756/125/GCF_000756125.1_ASM75612v1/GCF_000756125.1_ASM75612v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/GCF_000838265.1_ViralProj14158_feature_table.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/765/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_feature_table.txt.gz")

##### Upload files
for (i in features.tables.urls){
    filename <- tail(strsplit(i, "/")[[1]],1)
    system(paste("wget ", i, " -O ", file.path(dir.feature.table, filename), sep=""))
}

### to add kraken2 fractions we need to match NCBI ID to biological pathogen names
pathogen.names <- c("Bacillus_cereus", "Salmonella_enterica", "Salmonella_bongori", "Escherichia_coli", "Cutibacterium_acnes", 
    "Staphylococcus_aureus", "Staphylococcus_epidermidis", "Yersinia_pestis", "Ralstonia_insidiosa", "Ralstonia_pickettii", 
    "Stenotrophomonas_maltophilia", "Serratia_marcescens", "Acinetobacter_lwoffii", 
    "Acinetobacter_baumannii", "Acinetobacter_sp._NEB149", "Acinetobacter_johnsonii", "Pseudomonas_aeruginosa",
    "Klebsiella_pneumoniae", "Helicobacter_pylori", "Burkholderia_pseudomallei", 
    "Human_gammaherpesvirus_4", "Human_gammaherpesvirus_8", "Human_immunodeficiency_virus_1")

##### upload files with kraken2 bacteria and virus fractions
common.bacteria.file <- read.table(file.path(dir.figures, "Common.bacteria.kraken.dataset.txt"), sep="\t", header=TRUE)
common.virus.file <- read.table(file.path(dir.figures, "Common.virus.kraken.dataset.txt"), sep="\t", header=TRUE)
common.pathogen.file <- cbind(common.bacteria.file, common.virus.file)
for (i in c("I05.4080", "I04.3825", "I05.2671", "I04.1731")){
    rownames(common.pathogen.file) <- gsub(i, paste(strsplit(i, "\\.")[[1]][1], "-", strsplit(i, "\\.")[[1]][2], sep=""), rownames(common.pathogen.file))
}

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

for (i in 1:length(big.pathogen.stats.list)){
    sample.ids <- unlist(big.pathogen.stats.list[[i]]["Sample.ID"]) %in% a
    big.pathogen.stats.list[[i]] <- big.pathogen.stats.list[[i]][!sample.ids, ]
}

actual.pathogen.names <- c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Human gammaherpesvirus 8 (KSHV)", 
    "Human immunodeficiency virus 1 (HIV1)", "Salmonella enterica", "Ralstonia insidiosa", "Acinetobacter sp. NEB149", 
    "Stenotrophomonas maltophilia", "Bacillus cereus", 
    "Human gammaherpesvirus 4 (EBV)", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii")

species <- c("Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Virus", "Virus", "Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Virus", "Bacteria", "Bacteria", "Bacteria",
    "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria")

### subset error rates and average mapping quality for all pathogens
error.rate <- c(); ave.map.qual <- c(); pathogen.anno <- c(); species.anno <- c()
for (i in 1:length(big.pathogen.stats.list)){
    error.rate <- append(error.rate, big.pathogen.stats.list[[i]]$Error.rate)
    ave.map.qual <- append(ave.map.qual, big.pathogen.stats.list[[i]]$Ave.quality)
    pathogen.anno <- append(pathogen.anno, rep(actual.pathogen.names[i], length(big.pathogen.stats.list[[i]]$Error.rate)))
    species.anno <- append(species.anno, rep(species[i], length(big.pathogen.stats.list[[i]]$Error.rate)))
}

df <- data.frame(Error.rate = error.rate, Ave.Map.Qual = ave.map.qual, Pathogen = pathogen.anno, Species = species.anno)

axis.x.colors <- c("red", "blue"); names(axis.x.colors) <- unique(df$Species); df$color <- as.vector(axis.x.colors[df$Species])

df$Pathogen <- factor(x = df$Pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Human gammaherpesvirus 8 (KSHV)", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

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
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.05, color = "brown4", linetype="dashed") +
    xlab("Error rate") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("Error rate across pathogens")

png(paste(dir.final.figures, "/Error.rate.dist.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
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
# [1] 168  23
kraken2.fraction.df <- data.frame(do.call(cbind, kraken2.fraction.list), row.names = sample.ids)
dim(kraken2.fraction.df)
# [1] 168  23
wgcov.df <- data.frame(do.call(cbind, wgcov.list), row.names = sample.ids)
dim(wgcov.df)
# [1] 168  23
cdscov.df <- data.frame(do.call(cbind, cdscov.list), row.names = sample.ids)
dim(cdscov.df)
# [1] 168  23
 
### Heatmap: pathogen load of samples across pathogens
colnames(pathogen.load.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "firebrick1", "Virus" = "blue"))

h1_top <- HeatmapAnnotation(`Pathogen type` = species, col = species.colors)

viridis.colors <- colorRamp2(c(0, 1, 5), c("purple", "grey20", "yellow")) ### set color scheme

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

png(paste(dir.final.figures, "/Pathogen.load.heatmap.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
h1
dev.off()

### Heatmap: kraken2 fraction of samples across pathogens
colnames(kraken2.fraction.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "firebrick1", "Virus" = "blue"))

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

png(paste(dir.final.figures, "/Kraken2.fraction.heatmap.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
h2
dev.off()

### Heatmap: whole genome coverage of samples across pathogens
colnames(wgcov.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "firebrick1", "Virus" = "blue"))

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

png(paste(dir.final.figures, "/WG.Cov.heatmap.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
h3
dev.off()


### Heatmap: cds coverage of samples across pathogens
colnames(cdscov.df) <- actual.pathogen.names ### set proper pathogen names

species.colors <- list(`Pathogen type` = c("Bacteria" = "firebrick1", "Virus" = "blue"))

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

png(paste(dir.final.figures, "/CDS.Cov.heatmap.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
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

axis.x.colors <- c("red", "blue"); names(axis.x.colors) <- unique(df$Species); df$color <- as.vector(axis.x.colors[df$Species])

df$Pathogen <- factor(x = df$Pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Human gammaherpesvirus 8 (KSHV)", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

### plot wg cov distribution for all samples
wg.cov.p <- ggplot(df, aes(x = WG.cov, y = Pathogen, fill = Pathogen, color=Species)) + 
    geom_violin(trim = FALSE) +
    scale_fill_manual(values=magma(length(unique(df$Pathogen)))) +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.1, color = "brown4", linetype="dashed") +
    xlab("WG Coverage") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("WG Coverage across pathogens")

png(paste(dir.final.figures, "/WG.Cov.dist.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
wg.cov.p + guides(fill="none")
dev.off()

### plot cds cov distribution for all samples
cds.cov.p <- ggplot(df, aes(x = CDS.cov, y = Pathogen, fill = Pathogen, color=Species)) + 
    geom_violin(trim = FALSE) +
    scale_fill_manual(values=magma(length(unique(df$Pathogen)))) +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() + 
    geom_vline(xintercept = 0.3, color = "brown4", linetype="dashed") +
    xlab("CDS Coverage") + ylab("Pathogen") + labs(colour="Pathogen type") +
    ggtitle("CDS Coverage across pathogens")

png(paste(dir.final.figures, "/CDS.Cov.dist.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
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

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Human gammaherpesvirus 8 (KSHV)", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.36, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ Kraken2 fraction")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.kraken2_fraction.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

#### Pathogen load vs whole genome coverage
temp.df <- data.frame(R = as.vector(pl.wg.rho), p.val = as.vector(pl.wg.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Listeria monocytogenes", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica", "Staphylococcus haemolyticus", "Ralstonia mannitolilytica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Bacteroides thetaiotaomicron", "Human gammaherpesvirus 8 (KSHV)", "Human mastadenovirus C", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ WG Coverage")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.wg_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Pathogen load vs cds coverage
temp.df <- data.frame(R = as.vector(pl.cds.rho), p.val = as.vector(pl.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Listeria monocytogenes", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica", "Staphylococcus haemolyticus", "Ralstonia mannitolilytica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Bacteroides thetaiotaomicron", "Human gammaherpesvirus 8 (KSHV)", "Human mastadenovirus C", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("Pathogen load ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.pathogen_load.vs.cds_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Kraken2 fraction vs whole genome coverage
temp.df <- data.frame(R = as.vector(kf.wg.rho), p.val = as.vector(kf.wg.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Listeria monocytogenes", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica", "Staphylococcus haemolyticus", "Ralstonia mannitolilytica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Bacteroides thetaiotaomicron", "Human gammaherpesvirus 8 (KSHV)", "Human mastadenovirus C", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("Kraken2 fraction ~ WG Coverage")

png(paste(dir.final.figures, "/Line_plot.kraken2_fraction.vs.wg_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

#### Kraken2 fraction vs cds coverage
temp.df <- data.frame(R = as.vector(kf.cds.rho), p.val = as.vector(kf.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Listeria monocytogenes", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica", "Staphylococcus haemolyticus", "Ralstonia mannitolilytica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Bacteroides thetaiotaomicron", "Human gammaherpesvirus 8 (KSHV)", "Human mastadenovirus C", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("Kraken2 fraction ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.kraken2_fraction.vs.cds_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


#### Whole genome coverage vs cds coverage
temp.df <- data.frame(R = as.vector(wg.cds.rho), p.val = as.vector(wg.cds.pval), pathogen = actual.pathogen.names, 
    species = species)

temp.df$pathogen <- factor(x = temp.df$pathogen, levels = c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Listeria monocytogenes", 
    "Salmonella bongori", "Burkholderia pseudomallei", "Salmonella enterica", "Staphylococcus haemolyticus", "Ralstonia mannitolilytica",
    "Ralstonia insidiosa", "Acinetobacter sp. NEB149", "Stenotrophomonas maltophilia", "Bacillus cereus", "Helicobacter pylori", "Serratia marcescens", "Yersinia pestis",
    "Cutibacterium acnes", "Acinetobacter johnsonii", "Staphylococcus epidermidis", "Acinetobacter baumannii",
    "Klebsiella pneumoniae", "Ralstonia pickettii", "Acinetobacter lwoffii", "Bacteroides thetaiotaomicron", "Human gammaherpesvirus 8 (KSHV)", "Human mastadenovirus C", 
    "Human immunodeficiency virus 1 (HIV1)", "Human gammaherpesvirus 4 (EBV)"))

p <- ggplot(temp.df, aes(x = pathogen, y = R, color=species, group = 1)) +
    geom_line() +
    geom_point() +
    ylim(c(-0.3, 1)) +
    geom_hline(yintercept = 0.5, color = "red", linetype="dashed") +
    labs(x = "Pathogen", y = "Spearman correlation coefficient") +
    scale_color_manual(values = c("chocolate4", "blue")) + 
    theme_classic() +
    ggtitle("WG Coverage ~ CDS Coverage")

png(paste(dir.final.figures, "/Line_plot.kraken2_fraction.vs.cds_coverage.corr.all.pathogens.png", sep=""), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()


### Kraken2 fractions don't correlate well with rest of variables. Why?
### For pathogen with low correlation exclude kraken2 fraction

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

# Escherichia coli 166 2 
# Pseudomonas aeruginosa 166 2 
# Staphylococcus aureus 167 1 
# Salmonella bongori 167 1 
# Burkholderia pseudomallei 167 1 
# Human gammaherpesvirus 8 (KSHV) 133 35 
# Human immunodeficiency virus 1 (HIV1) 159 9 
# Salmonella enterica 167 1 
# Ralstonia insidiosa 162 6 
# Acinetobacter sp. NEB149 132 36 
# Stenotrophomonas maltophilia 163 5 
# Bacillus cereus 167 1 
# Human gammaherpesvirus 4 (EBV) 67 101 
# Helicobacter pylori 167 1 
# Serratia marcescens 166 2 
# Yersinia pestis 167 1 
# Cutibacterium acnes 162 6 
# Acinetobacter johnsonii 123 45 
# Staphylococcus epidermidis 167 1 
# Acinetobacter baumannii 145 23 
# Klebsiella pneumoniae 166 2 
# Ralstonia pickettii 167 1 
# Acinetobacter lwoffii 114 54 

### save the file with pathogen classification
saveRDS(big.pathogen.stats.list, file.path(dir.patho.detect, "Pathogen.sample.info.RData"))

### load pathogen file
big.pathogen.stats.list <- readRDS(file.path(dir.patho.detect, "Pathogen.sample.info.RData"))

### make barplot of samples infected and uninfected
temp.list <- list()
for (i in 1:length(big.pathogen.stats.list)){
    temp.list[[i]] <- big.pathogen.stats.list[[i]]$Classification_by_cover_only
    names(temp.list)[i] <- names(big.pathogen.stats.list)[i]
}

infection.status <- as.data.frame(do.call(cbind,temp.list))
colnames(infection.status) <- gsub(" ", ".", actual.pathogen.names)
rownames(infection.status) <- big.pathogen.stats.list[[1]]$Sample.ID

### make a table which sample is infected
infection.status <- infection.status %>% mutate(Acinetobacter = ifelse(Acinetobacter.baumannii == 1|
    Acinetobacter.johnsonii == 1| Acinetobacter.sp..NEB149 == 1| Acinetobacter.lwoffii == 1, 1, 0), .before = 2)

infection.status$Acinetobacter.baumannii <- NULL
infection.status$Acinetobacter.sp..NEB149 <- NULL
infection.status$Acinetobacter.lwoffii <- NULL
infection.status$Acinetobacter.johnsonii <- NULL

### save infection status file
save(infection.status, file = file.path(dir.patho.detect, "Infection.status.df.Robj"))

count_list <- list()
for (i in 1:ncol(infection.status)){
    count_list[[i]] <- table(infection.status[,i])
    names(count_list)[i] <- colnames(infection.status)[i]
}

counts <- do.call(rbind, count_list)
colnames(counts) <- c("Uninfected", "Infected")
rownames(counts) <- names(count_list)

counts <- melt(counts)
counts$value[counts$Var2 == "Uninfected"] <- counts$value[counts$Var2 == "Uninfected"] * -1
counts$Var1 <- gsub("\\.", " ", counts$Var1)

p1 <- ggplot(counts, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(position = "identity", stat="identity", colour="black") +
    geom_text(aes(label= abs(value), hjust=ifelse(sign(value)>0, -0.4, -0.1)), vjust=0.8) +
    scale_fill_manual(values = c("firebrick1", "blue")) +
    labs(x = "Pathogen", y = "Count", fill = "Status") +
    theme_classic() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    coord_flip()

png(paste(dir.final.figures, "/Infection.status.barplot.by.species.png", sep=""), res=200, unit="in", height=8, width=11)
p1
dev.off()

### load sample annotation file for huNSG mice samples
anno.file.dir <- "/home/cluster/abalay/scratch/moved_from_home/Pipelines_and_Files/Files_for_scripts/"
anno.file <- read_excel(paste0(anno.file.dir, "EBV-transformed_lymphoblast_cell_lines.xlsx"), 
    sheet = "huNSG mice samples")
anno.file <- anno.file[1:79,]

### load sample annotation file for Sammir's RNA-Seq sample
anno.file.1 <- read_excel(paste0(anno.file.dir, "EBV-transformed_lymphoblast_cell_lines.xlsx"), 
    sheet = "DLBCL ")
anno.file.1 <- anno.file.1[1:18,]

### load sample annotation file for DBGAP RNA-Seq sample
anno.file.2 <- read_excel(paste0(anno.file.dir, "EBV-transformed_lymphoblast_cell_lines.xlsx"), 
    sheet = "DBGAP RNA-SEQ DLBCL")

### get ids of samples with known and with no KSHV infection
hunsg.mice.ids <- paste(as.vector(unlist(anno.file["SRA-ID"]))[1:79], "_1_r1_unmapped_map", sep="")

### check if all ids present in list
all(hunsg.mice.ids %in% as.vector(unlist(big.pathogen.stats.list[['GCF_000838265']]["Sample.ID"])))
# [1] TRUE

### subset df with samle ids of known and no KSHV infection
temp.kshv.df <- big.pathogen.stats.list[['GCF_000838265']]
temp.kshv.df <- temp.kshv.df[as.vector(unlist(temp.kshv.df["Sample.ID"])) %in% hunsg.mice.ids, ]
dim(temp.kshv.df)
# [1] 79 10

### column with ground truth KSHV infection 
kshv.class <- rep(0, nrow(anno.file))
kshv.class[grep("KSHV", as.vector(unlist(anno.file["Infection agent"])))] = 1

anno.file["KSHV infection"] <- kshv.class

### arrange dataframe in order
temp.kshv.df <- temp.kshv.df %>% arrange(Sample.ID)

### check if sequence of samples matches
all(temp.kshv.df["Sample.ID"] == hunsg.mice.ids)
# [1] TRUE

### add ground truth classification for KSHV
temp.kshv.df["Ground Truth"] <- anno.file["KSHV infection"]

### create confusion matrix for kshv classification
cm <- confusionMatrix(factor(as.vector(unlist(temp.kshv.df["Classification_by_cover_only"]))), 
    factor(as.vector(unlist(temp.kshv.df["Ground Truth"]))))

TClass <- factor(c(0, 0, 1, 1))
PClass <- factor(c(0, 1, 0, 1))
cm.table <- as.numeric(cm$table)

df <- data.frame(TClass, PClass, cm.table)

stats.list <- list()
### plot confusion matrix
p1 <- ggplot(df, aes(x = TClass, y = PClass)) + 
    geom_tile(aes(fill = cm.table), colour = "white") +
    geom_richtext(label = sprintf("%1.0f", cm.table), fill = "firebrick2", text.colour = "black",
    vjust = 1, size = 10) +
    # scale_fill_gradient(low = viridis(1), high = viridis(2)) +
    scale_fill_viridis(option = "viridis") +
    theme_classic() +
    scale_x_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    scale_y_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    xlab("True class") + ylab("Predicted class") +
    theme(legend.title = element_blank(), axis.text = element_text(size = 10, angle=90), axis.title = element_text(size = 15), 
    legend.key.size = unit(2, 'cm'), plot.title = element_text(hjust = 0.5)) +
    ggtitle("KSHV")

png(paste(dir.final.figures, "/Confusion.matrix.KSHV.png", sep=""), res=200, unit="in", height=8, width=11)
p1
dev.off()

stats.list[[1]] <- cm



### get ids of samples with known and with no EBV infection
sammirs.sample.ids <- paste(as.vector(unlist(anno.file.1["Sample ID"])), "_unmapped_map", sep="")

### check if all ids present in list
all(sammirs.sample.ids %in% as.vector(unlist(big.pathogen.stats.list[['GCF_002402265']]["Sample.ID"])))
# [1] TRUE

### subset df with samle ids of known and no EBV infection
temp.ebv.df <- big.pathogen.stats.list[['GCF_002402265']]
temp.ebv.df <- temp.ebv.df[as.vector(unlist(temp.ebv.df["Sample.ID"])) %in% c(hunsg.mice.ids, sammirs.sample.ids), ]
dim(temp.ebv.df)
# [1] 97  10

### columns with ground truth EBV infection 
ebv.class <- rep(0, nrow(anno.file))
ebv.class[grep("EBV", as.vector(unlist(anno.file["Infection agent"])))] = 1

anno.file["EBV infection"] <- ebv.class

ebv.class.1 <- rep(0, nrow(anno.file.1))
ebv.class.1[grep("EBV", as.vector(unlist(anno.file.1["Infection agent"])))] = 1

anno.file.1["EBV infection"] <- ebv.class.1

### arrange dataframe in order
ids <- c(hunsg.mice.ids, sammirs.sample.ids)
gr.truth <- c(as.vector(unlist(anno.file["EBV infection"])), as.vector(unlist(anno.file.1["EBV infection"])))
temp.ebv.df <- temp.ebv.df[match(ids, as.vector(unlist(temp.ebv.df["Sample.ID"]))), ]
all(temp.ebv.df["Sample.ID"] == ids)
# [1] TRUE
temp.ebv.df["Ground Truth"] <- gr.truth

### crop HIV infected samples due to unknown EBV status
temp.ebv.df <- temp.ebv.df[1:89,]


### create confusion matrix for ebv classification
cm <- confusionMatrix(factor(as.vector(unlist(temp.ebv.df["Classification_by_cover_only"]))), 
    factor(as.vector(unlist(temp.ebv.df["Ground Truth"]))))

TClass <- factor(c(0, 0, 1, 1))
PClass <- factor(c(0, 1, 0, 1))
cm.table <- as.numeric(cm$table)

df <- data.frame(TClass, PClass, cm.table)

### plot confusion matrix
p2 <- ggplot(df, aes(x = TClass, y = PClass)) + 
    geom_tile(aes(fill = cm.table), colour = "white") +
    geom_richtext(label = sprintf("%1.0f", cm.table), fill = "firebrick2", text.colour = "black",
    vjust = 1, size = 10) +
    # scale_fill_gradient(low = viridis(1), high = viridis(2)) +
    scale_fill_viridis(option = "viridis") +
    theme_classic() +
    scale_x_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    scale_y_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    xlab("True class") + ylab("Predicted class") +
    theme(legend.title = element_blank(), axis.text = element_text(size = 10, angle=90), axis.title = element_text(size = 15), 
    legend.key.size = unit(2, 'cm'), plot.title = element_text(hjust = 0.5)) +
    ggtitle("EBV")

png(paste(dir.final.figures, "/Confusion.matrix.EBV.png", sep=""), res=200, unit="in", height=8, width=11)
p2
dev.off()

stats.list[[2]] <- cm


### get ids of samples with known and with no HIV infection
dbgap.sample.ids <- paste(as.vector(unlist(anno.file.2["SRA-ID"])), "_dbGaP-28434_1_r1_unmapped_map", sep="")
dbgap.sample.ids[dbgap.sample.ids == dbgap.sample.ids[2]] = paste(as.vector(unlist(anno.file.2["SRA-ID"]))[2], "_dbGaP-28434_unmapped_map", sep="")
dbgap.sample.ids[dbgap.sample.ids == dbgap.sample.ids[3]] = paste(as.vector(unlist(anno.file.2["SRA-ID"]))[3], "_dbGaP-28434_unmapped_map", sep="")
dbgap.sample.ids[dbgap.sample.ids == dbgap.sample.ids[4]] = paste(as.vector(unlist(anno.file.2["SRA-ID"]))[4], "_dbGaP-28434_unmapped_map", sep="")
dbgap.sample.ids[dbgap.sample.ids == dbgap.sample.ids[5]] = paste(as.vector(unlist(anno.file.2["SRA-ID"]))[5], "_dbGaP-28434_unmapped_map", sep="")
dbgap.sample.ids[dbgap.sample.ids == dbgap.sample.ids[10]] = paste(as.vector(unlist(anno.file.2["SRA-ID"]))[10], "_dbGaP-28434_unmapped_map", sep="")

### check if all ids present in list
all(dbgap.sample.ids %in% as.vector(unlist(big.pathogen.stats.list[['GCF_000864765']]["Sample.ID"])))
# [1] TRUE

### subset df with samle ids of known and no HIV infection
temp.hiv.df <- big.pathogen.stats.list[['GCF_000864765']]
temp.hiv.df <- temp.hiv.df[as.vector(unlist(temp.hiv.df["Sample.ID"])) %in% c(sammirs.sample.ids, dbgap.sample.ids), ]
dim(temp.hiv.df)
# [1] 89 10

### columns with ground truth HIV infection 
hiv.class <- rep(0, nrow(anno.file.1))
hiv.class[grep("HIV", as.vector(unlist(anno.file.1["Infection agent"])))] = 1

anno.file.1["HIV infection"] <- hiv.class

hiv.class.1 <- rep(0, nrow(anno.file.2))
hiv.class.1[grep("HIV", as.vector(unlist(anno.file.2["Infection agent"])))] = 1

anno.file.2["HIV infection"] <- hiv.class.1


### arrange dataframe in order
ids <- c(sammirs.sample.ids, dbgap.sample.ids)
gr.truth <- c(as.vector(unlist(anno.file.1["HIV infection"])), as.vector(unlist(anno.file.2["HIV infection"])))
temp.hiv.df <- temp.hiv.df[match(ids, as.vector(unlist(temp.hiv.df["Sample.ID"]))), ]
all(temp.hiv.df["Sample.ID"] == ids)
# [1] TRUE
temp.hiv.df["Ground Truth"] <- gr.truth

### create confusion matrix for hiv classification
cm <- confusionMatrix(factor(as.vector(unlist(temp.hiv.df["Classification_by_cover_only"]))), 
    factor(as.vector(unlist(temp.hiv.df["Ground Truth"]))))

TClass <- factor(c(0, 0, 1, 1))
PClass <- factor(c(0, 1, 0, 1))
cm.table <- as.numeric(cm$table)

df <- data.frame(TClass, PClass, cm.table)

### plot confusion matrix
p3 <- ggplot(df, aes(x = TClass, y = PClass)) + 
    geom_tile(aes(fill = cm.table), colour = "white") +
    geom_richtext(label = sprintf("%1.0f", cm.table), fill = "firebrick2", text.colour = "black",
    vjust = 1, size = 10) +
    # scale_fill_gradient(low = viridis(1), high = viridis(2)) +
    scale_fill_viridis(option = "viridis") +
    theme_classic() +
    scale_x_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    scale_y_discrete(breaks=c("0","1"),
        labels=c("Uninfected", "Infected")) +
    xlab("True class") + ylab("Predicted class") +
    theme(legend.title = element_blank(), axis.text = element_text(size = 10, angle=90), axis.title = element_text(size = 15), 
    legend.key.size = unit(2, 'cm'), plot.title = element_text(hjust = 0.5)) +
    ggtitle("HIV")

png(paste(dir.final.figures, "/Confusion.matrix.HIV.png", sep=""), res=200, unit="in", height=8, width=11)
p3
dev.off()

stats.list[[3]] <- cm

names(stats.list)[1:3] <- c("KSHV", "EBV", "HIV")

saveRDS(stats.list, file.path(dir.patho.detect, "Confusion.matrix.list.RData"))

























 



