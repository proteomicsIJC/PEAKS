#### PEAKS TO ProteinGroups
### Get the data and libraries----
library(dplyr)
library(Vennerable)

peaks <- read.csv2("./proteins_peaks.csv", sep = ",", header = T, check.names = F)
#-----

####Get the correct names to construct ProteinIDs----
# Only the accession before the last of the accessions in every row
# Split the column into a list of vector
g_acs <- strsplit(peaks$Accession, "\\|")
length(g_acs)

beforelast <- function(charlist){
  empty <- c()
  final_list <- list()
  
  for (j in 1:length(charlist)){
    
    for (i in 1:length(charlist[j])){
      befla <- charlist[[j]][length(charlist[[j]])-1]
      befla <- c(empty,befla)}

    final_list[[j]] <- befla
    the_only_good <- unlist(final_list)
  }
  
  return(the_only_good)
}

good_names <- beforelast(g_acs)
length(beforelast(g_acs)) == nrow(peaks)
length(unique(peaks$GoodAccession))

# Introduce the data to the dataset
peaks$GoodAccession <- good_names
peaks <- peaks %>%
  relocate(GoodAccession, .before = Accession)

head(peaks[,c("Accession","GoodAccession")])
#-----

### Collapse rows with the same information----
# (when redundant)
# Redundant columns 
# -10lgP Coverage (%), Coverage (%) C1_MTBE, Coverage (%) C2_EtOH, Coverage (%) T1_MTBE,
# Coverage (%) T2_EtOH, Coverage (%) G1_MTBE, Coverage (%) G2_EtOH, 
# "Area C1_MTBE", "Area C2_EtOH", "Area T1_MTBE", "Area T2_EtOH", "Area G1_MTBE", "Area G2_EtOH",
# "#Peptides", "#Unique", 
# "#Spec C1_MTBE","#Spec C2_EtOH","#Spec T1_MTBE","#Spec T2_EtOH","#Spec G1_MTBE","#Spec G2_EtOH"
# "PTM"

peaks$Protein_Group <- as.character(peaks$Protein_Group)

Redundant_peaks <- peaks %>% 
  group_by(Protein_Group) %>%
  #-10lgP
  summarise(ten_lgp = head(ten_lgp,1),
            #Coverages
            Coverage_per = head(Coverage_per,1),
            Coverage_C1_MTBE = head(Coverage_C1_MTBE,1),
            Coverage_C2_EtOH = head(Coverage_C2_EtOH,1),
            Coverage_T1_MTBE = head(Coverage_T1_MTBE,1),
            Coverage_T2_EtOH = head(Coverage_T2_EtOH,1),
            Coverage_G1_MTBE = head(Coverage_G1_MTBE,1),
            Coverage_G2_EtOH = head(Coverage_G2_EtOH,1),
            #Area
            Area_C1_MTBE = head(Area_C1_MTBE,1),
            Area_C2_EtOH = head(Area_C2_EtOH,1),
            Area_T1_MTBE = head(Area_T1_MTBE,1),
            Area_T2_EtOH = head(Area_T2_EtOH,1),
            Area_G1_MTBE = head(Area_G1_MTBE,1),
            Area_G2_EtOH = head(Area_G2_EtOH,1),
            #ns
            n_Peptides = head(n_Peptides,1),
            n_Unique = head(n_Unique,1),
            #specs
            n_Spec_C1_MTBE = head(n_Spec_C1_MTBE,1),
            n_Spec_C2_EtOH = head(n_Spec_C2_EtOH,1),
            n_Spec_T1_MTBE = head(n_Spec_T1_MTBE,1),
            n_Spec_T2_EtOH = head(n_Spec_T2_EtOH,1),
            n_Spec_G1_MTBE = head(n_Spec_G1_MTBE,1),
            n_Spec_G2_EtOH = head(n_Spec_G2_EtOH,1),
            #PTM
            PTM = head(PTM,1))
#----

### Concatenate rows with information to keep----
# (when non-redundant)
# Non-Redundant
# "Protein ID", "GoodAccession" "Accession", "Avg. Mass","Description"
to_concatenate <- peaks[,c("Protein_Group","Protein_ID","GoodAccession","Accession","Avg_Mass","Description")]

Concatenate_peaks <- to_concatenate %>%
  group_by(Protein_Group) %>%
  mutate(Protein_ID = paste(Protein_ID, collapse = ";"),
         GoodAccession = paste(GoodAccession, collapse = ";"),
         Accession = paste(Accession, collapse = ";"),
         Avg_Mass = paste(Avg_Mass, collapse = ";"),
         Description = paste(Description, collapse = ";"))

Concatenate_peaks <- Concatenate_peaks %>%
  distinct(Protein_Group, Protein_ID, GoodAccession, Accession, Avg_Mass, Description)
#----

### Change Peaks to ProteinGroup----
Peaks_PG <- merge(x = Concatenate_peaks, y = Redundant_peaks,
                   by = "Protein_Group")
#----