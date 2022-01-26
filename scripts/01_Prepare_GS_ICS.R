library(CytoML) 
library(flowCore) 
library(flowWorkspace) 
library(ggcyto)
library(here)
library(tidyverse)
library(readxl)
library(cytoqc) # Use to fix inconsistent channels

## Create output directories if needed ## 
if(!dir.exists(here::here("out/QC"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC")))
  dir.create(here::here("out/QC"), recursive = T)
}

if(!dir.exists(here::here("out/QC/Counts"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/Counts")))
  dir.create(here::here("out/QC/Counts"), recursive = T)
}

if(!dir.exists(here::here("out/QC/DMSO_Signal"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/QC/DMSO_Signal")))
  dir.create(here::here("out/QC/DMSO_Signal"), recursive = T)
}

if(!dir.exists(here::here("out/GatingSets"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/GatingSets")))
  dir.create(here::here("out/GatingSets"), recursive = T)
}

## Load workspaces and prepare GatingSets ##
xml_path_b1 <- here::here("data/20211018_HAARVIVAC_B1V1_JP.xml")
xml_path_b4 <- here::here("data/20211014_HAARVIVAC_B4V3_JP.xml")
xml_path_b5 <- here::here("data/20211013_HAARVIVAC_B5V1_JP.xml")
fcs_subfolder_b1 <- here::here("data/20210212_HAARVIVAC_FCS_B1/")
fcs_subfolder_b4 <- here::here("data/20210902_HAARVIVAC_FCS_B4/")
fcs_subfolder_b5 <- here::here("data/20211006_HAARVIVAC_FCS_B5/")
ws_b1 <- open_flowjo_xml(xml_path_b1)
ws_b4 <- open_flowjo_xml(xml_path_b4)
ws_b5 <- open_flowjo_xml(xml_path_b5)

# Batch 1 has different keywords than Batches 4 and 5
names(fj_ws_get_keywords(ws_b1, 117)) 
keywords2import_1 <- c("EXPERIMENT NAME",
                       "$DATE",
                       "SAMPLE ID",
                       "PATIENT ID",
                       "STIM",
                       "WELL ID",
                       "PLATE NAME") 

names(fj_ws_get_keywords(ws_b4, 117))
keywords2import_2 <- c("EXPERIMENT NAME",
                     "$DATE",
                     "SAMPLE ID",
                     "PATIENT ID",
                     "Stim",
                     "WELL ID",
                     "PLATE NAME",
                     "Timepoint")

sampleGroup <- "Samples"

# Passing Batch 1 to the parser flowjo_to_gatingset() gives the following... 
# Error in get_cytoset(x@pointer) : Found channel inconsistency across samples.
# 136668.fcs_148480 has the channel 'V655-A' that is not found in other samples!

# Need to fix inconsistent channels across FCS files
# Load Batch 1 FCS files
rawfiles <- list.files(fcs_subfolder_b1, ".fcs", full.names = TRUE)
cqc_data <- cqc_load_fcs(rawfiles)

# Check channels
check_res <- cqc_check(cqc_data, type = "channel")
check_res

# Pick reference to match
match_res <- cqc_match(check_res, ref = 2)
match_res

# Fix the problem
cqc_fix(match_res)

# Check channels again to make sure data is standardized
cqc_check_channel(cqc_data)

# Coerce data to cytoset 
cs <- cytoset(cqc_data)

# Need to import comp matrix to pass through flow_jo_gatingset, otherwise...
# Error in (function (ws, group_id, subset, execute, path, cytoset, backend_dir,  :
# compensation parameter 'B710-A' not found in cytoframe parameters!

# I exported the Batch 1 comp matrix from FlowJo, saved as a .csv file, and 
# removed the original comp matrix from the .csv file 
compmat <- read.csv(file = "data/compmat_B1V1_JP.csv", header = TRUE, skip = 2, check.names = FALSE)
comp <- compensation(compmat, compensationId = "comp1")

# Finish loading GatingSets
gs_b1 <- flowjo_to_gatingset(ws_b1,                                    
                             name=sampleGroup, 
                             cytoset=cs,
                             keywords=keywords2import_1,
                             compensation = comp,
                             extend_val=-10000)

gs_b4 <- flowjo_to_gatingset(ws_b4,                                    
                             name=sampleGroup, 
                             keywords=keywords2import_2,
                             path=fcs_subfolder_b4, 
                             extend_val=-10000)

gs_b5 <- flowjo_to_gatingset(ws_b5,                                    
                             name=sampleGroup, 
                             keywords=keywords2import_2,
                             path=fcs_subfolder_b5, 
                             extend_val=-10000)

# Check the gating trees
pop_lists <- lapply(gs_b1, gh_get_pop_paths)
unique(pop_lists)

pop_lists <- lapply(gs_b4, gh_get_pop_paths)
unique(pop_lists)

pop_lists <- lapply(gs_b5, gh_get_pop_paths)
unique(pop_lists)

# Remove channels from flow data that are not used by gates
gs_b1 <- gs_remove_redundant_channels(gs_b1) # drop SSC-H

gs_b4 <- gs_remove_redundant_channels(gs_b4) # drop SSC-H, V655-A, V570-A

gs_b5 <- gs_remove_redundant_channels(gs_b5) # drop SSC-H, V655-A, V570-A

# Add names to all channels
dput(unname(pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,2]))
markernames_b1 <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a", "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19", "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames_b1) <- pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,1]
markernames(gs_b1) <- markernames_b1
pData(parameters(gh_pop_get_data(gs_b1[[1]])))[,c(1,2)]

dput(unname(pData(parameters(gh_pop_get_data(gs_b4[[1]])))[,2]))
markernames_b4 <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a", "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19", "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames_b4) <- pData(parameters(gh_pop_get_data(gs_b4[[1]])))[,1]
markernames(gs_b4) <- markernames_b4
pData(parameters(gh_pop_get_data(gs_b4[[1]])))[,c(1,2)]

dput(unname(pData(parameters(gh_pop_get_data(gs_b5[[1]])))[,2]))
markernames_b5 <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a", "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19", "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames_b5) <- pData(parameters(gh_pop_get_data(gs_b5[[1]])))[,1]
markernames(gs_b5) <- markernames_b5
pData(parameters(gh_pop_get_data(gs_b5[[1]])))[,c(1,2)]

# Reformat Batch 1 pData to match Batches 4 and 5
# Rename and reformat "Stim" column
pData(gs_b1) <- rename(pData(gs_b1), Stim = STIM)
i1 <- grepl("Spike 1", pData(gs_b1)$Stim) | grepl("Spike1", pData(gs_b1)$Stim)
i2 <- grepl("Spike 2", pData(gs_b1)$Stim) 
pData(gs_b1)$Stim[i1] <- "S1"
pData(gs_b1)$Stim[i2] <- "S2"

# Rename "EXPERIMENT NAME" column
pData(gs_b1)[pData(gs_b1) == "20210212_HAARVIVAC_ICS"] <- "20210212 HAARVIVAC Batch 1"

# Reformat "PATIENT ID" column
pData(gs_b1)$'PATIENT ID' <-  str_replace(pData(gs_b1)$'PATIENT ID', "^.*-", "")

# Add "Timepoint" column
i3 <- grepl("PRE", pData(gs_b1)$`SAMPLE ID`) 
i4 <- grepl("POST", pData(gs_b1)$`SAMPLE ID`) 
pData(gs_b1)$Timepoint <- NULL
pData(gs_b1)$Timepoint[i3] <- "PRE"
pData(gs_b1)$Timepoint[i4] <- "POST"

# Reformat "SAMPLE ID" column 
pData(gs_b1)$'SAMPLE ID'[i3] <- paste0(pData(gs_b1)$'SAMPLE ID'[i3], "-1") %>%
  str_replace("PRE-", "")
pData(gs_b1)$'SAMPLE ID'[i4] <- paste0(pData(gs_b1)$'SAMPLE ID'[i4], "-2") %>%
  str_replace("POST-", "")

# Make sure nodes, pData, and markers are consistent among the three batches
setdiff(sort(gh_get_pop_paths(gs_b1)), sort(gh_get_pop_paths(gs_b4)))
setdiff(sort(gh_get_pop_paths(gs_b1)), sort(gh_get_pop_paths(gs_b5)))
all(sort(gh_get_pop_paths(gs_b1)) == sort(gh_get_pop_paths(gs_b4)))
all(sort(gh_get_pop_paths(gs_b1)) == sort(gh_get_pop_paths(gs_b5)))
all(markernames(gs_b1) == markernames(gs_b4))
all(markernames(gs_b1) == markernames(gs_b5))
all(colnames(pData(gs_b1)) == colnames(pData(gs_b4))) # FALSE
all(colnames(pData(gs_b1)) == colnames(pData(gs_b5))) # FALSE

# Reorder Batch 1 pData columns to match Batches 4 and 5
pData(gs_b1) <-  select(pData(gs_b1), colnames(pData(gs_b4)))

# Merge GatingSets from all batches
gs <- merge_list_to_gs(c(gs_b1, gs_b4, gs_b5))

# 3H and 44H were found out to be seropositive, so reassign them as convalescent
pData(gs)$`PATIENT ID`[pData(gs)$`PATIENT ID` == "3H"] <- "3C"
pData(gs)$`PATIENT ID`[pData(gs)$`PATIENT ID` == "44H"] <- "44C"
pData(gs)$`SAMPLE ID`[pData(gs)$`SAMPLE ID` == "3H-1"] <- "3C-1"
pData(gs)$`SAMPLE ID`[pData(gs)$`SAMPLE ID` == "3H-2"] <- "3C-2"
pData(gs)$`SAMPLE ID`[pData(gs)$`SAMPLE ID` == "44H-1"] <- "44C-1"
pData(gs)$`SAMPLE ID`[pData(gs)$`SAMPLE ID` == "44H-2"] <- "44C-2"

# Add column to pData indicating sample group (Naive vs. Conv and PRE vs. POST)
i5 <- grepl("H", pData(gs)$`PATIENT ID`) & grepl("PRE", pData(gs)$Timepoint)
i6 <- grepl("H", pData(gs)$`PATIENT ID`) & grepl("POST", pData(gs)$Timepoint)
i7 <- grepl("C", pData(gs)$`PATIENT ID`) & grepl("PRE", pData(gs)$Timepoint)
i8 <- grepl("C", pData(gs)$`PATIENT ID`) & grepl("POST", pData(gs)$Timepoint)
pData(gs)$Group <- NULL
pData(gs)$Group[i5] <- "Naive PRE"
pData(gs)$Group[i6] <- "Naive POST"
pData(gs)$Group[i7] <- "Conv PRE"
pData(gs)$Group[i8] <- "Conv POST"

## Grab the cytokine and memory gates from CD8+ and add under DN gate 
#gates_to_copy <- c("/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CCR7+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD45RA+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD107a+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD154+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IFNg+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL2+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL4_5_13+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL17a+",
#"/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/TNFa+")

#for(path in gates_to_copy) {
#  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=path),
#             parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/DN")
#}
#recompute(gs, "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/DN")

# Plot gating tree
png(here::here("out/QC/GatingTree.png"), width = 7, height = 5, units = "in", res = 300)
plot(gs, fontsize=15, bool=T)
dev.off()

# Save GatingSet 
save_gs(gs, here::here("out/GatingSets/HAARVIVAC_GatingSet"))

## Perform QC ##
# Load gating set if needed: 
#gs <- load_gs(here::here("out/GatingSets/HAARVIVAC_GatingSet"))

dput(gh_get_pop_paths(gs))

# Extract CD3, CD4, and CD8 counts and combine with phenotype data
cd3_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live" 
cd4_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+"
cd8_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+"

cd3_cd4_cd8_counts <- pData(gs) %>%
  rownames_to_column(var="sample") %>%
  left_join(gs_pop_get_stats(gs, nodes = c(cd3_path, cd4_path, cd8_path))) %>%
  pivot_wider(names_from = "pop", values_from = "count") %>%
  rename(CD3 = !!cd3_path,
         CD4 = !!cd4_path,
         CD8 = !!cd8_path)

# Plot CD3 Count
png(here::here("out/QC/Counts/CD3_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd3_path < 10000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 10000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Plot CD4 Count
png(here::here("out/QC/Counts/CD4_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd4_path < 3000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Plot CD8 Count
png(here::here("out/QC/Counts/CD8_Counts.png"), width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts %>% 
         mutate(Color = ifelse(cd8_path < 3000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD8)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Investigate the low count wells
low_count <- cd3_cd4_cd8_counts %>%
  dplyr::filter(CD3 < 10000 | CD4 < 3000 | CD8 < 3000) %>%
  select("SAMPLE ID", "Stim", "CD3", "CD4", "CD8") %>%
  arrange("SAMPLE ID")

low_count

# I will filter based on CD4 and CD8 count when constructing the COMPASSContainer

## Plot DMSO signal stratified by cohort ##
# Load gating set if needed: 
# gs <- load_gs(here::here("out/GatingSets/HAARVIVAC_GatingSet"))

# Get nodes of interest (include parent nodes and markers of interest)
dput(gh_get_pop_paths(gs))
nodes <- c("/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+",
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+",
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD107a+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD154+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IFNg+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL2+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL4_5_13+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL17a+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/TNFa+",
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD107a+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD154+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IFNg+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL2+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL4_5_13+", 
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL17a+",
           "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/TNFa+")
nodes_short <- str_replace(nodes, "\\/Time\\/CD3\\+\\/CD14-\\CD19\\-\\/Lymphocytes\\/Singlet\\/Live\\/", "")

# Get DMSO counts
dmso_freq <- subset(gs, Stim == "DMSO") %>%
  gs_pop_get_count_with_meta(subpopulations = nodes) %>%
  pivot_wider(names_from = Population, values_from = Count) %>%
  rename_at(vars(all_of(nodes)), ~ nodes_short) 

# Plot DMSO frequencies and perform Kruskal-Wallis test among cohorts
# Argument "pop" is the list of nodes of interest
plot_pop <- function(pop) {     
  parent <- sub("(.*)\\/.*", "\\1", pop)
  tmp_dat <- dmso_freq %>%
    mutate(prop = !!as.name(pop) / ParentCount)
  kw_p <- kruskal.test(prop ~ Group, data = tmp_dat)$p.value
  p.unadj.text <- sprintf("Kruskal-Wallis Test: p-unadj%s",
                          if_else(kw_p < 0.001, "<0.001", paste0("=", sub("0.", ".", round(kw_p, 3)))))
  
  ggplot(tmp_dat, aes(Group, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width = 0.15, height = 0, pch = 21, fill = "grey", alpha = 0.8) +
    labs(y = sprintf("%% %s of %s", sub(".*\\/(.*)", "\\1", pop), parent),
         caption = paste0("DMSO frequencies\n", p.unadj.text)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(color="black", size=15),
          axis.text.x = element_text(color="black", size=15),
          plot.title = element_blank(),
          plot.caption = element_text(size=12),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1.3, 0.2, 0, 0.2, "cm")) +
    scale_y_continuous(labels = function(x) paste0(x*100))
}

plot_pop(nodes_short[[4]]) # The first three nodes in the list are parent nodes

for(pop in nodes_short[4:length(nodes_short)]) {
  png(file=here::here(sprintf("out/QC/DMSO_Signal/%s_vs_Cohort.png", 
                              sub("\\/", "_", pop))), width=408, height=265, units = "px")
  print(plot_pop(pop))
  dev.off()
}
