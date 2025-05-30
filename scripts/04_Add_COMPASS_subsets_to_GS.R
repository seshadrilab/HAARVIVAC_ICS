library(here)
library(tidyverse)
library(flowWorkspace)
library(flowCore)

# The counts for the COMPASS subsets are only stored in COMPASSResult objects if they were discovered by COMPASS for that stim,
# so we have to manually add boolean gates for each subset to the GatingSets and then extract the count data later.

## Load data ##
gsPath <- here::here("out/GatingSets/HAARVIVAC_GatingSet")
gs <- load_gs(gsPath)

merged_cd4_compass_data <- readRDS("processed_data/Merged_CD4_COMPASS_Data.rds")
merged_cd8_compass_data <- readRDS("processed_data/Merged_CD8_COMPASS_Data.rds")

## Add the CD4+ COMPASS boolean cytokine subset gates to the GatingSet ##
mapMarkers <- c("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
cd4NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd4NodeMarkerMap) <- paste0("CD4+", "/", c("IL2+", "IL4_5_13+", "IFNg+", "TNFa+", "IL17a+", "CD154+", "CD107a+"))

cd4_cats_mod <- as.data.frame(merged_cd4_compass_data$catsMerged) %>% 
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>% 
  dplyr::rename_at(vars(cd4NodeMarkerMap), ~ names(cd4NodeMarkerMap))
cd4_booleanSubsets <- cd4_cats_mod %>% 
  rowwise() %>% 
  do(booleanSubset = paste(paste0(., colnames(cd4_cats_mod)), collapse="&")) %>% 
  ungroup() %>% 
  dplyr::pull(booleanSubset) %>% 
  unlist()
names(cd4_booleanSubsets) <- paste0("CD4_", gsub("CD4\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd4_booleanSubsets))))
for(booleanSubsetName in names(cd4_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd4_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "CD4+", name=booleanSubsetName))
}

dput(names(cd4_booleanSubsets))

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(cd4_booleanSubsets), collapse="|"))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+", name = "CD4_COMPASS_Subsets")

## Add the CD8 COMPASS boolean cytokine subset gates to the GatingSet ##
mapMarkers <- c("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
cd8NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd8NodeMarkerMap) <- paste0("CD8+", "/", c("IL2+", "IL4_5_13+", "IFNg+", "TNFa+", "IL17a+", "CD154+", "CD107a+"))

cd8_cats_mod <- as.data.frame(merged_cd8_compass_data$catsMerged) %>%
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>%
  dplyr::rename_at(vars(cd8NodeMarkerMap), ~ names(cd8NodeMarkerMap))
cd8_booleanSubsets <- cd8_cats_mod %>%
  rowwise() %>%
  do(booleanSubset = paste(paste0(., colnames(cd8_cats_mod)), collapse="&")) %>%
  ungroup() %>%
  dplyr::pull(booleanSubset) %>%
  unlist()
names(cd8_booleanSubsets) <- paste0("CD8_", gsub("CD8\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd8_booleanSubsets))))
for(booleanSubsetName in names(cd8_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd8_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+", name=booleanSubsetName))
}
dput(names(cd8_booleanSubsets))

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(cd8_booleanSubsets), collapse="|"))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+", name = "CD8_COMPASS_Subsets")

## Define the memory subpopulations under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd4_ccr7_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CCR7+"
cd4_cd45ra_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD45RA+"
cd8_ccr7_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CCR7+"
cd8_cd45ra_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD45RA+" 

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets", name = "Naive")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets", name = "TCM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets", name = "TEMRA")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets", name = "TEM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets", name = "Naive")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets", name = "TCM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets", name = "TEMRA")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets", name = "TEM")

## Add activation gates (HLADR+CD38+, HLADR+, and CD38+) under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets ##
cd3_hladr_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/HLADR+"
cd3_cd38_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD38+"

for (parent_path in c("/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets",
                      "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets")) {
  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd3_hladr_path,
                                                           "&", cd3_cd38_path))))),
             parent = parent_path, name = "HLADR+CD38+")
}

for (parent_path in c("/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets",
                      "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets")) {
  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd3_cd38_path),
             parent = parent_path, name = "CD38+")
}

for (parent_path in c("/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets",
                      "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets")) {
  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y = cd3_hladr_path),
             parent = parent_path, name = "HLADR+")
}

## Recompute the gating set with the new gates ##
flowWorkspace::recompute(gs)

# Counts from CD4_COMPASS_Subsets (excluding SEB)
cd4_compass_subsets_parentGate <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_COMPASS_Subsets"

cd4_pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = cd4_compass_subsets_parentGate), by = c("rowname" = "name")) %>%
  dplyr::rename(CD4_COMPASS_Subsets = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD4_COMPASS_Subsets) %>%
  dplyr::filter(Stim != "SEB")

sum(cd4_pop_counts$CD4_COMPASS_Subsets) 

# Counts from the individual CD4 boolean subset gates (excluding SEB)
cd4_indiv_pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = gs_get_pop_paths(gs)[32:49]), by = c("rowname" = "name")) %>%
  dplyr::rename(CD4_COMPASS_Subsets = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD4_COMPASS_Subsets) %>%
  dplyr::filter(Stim != "SEB")

sum(cd4_indiv_pop_counts$CD4_COMPASS_Subsets) 

# Counts from CD8_COMPASS_Subsets (excluding SEB)
cd8_compass_subsets_parentGate <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD8_COMPASS_Subsets"

cd8_pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = cd8_compass_subsets_parentGate), by = c("rowname" = "name")) %>%
  dplyr::rename(CD8_COMPASS_Subsets = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD8_COMPASS_Subsets) %>%
  dplyr::filter(Stim != "SEB")

sum(cd8_pop_counts$CD8_COMPASS_Subsets)

# Counts from the individual CD8 boolean subset gates (excluding SEB)
cd8_indiv_pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = gs_get_pop_paths(gs)[51:61]), by = c("rowname" = "name")) %>%
  dplyr::rename(CD8_COMPASS_Subsets = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD8_COMPASS_Subsets) %>%
  dplyr::filter(Stim != "SEB")

sum(cd8_indiv_pop_counts$CD8_COMPASS_Subsets) 

# Plot the gating tree
png(here::here("out/GatingTree_with_COMPASS_Subsets.png"), width = 7, height = 5, units = "in", res = 300)
plot(gs, bool = T, fontsize = 10)
dev.off()

# Save the gating set 
save_gs(gs, here::here("out/GatingSets/HAARVIVAC_GatingSet_with_COMPASS_Subsets"))
