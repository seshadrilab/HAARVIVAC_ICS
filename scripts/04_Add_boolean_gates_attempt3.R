library(here)
library(tidyverse)
library(flowWorkspace)
library(stringi)

## Load data ##
gsPath <- here::here("out/GatingSets/HAARVIVAC_GatingSet")
gs <- load_gs(gsPath)

## Add CD4 boolean gates ##
# Look at the nodes in the gating hierarchy
gh_get_pop_paths(gs) # 31 nodes total

# List of partial gating paths (these are existing gates)
cd4_paths <- c("CD4+/IL2+", "CD4+/IL4_5_13+", "CD4+/IFNg+", "CD4+/TNFa+",
             "CD4+/IL17a+", "CD4+/CD154+", "CD4+/CD107a+")

# Get all possible boolean subsets
bool_list <- do.call(c, lapply(seq_along(cd4_paths),
                               combn, x = cd4_paths,
                               simplify = FALSE))

# Drop the single marker subsets (e.g., CD4+/IL2+)
i <- lengths(bool_list)
single <- i == 1
bool_list_polyf <- bool_list[!single]

# Concatenate the boolean subsets because they're a list of lists right now
bool_list_polyf_cat <- stri_join_list(bool_list_polyf, sep = "&", collapse = NULL)

# Reformat the subsets in proper boolean logic (e.g., CD4+/IL2+&CD4+/IFNg+&!CD4+/TNFa+)
bool_list_polyf_cat_long <- vector("list", 120)

for(i in seq_along(bool_list_polyf_cat)){
  out <- bool_list_polyf_cat[i]
  for(k in seq_along(cd4_paths)){
    new_subset <- out
    # Add logical negation (!) for every marker not part of the subset
    if(isFALSE(grepl(cd4_paths[k], new_subset, fixed = TRUE))){
      out <- paste(new_subset, cd4_paths[k], sep = "&!")
    }
  }
  bool_list_polyf_cat_long[i] <- out
}

bool_list_polyf_cat_long <- unlist(bool_list_polyf_cat_long)

# Assign shortened names to each boolean subset in the list
names(bool_list_polyf_cat_long) <- bool_list_polyf_cat

# Add individual gates
for(i in seq_along(bool_list_polyf_cat_long)) {
  gs_pop_add(gs,  eval(substitute(booleanFilter(v), list(v = as.symbol(bool_list_polyf_cat_long[[i]])))),
             parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+")
}

# Add union gate
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(bool_list_polyf_cat_long), collapse="|"))))),
           parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+", name = "CD4_PolyF_Subsets")

# Check gating tree
gh_get_pop_paths(gs) # 152 nodes total

## Add activation gate (HLADR+CD38+) under CD4_PolyF_Subsets ##
cd4_polyf_parent <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD4_PolyF_Subsets"

hladr_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/HLADR+"
cd38_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD38+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", hladr_path,
                                                         "&", cd38_path))))),
           parent = cd4_polyf_parent, name = "HLADR+CD38+")

## Add memory gates under CD4_PolyF_Subsets ##
cd4_ccr7_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CCR7+"
cd4_cd45ra_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD45RA+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = cd4_polyf_parent, name = "Naive")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = cd4_polyf_parent, name = "TCM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = cd4_polyf_parent, name = "TEMRA")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = cd4_polyf_parent, name = "TEM")

## Recompute GatingSet ##
flowWorkspace::recompute(gs, "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+")

## Compare counts between the sum of the individual boolean gates vs CD4_PolyF_Subsets ##
# Counts from the individual boolean subset gates (excluding SEB)
indiv_pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = gs_get_pop_paths(gs)[32:151]), by = c("rowname" = "name")) %>%
  dplyr::rename(CD4_PolyF = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD4_PolyF) %>%
  dplyr::filter(Stim != "SEB")

sum(indiv_pop_counts$CD4_PolyF) # 18818

# Counts from CD4_PolyF_Subsets gate (excluding SEB)
pop_counts <- pData(gs) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(gs_pop_get_count_fast(gs, subpopulations = cd4_polyf_parent), by = c("rowname" = "name")) %>%
  dplyr::rename(CD4_PolyF = Count) %>%
  dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Group, CD4_PolyF) %>%
  dplyr::filter(Stim != "SEB")

sum(pop_counts$CD4_PolyF) # 3

# Something is wrong, the counts are really different and they keep changing

# ## Plot gating tree ##
# png(here::here("out/GatingTree_with_Bool_Subsets.png"), width = 7, height = 5, units = "in", res = 300)
# plot(gs, bool = T, fontsize = 10)
# dev.off()
# 
# ## Save GatingSet ##
# save_gs(gs, here::here("out/GatingSets/HAARVIVAC_GatingSet_with_bool"))
#gs <- load_gs(here::here("out/GatingSets/HAARVIVAC_GatingSet_with_bool"))
