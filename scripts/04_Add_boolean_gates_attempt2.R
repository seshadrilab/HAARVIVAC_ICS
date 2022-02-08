library(here)
library(tidyverse)
library(flowWorkspace)
library(openCyto)
library(stringr)

## Load data ##
gsPath <- here::here("out/GatingSets/HAARVIVAC_GatingSet")
gs <- load_gs(gsPath)

## Add CD4 boolean gates ##
gh_get_pop_paths(gs) # 31 nodes total

# List of partial gating paths
cd4_paths <- c("CD4+/IL2+", "CD4+/IL4_5_13+", "CD4+/IFNg+", "CD4+/TNFa+",
               "CD4+/IL17a+", "CD4+/CD154+", "CD4+/CD107a+")

# Input for gating_args in gs_add_gating_method()
cd4_paths_cat <- paste(cd4_paths, collapse = ":")

# Add boolean combinations of existing gates
gs_add_gating_method(gs, gating_method = "polyFunctions",
                     parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+",
                     gating_args = cd4_paths_cat)

# Remove non-polyfunctional (<2 markers) boolean gates
gh <- gh_get_pop_paths(gs) 

for(i in seq_along(gh)){
  if(str_count(gh[i], "\\!") >= 6){
    gs_pop_remove(gs, gh[i])
  }
}

# Check the gating trees
pop_lists <- lapply(gs, gh_get_pop_paths)
unique(pop_lists) # 151 nodes total

# Check to make sure we got rid of all non-polyfunctional boolean gates
str_count(gh_get_pop_paths(gs), "\\!") 

# # Grab the boolean gates and define them all under "CD4_PolyF_Subsets"
# bool_paths <- grep(":", gh_get_pop_paths(gs), value = TRUE)
# bool_paths_short <- str_replace(bool_paths, "\\/Time\\/CD3\\+\\/CD14-\\CD19\\-\\/Lymphocytes\\/Singlet\\/Live\\/CD4\\+\\/", "")
# bool_paths_cat <- paste(bool_paths_short, collapse="|") # | means "or"
# 
# gs_pop_add(gs, eval(substitute(booleanFilter(v), list(v = as.symbol(bool_paths_cat)))),
#            parent = "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+",
#            name = "CD4_PolyF_Subsets")
# 
# #Error in as.symbol(bool_paths_cat) : 
# #  variable names are limited to 10000 bytes
# 
# # Need to find a different way to combine the boolean gates

## Add activation gate (HLADR+CD38+) under boolean gates ##
cd3_hladr_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/HLADR+"
cd3_cd38_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD38+"

bool_paths <- grep(":", gh_get_pop_paths(gs), value = TRUE)

for(i in seq_along(bool_paths)){
  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd3_hladr_path, "&", cd3_cd38_path))))),
             parent = bool_paths[i],
             name = "HLADR+CD38+")
}

## Add memory gates under boolean gates ##
cd4_ccr7_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CCR7+"
cd4_cd45ra_path <- "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD45RA+"

for(i in seq_along(bool_paths)){
  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                           "&", cd4_cd45ra_path))))),
             parent = bool_paths[i], name = "Naive")

  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                           "&!", cd4_cd45ra_path))))),
             parent = bool_paths[i], name = "TCM")

  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                           "&", cd4_cd45ra_path))))),
             parent = bool_paths[i], name = "TEMRA")

  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                           "&!", cd4_cd45ra_path))))),
             parent = bool_paths[i], name = "TEM")
}

## Recompute the GatingSet with the new gates ##
recompute(gs, "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+")

# Plot the gating tree
png(here::here("out/GatingTree_with_bool.png"), width = 7, height = 7, units = "in", res = 300)
plot(gs, bool = T, fontsize = 10)
dev.off()

# Save GatingSet
save_gs(gs, here::here("out/GatingSets/HAARVIVAC_GatingSet_with_bool"))
