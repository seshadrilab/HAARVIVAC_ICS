#!/usr/bin/Rscript

# 1) Ideally, run this script from your terminal, i.e., NOT in RStudio, to take full advantage of parallel computation.
# The rest of the pipeline can and should be run in RStudio.
# 
# One way to do this is to open your terminal and navigate to the top level of the project folder, start up R in interactive mode,
# and paste in the contents of this script.
# 
# Or run: R CMD scripts/02_Run_Compass.R
# If you are using Unix/macOS, multicore processing (efficient forked/shared memory) only works when run from the terminal.
# Windows does not support forking.
# See: future::supportsMulticore()
# 
# 2) It's also a good idea to log the run: type "script out/CompassOutput/COMPASS_log.txt" in your terminal prior to running R.
# The COMPASS stdout and stderr will then get printed into COMPASS_log at the end of the furrr::future_pmap loop.
# Then after you quit R, type "exit" to close the log.

library(furrr)
library(COMPASS)
library(grid)
library(flowWorkspace)

source(here::here("scripts/Helper_Functions.R"))

date <- 20211020

stims_for_compass_runs <- c("S1", "S2", "NCAP")
parent_nodes_for_compass_runs <- c("CD4+", "CD8+")
seeds_for_compass_runs <- as.list(date:(date + length(stims_for_compass_runs)*length(parent_nodes_for_compass_runs) - 1))

stims_for_compass_runs_rep <- rep(stims_for_compass_runs, each = length(parent_nodes_for_compass_runs))
parent_nodes_for_compass_runs_rep <- rep(parent_nodes_for_compass_runs, times = length(stims_for_compass_runs))

gsPath <- here::here("out/GatingSets/HAARVIVAC_GatingSet/")
gs <- load_gs(gsPath)


# mapMarkers contains output of markernames(gs)
mapMarkers <- list("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a") 

future::supportsMulticore() # Run in terminal to get TRUE
# If you run this script in RStudio, the next line throws the following warning:
# "Warning message:
# [ONE-TIME WARNING] Forked processing ('multicore') is disabled in future (>= 1.13.0) when running R from RStudio,
# because it is considered unstable. Because of this, plan("multicore") will fall back to plan("sequential"),
# and plan("multiprocess") will fall back to plan("multisession") - not plan("multicore") as in the past.
# For more details, how to control forked processing or not, and how to silence this warning in future R sessions,
# see ?future::supportsMulticore "
future::plan(multicore(workers = max(1, availableCores() - 2)))

system.time({
  out <- furrr::future_pmap(.l = list(stims_for_compass_runs_rep,
                                      parent_nodes_for_compass_runs_rep,
                                      seeds_for_compass_runs),
                            .f = function(currentStim, parent, currentSeed) {
                              
                              o <- tryCatch( {
                                gsSub <- subset(gs, Stim %in% c("DMSO", currentStim))
                                
                                
                                
                                currentNodeMarkerMap <- mapMarkers
                                # currentNodeMarkerMap names are gating tree paths
                                names(currentNodeMarkerMap) <- paste0(parent, "/", c("IL2+", "IL4_5_13+", "IFNg+", "TNFa+", "IL17a+", "CD154+", "CD107a+"))
                                outDir <- here::here(sprintf("out/CompassOutput/%s/%s", parent, gsub(" ", "_", currentStim)))
                                if(!dir.exists(outDir)) {
                                  dir.create(outDir, recursive = T)
                                }
                                
                                runCompassOnce(gs=gsSub,
                                               seed=currentSeed,
                                               outDir=outDir,
                                               parentNode=parent,
                                               nodeMarkerMap=currentNodeMarkerMap,
                                               uniqueIdentifier="SAMPLE ID",
                                               treatmentCol="Stim",
                                               currentTreatment=currentStim,
                                               currentControl="DMSO",
                                               stratifyBy=NULL, 
                                               iter=40000,
                                               eventCountFilterThreshold=3000,
                                               textForRunOutputId=paste0(parent, "_", gsub(" ", "_", currentStim)))
                                gc()
                              }, error = function(e) { print(e) })
                              o
                            },
                            # Progress bar reflects how many COMPASS runs have completed
                            .progress = T,
                            .options = furrr_options(seed = T))
})