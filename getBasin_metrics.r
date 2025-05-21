library(OCNet)
library(igraph)
library(tidygraph)
library(gridExtra)
library(Matrix)


####
###############
# R script to generate a square basin with area and number of nodes matching a real-word basin.
# In the terminal run the script using: Rscript getBasin_metrics.r basin_id trial
# basin_id is the basin index in the metadata file


harmonic_centrality_manual <- function() {
  centrality_manual('dist_sp', FUN = netrankr::dist_inv)
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

harmonic_random_walk <- function() {
  centrality_manual('dist_rwalk', FUN = netrankr::dist_inv)
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

##################
args = commandArgs(trailingOnly=TRUE)
basin_id = as.numeric(args[1])
trial = as.numeric(args[2])
metadata = read.csv("./basins_info.csv")
river_geo = metadata[basin_id,]
river_geo
print(basin_id)
print(trial)
print(river_geo)
##################################################
set.seed(trial+basin_id*1000)	
### args
req_area=river_geo['area'][1,1] # meters
req_asp=1
req_nodes=river_geo['nodes'][1,1]
river_basin=river_geo['MAIN_BAS'][1,1]
######################
DIM=150 #as.integer(5*sqrt(req_nodes))
cellsize=sqrt(req_area)/DIM  # meters
######################
dimX=as.integer(DIM*sqrt(req_asp))
dimY=as.integer(DIM/sqrt(req_asp))
area=dimX*dimY*cellsize^2
OCN<-create_OCN(dimX,dimY, nOutlet = 1, outletSide = "S",
	outletPos = round(dimX/2), periodicBoundaries = FALSE,
	typeInitialState = "V", flowDirStart = NULL, expEnergy = 0.5,
	cellsize, xllcorner = 0.5 * cellsize, yllcorner = 0.5 *
	cellsize, nIter = 50 * dimX * dimY, nUpdates = 100,
	initialNoCoolingPhase = 0.1, coolingRate = 0.5,
	showIntermediatePlots = FALSE, thrADraw = 0.002 * dimX * dimY *
	cellsize^2, easyDraw = NULL, saveEnergy = TRUE, saveExitFlag = FALSE,
	saveN8 = FALSE, saveN4 = FALSE, displayUpdates = 1)

OCN0=OCN
OCN <- landscape_OCN(OCN0, slope0 = 0.01)

thr <- find_area_threshold_OCN(OCN)
# find index corresponding to thr$Nnodes ~= 20
indThr <- which(abs(thr$nNodesAG - req_nodes) == min(abs(thr$nNodesAG - req_nodes)))
indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
thr <- thr$thrValues[indThr] # corresponding threshold area

OCN = aggregate_OCN(OCN, thrA = thr)
OCN <- rivergeometry_OCN(OCN) 
A = as.matrix(OCN$AG$W)
nodes = length(OCN$AG$X)


print("Getting metrics...")
###################
#UNDIRECTED
##########################################
g = graph_from_adjacency_matrix(A,mode="max")
#saveRDS(OCN,file=paste0("./final_trial/OCN_D",as.character(DIM),"_",as.character(basin_id),"_HS_trial_",as.character(trial),".rds"))
#
DG = mean(degree(g,normalized = FALSE)) # no variation, only for small
#
########################################
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = harmonic_centrality_manual()))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HC = mean(m0) # 
#
BW = mean(betweenness(g, normalized = FALSE)) 
#
CL = mean(closeness(g)) #shortest distance
#
EI = mean(eigen_centrality(g)$value)
datum=c()
datum=rbind(datum,c(river_basin,area,nodes,req_nodes,dimX/dimY,DG,HC,BW,CL,EI))
datum=data.frame(datum)
colnames(datum)=c("BasinName","area","nodes","req_nodes","asp","DG","HC","BW","CL","EI")

saveRDS(datum,file=paste0("./metrics_undirected_D",as.character(DIM),"_",as.character(basin_id),"_HS_ratio_",as.character(req_asp),"_",as.character(trial),".rds"))
print(datum)
