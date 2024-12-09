library(OCNet)
library(igraph)
library(tidygraph)
library(gridExtra)
library(Matrix)


harmonic_centrality_manual <- function() {
  centrality_manual('dist_sp', FUN = netrankr::dist_inv)
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

harmonic_random_walk <- function() {
  centrality_manual('dist_rwalk', FUN = netrankr::dist_inv)
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

exponential_harmonic_centrality <- function(ell) {
  centrality_manual('dist_sp', FUN = function(x) exp(-x/ell))
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

exponential_random_walk <- function(ell) {
  centrality_manual('dist_rwalk', FUN = function(x) exp(-x/ell))
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}


betweenness_communicability_manual<- function(g) {
	depend_exp <- function(g) {
		A <- igraph::get.adjacency(g, "both", sparse = FALSE)
		eigen_A <- eigen(A)
		n <- nrow(A)
		expA <- expm(A)
		C <- (n - 1)^2 - (n - 1)
		combet <- matrix(0, n, n)
		for (i in 1:n) {
			E <- matrix(0, n, n)
			E[which(A[, i] == 1), i] <- -1
			E[i, which(A[i, ] == 1)] <- -1
			E <- A + E
			eigen_E <- eigen(E)
			expE <- expm(E)
			expE <- (expA - expE) / expA
			expE[i, ] <- 0
			expE[, i] <- 0
			diag(expE) <- 0
			combet[i, ] <-  rowSums(expE)/n
		}
		return(combet)
	}
	sum(depend_exp(g))
}

##################
data = read.csv("./rivers10to1000_sample.csv")
DIM=100
cell_area = 1000*sqrt(data$area)/DIM
##################
#Filter IDS
#IDs = which(nodes_data>=10 & nodes_data<=1000)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#i=as.numeric(args[1]) #USE THIS FOR EXTERNAL CALLS
i=5


set.seed(i)	
dimX=DIM
dimY=DIM
cellsize = cell_area[i]
N_nodes = data$nodes[i]
river_name = data$river.name[i]
river_area = data$area[i]


print(river_name)
print(N_nodes)
OCN<-create_OCN(dimX,dimY, nOutlet = 1, outletSide = "S",
	outletPos = round(dimX/2), periodicBoundaries = FALSE,
	typeInitialState = "I", flowDirStart = NULL, expEnergy = 0.5,
	cellsize, xllcorner = 0.5 * cellsize, yllcorner = 0.5 *
	cellsize, nIter = 50 * dimX * dimY, nUpdates = 100,
	initialNoCoolingPhase = 0, coolingRate = 0.25,
	showIntermediatePlots = FALSE, thrADraw = 0.002 * dimX * dimY *
	cellsize^2, easyDraw = NULL, saveEnergy = FALSE, saveExitFlag = FALSE,
	saveN8 = FALSE, saveN4 = FALSE, displayUpdates = 1)

OCN0=OCN
saveRDS(OCN,file=paste0("./data/OCN_D",as.character(DIM),"_",as.character(i),".rds"))
OCN <- landscape_OCN(OCN0, slope0 = 0.01)

thr <- find_area_threshold_OCN(OCN)
# find index corresponding to thr$Nnodes ~= 20
indThr <- which(abs(thr$nNodesAG - N_nodes) == min(abs(thr$nNodesAG - N_nodes)))
indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
thr <- thr$thrValues[indThr] # corresponding threshold area

OCN = aggregate_OCN(OCN, thrA = thr)
OCN <- rivergeometry_OCN(OCN) 
A = as.matrix(OCN$AG$W)

diag(A) = 0 #remove loops
sim_nodes = dim(A)[1]
print(sim_nodes)

print("Getting metrics...")
###################
#UNDIRECTED
##########################################
g = graph_from_adjacency_matrix(A,mode="max")
saveRDS(g,file=paste0("./data/UNDgraphD",as.character(DIM),"_",as.character(i),".rds"))
#
DG = mean(degree(g,normalized = FALSE)) # no variation, only for small
#
########################################
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = harmonic_centrality_manual()))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HC = mean(m0) # 
#
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = exponential_harmonic_centrality(5)))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HCE5 = mean(m0) # 
#
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = exponential_harmonic_centrality(10)))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HCE10 = mean(m0) # 
#
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = harmonic_random_walk()))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HCR = mean(m0) #
####
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = exponential_random_walk(5)))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HCRE5 = mean(m0) #
#####
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = exponential_random_walk(10)))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
HCRE10 = mean(m0) #
#####
####################################
BW = mean(betweenness(g, normalized = FALSE)) 
#
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = centrality_betweenness_current()))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
BWR = mean(m0) 
#
# g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = centrality_betweenness_communicability()))
# m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
BWC = betweenness_communicability_manual(g)
#
####################################
CL = mean(closeness(g)) #shortest distance
#
g0 = (as_tbl_graph(g) %>% activate(nodes) %>% mutate(importance = centrality_random_walk()))
m0 = (activate(g0, nodes) %>% as_tibble())[[1]]
RW = mean(m0) # clossness but using random-walk distances
#############
	#
CO = vertex_connectivity(g) # always 1
#
EI = mean(eigen_centrality(g)$value)

datum=c()
datum=rbind(datum,c(river_name,river_area,N_nodes,sim_nodes,DG,HC,HCE5,HCE10,HCR,HCRE5,HCRE10,BW,BWR,BWC,CL,RW,CO,EI))
datum=data.frame(datum)
colnames(datum)=c("name","area","req_nodes","nodes","DG","HC","HCE5","HCE10","HCR","HCRE5","HCRE10","BW","BWR","BWC","CL","RW","CO","EI")

saveRDS(datum,file=paste0("./data/metrics_undirected_D",as.character(DIM),"_",as.character(i),".rds"))


##################
##################
### DIRECTED
##################
##################

g = graph_from_adjacency_matrix(A,mode="directed")
saveRDS(g,file=paste0("./data/DIRgraphD",as.character(DIM),"_",as.character(i),".rds"))

#
DG = mean(degree(g, normalized = FALSE, mode = "in")) # no variation, only for small; mode "in" for "in-degree" 
#
########################################
HC = mean(harmonic_centrality(g, normalized = FALSE, mode = "out")) # “out” follows paths along the edge directions only, “in” traverses the edges in reverse, while “all” ignores edge directions
#
####################################
BW = mean(betweenness(g, directed = TRUE, normalized = FALSE)) 
#
####################################
CL = mean(closeness(g, mode = "in")) #shortest distance, mode "in" measures the paths to a vertex
#
CO = vertex_connectivity(g) # always 1
#
EI = mean(eigen_centrality(g, directed = TRUE)$value) # In the directed case, the centrality of a vertex is proportional to the sum of the centralities of vertices pointing to it.

datum=c()
datum=rbind(datum,c(river_name,river_area,sim_nodes,N_nodes,DG,HC,BW,CL,RW,CO,EI))
datum=data.frame(datum)
colnames(datum)=c("name","area","req_nodes","nodes","DG","HC","BW","CL","RW","CO","EI")

saveRDS(datum,file=paste0("./data/metrics_directed_D",as.character(DIM),"_",as.character(i),".rds"))
