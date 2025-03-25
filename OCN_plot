library(OCNet)
library(igraph)
library(tidygraph)
library(gridExtra)
library(centiserve)
library(expm)

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

# betweenness_centrality_exponential <- function(a) {
#   centrality_manual('depend_sp', FUN = function(x) exp(-x/a))
# #centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
# }

betweenness_centrality_manual<- function() {
  centrality_manual('depend_sp')
#centrality_manual('dist_sp', FUN = netrankr::dist_inv) #normal harmonic centrality
}

betweenness_communicability_manual<- function() {
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
data = read.csv("./metrics_undirected_real.csv")
data$river_name
DIM=150
cell_area = 1000*sqrt(data$river_area)/DIM
##################
#Filter IDS
#IDs = which(nodes_data>=10 & nodes_data<=1000)
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#i=as.numeric(args[1])
i=918
data[i,]
set.seed(i)	
dimX=DIM/2
dimY=DIM
cellsize = cell_area[i]
N_nodes = data$N_nodes[i]
river_name = data$river_name[i]
river_area = data$river_area[i]
print(i)
print(river_name)

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

OCN <- landscape_OCN(OCN0, slope0 = 0.01)

thr <- find_area_threshold_OCN(OCN)
# find index corresponding to thr$Nnodes ~= 20
indThr <- which(abs(thr$nNodesAG - N_nodes) == min(abs(thr$nNodesAG - N_nodes)))
indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
thr <- thr$thrValues[indThr] # corresponding threshold area

OCN = aggregate_OCN(OCN, thrA = thr)
OCN <- rivergeometry_OCN(OCN) 
A = as.matrix(OCN$AG$W)


draw_subcatchments_OCN(OCN,riverColor = '#2171b5',drawRiver = TRUE,min_lwd=5, max_lwd=5,colPalette = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b'))
draw_simple_OCN(OCN, thrADraw = thr,easyDraw = TRUE)
# g <- OCN_to_igraph(OCN, level = "AG",layout = matrix(c(OCN$AG$X,OCN$AG$Y))
# plot.igraph(g)
