#' Calculate the consensus matrix of the network
#'
#' The function takes as input the results from the consensus clustering algorithm and calculates
#' the consensus matrix
#'
#' @param cl the results from the consensus clustering algorithm
#' @return The corresponding consensus matrix
#' @keywords internal
#' @noRd
getConsensusMatrix <- function(cluster_results){

  ##set up parameters
  num_iters <- length(cluster_results)
  num_features <- length(cluster_results[[1]]$membership)
  consensus_matrix <- matrix(0, num_features, num_features)

  #check consensus
  for (feat1 in 1:(num_features-1)){
    for (feat2 in (feat1+1):num_features){
      for (algo_iter in 1:num_iters){
        cluster_iteration = cluster_results[[algo_iter]]$membership
        consensus_matrix[feat1,feat2] <- consensus_matrix[feat1,feat2] + (length(unique(cluster_iteration[c(feat1,feat2)])) == 1)
      }
    }
  }

  #return conensus matrix
  consensus_matrix <- consensus_matrix + t(consensus_matrix)
  consensus_matrix <- consensus_matrix/num_iters
  diag(consensus_matrix) <- rep(1, num_features)
  return(consensus_matrix)
}

#' Cluster an adjacency graph to identify metabolic modules
#'
#' This function takes as input an igraph graph object made from an adjacency matrix performs seven clustering
#' algorithms from the \code{\link{igraph}} package:
#'
#' 1. \code{\link{igraph::cluster_edge_betweenness}}
#'
#' 2. \code{\link{igraph::cluster_fast_greedy}}
#'
#' 3. \code{\link{igraph::cluster_infomap}}
#'
#' 4. \code{\link{igraph::cluster_label_prop}}
#'
#' 5. \code{\link{igraph::cluster_louvain}}
#'
#' 6. \code{\link{igraph::cluster_walktrap}}
#'
#' 7. \code{\link{igraph::cluster_leading_eigen}}
#'
#' and outputs a list of the results
#'
#' @param adjacency_graph An igraph graph object created from an adjacency matrix
#' @param graph_weights Edge weights to be used during clustering ie. the partial correlations for each feature-feature
#' interaction
#'
#' @returns A list containing the clustering results from each of the seven algorithms
#' @keywords internal
#' @noRd
ensembl_cluster <- function(adjacency_graph,
                            graph_weights = NULL){


  ##initiate list
  clustering_results <- vector("list", length = 7)

  ##perform clustering
  clustering_results[[1]] <- cluster_edge_betweenness(adjacency_graph, weights = graph_weights)
  clustering_results[[2]] <- cluster_fast_greedy(adjacency_graph, weights = graph_weights)
  clustering_results[[3]] <- cluster_infomap(adjacency_graph, e.weights = graph_weights)
  clustering_results[[4]] <- cluster_label_prop(adjacency_graph, weights = graph_weights)
  clustering_results[[5]] <- cluster_louvain(adjacency_graph, weights = graph_weights)
  clustering_results[[6]] <- cluster_walktrap(adjacency_graph, weights = graph_weights)
  clustering_results[[7]] <- tryCatch(cluster_leading_eigen(adjacency_graph, weights = graph_weights),
                                      error = function(some_error){
                                        message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
                                        message('This is a known issue with a dependency and will not affect your results')
                                        return(NA)
                                      })

  return(clustering_results)
}
#' Perform consensus clustering
#'
#' This function clusters the biological networks constructed using \code{\link{getNetworks}} using a consensus clustering
#' approach described in Ma et al. *(Please see references for more details)*
#'
#' Seven clustering algorithms from the \code{\link{igraph}} package:
#' 1. \code{\link{igraph::cluster_edge_betweenness}}
#'
#' 2. \code{\link{igraph::cluster_fast_greedy}}
#'
#' 3. \code{\link{igraph::cluster_infomap}}
#'
#' 4. \code{\link{igraph::cluster_label_prop}}
#'
#' 5. \code{\link{igraph::cluster_louvain}}
#'
#' 6. \code{\link{igraph::cluster_walktrap}}
#'
#' 7. \code{\link{igraph::cluster_leading_eigen}}
#'
#' are performed iteratively on the adjacency matrix constructed using \code{\link{getNetworks}} until a consensus
#' is reached on resulting subnetwork membership, or the specified max_iterations is reached.
#'
#' @param adjacency_graph An adjacency matrix of the determined network
#' @param tau The consensus probability threshold for agreement among clustering runs
#' @param max_iterations Maximum number of iterations to perform trying to reach consensus.
#' @return Sub-network determinations for the nodes within the input network
#'
#' @import igraph
#' @keywords internal
#' @noRd
run_consensus_cluster <- function(adjacency_graph, tau=0.5, max_iterations = 5){

  ##cluster the adjacency graph
  clustering_results <- ensembl_cluster(adjacency_graph, graph_weights = NULL)

  ##get consensus matrix
  consensus_matrix <- getConsensusMatrix(clustering_results[!(is.na(clustering_results))])

  ##set iter and start loop
  iter <- 0
  for(x in 1:max_iterations){

    ##stop iterations if at max
    if(length(table(consensus_matrix)) < 3){

      message(paste0("Consensus was reached in: ", iter, " iterations!"))
      break
    }

    ##get thresholded consensus matrix
    diag(consensus_matrix) <- 0
    threshold_consensus_matrix <- consensus_matrix * (consensus_matrix > tau)

    ##great graph from consensus matrix
    threshold_consensus_graph <- graph.adjacency(threshold_consensus_matrix,
                                                 mode="undirected",
                                                 weighted = TRUE)

    ##run clustering algorithms
    final_consensus_cluster <- ensembl_cluster(threshold_consensus_graph,
                                               graph_weights = E(threshold_consensus_graph)$weight)


    #get new consensus matrix
    consensus_matrix <- getConsensusMatrix(final_consensus_cluster[!(is.na(final_consensus_cluster))])

    #add to iter
    iter <- iter + 1
  }

  new.order <- order(final_consensus_cluster[[1]]$membership)
  return(list(final_consensus_cluster = final_consensus_cluster[[1]]$membership,
              consensus_matrix = consensus_matrix,
              order = new.order,iter = iter))
}
