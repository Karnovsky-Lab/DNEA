
#' getConsensusMatrix will calculate the consensus matrix of the network
#'
#' The function takes as input the results from the consensus clustering algorithm and calculates
#' the consensus matrix
#'
#' @param cl the results from the consensus clustering algorithm
#' @return The corresponding consensus matrix
#' @keywords internal
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

#' Performs consensus clustering
#'
#' This function will take as input an adjacency matrix graph from the determined networks and perform
#' consensus clustering using one of several methods to choose from. The output results in sub-network
#' classification for the nodes within the network.
#'
#' @param graph An adjacency matrix of the determined network
#' @param K The number of iterations for clustering. This parameter is not necessary for default "ensemble"
#' @param tau The consensus probabilty threshold for agreement among clustering runs
#' @param method The consensus clustering method to be used. The options are as follows:\n
#'        *1. "ensemble" - *indicates that all seven of the available clustering methods
#'        (cluster_edge_betweenness, cluster_fast_greedy, cluster_infomap, cluster_label_prop,
#'        cluster_leading_eigen, cluster_louvain, cluster_walktrap) should be used.\n
#'        *2. "lpm" - *utilizes the cluster_label_prop, cluster_infomap, and cluster_walktrap methods.\n
#'        *3. "walktrap" - *utilizes only cluster_walktrap.\n
#'        *4. "infomap" - *utilizes only infomap.\n
#' @param num_iterations The number of clustering iterations to perform - this parameter not relevant
#'        for the "ensemble" method. Default is 10 iterations.
#' @param maxIter Maximum number of iterations to perform trying to reach consensus.
#'
#' @return Sub-network determinations for the nodes within the input network
#'
#' @import igraph
#' @import furrr
#' @noRd
run_consensus_cluster <- function(adjacency_graph, num_iterations=10, tau=0.5,
                                  method= c("ensemble", "lpm", "infomap", "walktrap"),
                                  maxIter=5){

  method <- match.arg(method)
  if (method=="ensemble"){
    clustering_results <- list()
    set.seed(1)

    clustering_results[[1]] <- cluster_edge_betweenness(adjacency_graph, weights = NULL)
    clustering_results[[2]] <- cluster_fast_greedy(adjacency_graph, weights = NULL)
    clustering_results[[3]] <- cluster_infomap(adjacency_graph, e.weights = NULL)
    clustering_results[[4]] <- cluster_label_prop(adjacency_graph, weights = NULL)
    clustering_results[[5]] <- cluster_louvain(adjacency_graph, weights = NULL)
    clustering_results[[6]] <- cluster_walktrap(adjacency_graph, weights = NULL)
    clustering_results[[7]] <- tryCatch(cluster_leading_eigen(adjacency_graph, weights = NULL),
                                        error = function(some_error){
                                          message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
                                          message('This is a known issue with a dependency and will not affect your results')
                                          return(NA)
                                        })
  } else {

    clustering_results <- vector("list", num_iterations)
    for (algo_iter in 1:num_iterations){
      set.seed(algo_iter)
      if (method=="lpm"){
        clustering_results[[algo_iter]] <- cluster_label_prop(adjacency_graph, weights=NULL)
      } else if (method=="infomap"){
        clustering_results[[algo_iter]] <- cluster_infomap(adjacency_graph, e.weights=NULL)
      } else if (method=="walktrap"){
        clustering_results[[algo_iter]] <- cluster_walktrap(adjacency_graph, weights=NULL, steps=algo_iter)
      }
    }
  }

  consensus_matrix <- getConsensusMatrix(clustering_results[!(is.na(clustering_results))])
  iter <- 0
  while(length(table(consensus_matrix)) > 2 && iter < maxIter){
    diag(consensus_matrix) <- 0
    threshold_consensus_matrix <- consensus_matrix * (consensus_matrix > tau)

    threshold_consensus_graph <- graph.adjacency(threshold_consensus_matrix, mode="undirected", weighted = TRUE)
    if (method=="ensemble"){
      final_consensus_cluster <- list()
      set.seed(iter)
      final_consensus_cluster[[1]] <- cluster_edge_betweenness(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[2]] <- cluster_fast_greedy(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[3]] <- cluster_infomap(threshold_consensus_graph, e.weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[4]] <- cluster_label_prop(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[5]] <- cluster_louvain(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[6]] <- cluster_walktrap(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
      final_consensus_cluster[[7]] <- tryCatch(cluster_leading_eigen(threshold_consensus_graph,
                                                 weights = E(threshold_consensus_graph)$weight),
                           error = function(some_error){
                             message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
                             message('This is a known issue with a dependency and will not affect your results')
                             return(NA)
                           })
    } else {
      final_consensus_cluster <- vector("list",num_iterations)
      for (k in 1:num_iterations){
        set.seed(k)
        if (method=="lpm"){
          final_consensus_cluster[[k]] <- cluster_label_prop(threshold_consensus_graph, weights = E(threshold_consensus_graph)$weight)
        } else if (method=="infomap"){
          final_consensus_cluster[[k]] <- cluster_infomap(threshold_consensus_graph, e.weights = E(threshold_consensus_graph)$weight)
        } else if (method=="walktrap"){
          final_consensus_cluster[[k]] <- cluster_walktrap(threshold_consensus_graph, weights = NULL, steps=k)
        }
      }
    }
    consensus_matrix <- getConsensusMatrix(final_consensus_cluster[!(is.na(final_consensus_cluster))])
    iter <- iter + 1
  }
  new.order <- order(final_consensus_cluster[[1]]$membership)
  return(list(final_consensus_cluster = final_consensus_cluster[[1]]$membership,
              consensus_matrix = consensus_matrix,
              order = new.order,iter = iter))
}

# run_consensus_cluster <- function(adjacency_graph, K=10, tau=0.5,
#                                   method=c("infomap","lpm","walktrap","ensemble"),
#                                   maxIter=5,
#                                   runParallel = FALSE,
#                                   nCores = nCores,
#                                   main.seed = 101){
#   method <- match.arg(method)
#   if (method=="ensemble"){
#     #clustering_results <- list()
#     set.seed(main.seed)
#     if(runParallel){
#
#       plan('multisession', workers = nCores)
#
#       clustering_results <- furrr::future_invoke_map(.f = c('cluster_edge_betweenness',
#                                                             'cluster_fast_greedy',
#                                                             'cluster_infomap',
#                                                             'cluster_label_prop',
#                                                             'catch_cluster_leading_eigen',
#                                                             'cluster_louvain',
#                                                             'cluster_walktrap'),
#                                                      .x = list(adjacency_graph),
#                                                      .progress = TRUE,
#                                                      .options = furrr_options(packages = 'igraph'))
#
#     } else{
#
#       clustering_results <- furrr::future_invoke_map(.f = c('cluster_edge_betweenness',
#                                                             'cluster_fast_greedy',
#                                                             'cluster_infomap',
#                                                             'cluster_label_prop',
#                                                             'catch_cluster_leading_eigen',
#                                                             'cluster_louvain',
#                                                             'cluster_walktrap'),
#                                                      .x = list(adjacency_graph),
#                                                      .progress = TRUE,
#                                                      .options = furrr_options(packages = 'igraph'))
#       # clustering_results[[1]] <- cluster_edge_betweenness(graph, weights = NULL)
#       # clustering_results[[2]] <- cluster_fast_greedy(graph, weights = NULL)
#       # clustering_results[[3]] <- cluster_infomap(graph, e.weights = NULL)
#       # clustering_results[[4]] <- cluster_label_prop(graph, weights = NULL)
#       # clustering_results[[5]] <- catch_cluster_leading_eigen(graph, weights = NULL)
#       # clustering_results[[6]] <- cluster_louvain(graph, weights = NULL)
#       # clustering_results[[7]] <- cluster_walktrap(graph, weights = NULL)
#     }
#
#   } else {
#     clustering_results <- vector("list", K)
#     for (k in 1:K){
#       set.seed(k)
#       if (method=="lpm"){
#         clustering_results[[k]] <- cluster_label_prop(graph, weights = NULL)
#       } else if (method=="infomap"){
#         clustering_results[[k]] <- cluster_infomap(graph, e.weights = NULL)
#       } else if (method=="walktrap"){
#         clustering_results[[k]] <- cluster_walktrap(graph, weights = NULL, steps=k)
#       }
#     }
#   }
#   D <- getConsensusMatrix(clustering_results)
#   iter <- 0
#   while(length(table(D))>2 && iter<maxIter){
#     diag(D) <- 0
#     thresholded_D <- D*(D>tau)
#
#     Dgraph <- graph.adjacency(thresholded_D, mode="undirected", weighted = TRUE)
#     if (method=="ensemble"){
#       dcl <- list()
#       set.seed(iter)
#       dcl[[1]] <- cluster_edge_betweenness(Dgraph, weights = E(Dgraph)$weight)
#       dcl[[2]] <- cluster_fast_greedy(Dgraph, weights = E(Dgraph)$weight)
#       dcl[[3]] <- cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
#       dcl[[4]] <- cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
#       dcl[[5]] <- catch_cluster_leading_eigen(Dgraph, weights = E(Dgraph)$weight)
#       dcl[[6]] <- cluster_louvain(Dgraph, weights = E(Dgraph)$weight)
#       dcl[[7]] <- cluster_walktrap(Dgraph, weights = E(Dgraph)$weight)
#     } else {
#       dcl <- vector("list",K)
#       for (k in 1:K){
#         set.seed(k)
#         if (method=="lpm"){
#           dcl[[k]] <- cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
#         } else if (method=="infomap"){
#           dcl[[k]] <- cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
#         } else if (method=="walktrap"){
#           dcl[[k]] <- cluster_walktrap(Dgraph, weights = NULL, steps=k)
#         }
#       }
#     }
#     D <- getConsensusMatrix(dcl)
#     iter <- iter + 1
#   }
#   new.order <- order(dcl[[1]]$membership)
#   return(list(dcl=dcl[[1]]$membership,D=D,order=new.order,iter=iter))
# }
