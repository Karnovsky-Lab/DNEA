#'
#'
#'
#' Add additional metadata to the DNEAresults object
#'
#' This function will take additional metadata and add it to the specified dataframe in the metadata
#' slot. \strong{\emph{NOTE:}} The rownames of the new metadata must match the order of the input sample names or feature names,
#' respectively
#'
#' @param object A DNEAresults object
#' @param type sample or feature metadata
#' @param metadata a data.frame containing metadata to add
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{featureNames}}, \code{\link{sampleNames}},
#'
#' @return A DNEAresults object with the specified additions
#'
#' @examples
#' #import example data
#' data(TEDDYresults)
#'
#' #create sample metadata file
#' new_metadat <- data.frame(new_group = c(rep("group1", 50), rep('group2', 50)),
#'                           row.names = sampleNames(TEDDYresults))
#'
#' #add new metadata to DNEAresults object
#' TEDDYresults <- includeMetadata(object = TEDDYresults, type = "sample", metadata = new_metadat)
#'
#' @export
includeMetadata <- function(object, type = c('sample', 'feature'), metadata){

  type = match.arg(type)
  if(type == 'sample'){
    if(all(sampleNames(object) == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){


        object@metadata[["samples"]][[colnames(metadata)[i]]] <- metadata[, i]
      }
    } else{

      stop('new metadata order does not match sample order in DNEAobject')

    }
  } else{
    if(all(featureNames(object) == rownames(metadata)) |
       all(object@metadata[["features"]]$clean_feature_names == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["features"]][[metadata[, i]]] <- metadata[, i]
      }
    } else{
      stop('new metadata order does not match feature order in DNEAobject')
    }
  }

  return(object)
}

#' Save network information for input to Cytoscape
#'
#' This function will save the node and edge information as .csv files in the working directory.
#' The files are already formatted for input into Cytoscape.
#'
#' @param object A DNEAresults object
#' @param file_path The filepath to save the node and edge lists to. If **NULL**, the files will be saved to the working
#' directory
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{edgeList}}, \code{\link{nodeList}}
#'
#' @return Two .csv files, one for the node list and one for the edge list, saved to the specified file path.
#' @examples
#' #import example data
#' data(TEDDYresults)
#'
#' #save node and edge list for input to cytoscape
#' \donttest{getNetworkFiles(TEDDYresults)}
#'
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object, file_path=NULL){

  if(missing(file_path)){
    file_path <- paste0(getwd(), "/")
  }
  #save node list
  write.csv(object@node_list, paste0(file_path, object@project_name,'_nodelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  #save edge list
  write.csv(object@edge_list, paste0(file_path, object@project_name,'_edgelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  return(object)
}

#' Visualize the biological networks identified via DNEA
#'
#' The function plots the condition networks as well as the specified subnetworks identified via DNEA
#'
#' @param object A DNEAresults object
#' @param type "group_networks" will plot the case, control, and total network. "subnetworks" will plot the subnetwork
#' specified by teh subnetwork parameter
#' @param subnetwork The subnetwork to plot
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{getNetworks}}, \code{\link{clusterNet}}
#'
#'
#' @returns a plot of the specified network
#'
#' @examples
#' #import example data
#' data(TEDDYresults)

#'
#' #plot the networks
#' plotNetworks(object = TEDDYresults, type = "group_networks")
#' plotNetworks(object = TEDDYresults, type = "subnetworks", subnetwork = 1)
#'
#' @import igraph
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#' @importFrom autoimage reset.par
#' @export
plotNetworks <- function(object,
                         type = c("group_networks", "subnetworks"),
                         subnetwork){

  #get type argument
  type <- match.arg(type)

  #grab network graph
  network_graph <- adjacencyGraph(object, "joint_graph")

  #grab node list
  node_list <- nodeList(object)

  if(type == "group_networks"){

    #grab necessary info
    edge_list <- edgeList(object)
    group1_nodes <- unique(c(edge_list$Metabolite.A[edge_list$edge == networkGroups(object)[[1]]],
                             edge_list$Metabolite.B[edge_list$edge == networkGroups(object)[[1]]]))

    group2_nodes <- unique(c(edge_list$Metabolite.A[edge_list$edge == networkGroups(object)[[2]]],
                             edge_list$Metabolite.B[edge_list$edge == networkGroups(object)[[2]]]))

    #graph for control network
    cluster_c1 <- induced.subgraph(network_graph, V(network_graph)$name[match(group1_nodes, V(network_graph)$name)])

    #graph for total network
    cluster_c3 <- induced.subgraph(network_graph, V(network_graph)$name)

    #graph for case network
    cluster_c2 <- induced.subgraph(network_graph, V(network_graph)$name[match(group2_nodes, V(network_graph)$name)])

    ##plot
    #set layout
    par(mfrow = c(1,3))

    #control network
    plot(cluster_c1, vertex.label = V(cluster_c1)$name, vertex.label.cex = 1,
         layout = layout.fruchterman.reingold(cluster_c1),
         main = networkGroups(object)[[1]])

    #total network
    plot(cluster_c3, vertex.label = V(cluster_c3)$name, vertex.label.cex = 1,
         layout = layout.fruchterman.reingold(cluster_c3),
         main = "Total Network")

    #case network
    plot(cluster_c2, vertex.label = V(cluster_c2)$name, vertex.label.cex = 1,
         layout = layout.fruchterman.reingold(cluster_c2),
         main = networkGroups(object)[[2]])

    reset.par()
  }else if(type == "subnetworks"){

    #check that subnetwork given is relevant
    if(all(is.na(match(subnetwork, node_list$membership)))) stop("The subnetwork specified does not exist!\nPlease specify a value contained in the membership column of the node list")
    cluster_c <- induced.subgraph(network_graph, V(network_graph)$name[node_list$membership == subnetwork])

    plot(cluster_c, vertex.label = V(cluster_c)$name, vertex.label.cex = 1,
         layout = layout.fruchterman.reingold(cluster_c),
         main = paste0("subnetwork: ", subnetwork))
  }
}
filterNetworks.DNEAresults <- function(data,
                                       pcor,
                                       top_percent_edges){

  ##grab adjacency matrices
  weighted_adjacency_matrices <- adjacencyMatrix(data, weighted = TRUE)

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length = length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  #can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){

    #filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      #grab quantile value
      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs = (1 - top_percent_edges))

      #filter
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    #must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  #store the adjacency matrices in DNEAresults object
  adjacencyMatrix(x = data, weighted = TRUE) <- weighted_adjacency_matrices
  adjacencyMatrix(x = data, weighted = FALSE) <- unweighted_adjacency_matrices

  ##update edge list
  #initiate output dataframe
  pairs <- combn(as.character(featureNames(data)), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names = FALSE)

  #concatenate results into dataframe
  edge_list[,1:2] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs)) #non-edge
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 == 0] <- names(weighted_adjacency_matrices)[[1]] #control
  edge_list$edge[edge_list$pcor.0 == 0 & edge_list$pcor.1 != 0] <- names(weighted_adjacency_matrices)[[2]] #case
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 != 0] <- "Both" #Both
  edge_list <- edge_list[!is.na(edge_list$edge),] #remove non-edges
  rownames(edge_list) <- NULL

  #replace edgelist
  edgeList(data) <- edge_list


  #print message for total edges
  message(paste0(names(unweighted_adjacency_matrices)[[1]]," network specific edges: ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
  message(paste0(names(unweighted_adjacency_matrices)[[2]]," network specific edges: ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)
  message(rep("-", 35), appendLF = TRUE)
  message(paste0("Number of edges shared by both networks: ", sum(edge_list$edge == "Both")))
  message(paste0("Total number of edges in dataset: ", nrow(edge_list)))

  return(data)
}
#' Filter the adjacency matrices to only the edges that meet the filter conditions
#'
#' @description
#' This function takes as input a DNEAresults object and allows the user to filter the network edges by one of two methods:
#'
#' \enumerate{
#' \item \strong{Partial Correlation} - The networks can be filtered to only include edges greater than or equal to a specified partial correlation (pcor) value.
#' \item \strong{Top \emph{X%} of edges} - The networks can be filtered to only include the strongest \emph{X%} of edges by pcor value.}
#'
#' Filtering is performed on the case and control adjacency matrices separately. \cr
#'
#' @param data A DNEAresults object or list of adjacency matrices
#' @param pcor A pcor value of which to threshold the adjacency matrices. Edges with pcor values <= to this value will be removed.
#' @param top_percent_edges A value between 0-1 that corresponds to the top x% edges to keep in the networks
#' ie. top_percent_edges = 0.1 will keep only the top 10% strongest edges in the networks.
#' @author Christopher Patsalis
#' @seealso \code{\link{getNetworks}}, \code{\link{adjacencyMatrix}}
#' @returns The input object after filtering the egdes in the network according to the specified parameters
#'
#' @importFrom stats quantile
#' @rdname filterNetworks
#' @export
setMethod("filterNetworks", signature(data = "DNEAresults"), filterNetworks.DNEAresults)

filterNetworks.list <- function(data, pcor, top_percent_edges){

  ##rename input variable
  weighted_adjacency_matrices <- data

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length = length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  #can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){

    #filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      #grab quantile value
      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs = (1 - top_percent_edges))

      #filter
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    #must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  ##print message for total edges
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[1]],": ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[2]],": ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)

  return(list(weighted = weighted_adjacency_matrices, unweighted = unweighted_adjacency_matrices))
}
#'
#' @rdname filterNetworks
#' @keywords internal
setMethod("filterNetworks", signature(data = "list"), filterNetworks.list)
