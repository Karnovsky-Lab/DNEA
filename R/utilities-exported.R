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
#' @returns A DNEAresults object with the specified additions
#'
#' @examples
#' #import example data
#' data(dnw)
#' data(T1Dmeta)
#'
#' #make sure metadata has same sample order as DNEAresults object
#' T1Dmeta <- T1Dmeta[sampleNames(dnw), ]
#'
#' #add new metadata to DNEAresults object
#' dnw <- includeMetadata(object = dnw, type = "sample", metadata = T1Dmeta)
#'
#' @export
includeMetadata <- function(object, type = c('sample', 'feature'), metadata){

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')
  if(!any(inherits(metadata, "matrix") | inherits(metadata, "data.frame"))) stop('the input metadata should be of class "matrix" or "data.frame"!')

  type <- match.arg(type)
  if(type == 'sample'){
    if(all(sampleNames(object) == rownames(metadata))){
      for(i in seq(1, length(colnames(metadata)))){


        object@metadata[["samples"]][[colnames(metadata)[i]]] <- metadata[, i]
      }
    } else{

      stop('new metadata order does not match sample order in DNEAobject')

    }
  } else{
    if(all(featureNames(object) == rownames(metadata)) |
       all(object@metadata[["features"]]$clean_feature_names == rownames(metadata))){
      for(i in seq(1, length(colnames(metadata)))){

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["features"]][[metadata[, i]]] <- metadata[, i]
      }
    } else{
      stop('new metadata order does not match feature order in DNEAobject')
    }
  }

  return(object)
}

#' Include custom normalized data in the DNEAresults object
#'
#' This function allows the user to input custom-normalized data into the DNEAresults object for use in DNEA analysis.
#'
#' @param object A \code{DNEAresults} object
#' @param data An \emph{m x n} numeric matrix of custom-normalized expression expression data
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{DNEAresults},
#'
#' @returns A \code{DNEAresults} object with the added expression data in the @@assays slot
#'
#' @examples
#' #import example data
#' data(TEDDY)
#' data(dnw)
#' data(T1Dmeta)
#'
#' #transpose TEDDY data
#' TEDDY <- t(TEDDY)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[rownames(TEDDY),]
#'
#' dat <- list('DM:control' = TEDDY[T1Dmeta$group == "DM:control",],
#'             'DM:case' = TEDDY[T1Dmeta$group == "DM:case",])
#'
#' #log-transform and median center the expression data without scaling
#' newdat <- NULL
#' for(cond in dat){
#'  for(i in seq(1, ncol(cond))){
#'   my_median = median(cond[, i], na.rm = TRUE)
#'    my_range = range(cond[, i], na.rm = TRUE)
#'     scale_factor = max(abs(my_range-my_median))
#'     cond[, i] <- (cond[, i] - my_median) / scale_factor
#'   }
#'   newdat <- rbind(newdat, cond)
#' }
#'
#' #reorder to match dnw
#' newdat <- newdat[sampleNames(dnw), featureNames(dnw)]
#' newdat <- t(newdat)
#'
#' #add data
#' dnw <- addExpressionData(object = dnw, data = newdat)
#'
#' @rdname addExpressionData
#' @export
addExpressionData <- function(object, data){

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')
  if(!inherits(data, "matrix")) stop('the input metadata should be of class "matrix"!')
  if(!all(rownames(data) == rownames(expressionData(object, normalized = TRUE)))) stop("The feature order of new data does not match the feature order in the DNEAresults object!")
  if(!all(colnames(data) == colnames(expressionData(object, normalized = TRUE)))) stop("The sample order of new data does not match the sample order in the DNEAresults object!")
  if(!is.numeric(data)) stop("The new data should be a numeric matrix!")

  object@assays[["DNEA_scaled_data"]] <- expressionData(object, normalized = TRUE)
  object@assays[["scaled_expression_data"]] <- data

  validObject(object)
  return(object)
}


#' Save network information for input to Cytoscape
#'
#' This function will save the node and edge information as .csv files in the working directory.
#' The files are already formatted for input into Cytoscape.
#'
#' @param object A \code{DNEAresults} object
#' @param file_path The filepath to save the node and edge lists to. If **NULL**, the files will be saved to the working
#' directory
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{edgeList}}, \code{\link{nodeList}}
#'
#' @returns Two .csv files, one for the node list and one for the edge list, saved to the specified file path
#' @examples
#' #import example data
#' data(dnw)
#'
#' #save node and edge list for input to cytoscape
#'
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object, file_path = NULL){

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')

  if(missing(file_path)){
    file_path <- paste0(getwd(), "/")
  }else{
    if(!is.character(file_path)) stop("file_path should be a character string corresponding to the directory path in which to save the network files")
  }
  #save node list
  write.csv(nodeList(object), paste0(file_path, object@project_name,'_nodelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  #save edge list
  write.csv(edgeList(object), paste0(file_path, object@project_name,'_edgelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  return(object)
}

#' Visualize the biological networks identified via DNEA
#'
#' This function plots the total network, condition networks, or subnetworks as specified by the user. Purple nodes are
#' differential features, green indicates edges specific to group 1, and red indicates edges specific to group 2.
#'
#' @param object A \code{DNEAresults} object
#' @param type There are two possible arguments to \strong{type}: \emph{"group_networks"} specifies the whole network or condition networks.
#' \emph{"subnetworks"} specifies that one of the sub networks should be plot. Additional input via the \strong{subtype} parameter is required
#' @param subtype There are several possible arguments to \strong{subtype}. If \emph{type == "group_networks"}, \strong{subtype}
#' can be \emph{"All"} to plot the whole network (ie. both conditions in the data returned by \code{\link{networkGroups}}), or one of the
#' condition network names to plot the network corresponding to that condition. If \emph{type == "subnetworks"}, \strong{subtype}
#' should be a single-value numeric vector corresponding to the sub network to plot.
#' @param layout_func The layout in which to plot the specified network. Please see \code{\link[igraph]{plot.igraph}} for more information
#' @param main A character string to use as the plot title
#' @param node_size The size of the nodes in the plot. The default is 15. Please see \emph{vertex.size} parameter in
#' \code{\link[igraph]{tkplot}} for more details
#' @param edge_width The width of the edges in the plot. The default is 1. Please see \emph{width} parameter in
#' \code{\link[igraph]{tkplot}} for more details
#' @param label_size The size of the node labels in the plot. The default is 1. Please see \emph{label.size} in
#' \code{\link[igraph]{tkplot}} for more details.
#' @param label_font Specifies the font type to use in the plot. 1 is normal font, 2 is bold-type, 3 is italic-type,
#' 4 is bold- and italic-type. Please see the \emph{label.font} parameter in \code{\link[igraph]{plot.igraph}} for more details
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{getNetworks}}, \code{\link{clusterNet}}, \code{\link{networkGroups}}
#'
#'
#' @returns A plot of the specified network
#'
#' @examples
#' #import example data
#' data(dnw)

#'
#' #plot the networks
#' plotNetworks(object = dnw, type = "group_networks", subtype = "All")
#' plotNetworks(object = dnw, type = "subnetworks", subtype = 1)
#'
#' @import igraph
#' @export
plotNetworks <- function(object,
                         type = c("group_networks", "subnetworks"),
                         subtype = "All",
                         layout_func,
                         main = "",
                         node_size = 15,
                         edge_width = 1,
                         label_size = 1,
                         label_font = 1){

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')

  #get type argument
  type <- match.arg(type)

  #grab network graph
  network_graph <- adjacencyGraph(object, "joint_graph")

  #grab node list
  node_list <- nodeList(object)

  if(type == "group_networks"){

    #grab edge list
    edge_list <- edgeList(object)

    if(subtype == "All"){

      #graph for total network
      subtype_network <- induced.subgraph(network_graph, V(network_graph)$name)
    }else{

      #check the subtype exists
      if(!(subtype %in% c(networkGroups(object), "Both"))) stop("The subtype provided does not match a network within the data!")

      #grab metabolite names present in specified subtype
      subgroup_nodes <- unique(c(edge_list$Metabolite.A[edge_list$edge == subtype | edge_list$edge == "Both"],
                                 edge_list$Metabolite.B[edge_list$edge == subtype | edge_list$edge == "Both"]))

      #threshold network_graph network
      subtype_network <- induced.subgraph(network_graph, V(network_graph)$name[match(subgroup_nodes, V(network_graph)$name)])
    }

  }else if(type == "subnetworks"){

    #check that subnetwork given is relevant
    if(all(is.na(match(subtype, node_list$membership)))) stop("The subnetwork specified does not exist!\nPlease specify a value contained in the membership column of the node list")

    #threshold network_graph network
    subtype_network <- induced.subgraph(network_graph, V(network_graph)$name[node_list$membership == subtype])
  }

  #set graph layout
  if(missing(layout_func)){

    graph_layout <- layout.fruchterman.reingold(graph = subtype_network)
  }else{

    graph_layout <- layout_func(graph = subtype_network)
  }

  #plot the specified network
  plot(subtype_network, main = main, layout = graph_layout, vertex.size = node_size, width = edge_width,
       vertex.label = V(subtype_network)$name, vertex.label.cex = label_size, vertex.label.font = label_font)

}

#'
#'
filterNetworks.DNEAresults <- function(data,
                                       pcor,
                                       top_percent_edges){

  ##test for proper input
  if(!inherits(data, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')

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

    #check for valid input
    if(pcor <=0 | pcor >= 1) stop("the specified partial correlation threshold should be between 0 and 1!")

    #filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #check for valid input
    if(top_percent_edges <=0 | top_percent_edges >= 1) stop("the specified % threshold should be between 0 and 1!")

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
  edge_list[,c(1,2)] <- do.call(rbind, pairs)
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


  ##return to console messages for total edges
  #control network
  message(names(unweighted_adjacency_matrices)[[1]], " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[1]])/2) - sum(edge_list$edge == "Both"),
          appendLF = TRUE)

  #case network
  message(names(unweighted_adjacency_matrices)[[2]],
          " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[2]])/2) - sum(edge_list$edge == "Both"),
          appendLF = TRUE)

  #shared edges
  message(rep("-", 35), appendLF = TRUE)
  message("Number of edges shared by both networks: ",
          sum(edge_list$edge == "Both"), appendLF = TRUE)
  message("Total number of edges in dataset: ",
          nrow(edge_list), appendLF = TRUE)

  return(data)
}
#' Filter the adjacency matrices to only the edges that meet the filter conditions
#'
#' @description
#' This function takes as input a \code{DNEAresults} object and allows the user to filter the network edges by one of two methods:
#'
#' \enumerate{
#' \item \strong{Partial Correlation} - The networks can be filtered to only include edges greater than or equal to a specified partial correlation (pcor) value.
#' \item \strong{Top \emph{X%} of edges} - The networks can be filtered to only include the strongest \emph{X%} of edges by pcor value.}
#'
#' Filtering is performed on the case and control adjacency matrices separately. \cr
#'
#' @param data A \code{DNEAresults} object or list of adjacency matrices
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

  ##create edge list
  #initiate output dataframe
  pairs <- combn(as.character(colnames(weighted_adjacency_matrices[[1]])), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names = FALSE)

  #concatenate results into dataframe
  edge_list[,c(1,2)] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs)) #non-edge
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 == 0] <- names(weighted_adjacency_matrices)[[1]] #control
  edge_list$edge[edge_list$pcor.0 == 0 & edge_list$pcor.1 != 0] <- names(weighted_adjacency_matrices)[[2]] #case
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 != 0] <- "Both" #Both
  edge_list <- edge_list[!is.na(edge_list$edge),] #remove non-edges
  rownames(edge_list) <- NULL

  ##return to console messages for total edges
  #control network
  message(names(unweighted_adjacency_matrices)[[1]], " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[1]])/2) - sum(edge_list$edge == "Both"),
          appendLF = TRUE)

  #case network
  message(names(unweighted_adjacency_matrices)[[2]],
          " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[2]])/2) - sum(edge_list$edge == "Both"),
          appendLF = TRUE)

  #shared edges
  message(rep("-", 35), appendLF = TRUE)
  message("Number of edges shared by both networks: ",
          sum(edge_list$edge == "Both"), appendLF = TRUE)
  message("Total number of edges in dataset: ",
          nrow(edge_list), appendLF = TRUE)

  return(list(weighted = weighted_adjacency_matrices, unweighted = unweighted_adjacency_matrices,
              edge_list = edge_list))
}
#'
#' @rdname filterNetworks
#' @keywords internal
setMethod("filterNetworks", signature(data = "list"), filterNetworks.list)
