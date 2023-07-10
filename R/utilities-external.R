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
#' @param metadata a dataframe containing metadata to add
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{featureNames}}, \code{\link{sampleNames}},
#'
#' @return The same object as input with specified additions
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #create sample metadata file
#' new_metadat <- data.frame(new_group = c(rep("group1", 72),
#'                           rep("group2", 72)),
#'                           row.names(sampleNames(DNEA)))
#'
#' #add new metadata to DNEAresults object
#' DNEA <- includeMetadata(object = DNEA, type = "sample", metadata = new_metadat)
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
#' @return The same object as input and saves the node and edge information as .csv files in
#'          the working directory
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 5, BPPARAM = bpparam())
#'
#' #construct the networks
#' DNEA <- getNetworks(object = DNEA)
#'
#' #identify metabolic modules via consensus clustering
#' DNEA <- runConsensusCluster(object = DNEA)
#'
#' #perform pathway enrichment analysis using netGSA
#' DNEA <- runNetGSA(object = DNEA)
#'
#'
#' #save node and edge list for input to cytoscape
#' getNetworkFiles(DNEA)
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
#' @seealso \code{\link{getNetworks}}, \code{\link{runConsensusCluster}}
#'
#'
#' @returns a plot of the specified network
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 5, BPPARAM = bpparam())
#'
#' #construct the networks
#' DNEA <- getNetworks(object = DNEA)
#'
#' #identify metabolic modules via consensus clustering
#' DNEA <- runConsensusCluster(object = DNEA)
#'
#' #plot the networks
#' plotNetworks(object = DNEA, type = "group_networks")
#' plotNetworks(object = DNEA, type = "subnetworks", subnetwork = 1)
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
