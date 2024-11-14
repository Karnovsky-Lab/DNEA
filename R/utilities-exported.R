
#' Add additional metadata to the DNEAobj object
#'
#' This function will take additional metadata and add it to the specified
#' data frame in the metadata slot. \strong{\emph{NOTE:}} The row names
#' of the new metadata must match the order of the input sample names
#' or feature names, respectively.
#'
#' @param object A \code{\link{DNEAobj}} object.
#'
#' @param type A character string corresponding to the type of metadata
#' being included. Can be either "samples" or "features"
#'
#' @param metadata a data frame containing metadata to add. The row
#' names should be either the sample names or feature names, respectively
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{featureNames}},\code{\link{sampleNames}},
#' \code{\link{metaData}}
#'
#' @returns A \code{\link{DNEAobj}} object with the specified additions.
#'
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEAobj
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' #make sure metadata has same sample order as DNEAobj object
#' T1Dmeta <- T1Dmeta[sampleNames(dnw), ]
#'
#' #add new metadata to DNEAobj object
#' dnw <- includeMetadata(object=dnw, type="samples", metadata=T1Dmeta)
#'
#' @export
includeMetadata <- function(object,
                            type=c('samples', 'features'),
                            metadata){

  ##test for proper input
  if(!type %in% c('samples', 'features')) stop('type can be "samples" or "features"')
  if(!inherits(object, "DNEAobj")) stop('the input object should be of class "DNEAobj"!')
  if(!any(inherits(metadata, "matrix") | inherits(metadata, "data.frame"))) {
    stop('the input metadata should be of class "matrix" or "data.frame"!')
  }
  if(is.null(colnames(metadata))) stop("dat must have column names!")
  if(is.null(rownames(metadata))) stop("dat must have row names!")

  type <- match.arg(type)
  if(vector_compare(rownames(metaData(x=object, type=type)),
                    rownames(metadata))){
    new_metadata <- cbind(metaData(object, type=type),
                          metadata)
    metaData(object, type=type) <- new_metadata
  }else{
    stop('new metadata order does not match ', type, ' order!')
  }
  return(object)
}

#' Include custom normalized data in the DNEAobj object
#'
#' This function allows the user to input custom-normalized data into the
#' \code{\link{DNEAobj}} object for use in DNEA analysis.
#'
#' @param object A \code{\link{DNEAobj}} object.
#' @param dat A list of \emph{m x n} numeric matrices of
#' custom-normalized expression data, one matrix for each
#' experimental condition. The list elements should be
#' labeled for their respective condition. These should
#' match the labels returned by
#' \code{\link{networkGroups}}.
#'
#' @param assay_name A character string corresponding
#' to the name the new data will be stored under in the
#' assays slot of the \code{\link{DNEAobj}}.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{DNEAobj}},
#'
#' @returns A \code{\link{DNEAobj}} object with the added
#' expression data in the @@assays slot
#'
#' @examples
#' #load example data
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEAobj
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                         group_labels = group_labels)
#'
#' #transpose TEDDY data
#' TEDDY <- t(log(TEDDY))
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[rownames(TEDDY),]
#'
#' dat <- list()
#' for(cond in networkGroups(dnw)){
#'   dat[[cond]] <- TEDDY[names(group_labels)[group_labels == cond],]
#' }
#'
#'
#' #log-transform and median center the expression data without scaling
#' newdat <- list()
#' for(cond in networkGroups(dnw)){
#'
#'   group_dat <- dat[[cond]]
#'   for(i in seq(1, ncol(group_dat))){
#'    metab_median=median(group_dat[, i], na.rm=TRUE)
#'     metab_range=range(group_dat[, i], na.rm=TRUE)
#'     scale_factor=max(abs(metab_range - metab_median))
#'     group_dat[, i] <- (group_dat[, i] - metab_median) / scale_factor
#'
#'     rm(metab_median, metab_range, scale_factor)
#'   }
#'
#'   group_dat <- t(group_dat)
#'   newdat <- append(newdat, list(group_dat))
#'
#'   rm(i, group_dat)
#' }
#'
#' #add names
#' names(newdat) <- names(dat)
#'
#' #add data
#' dnw <- addExpressionData(object=dnw,
#'                          dat=newdat,
#'                          assay_name="median_scaled_data")
#'
#' @rdname addExpressionData
#' @export
addExpressionData <- function(object,
                              dat,
                              assay_name){
  if(!inherits(object, "DNEAobj")) stop('the input object should be of class "DNEAobj"!')
  feature_names <- featureNames(object)

  for(i in seq(length(dat))){
    tmp_samples <- names(networkGroupIDs(object)[networkGroupIDs(object)==networkGroups(object)[i]])
    table_metadata_check(dat[[i]],
                         sample_names = tmp_samples,
                         feature_names = feature_names)
  }

  assays(object)[[assay_name]] <- dat

  validObject(object)
  return(object)
}


#' Save network information to .csv files
#'
#' This function will save the node and edge information as .csv files
#' to the specified directory. The files are already formatted for
#' input into Cytoscape.
#'
#' @param object A \code{\link{DNEAobj}} object.
#' @param file_path The file path to save the node and edge lists
#' to. If **NULL**, the files will be saved to the working directory.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{edgeList}},\code{\link{nodeList}}
#'
#' @returns Two .csv files, one for the node list and one for
#' the edge list, saved to the specified file path
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #filepath wherein to save the networks files
#' filepath <- tempdir()
#' #save node and edge list for input to cytoscape
#' getNetworkFiles(dnw, file_path=filepath)
#'
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object,
                            file_path=getwd()){

  ##test for proper input
  if(!inherits(object, "DNEAobj")) {
    stop('the input object should be of class "DNEAobj"!')
  }
  if(!is.character(file_path)) {
    stop('"file_path" should be a character string file path!')
  }
  #save node list
  write.csv(nodeList(object), paste0(file_path, "/", object@project_name,'_nodelist_',Sys.Date(),'.csv'),
            row.names=FALSE, quote=FALSE)

  #save edge list
  write.csv(edgeList(object), paste0(file_path, "/", object@project_name,'_edgelist_',Sys.Date(),'.csv'),
            row.names=FALSE, quote=FALSE)

  return(object)
}
################################################################################
#' Visualize the biological networks
#'
#' This function plots the total network, condition networks, or
#' sub networks as specified by the user. Purple nodes are differential
#' features, green indicates edges specific to group 1, and red indicates
#' edges specific to group 2.
#'
#' @param object A \code{\link{DNEAobj}} object.
#'
#' @param type There are two possible arguments to \strong{type}:
#' \emph{"group_networks"} specifies the whole network or condition networks.
#' \emph{"sub_networks"} specifies that one of the sub networks should
#' be plotted.
#' *Additional input via the \strong{subtype} parameter is required.*
#'
#' @param subtype There are several possible arguments to \strong{subtype}.
#' If \emph{type == "group_networks"}, \strong{subtype} can be \emph{"All"}
#' to plot the whole network (ie. both conditions in the data returned by
#' \code{\link{networkGroups}}), or one of the condition network names to
#' plot the network corresponding to that condition. If
#' \emph{type == "sub_networks"}, \strong{subtype} should be a single-value
#' numeric vector corresponding to the sub network to plot.
#'
#' @param layout_func The layout in which to plot the specified network.
#' Please see \code{\link[igraph]{plot.igraph}} for more information.
#'
#' @param main A character string to use as the plot title.
#'
#' @param node_size The size of the nodes in the plot. The default is 15.
#' Please see \emph{vertex.size} parameter in
#' \code{\link[igraph]{tkplot}} for more details.
#'
#' @param edge_width The width of the edges in the plot. The default is 1.
#' Please see \emph{width} parameter in
#' \code{\link[igraph]{tkplot}} for more details.
#'
#' @param label_size The size of the node labels in the plot.
#' The default is 1. Please see \emph{label.size} in
#' \code{\link[igraph]{tkplot}} for more details.
#'
#' @param label_font Specifies the font type to use in the plot.
#' 1 is normal font, 2 is bold-type, 3 is italic-type, 4 is bold- and
#' italic-type. Please see the \emph{label.font} parameter in
#' \code{\link[igraph]{plot.igraph}} for more details.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{getNetworks}},\code{\link{clusterNet}},
#' \code{\link{networkGroups}}
#'
#'
#' @returns A plot of the specified network
#'
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #plot the networks
#' plotNetworks(object=dnw, type="group_networks", subtype="All")
#' plotNetworks(object=dnw, type="sub_networks", subtype=1)
#'
#' @import igraph
#' @export
plotNetworks <- function(object,
                         type=c("group_networks", "sub_networks"),
                         subtype="All",
                         layout_func,
                         main="",
                         node_size=15,
                         edge_width=1,
                         label_size=1,
                         label_font=1){

  ##test for proper input
  if(!inherits(object, "DNEAobj")) {
    stop('the input object should be of class "DNEAobj"!')
  }

  #set up parameters
  type <- match.arg(type)
  network_graph <- adjacencyGraph(object, "joint_graph")
  node_list <- nodeList(object)

  if(type == "group_networks"){

    edge_list <- edgeList(object)
    if(subtype == "All"){

      subtype_network <- induced_subgraph(network_graph, V(network_graph)$name)
    }else{

      ##check the subtype exists
      if(!(subtype %in% c(networkGroups(object), "Both"))) {
        stop("The subtype provided does not match a network within the data!")
      }
      ##grab metabolite names present in specified subtype
      subgroup_nodes <- unique(c(edge_list$Metabolite.A[edge_list$edge == subtype | edge_list$edge == "Both"],
                                 edge_list$Metabolite.B[edge_list$edge == subtype | edge_list$edge == "Both"]))
      ##threshold network_graph network
      subtype_network <- induced_subgraph(network_graph, V(network_graph)$name[match(subgroup_nodes, V(network_graph)$name)])
    }

  }else if(type == "sub_networks"){

    ##check that subnetwork given is relevant
    if(all(is.na(match(subtype, node_list$membership)))) {
      stop("The subnetwork specified does not exist!",
      "\nPlease specify a value contained in the membership column of the node list")
    }
    ##threshold network_graph network
    subtype_network <- induced_subgraph(network_graph, V(network_graph)$name[node_list$membership == subtype])
  }
  ##set graph layout
  if(missing(layout_func)){

    graph_layout <- layout_with_fr(graph=subtype_network)
  }else{

    graph_layout <- layout_func(graph=subtype_network)
  }

  ##plot the specified network
  plot(subtype_network, main=main, layout=graph_layout,
       vertex.size=node_size, width=edge_width,
       vertex.label=V(subtype_network)$name,
       vertex.label.cex=label_size,
       vertex.label.font=label_font)

}

filterNetworks.DNEAobj <- function(data,
                                   pcor,
                                   top_percent_edges){

  ##test for proper input
  if(!inherits(data, "DNEAobj")) {
    stop('the input object should be of class "DNEAobj"!')
  }

  ##grab adjacency matrices
  weighted_adjacency_matrices <- adjacencyMatrix(data, weighted=TRUE)

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length=length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  ##can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){
    if(pcor <=0 | pcor >= 1) {
      stop("the specified partial correlation threshold should be between 0 and 1!")
    }

    ##filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    ##check for valid input
    if(top_percent_edges <=0 | top_percent_edges >= 1) {
      stop("the specified % threshold should be between 0 and 1!")
    }

    ##filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs=(1 - top_percent_edges))
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    ##must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  ##store the adjacency matrices in DNEAobj object
  adjacencyMatrix(x=data, weighted=TRUE) <- weighted_adjacency_matrices
  adjacencyMatrix(x=data, weighted=FALSE) <- unweighted_adjacency_matrices

  ##update edge list and create output
  pairs <- combn(as.character(featureNames(data)), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names=FALSE)

  ##concatenate results into dataframe
  edge_list[,c(1,2)] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs))
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 == 0] <- names(weighted_adjacency_matrices)[[1]]
  edge_list$edge[edge_list$pcor.0 == 0 & edge_list$pcor.1 != 0] <- names(weighted_adjacency_matrices)[[2]]
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 != 0] <- "Both"
  edge_list <- edge_list[!is.na(edge_list$edge),]
  rownames(edge_list) <- NULL

  ##replace edgelist
  edgeList(data) <- edge_list


  ##inform user of remaining edges
  ##control network
  message(names(unweighted_adjacency_matrices)[[1]], " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[1]])/2) - sum(edge_list$edge == "Both"))

  ##case network
  message(names(unweighted_adjacency_matrices)[[2]],
          " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[2]])/2) - sum(edge_list$edge == "Both"))
  ##shared edges
  message(rep("-", 35))
  message("Number of edges shared by both networks: ",
          sum(edge_list$edge == "Both"))
  message("Total number of edges in dataset: ",
          nrow(edge_list))

  return(data)
}

#' Filter the adjacency matrices to only the edges that meet the
#' filter conditions
#'
#' @description
#' This function takes as input a \code{\link{DNEAobj}} object and
#' allows the user to filter the network edges by one of two methods:
#'
#' \enumerate{
#' \item \strong{Partial Correlation} - The networks can be filtered to
#' only include edges greater than or equal to a specified partial
#' correlation (pcor) value.
#' \item \strong{Top \emph{X%} of edges} - The networks can be filtered
#' to only include the strongest \emph{X%} of edges determined
#' by their partial correlation values.}
#'
#' Filtering is performed on the case and control adjacency matrices
#' separately. \cr
#'
#' @param data A \code{\link{DNEAobj}} object.
#'
#' @param pcor A partial correlation value of which to threshold
#' the adjacency matrices. Edges with pcor values <= to this
#' value will be removed.
#' @param top_percent_edges A value between 0-1 that corresponds to
#' the top x% edges to keep in each network, respectively
#' (i.e. top_percent_edges = 0.1 will keep only the top 10%
#' strongest edges in the networks).
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{getNetworks}},\code{\link{adjacencyMatrix}}
#'
#' @returns The input object after filtering the egdes in the
#' network according to the specified parameters.
#'
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #filter the networks by a correlation threshold of 0.166
#' dnw <- filterNetworks(dnw, pcor=0.166)
#'
#' #filter networks for the top 40% strongest correlations
#' dnw <- filterNetworks(dnw, top_percent_edges=0.4)
#'
#' @importFrom stats quantile
#' @rdname filterNetworks-methods
#' @aliases filterNetworks
#' @export
setMethod("filterNetworks", signature(data="DNEAobj"), filterNetworks.DNEAobj)

filterNetworks.list <- function(data,
                                pcor,
                                top_percent_edges){

  ##rename input variable
  weighted_adjacency_matrices <- data

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length=length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  #can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){

    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs=(1 - top_percent_edges))
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    ##must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  ##create edge list and output
  pairs <- combn(as.character(colnames(weighted_adjacency_matrices[[1]])), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names=FALSE)

  ##concatenate results into dataframe
  edge_list[,c(1,2)] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs)) #non-edge
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 == 0] <- names(weighted_adjacency_matrices)[[1]] #control
  edge_list$edge[edge_list$pcor.0 == 0 & edge_list$pcor.1 != 0] <- names(weighted_adjacency_matrices)[[2]] #case
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 != 0] <- "Both" #Both
  edge_list <- edge_list[!is.na(edge_list$edge),] #remove non-edges
  rownames(edge_list) <- NULL

  ##inform user of remaining edges
  ##control network
  message(names(unweighted_adjacency_matrices)[[1]], " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[1]])/2) - sum(edge_list$edge == "Both"))

  ##case network
  message(names(unweighted_adjacency_matrices)[[2]],
          " network specific edges: ",
          (sum(unweighted_adjacency_matrices[[2]])/2) - sum(edge_list$edge == "Both"))

  ##shared edges
  message(rep("-", 35))
  message("Number of edges shared by both networks: ",
          sum(edge_list$edge == "Both"))
  message("Total number of edges in dataset: ",
          nrow(edge_list))

  return(list(weighted=weighted_adjacency_matrices,
              unweighted=unweighted_adjacency_matrices,
              edge_list=edge_list))
}
#'
#' @rdname filterNetworks-methods
#' @aliases filterNetworks
setMethod("filterNetworks", signature(data="list"), filterNetworks.list)
