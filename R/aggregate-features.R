#' Collapse correlated features into a single group
#'
#' @description
#' This function takes as input a DNEAobj object and collapses highly
#' correlated features within the non-normalized, non-transformed data
#' using one of three methods:
#'
#'   \enumerate{
#'   \item \strong{correlation-based}
#'   \item \strong{knowledge-based}
#'   \item \strong{hybrid}}
#'
#' More info about the different approaches can be found in the
#' \strong{\emph{Details}} section. Highly correlated groups of features
#' are collapsed by taking the mean expression of all features.
#'
#' \strong{\emph{NOTE:}} This method was developed using non-normalized,
#' non-transformed data and this fact is critical for feature collapsing
#' since the mean expression value is used for each group. Normalized data
#' may alter the results of collapsing.
#'
#' @param object A \code{\link{DNEAobj}} object
#' @param method A character string that dictates the collapsing method
#' to use. The available methods are: "correlation", "knowledge", or "hybrid"
#' @param correlation_threshold A threshold wherein features correlated
#' above correlation_threshold will be collapsed into one. This parameter
#' is only necessary for the correlation and hybrid methods
#' @param feature_groups A data.frame containing group information for
#' the collapsing algorithm indicated by the "knowledge" and
#' "hybrid" methods
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}},
#' \code{\link{stabilitySelection}}
#'
#' @details
#' Due to the computational complexity of the DNEAdev algorithm,
#' the processing time for a given dataset increases dramatically as the
#' number of features increases. The ability to process each replicate
#' performed in  \code{\link{stabilitySelection}} in parallel helps
#' circumvent this issue, however, a user may still be constrained by
#' the resources available to them (ie. a limited number of cpu cores
#' or memory). Collapsing related features into groups is another
#' method by which the user can reduce the complexity of the
#' analysis, and as a result decrease the necessary resources.  \cr
#'
#' In a related scenario, you may also have many
#' highly-correlated features of the same class of compounds (ie.
#' fatty acids, carnitines, etc.), and network analysis at the
#' resolution of these individual features is not important.
#' Collapsing features would decrease the computing time without
#' losing critical information to the analysis (Please see the
#' \strong{\emph{Details}} section of
#' \code{\link{createDNEAobject}} for more information about
#' the motivation behind collapsing highly correlated features). \cr
#'
#' Ultimately, this function allows the user to reduce the complexity
#' of the dataset and reduce the computational power necessary for the
#' analysis and/or improve the quality of the results. The most appropriate
#' method to use when collapsing data is dependent on the dataset and
#' prior information known about the features. The following text
#' explains more about each method and the best use cases:
#'
#'   \enumerate{
#'   \item \strong{correlation-based - } The user specifies a correlation
#'   threshold wherein features with a higher pearson correlation value
#'   than the threshold are collapsed into one group. This approach is
#'   best when no prior information about the features are known.
#'
#'
#'   \item \strong{knowledge-based - } The user specifies feature groups
#'   based on a priori information (ie. all of the carnitines in a dataset
#'   are specified as a single group) and the features within each group
#'   are collapsed into one feature. This approach is best in
#'   experiments where the dataset contains many highly similar compounds,
#'   like fatty acids, carnitines, ceramides, etc.
#'
#'
#'   \item \strong{hybrid - } The user specifies both a correlation
#'   threshold, like in the correlation-based approach, and feature
#'   groups based on \emph{a priori} information similar to the
#'   knowledge-based approach. The features within each user-specified
#'   group that have a higher pearson correlation than the provided
#'   threshold are collapsed into one group. This approach is best
#'   in experiments where the dataset contains many compounds of a
#'   similar class, but the user is unsure how correlated the
#'   features of said class will be. This method prevents poorly
#'   correlated or uncorrelated features from being collapsed
#'   into a single feature.}
#'
#' @returns A \code{collapsed_DNEAobj} object.
#'
#' @examples
#' #import example data
#' data(dnw)
#'
#' #simulate group labels
#' TEDDY_groups <- data.frame(features = rownames(expressionData(x = dnw,
#'                                                assay = "input_data")),
#'                            groups = rownames(expressionData(x = dnw,
#'                                              assay = "input_data")),
#'                            row.names = rownames(expressionData(x = dnw,
#'                                                 assay = "input_data")))
#'
#' TEDDY_groups$groups[TEDDY_groups$groups %in% c("isoleucine",
#'                                                "leucine",
#'                                                "valine")] <- "BCAAs"
#' TEDDY_groups$groups[grep("acid", TEDDY_groups$groups)] <- "fatty_acids"
#'
#'
#' collapsed_TEDDY <- aggregateFeatures(object = dnw,
#'                                   method = "hybrid",
#'                                   correlation_threshold = 0.7,
#'                                   feature_groups = TEDDY_groups)
#'
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @export
aggregateFeatures <- function(object,
                              method = c("correlation",
                                         "knowledge",
                                         "hybrid"),
                              correlation_threshold = NULL,
                              feature_groups = NULL){

  ##set method
  method <- match.arg(method)

  ##bind local variables to function
  feature_membership <- NULL

  ##Check to see if there is non-normalized, non-transformed expression data
  ##provided - Feature reduction requires this
  if(is.null(expressionData(x = object, assay = "input_data"))){
    stop("AGGREGATION MUST BE DONE ON RAW PEAK INTENSITIES/CONCENTRATIONS!",
         '\nTo proceed, please insert un-scaled expression data into the DNEAobject')
  }

  #check that rownames of feature_groups match feature order
  if(!(all(rownames(expressionData(x = object, assay = "input_data")) ==
           rownames(feature_groups)))){
    stop("The feature order of feature_groups does not match the expression data!")
  }

  #warning if correlation_threshold was provided for knowledge-based collapsing
  if(method == "knowledge" & !is.null(correlation_threshold)){
    warning("correlation_threshold is only used in the correlation-based and hybrid node collapsing methods...\n",
            "...The threshold will be ignored and nodes will be collapsed based on the provided feature_groups.\n")
  }

  ##create dataframe input for node collapsing algorithm
  collapse_dat <- data.frame(samples = sampleNames(object),
                             groups = networkGroupIDs(object),
                             t(expressionData(x = object, assay = "input_data")))

  ##collapse features based on specified approach
  if (method == "correlation") {
    res <- collapseNodes_cor(dat = collapse_dat,
                             correlation_threshold = correlation_threshold)
  }else if (method == "knowledge") {
    res <- collapseNodes_knowledge(dat = collapse_dat,
                                   feature_groups = feature_groups)
  }else if (method == "hybrid") {
    res <- collapseNodes_hybrid(dat = collapse_dat,
                                correlation_threshold = correlation_threshold,
                                feature_groups = feature_groups)
  }

  message("The raw peak intensity data was used for aggregation\n",
          "The aggregated log-scaled data is in the @assay slot\n",
          "(The orginal DNEAobj can be found in the @original_experiment slot)")

  #new input
  new_dat <- t(res[["collapsed_data"]][,-c(1,2)])
  new_group_labels <- res[["collapsed_data"]][["groups"]]
  names(new_group_labels) <- res[["collapsed_data"]][["samples"]]

  #initialize new collapsed_DNEAobj object
  reduced_object <- createDNEAobject(project_name = projectName(object),
                                     expression_data = new_dat,
                                     group_labels = new_group_labels)

  collapsed_object <- new("collapsed_DNEAobj",
                          project_name = projectName(reduced_object),
                          assays =  list(input_data = expressionData(x = reduced_object, assay = "input_data"),
                                         log_input_data = expressionData(x = reduced_object, assay = "log_input_data"),
                                         scaled_expression_data = expressionData(x = reduced_object, assay = "scaled_expression_data")),
                          metadata = list(samples = data.frame(samples = sampleNames(reduced_object),
                                                               conditions = networkGroupIDs(reduced_object),
                                                               row.names = sampleNames(reduced_object)),
                                          features = data.frame(feature_names = featureNames(reduced_object, original = TRUE),
                                                                clean_feature_names = featureNames(reduced_object, original = FALSE),
                                                                row.names = featureNames(reduced_object, original = TRUE)),
                                          network_group_IDs = networkGroupIDs(reduced_object),
                                          network_groups = networkGroups(reduced_object)),
                          hyperparameter = list(BIC_scores = NULL, optimized_lambda = NULL, tested_lambda_values = NULL),
                          adjacency_matrix = list(weighted_adjacency = NULL, unweighted_adjacency = NULL),
                          stable_networks = list(selection_results = NULL, selection_probabilities = NULL),
                          original_experiment = reduced_object,
                          feature_membership = res[["feature_membership"]])

  #add diagnostics
  datasetSummary(collapsed_object) <- new("DNEAinputSummary",
                                          num_samples = numSamples(reduced_object),
                                          num_features = numFeatures(reduced_object),
                                          diagnostic_values = diagnostics(reduced_object))

  ##add node list
  nodeList(collapsed_object) <- nodeList(reduced_object)
  return(collapsed_object)

}

################################################################################
# main functions #

#' correlation-based node-collapsing
#'
#' This function calculates a pearson correlation matrix for the features
#' within the input dataset and then clusters and collapses features more
#' highly-correlated than the user-specified cut-off. The mean expression of
#' the features in each group is calculated as the new expression value for
#' said group. *NOTE:* This method was developed using non-normalized,
#' non-transformed data and this fact is critical for feature collapsing
#' since the mean expression value is used for each group. Normalized data
#' may alter the results of collapsing.
#'
#' @param dat non-normalized, non-transformed expression data
#' @param correlation_threshold correlation threshold wherein
#' features with a higher pearson correlation will be collapsed
#' into a single node.
#'
#' @author Gayatri Iyer
#'
#' @seealso
#' \code{\link{aggregateFeatures}}
#'
#' @returns The collapsed expression data and a group table
#' indicating which features were collapsed into a given group.
#'
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
#' @noRd
collapseNodes_cor <- function(dat, correlation_threshold = 0.9) {

  ##bind local variables to function
  feature_membership <- NULL

  ##set-up parameters and split data
  num_features <- ncol(dat)-2
  sample_groups <- unique(dat[,2])
  metabs <- colnames(dat)[-c(1,2)]

  input_dat <- list()
  input_dat[[1]] <- dat[dat[,2] == sample_groups[1],]
  input_dat[[2]] <- dat[dat[,2] == sample_groups[2],]

  #calculate pearson correlations and cluster
  cor_mat <- list()
  clust <- list()
  clust_groups <- list()
  for (a in seq(1, length(input_dat))) {
    cor_mat[[a]] <- cor(as.matrix(input_dat[[a]][,-c(1, 2)]),
                        use = "pairwise.complete.obs",
                        method = "pearson")
    clust[[a]] <- hclust(as.dist(1-abs(cor_mat[[a]])))
    clust_groups[[a]] <- cutree(clust[[a]], h = 1-correlation_threshold)
  }

  #create adjacency matrix
  group_matrix <- matrix(0, num_features, num_features)
  colnames(group_matrix) <- rownames(group_matrix) <- metabs
  for(i in seq(1, (num_features-1))){
    for(j in seq((i+1), num_features)){
      for(k in seq(1, length(clust_groups))){
        temp_group <- clust_groups[[k]]
        group_matrix[i,j] <- group_matrix[i,j] + (length(unique(temp_group[c(i,j)])) == 1)
      }
    }
  }


  #create adjacency graph
  adjacency_graph <- graph.adjacency(group_matrix, mode = "max", weighted = TRUE)

  #edge list from adjacency graph
  edge_list <- cbind.data.frame(get.edgelist(adjacency_graph),
                                data.frame(weights = E(adjacency_graph)$weight))


  #if no correlation in either condition send error
  if(dim(edge_list)[1] == 0){

    stop("Features in teh data are not sufficientl correlated for aggregating...",
         "\nTry using a lower correlation to aggregate features")
  }

  #filter to features clustering in both sample conditions only
  filtered_edge_list <- edge_list[edge_list[, 3] == 2,c(1, 2)]

  #if no correlation in shared by both conditions send error
  if(dim(filtered_edge_list)[1] == 0){

    stop("Features in teh data are not sufficientl correlated for aggregating...",
         "\nTry using a lower correlation to aggregate features")
  } else {

    #new graph from filtered edge list
    filtered_adjacency_graph <- graph_from_edgelist(as.matrix(filtered_edge_list[,c(1, 2)]))
    graph_components <- components(filtered_adjacency_graph)$membership

    #combine new groups with independet features
    final_membership <- data.frame(feature = names(graph_components),
                                   feature_membership=graph_components)
    final_membership$feature_membership <- paste0("group",final_membership$feature_membership)
    independent_features <- data.frame(feature = metabs[!(metabs %in% names(graph_components))],
                                    feature_membership=metabs[!(metabs %in% names(graph_components))])
    rownames(independent_features) <- independent_features$feature_membership
    final_membership <- rbind.data.frame(final_membership, independent_features)

    newdat <- NULL
    for (a in seq(1, length(input_dat))) {

      #reformat input data by condition
      dat_by_cond <- cbind.data.frame(data.frame(feature = colnames(input_dat[[a]])[-c(1,2)]),
                                      data.frame(t(input_dat[[a]][,-c(1,2)])))

      #merge with group membership info
      dat_by_cond <- merge(final_membership, dat_by_cond, by = "feature", all.y = TRUE)

      #collapse expression data within groups to singular value per condition
      dat_by_cond <- as.data.frame(dat_by_cond[,-1] %>% group_by(feature_membership) %>%
                                     summarise(across(everything(), mean)))

      #reformat result
      rownames(dat_by_cond) <- dat_by_cond$feature_membership
      dat_by_cond <- t(dat_by_cond[,-1])
      dat_by_cond <- dat_by_cond[unlist(input_dat[[a]][1]), ]
      dat_by_cond <- cbind.data.frame(input_dat[[a]][,c(1, 2)], dat_by_cond)

      #bind with data output
      newdat <- rbind(newdat, dat_by_cond)
    }


    return(list(feature_membership = final_membership,
                collapsed_data = newdat))
  }
}

#' Node-collapsing based on user-supplied metabolite groups
#'
#' This function collapses related features based on user-specified
#' groups via the feature_groups input parameter. The mean expression value
#' for each group is used as the new value. *NOTE:* This method was developed
#' using non-normalized, non-transformed data and this fact is critical for
#' feature collapsing since the mean expression value is used for each group.
#' Normalized data may alter the results of collapsing.
#'
#' @param dat unscaled expression data
#' @param feature_groups metabolite groupings
#'
#' @author Gayatri Iyer
#'
#' @seealso
#' \code{\link{aggregateFeatures}}
#'
#' @returns The collapsed expression data and a group table indicating
#' which features were collapsed into a given group.
#'
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
#' @noRd
collapseNodes_knowledge <- function (dat,
                                     feature_groups) {

  metab_group <- NULL
  colnames(feature_groups) <- c("metabolite", "metab_group")
  sample_groups <- unique(dat[,2])

  input_dat <- list()
  input_dat[[1]] <- dat[dat[,2] == sample_groups[1],]
  input_dat[[2]] <- dat[dat[,2] == sample_groups[2],]

  newdat <- NULL
  for (a in seq(1, length(input_dat))) {

    #merge expression data with user-specified groups
    dat_by_cond <- cbind.data.frame(data.frame(metabolite = colnames(input_dat[[a]])[-c(1,2)]),
                                    as.data.frame(t(input_dat[[a]][,-c(1,2)])))
    dat_by_cond <- merge(feature_groups, dat_by_cond,
                         by = "metabolite")

    #collapse groups by mean
    dat_by_cond <- as.data.frame(dat_by_cond[,-1] %>% group_by(metab_group) %>%
                                   summarise(across(everything(), mean)))

    #reformat and merge with sample/sample condition info
    rownames(dat_by_cond) <- dat_by_cond$metab_group
    dat_by_cond <- t(dat_by_cond[,-1])
    dat_by_cond <- dat_by_cond[unlist(input_dat[[a]][1]), ]
    dat_by_cond <- cbind.data.frame(input_dat[[a]][,c(1, 2)], dat_by_cond)

    #bind to output
    newdat <- rbind(newdat, dat_by_cond)
  }

  return(list(feature_membership = feature_groups,
              collapsed_data = newdat))
}

#' Node-collapsing based on correlations and user-supplied metabolite groups
#'
#' This function is a hybrid of \code{\link{collapseNodes_cor}} and
#' \code{\link{collapseNodes_knowlege}}. It initially groups metabolites
#' by user-specified groups via the feature_groups parameter. Then
#' \code{\link{collapseNodes_cor}} is called on each group to identify
#' sub-groups via the correlation based collapsing approach. More details on
#' this method can be found at \code{\link{collapseNodes_cor}}. The mean of
#' the features in each group is calculated as the new expression value for said
#' group. *NOTE:* This method was developed using non-normalized,
#' non-transformed data and this fact is critical for feature collapsing since
#' the mean expression value is used for each group. Normalized data may
#' alter the results of collapsing.
#'
#' @param dat unscaled expression data
#' @param feature_groups metabolite groupings
#' @param correlation_threshold correlation threshold
#'
#' @author Gayatri Iyer
#'
#' @seealso
#' \code{\link{collapseNodes_cor}},\code{\link{collapseNodes_knowlege}},
#' \code{\link{aggregateFeatures}}
#'
#' @returns The collapsed expression data and a group table indicating
#' which features were collapsed into a given group.
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
#' @noRd
collapseNodes_hybrid <- function (dat,
                                  feature_groups,
                                  correlation_threshold = 0.9) {

  ##bind local variables to function
  feature_membership <- NULL

  ##set parameters
  colnames(feature_groups) <- c("metabolite", "metab_group")
  final_membership <- list()
  newdat <- list()

  for (feature_group in unique(feature_groups$metab_group)) {

    #filter data to knowledge-based group
    feature_group_membership <- feature_groups[feature_groups$metab_group == feature_group,]
    feature_group_data <- select(dat, c(1, 2, match(feature_group_membership$metabolite, names(dat))))
    if(dim(feature_group_membership)[1] == 1) {

      #features that are independent stay the same
      final_membership[[feature_group]] <- data.frame(feature=feature_group_membership$metab_group,
                                                      feature_membership=feature_group_membership$metab_group,
                                                      row.names=feature_group_membership$metabolite)
      newdat[[feature_group]] <- feature_group_data
    } else {

      #correlation-based collapsed within knowledge-based groups
      correlation_based_collapse <- tryCatch(expr={
        collapseNodes_cor(feature_group_data,
                          correlation_threshold = correlation_threshold)
      }, error=function(e){

        stop("Features in ", feature_group, " are not sufficiently correlated...",
             "\nTry using a lower coefficient to aggregate features")
      })

      #add new data to final_membership and rename groups
      final_membership[[feature_group]] <- data.frame(correlation_based_collapse$feature_membership)
      new_membership_groups <- str_detect(final_membership[[feature_group]]$feature_membership, "group")
      final_membership[[feature_group]]$feature_membership[new_membership_groups] <- paste0("knowledge: ", feature_group, "-correlation: ",
                                                                                            final_membership[[feature_group]]$feature_membership[new_membership_groups])

      #add new data to output expression data and rename groups
      newdat[[feature_group]] <- correlation_based_collapse$collapsed_data
      new_dat_groups <- str_detect(colnames(newdat[[feature_group]]), "group[[:digit:]]")
      colnames(newdat[[feature_group]])[new_dat_groups] <- paste0("knowledge: ", feature_group, "-correlation: ",
                                                        colnames(newdat[[feature_group]])[new_dat_groups])

    }
  }

  #format final_membership for output
  names(final_membership) <- NULL
  final_membership <- do.call("rbind", final_membership)

  #format newdat for output
  names(newdat) <- NULL
  newdat <- do.call("cbind", lapply(newdat, function(x) x[-c(1,2)]))
  newdat <- cbind.data.frame(dat[,c(1, 2)], newdat)

  return(list(feature_membership = final_membership,
              collapsed_data = newdat))
}
