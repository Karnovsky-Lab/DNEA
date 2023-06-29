
#' reduceFeatures() will collapse correlated features in the input non-normalized expression data
#'
#' reduceFeatures collapses correlated features in the input non-normalized input expression data using one of 3 methods:
#'
#' **1. correlation-based**
#'
#' **2. knowledge-based**
#'
#' **3. hybrid**
#'
#' Identified groups are collapsed by taking the mean expression of all features. More info about the different approaches
#' can be found in the **details** section.
#'
#' @details
#' The three collapsing methods have their strengths dependent on the dataset and prior information known about the features.
#' The following text explains more about each method and the best use cases:
#'
#'
#' **1. correlation-based:** The user specifies a correlation threshold wherein metabolites with a higher pearson correlation value will be
#' than the threshold are collapsed into one group. This approach is best when no prior information about the features is known.
#'
#'
#' **2. knowledge-based:** The user specifies feature groups based on a priori information (ie. all of the carnitines in a dataset
#' are specified as a single group) and the features within each group are collapsed into one feature. This approach is best in
#' experiments where the dataset contains many highly similar compounds, like Fatty acids, carnitines, ceramides, etc.
#'
#'
#' **3. hybrid:** The user specifies both a correlation threshold like in the correlation-based approach and feature groups
#'  based on a priori information similar to the knowledge-based approach. The features within each user-specified group that
#'  have a higher pearson correlation than the provided threshold are collapsed into one group. This approach is best in
#' experiments where the dataset contains many highly similar compounds and prevents poorly correlated features from being
#' collapsed into a single group.
#'
#'
#'
#' @param object A DNEAresults object
#' @param method A parameter that dictates the collapsing method to use. The options are
#'        as follows:
#'            1. correlation
#'            2. knowledge
#'            3. hybrid
#' @param correlation_threshold A threshold wherein features correlated above correlation_threshold
#'        will be grouped into one. This parameter is only necessary for the correlation and hybrid
#'        methods
#' @param metabolite_groups A dataframe containing group information for the collapsing algorithm
#'        indicated by the "knowledge" and "hybrid" methods
#' @returns A collapsed_DNEAresults object
#'
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @export
reduceFeatures <- function(object,
                           method = c("correlation",
                                      "knowledge",
                                      "hybrid"),
                           correlation_threshold = 0.9,
                           metabolite_groups = NULL){

  ################################################################
  #**Check that input data is correct and initialize parameters**#
  ################################################################

  #set method
  method <- match.arg(method)

  #Check to see if there is raw expression data provided - Feature reduction must be done on raw expression data
  if(is.null(expressionData(object, type = "input"))) stop(paste0('\n','FEATURE REDUCTION MUST BE DONE ON RAW EXPRESSION DATA!','\n',
                                                  'To proceed, please insert un-scaled expression data into the DNEAobject using the Expression(x)<- function', '\n'))

  #check to see if data looks scaled
  feature_mean <- colMeans(expressionData(object, type = "input"))
  if(all(feature_mean < 0.05 & feature_mean > -0.05)) warning(paste0("Data in expression_data assay looks to be scaled...",
                                                                     " Feature reduction must be done on raw data!"))

  #create dataframe input for node collapsing algorithm
  collapse_dat <- data.frame(samples = sampleNames(object),
                             groups = networkGroupIDs(object),
                             expressionData(object, type = "input"))
  #######################################
  #**perform node collapsing algorithm**#
  #######################################

  if (method == "correlation") {
    res <- collapseNodes_cor(dat = collapse_dat,
                             correlation_threshold = correlation_threshold)
  }else if (method == "knowedge") {
    res <- collapseNodes_knowledge(dat = collapse_dat,
                                   metabolite_groups = metabolite_groups)
  }else if (method == "hybrid") {
    res <- collapseNodes_hybrid(dat = collapse_dat,
                                correlation_threshold = correlation_threshold,
                                metabolite_groups = metabolite_groups)
  }


  message(paste0('The un-normalized data from the Expression Assay was used for feature reduction.','\n',
                                                              'The data in the NormalExpression Assay was replaced with log-scaled collapsed data.', '\n',
                                                              'If you prefer another normalization method replace this data prior to proceeding!', '\n\n',
                                                              '(orginal DNEAobject can be found in the original_experiment slot)', '\n'))

  ###################################################
  #**use new data to initialize reduced DNEAobject**#
  ###################################################

  #convert reduced data to numeric
  res[["collapsed_data"]] <-cbind.data.frame(res[["collapsed_data"]][,c(1,2)],
                                             apply(res[["collapsed_data"]][,-c(1,2)], 2, as.numeric))[,-1]

  #initialize new collapsed_DNEAresults object
  reduced_object <- createDNEAobject(project_name = projectName(object),
                                     expression_data = res[["collapsed_data"]],
                                     control = networkGroups(object)[[1]],
                                     case = networkGroups(object)[[2]])

  collapsed_object <- new("collapsed_DNEAresults",
                          project_name = projectName(reduced_object),
                          assays =  list(expression_data = expressionData(reduced_object, type = "input"),
                                         scaled_expression_data = expressionData(reduced_object, type = "normalized")),
                          metadata = list(samples = data.frame(samples = sampleNames(reduced_object),
                                                               conditions = networkGroupIDs(reduced_object),
                                                               row.names = sampleNames(reduced_object)),
                                          features = data.frame(feature_names = featureNames(reduced_object, original = TRUE),
                                                                clean_feature_names = featureNames(reduced_object, original = FALSE),
                                                                row.names = featureNames(reduced_object, original = TRUE)),
                                          network_group_IDs = networkGroupIDs(reduced_object),
                                          network_groups = networkGroups(reduced_object)),
                          joint_graph = make_empty_graph(n = numFeatures(reduced_object),
                                                         directed = TRUE),
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

# correlation-based node-collapsing
#' collapseNodes_cor will collapse nodes based on correlation
#'
#' more info about correlation
#' @param dat unscaled expression data
#' @param correlation_threshold correlation threshold
#'
#' @return collapsed data
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
collapseNodes_cor <- function(dat, correlation_threshold = 0.9) {


  ##set-up parameters and split data
  num_features <- ncol(dat)-2
  sample_groups <- unique(dat[,2])
  metabs <- colnames(dat)[-c(1:2)]

  input_dat <- list()
  input_dat[[1]] <- dat[dat[,2] == sample_groups[1],]
  input_dat[[2]] <- dat[dat[,2] == sample_groups[2],]

  #calculate pearson correlations and cluster
  cor_mat <- list()
  clust <- list()
  clust_groups <- list()
  for (a in 1:length(input_dat)) {
    cor_mat[[a]] <- cor(as.matrix(input_dat[[a]][,-c(1:2)]),
                        use = "pairwise.complete.obs",
                        method = "pearson")
    clust[[a]] <- hclust(as.dist(1-abs(cor_mat[[a]])))
    clust_groups[[a]] <- cutree(clust[[a]], h = 1-correlation_threshold)
  }

  #create adjacency matrix
  group_matrix <- matrix(0, num_features, num_features)
  colnames(group_matrix) <- rownames(group_matrix) <- metabs
  for (i in 1:(num_features-1)){
    for (j in (i+1):num_features){
      for (k in 1:length(clust_groups)){
        temp_group = clust_groups[[k]]
        group_matrix[i,j] <- group_matrix[i,j] + (length(unique(temp_group[c(i,j)]))==1)
      }
    }
  }


  #create adjacency graph
  adjacency_graph <- graph.adjacency(group_matrix, mode = "undirected", weighted = TRUE)

  #edge list from adjacency graph
  edge_list <- cbind.data.frame(get.edgelist(adjacency_graph), data.frame(weights = E(adjacency_graph)$weight))


  #if no correlation in either condition send error
  if(dim(edge_list)[1] == 0){

    stop("Metabolites in the data are not sufficiently correlated for collapsing...\n Try using a lower correlation coefficient threshold to collapse metabolites.")
  }

  #filter to features clustering in both sample conditions only
  filtered_edge_list <- edge_list[edge_list[, 3] == 2,c(1:2)]

  #if no correlation in shared by both conditions send error
  if(dim(filtered_edge_list)[1] == 0){

    stop("Metabolites in the data are not sufficiently correlated for collapsing...\n Try using a lower correlation coefficient threshold to collapse metabolites.")
  } else {

    #new graph from filtered edge list
    filtered_adjacency_graph <- graph_from_edgelist(as.matrix(filtered_edge_list[,c(1:2)]))
    graph_components <- components(filtered_adjacency_graph)$membership

    #combine new groups with independet features
    final_membership <- data.frame(feature = names(graph_components), feature_membership=graph_components)
    final_membership$feature_membership <- paste0("group",
                                                  final_membership$feature_membership)
    independent_features <- data.frame(feature = metabs[!(metabs %in% names(graph_components))],
                                    feature_membership=metabs[!(metabs %in% names(graph_components))])
    rownames(independent_features) <- independent_features$feature_membership
    final_membership <- rbind.data.frame(final_membership, independent_features)

    newdat <- NULL
    for (a in 1:length(input_dat)) {

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
      dat_by_cond <- cbind.data.frame(input_dat[[a]][,c(1:2)], dat_by_cond)

      #bind with data output
      newdat <- rbind(newdat, dat_by_cond)
    }


    return(list(feature_membership = final_membership,
                collapsed_data = newdat))
  }
}


#' collapseNodes_knowledge node-collapsing based on user-supplied metabolite groups
#'
#' more info about correlation
#'
#' @param dat unscaled expression data
#' @param metabolite_groups metabolite groupings
#'
#' @return collapsed data
#'
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
collapseNodes_knowledge <- function (dat,
                                     metabolite_groups) {

  metab_group <- NULL
  colnames(metabolite_groups) <- c("metabolite", "metab_group")
  sample_groups <- unique(dat[,2])

  input_dat <- list()
  input_dat[[1]] <- dat[dat[,2] == sample_groups[1],]
  input_dat[[2]] <- dat[dat[,2] == sample_groups[2],]

  newdat <- NULL
  for (a in 1:length(input_dat)) {

    #merge expression data with user-specified groups
    dat_by_cond <- cbind.data.frame(data.frame(metabolite = colnames(input_dat[[a]])[-c(1,2)]),
                                    as.data.frame(t(input_dat[[a]][,-c(1,2)])))
    dat_by_cond <- merge(metabolite_groups, dat_by_cond,
                         by = "metabolite")

    #collapse groups by mean
    dat_by_cond <- as.data.frame(dat_by_cond[,-1] %>% group_by(metab_group) %>%
                                   summarise(across(everything(), mean)))

    #reformat and merge with sample/sample condition info
    rownames(dat_by_cond) <- dat_by_cond$metab_group
    dat_by_cond <- t(dat_by_cond[,-1])
    dat_by_cond <- dat_by_cond[unlist(input_dat[[a]][1]), ]
    dat_by_cond <- cbind.data.frame(input_dat[[a]][,c(1:2)],dat_by_cond)

    #bind to output
    newdat <- rbind(newdat, dat_by_cond)
  }

  return(list(final_membership=metabolite_groups,
              collapsed_data = do.call("rbind",newdat)))
}




#' collapseNodes_hybrid node-collapsing based on correlations and user-supplied metabolite groups
#'
#' more info about correlation
#'
#' @param dat unscaled expression data
#' @param metabolite_groups metabolite groupings
#' @param correlation_threshold correlation threshold
#'
#' @return collapsed data
#' @import igraph
#' @importFrom dplyr %>% summarise across everything group_by select
#' @importFrom stringr str_detect
#' @importFrom stats hclust cutree as.dist
#' @keywords internal
collapseNodes_hybrid <- function (dat,
                                  metabolite_groups,
                                  correlation_threshold = 0.9) {

  colnames(metabolite_groups) <- c("metabolite", "metab_group")
  final_membership <- list()
  newdat <- list()

  for (feature_group in unique(metabolite_groups$metab_group)) {

    #filter data to knowledge-based group
    feature_group_membership <- metabolite_groups[metabolite_groups$metab_group == feature_group,]
    feature_group_data <- select(dat, c(1, 2, match(feature_group_membership$metabolite, names(dat))))
    if(dim(feature_group_membership)[1] == 1) {

      #features that are independent stay the same
      final_membership[[feature_group]] <- data.frame(feature=feature_group_membership$metab_group,
                                                      feature_membership=feature_group_membership$metab_group,
                                                      row.names=feature_group_membership$metabolite)
      newdat[[feature_group]] <- feature_group_data
    } else {

      #correlation-based collapsed within knowledge-based groups
      correlation_based_collapse <- tryCatch(expr = {
        collapseNodes_cor(feature_group_data,
                          correlation_threshold = correlation_threshold)
      }, error = function(e){
        stop(paste0("Features in ", feature_group, " are not sufficiently correlated for collapsing...\n Try using a lower correlation coefficient threshold to collapse features or maintain them as independent features"))
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
  newdat <- cbind.data.frame(dat[,c(1:2)], newdat)

  return(list(final_membership=final_membership,
              collapsed_data=newdat))
}
