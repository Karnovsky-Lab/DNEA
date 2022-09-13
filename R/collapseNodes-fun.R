# wrapper function
#'ReduceFeatures
#'@export
#'@import igraph
#'@import dplyr
#'@import stringr
ReduceFeatures <- function(object,
                          method = c("correlation",
                                     "knowledge",
                                     "hybrid"),
                          corr.threshold = 0.9,
                          metabolite.groups = NULL){

  #Check to see if there is raw expression data provided - Feature reduction must be done on raw expression data
  if(is.null(Expression(object))) stop(paste0('\n','FEATURE REDUCTION MUST BE DONE ON RAW EXPRESSION DATA!','\n',
                                              'To proceed, please insert un-scaled expression data into the DNEAobject using the Expression(x)<- function', '\n'))

  feature_mean_expression <- apply(Expression(object),2,mean)
  feature_sd_expression <- apply(Expression(object), 2, sd)
  if(all(feature_mean_expression < 0.05 & feature_mean_expression > -0.05) |
         all(feature_sd_expression < 1.05 & feature_sd_expression > 0.95)) warning("Data in Expression Assay looks to be normalized. Feature reduction must be done on raw data!")
  if(is.null(Expression(object))) stop("Feature reduction must be done on raw data")

  dat <- data.frame(sample = sampleNames(object), group = condition(object))
  dat <-cbind.data.frame(dat, Expression(object))
  rownames(dat) <- NULL
  control <- object@Dataset_summary$condition_levels[[1]]
  case <- object@Dataset_summary$condition_levels[[2]]

  uncollapsed_data <- list()
  uncollapsed_data[['Assays']] <- object@Assays
  uncollapsed_data[['Metadata']] <- object@Metadata
  uncollapsed_data[['Dataset_summary']] <- object@Dataset_summary
  uncollapsed_data[['Nodes']] <- object@Nodes

  object@Assays[['Expression']] <- NULL
  object@Assays[['NormalExpression']] <- NULL

  object@Metadata[['Samples']] <-NULL
  object@Metadata[['Features']] <- NULL
  object@Metadata[['Samples']] <- NULL
  object@Metadata[['clean_Feature_Names']] <- NULL
  object@Metadata[['condition']] <- NULL

  object@Nodes <- data.frame()
  object@Dataset_summary <- list()

  if (method == "correlation") {
    res <- collapseNodes_cor(dat = dat,
                      corr.threshold = corr.threshold)
  }else if (method == "knowedge") {
    res <- collapseNodes_knowledge(dat = dat,
                            metabolite.groups = metabolite.groups)
  }else if (method == "hybrid") {
    res <- collapseNodes_hybrid(dat = dat,
                         metabolite.groups = metabolite.groups,
                         corr.threshold = corr.threshold)
  }
  if(!(is.null(NormalExpression(object)))) warning(paste0('The un-normalized data from the Expression Assay was used for feature reduction.','\n',
                                                          'The data in the NormalExpression Assay was replaced with log-scaled collapsed data.', '\n',
                                                          'If you prefer another normalization method replace this data prior to proceeding!', '\n\n',
                                                          '(orginal data can be accessed with unCollapsedData function)', '\n'))
  res[['new.data']] <-cbind.data.frame(res[['new.data']][,2],apply(res[['new.data']][,-c(1,2)], 2, as.numeric))
  rownames(res[['new.data']]) <- uncollapsed_data[["Metadata"]][["Samples"]]
  restructured_data <- restructure_input_data(Expression = res[['new.data']], control = control, case = case)
  restructured_data[[2]][["group_membership"]] <- res[['final.membership']]
  restructured_data[[2]][["uncollapsed_data"]] <- uncollapsed_data
  object@Assays <- restructured_data[[1]]
  object@Metadata <-restructured_data[[2]]

  #Perform diagnostics on collapsed dataset
  diagnostic_values <- dataDiagnostics(object)
  object@Dataset_summary <- diagnostic_values[[1]]
  object@Nodes <- diagnostic_values[[2]]


  return(object)

}


################################################################################
# main functions #

# correlation-based node-collapsing
collapseNodes_cor <- function(dat, corr.threshold = 0.9) {
  n <- nrow(dat)
  p <- ncol(dat)-2
  sample.groups <- unique(dat[,2])
  metabs <- colnames(dat)[-c(1:2)]

  x <- list()
  x[[1]] <- dat[dat[,2] == sample.groups[1],]
  x[[2]] <- dat[dat[,2] == sample.groups[2],]

  cor.mat <- list()
  clust <- list()
  mycl <- list()
  for (a in 1:length(x)) {
    cor.mat[[a]] <- cor(x[[a]][,-c(1:2)],
                        use = "pairwise.complete.obs",
                        method = "pearson")
    clust[[a]] <- hclust(as.dist(1-abs(cor.mat[[a]])))
    mycl[[a]] <- cutree(clust[[a]], h = 1-corr.threshold)
  }

  D <- matrix(0,p,p)
  colnames(D) <- rownames(D) <- metabs
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      for (k in 1:length(mycl)){
        tmp = mycl[[k]]
        D[i,j] <- D[i,j] + (length(unique(tmp[c(i,j)]))==1)
      }
    }
  }

  D.graph <- graph.adjacency(D, mode = "undirected", weighted = T)
  D.edgelist <- cbind.data.frame(get.edgelist(D.graph), E(D.graph)$weight)
  D.edgelist2 <- D.edgelist[D.edgelist$`E(D.graph)$weight` == 2,c(1:2)]

  if (dim(D.edgelist2)[1] == 0) {
    warning("Metabolites in the data are not sufficiently correlated for collapsing.\n Returning original data without collapsing...\n Try using a lower correlation coefficient threshold to collapse metabolites.")

    final.membership <- data.frame(feature.membership=colnames(dat)[-c(1:2)])
    rownames(final.membership) <- final.membership$feature.membership

    return(list(final.membership = final.membership, new.data = dat))

  } else {
    D.graph2 <- graph_from_edgelist(as.matrix(D.edgelist2[,c(1:2)]))
    D.components <- components(D.graph2)$membership

    final.membership <- data.frame(feature.membership=D.components)
    final.membership$feature.membership <- paste0("FeatureGroup_",
                                                  final.membership$feature.membership)
    final.membership2 <- data.frame(feature.membership=metabs[!(metabs %in% names(D.components))])
    rownames(final.membership2) <- final.membership2$feature.membership
    final.membership <- rbind.data.frame(final.membership, final.membership2)

    newdat <- list()
    for (a in 1:length(x)) {
      newdat[[a]] <- t(x[[a]][,-c(1:2)])
      newdat[[a]] <- merge(final.membership, newdat[[a]],
                           by.x = "row.names", by.y = "row.names")
      newdat[[a]] <- newdat[[a]] %>% group_by(feature.membership) %>%
        summarise(across(everything(), mean))
      rownames(newdat[[a]]) <- newdat[[a]]$feature.membership
      newdat[[a]] <- t(newdat[[a]][,-2])
      colnames(newdat[[a]]) <- newdat[[a]][1,]
      newdat[[a]] <- newdat[[a]][-1,]
      newdat[[a]] <- cbind.data.frame(x[[a]][,c(1:2)],newdat[[a]])
    }

    return(list(final.membership = final.membership,
                new.data = do.call("rbind",newdat)))
  }
}


# node-collapsing based on user-supplied metabolite groups
collapseNodes_knowledge <- function (dat,
                                     metabolite.groups) {

  colnames(metabolite.groups) <- c("metabolite", "metab_group")
  sample.groups <- unique(dat[,2])

  x <- list()
  x[[1]] <- dat[dat[,2] == sample.groups[1],]
  x[[2]] <- dat[dat[,2] == sample.groups[2],]

  newdat <- list()
  for (a in 1:length(x)) {
    newdat[[a]] <- t(x[[a]][,-c(1:2)])
    newdat[[a]] <- merge(metabolite.groups, newdat[[a]],
                         by.x = "metabolite", by.y = "row.names")
    newdat[[a]] <- newdat[[a]] %>% group_by(metab_group) %>%
      summarise(across(everything(), mean))
    rownames(newdat[[a]]) <- newdat[[a]]$metab_group
    newdat[[a]] <- t(newdat[[a]][,-2])
    colnames(newdat[[a]]) <- newdat[[a]][1,]
    newdat[[a]] <- newdat[[a]][-1,]
    newdat[[a]] <- cbind.data.frame(x[[a]][,c(1:2)],newdat[[a]])
  }

  return(new.data = do.call("rbind",newdat))
}



# node-collapsing based on correlations and user-supplied metabolite groups
collapseNodes_hybrid <- function (dat,
                                  metabolite.groups,
                                  corr.threshold = 0.9) {

  colnames(metabolite.groups) <- c("metabolite", "metab_group")
  final.membership <- list()
  new.dat <- list()

  for (z in unique(metabolite.groups$metab_group)) {
    tmp1 <- metabolite.groups[metabolite.groups$metab_group==z,]
    tmp2 <- select(dat, c(1,2,which(names(dat) %in% tmp1$metabolite)))
    if( dim(tmp1)[1] == 1 ) {
      final.membership[[z]] <- data.frame(feature.membership=tmp1$metab_group)
      rownames(final.membership[[z]]) <- tmp1$metabolite
      new.dat[[z]] <- tmp2
    } else {
      final.membership[[z]] <- collapseNodes_cor(tmp2,
                                                 corr.threshold = corr.threshold)$final.membership
      tmp3 <- str_detect(final.membership[[z]]$feature.membership, "FeatureGroup_")
      final.membership[[z]]$feature.membership[tmp3] <- paste0(z, "_",
                                                               final.membership[[z]]$feature.membership[tmp3])

      new.dat[[z]] <- collapseNodes_cor(tmp2,
                                        corr.threshold = corr.threshold)$new.data
      tmp4 <- str_detect(colnames(new.dat[[z]]), "FeatureGroup_")
      colnames(new.dat[[z]])[tmp4] <- paste0(z, "_", colnames(new.dat[[z]])[tmp4])

    }
  }

  names(final.membership) <- NULL
  names(new.dat) <- NULL

  final.membership <- do.call("rbind", final.membership)
  new.dat <- do.call("cbind", lapply(new.dat, function(x) x[-c(1,2)]))
  new.dat <- cbind.data.frame(dat[,c(1:2)], new.dat)

  return(list(final.membership=final.membership,
              new.data=new.dat))
}
