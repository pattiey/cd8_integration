qq <- function(...) {
  sapply(match.call()[-1], deparse)
}

side_by_side_cluster_plot <- function(data, reduction="pca", group.by1, group.by2, clusterColors){
  plot1 <- DimPlot(data.seurat,
      reduction = reduction,
      group.by = group.by1,
      cols=clusterColors) +
    ggtitle(group.by1, reduction) +
    theme(aspect.ratio = 1, legend.position="right")
  plot2 <- DimPlot(data.seurat,
      reduction = reduction,
      group.by = group.by2,
      cols=brewer.pal(4, "Dark2")) +
    ggtitle(paste(group.by2, reduction)) +
    theme(aspect.ratio = 1, legend.position="right")
  return(plot1 + plot2)
}

subset_setup  <- function(data, idents, seed = 123) {
  subset <- subset(x = data, idents = idents)
  subset <- SCTransform(subset, variable.features.n = 10000, verbose = FALSE, seed.use = seed)
  DefaultAssay(object = subset) = "SCT"
  return(subset)
}

dim_reductions <- function(data, seed = 123) {
  set.seed(seed)
  data <- RunPCA(data, ndims.print = 1:5, nfeatures.print = 10)
  data <- RunUMAP(data, reduction = "pca", dims = 1:10, seed.use = seed)
  data <- RunTSNE(data, reduction = "pca", dims = 1:10, seed.use = seed, perplexity = 30)
  return(data)
}

subset_correlations  <- function(data, genes_of_interest){
  matrix <- data@assays$SCT@scale.data
  matrix_mod <- as.matrix(matrix)
  matrix_mod <- matrix_mod[rowSums(matrix_mod) != 0,]
  genes <- matrix_mod[genes_of_interest,]
  print(Heatmap(cor(t(genes))))
  correlations <- apply(genes, 1,
                              function(x) {apply(matrix_mod, 1, function(y) {cor(x,y)})}) %>%
                        as_tibble() %>%
                        mutate(gene = row.names(matrix_mod)) %>%
                        select(gene, everything())
  return(correlations)
}

subset_highly_correlated <- function(correlations, genes_of_interest, n_to_pull=20){
  highly_correlated <- list()
  for(gene_name in genes_of_interest){
    correlated_genes <- correlations %>%
      slice_max(n = n_to_pull, order_by = abs(correlations %>% select(one_of(gene_name)))) %>%
      pull(gene)

    highly_correlated <- append(highly_correlated, correlated_genes)
  }
  highly_correlated <- correlations %>%
    filter(gene %in% highly_correlated)
  gene_names <- highly_correlated %>% pull(gene)
  highly_correlated <- highly_correlated %>%
    select(-gene) %>%
    as.matrix()
  row.names(highly_correlated) <- gene_names
  return(highly_correlated)
}

heatmap_correlated_genes <- function(data, highly_correlated){
  gene_names <- row.names(highly_correlated)
  matrix <- data@assays$SCT@scale.data
  matrix_mod <- as.matrix(matrix)
  matrix_mod <- matrix_mod[rowSums(matrix_mod) != 0,]
  return(Heatmap(cor(t(matrix_mod[gene_names,])), show_row_names = FALSE, show_column_names = FALSE))
}
