################################################################################
#' Sankey plot of \code{\link{phyloseq-class}} object
#' @description
#' `r lifecycle::badge("maturing")`
#' @inheritParams clean_pq
#' @param fact Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param taxa a vector of taxonomic rank to plot
#' @param add_nb_seq Represent the number of sequences or the
#'   number of OTUs (add_nb_seq = FALSE). Note that plotting the number of
#'   sequences is slower.
#' @param min_prop_tax (default: 0) The minimum proportion for taxon to be
#'  plotted. EXPERIMENTAL. For the moment each links below the min.prop.
#'  tax is discard from the sankey network resulting in sometimes weird plot.
#' @param tax2remove  a vector of taxonomic groups to remove from the analysis
#'   (e.g. \code{c('Incertae sedis', 'unidentified')})
#' @param units  character string describing physical units (if any) for Value
#' @param symbol2sub (default: c('\\.', '-')) vector of symbol to delete in
#'   the taxonomy
#' @param ... Additional arguments passed on to
#'   \code{\link[networkD3]{sankeyNetwork}}
#'
#' @examples
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' sankey_pq(GP, fact = "SampleType")
#' \donttest{
#' sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01)
#' sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01, add_nb_seq = TRUE)
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link[networkD3]{sankeyNetwork}} plot representing the
#'  taxonomic distribution of OTUs or sequences. If \code{fact} is set,
#'  represent the distribution of the last taxonomic level in the modalities
#'  of \code{fact}
#'
#' @export
#' @seealso \code{\link[networkD3]{sankeyNetwork}}

sankey_pq <-
  function(physeq = NULL,
           fact = NULL,
           taxa = 1:4,
           add_nb_seq = FALSE,
           min_prop_tax = 0,
           tax2remove = NULL,
           units = NULL,
           symbol2sub = c("\\.", "-"),
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }
    
    if (!physeq@otu_table@taxa_are_rows) {
      otu_tab <- t(physeq@otu_table)
    } else {
      otu_tab <- physeq@otu_table
    }
    
    if (!add_nb_seq) {
      otu_tab[otu_tab > 0] <- 1
      mat_interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames(mat) <- c("Var1", "Var2", "value")
      for (i in 1:(length(taxa) - 1)) {
        res_interm <-
          table(physeq@tax_table[, taxa[i]], physeq@tax_table[, taxa[i + 1]])
        mat_interm <- reshape2::melt(res_interm)
        mat_interm <- mat_interm[mat_interm[, 3] > 0, ]
        mat <- rbind(mat, mat_interm)
      }
    } else if (add_nb_seq) {
      mat_interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames(mat) <- c("Var1", "Var2", "value")
      tax_table_interm <-
        physeq@tax_table[rep(seq(1, ntaxa(physeq)),
                             times = taxa_sums(physeq)
        )]
      
      for (i in 1:(length(taxa) - 1)) {
        res_interm <-
          table(tax_table_interm[, taxa[i]], tax_table_interm[, taxa[i + 1]])
        mat_interm <- reshape2::melt(res_interm)
        mat_interm <- mat_interm[mat_interm[, 3] > 0, ]
        mat <- rbind(mat, mat_interm)
      }
    }
    
    if (!is.null(fact)) {
      net_matrix2links <- function(m = NULL) {
        res <- matrix(ncol = 3)
        for (i in seq_len(dim(m)[1])) {
          for (j in seq_len(dim(m)[2])) {
            if (m[i, j] > 0) {
              res <- rbind(res, c(rownames(m)[i], colnames(m)[j], m[i, j]))
            }
          }
        }
        return(res)
      }
      
      mat_interm <-
        apply(otu_tab, 1, function(x) {
          tapply(
            x, physeq@sam_data[, fact],
            sum
          )
        })
      
      if (!add_nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x) {
            tapply(x, physeq@tax_table[
              ,
              taxa[length(taxa)]
            ], function(x) {
              sum(x > 0)
            })
          })
      } else if (add_nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x) {
            tapply(x, physeq@tax_table[
              ,
              taxa[length(taxa)]
            ], sum)
          })
      }
      
      samp_links <- net_matrix2links(mat_interm)
      samp_links[, 2] <- toupper(samp_links[, 2])
      colnames(samp_links) <- colnames(mat)
      mat <- rbind(mat, samp_links)
    }
    
    mat <- as.data.frame(mat[rowSums(is.na(mat)) == 0, ])
    mat[, 3] <- as.numeric(as.vector(mat[, 3]))
    mat <- mat[rowSums(is.na(mat)) == 0, ]
    
    
    if (!is.null(tax2remove)) {
      mat <- mat[!mat[, 1] %in% tax2remove, ]
      mat <- mat[!mat[, 2] %in% tax2remove, ]
    }
    
    if (min_prop_tax != 0) {
      min_nb_tax <- min_prop_tax * sum(mat[, 3]) / length(taxa)
      mat <- mat[mat[, 3] >= min_nb_tax, ]
    }
    
    for (i in seq_len(length(symbol2sub))) {
      mat <- apply(mat, 2, function(x) {
        gsub(symbol2sub[i], "", x)
      })
    }
    
    tax_sank <- list()
    names_nodes <-
      unique(c(as.vector(mat[, 1]), as.vector(mat[, 2])))
    names_nodes <- names_nodes[!is.na(names_nodes)]
    tax_sank$nodes <-
      data.frame((seq_len(length(names_nodes))) - 1, names_nodes)
    names(tax_sank$nodes) <- c("code", "name")
    mat2 <- mat
    for (i in seq_len(nrow(tax_sank$nodes))) {
      mat2[, 1] <-
        gsub(
          paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
          tax_sank$nodes[
            i,
            1
          ],
          mat2[, 1]
        )
      mat2[, 2] <-
        gsub(
          paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
          tax_sank$nodes[
            i,
            1
          ],
          mat2[, 2]
        )
    }
    
    tax_sank$links <- apply(mat2, 2, as.numeric)
    tax_sank$links <-
      data.frame(tax_sank$links[rowSums(is.na(tax_sank$links)) == 0, ])
    tax_sank$nodes <-
      as.data.frame(as.character(tax_sank$nodes[, 2]))
    names(tax_sank$nodes) <- "name"
    names(tax_sank$links) <- c("source", "target", "value")
    if (is.null(units)) {
      if (!add_nb_seq) {
        units <- "OTUs"
      } else if (add_nb_seq) {
        units <- "Sequences"
      }
    }
    networkD3::sankeyNetwork(
      Links = tax_sank$links,
      Nodes = tax_sank$nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      units = units,
      ...
    )
  }
################################################################################
