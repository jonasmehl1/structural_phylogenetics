#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# This script:
#   1) Reads a reference tree (e.g., your best IMD tree).
#   2) Reads all bootstrap replicate trees from a folder.
#   3) Computes bootstrap support (Felsenstein style) for each node
#      in the reference tree.
#   4) Outputs a single Newick file with the support in node labels.
#
# Usage example:
#   ./compute_imd_bs.R \
#       --ref_tree best_imd_tree.nwk \
#       --rep_dir   path/to/replicates/ \
#       --out       imd_bs_tree.nwk
# ------------------------------------------------------------------

library(phangorn)
library(argparser)

# 1. Parse command-line arguments
p <- arg_parser("Compute bootstrap support from multiple tree files in a directory.")
p <- add_argument(p, "--ref_tree", help="Reference tree file (Newick)",     type="character")
p <- add_argument(p, "--rep_dir",  help="Directory containing bootstrap replicate trees", type="character")
p <- add_argument(p, "--out",      help="Output tree with bootstrap support", type="character")
args <- parse_args(p)

# 2. Define a helper function to compute bootstrap support
compute_bs_support <- function(topology, replicates) {
  # 'prop.clades' returns the fraction of replicates that contain each clade in 'topology'
  bs_support <- prop.clades(topology, replicates)
  
  # If any node doesn't appear in any replicate, we get NA -> replace with 0
  bs_support[is.na(bs_support)] <- 0
  
  # Assign those bootstrap values to the node labels in 'topology'
  topology$node.label <- bs_support
  
  # (Optionally) remove label from the root node if desired
  # topology$node.label[1] <- ""
  
  return(topology)
}

# 3. Read the reference tree
ref_topology <- read.tree(args$ref_tree)

# 4. Gather all replicate tree files from the specified folder
rep_files <- list.files(
  path = args$rep_dir,
  pattern = "\\.nwk$|\\.tre$|\\.tree$",  # or another pattern matching your replicate tree suffix
  full.names = TRUE
)

if (length(rep_files) == 0) {
  stop("No tree files found in the replicate directory: ", args$rep_dir)
}

# 5. Read in all replicate trees, concatenating them into a single multiPhylo object
#    Each file may contain 1 or more trees. If each file has only 1 tree, it's straightforward.
#    If each file has multiple trees, read.tree() will return a multiPhylo, so we combine them all.
all_reps <- vector("list", length(rep_files))
for (i in seq_along(rep_files)) {
  # read.tree can return a phylo or multiPhylo
  tmp <- read.tree(rep_files[i])
  if (class(tmp) == "phylo") {
    # Single tree in this file -> convert to multiPhylo
    tmp <- c(tmp)
  }
  all_reps[[i]] <- tmp
}
# Combine list of multiPhylo objects into one big multiPhylo:
all_replicate_trees <- do.call(c, all_reps)

message("Total number of replicate trees loaded: ", length(all_replicate_trees))

# 6. Compute the bootstrap support on the reference tree from these replicates
bs_tree <- compute_bs_support(ref_topology, all_replicate_trees)

# 7. Write out the final tree with node labels = bootstrap supports
write.tree(bs_tree, file = args$out)

message("Wrote annotated tree with bootstrap support to: ", args$out)

