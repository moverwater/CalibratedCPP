#!/usr/bin/env Rscript
# Convert the IQ-TREE ML tree (genetic-distance branch lengths) into an ultrametric
# chronogram (time-based branch lengths, tips aligned at present) using ape::chronos(),
# supplying the same 13 fossil calibration bounds used throughout this project's BEAST
# analyses (from primates.lphy) as calibration constraints, rather than just rescaling
# relative to one arbitrary point.
#
# Usage: Rscript make_ultrametric.R [input.treefile] [output.nwk]

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
input_path  <- if (length(args) >= 1) args[1] else "primates_iqtree.treefile"
output_path <- if (length(args) >= 2) args[2] else "primates_iqtree_startingTree.nwk"

tree <- read.tree(input_path)
cat("Loaded tree:", length(tree$tip.label), "tips\n")

# IQ-TREE's -o flag only repositions the named outgroup for display; the underlying tree
# is still an unrooted bifurcation with a basal trifurcation (BigClade, Mus_musculus,
# Tupaia_chinensis) rather than the properly nested (Mus_musculus,(Tupaia_chinensis,...))
# structure our calibration scheme assumes. Without this, getMRCA() for "all tips except
# Mus_musculus" resolves to the same node as the full-tree root, corrupting the
# Euarchontoglires/Euarchonta calibrations. Explicitly re-root and force full bifurcation.
tree <- root(tree, outgroup = "Mus_musculus", resolve.root = TRUE)
cat("Re-rooted on Mus_musculus with forced bifurcation\n")

# Same 13 calibration clades (taxa + bounds, in Ma) as calibratedcpp-beast/examples/primates/primates.lphy
calibrations <- list(
  list(taxa = tree$tip.label, lower = 65.79, upper = 125.816),  # Euarchontoglires (root)
  list(taxa = setdiff(tree$tip.label, "Mus_musculus"), lower = 65.79, upper = 125.816),  # Euarchonta
  list(taxa = c("Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli",
                "Carlito_syrichta", "Cebus_capucinus_imitator", "Saimiri_boliviensis",
                "Aotus_nancymaae", "Callithrix_jacchus", "Nomascus_leucogenys", "Pongo_abelii",
                "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes",
                "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana", "Chlorocebus_sabaeus", "Macaca_nemestrina",
                "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys", "Mandrillus_leucophaeus",
                "Papio_anubis", "Theropithecus_gelada"), lower = 55.935, upper = 66.095),  # Primates
  list(taxa = c("Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli"),
       lower = 18.5, upper = 55.8),  # Lorisiformes
  list(taxa = c("Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana", "Chlorocebus_sabaeus", "Macaca_nemestrina",
                "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys", "Mandrillus_leucophaeus",
                "Papio_anubis", "Theropithecus_gelada"), lower = 12.47, upper = 25.235),  # Cercopithecidae
  list(taxa = c("Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana"), lower = 8.125, upper = 15.0),  # Colobinae
  list(taxa = c("Chlorocebus_sabaeus", "Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta",
                "Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"),
       lower = 6.5, upper = 15.0),  # Cercopithecinae
  list(taxa = c("Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys",
                "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"),
       lower = 5.33, upper = 12.51),  # Papionini
  list(taxa = c("Nomascus_leucogenys", "Pongo_abelii", "Gorilla_gorilla", "Homo_sapiens",
                "Pan_paniscus", "Pan_troglodytes"), lower = 13.4, upper = 25.235),  # Hominoidea
  list(taxa = c("Pongo_abelii", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"),
       lower = 12.3, upper = 25.235),  # Hominidae
  list(taxa = c("Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"), lower = 4.631, upper = 15.0),  # Homo+Pan
  list(taxa = c("Cebus_capucinus_imitator", "Saimiri_boliviensis", "Aotus_nancymaae", "Callithrix_jacchus"),
       lower = 13.183, upper = 34.5),  # Callitrichidae+Cebidae
  list(taxa = c("Cebus_capucinus_imitator", "Saimiri_boliviensis"), lower = 13.032, upper = 34.5)  # Cebidae
)

calib_rows <- do.call(rbind, lapply(calibrations, function(cal) {
  node <- if (length(cal$taxa) == length(tree$tip.label)) {
    Ntip(tree) + 1  # root node
  } else {
    getMRCA(tree, cal$taxa)
  }
  data.frame(node = node, age.min = cal$lower, age.max = cal$upper, soft.bounds = FALSE)
}))

cat("Calibration nodes resolved:\n")
print(calib_rows)

chrono_tree <- chronos(tree, calibration = calib_rows, model = "relaxed")

write.tree(chrono_tree, output_path)
cat("Wrote ultrametric tree to", output_path, "\n")
cat("Root age:", max(node.depth.edgelength(chrono_tree)), "\n")
