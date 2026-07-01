## Graph the recursive *strong* dependencies of lme4 (Depends/Imports/LinkingTo),
## distinguishing base R packages, Recommended packages, and ordinary CRAN
## dependencies. Requires the 'igraph' package and (optionally) the
## Graphviz 'dot' command-line tool for nicer layout.

options(repos = c(CRAN = "https://cloud.r-project.org"))

db <- tools::CRAN_package_db()
strong <- c("Depends", "Imports", "LinkingTo")
root <- "lme4"

## lme4's own direct strong deps come from the local DESCRIPTION, not CRAN,
## since the CRAN release may lag behind the working tree (e.g. it may still
## list a dependency, such as rlang, that has since been removed locally)
parse_field <- function(x) {
  if (is.na(x)) return(character(0))
  trimws(gsub("\\(.*\\)", "", strsplit(x, ",")[[1]]))
}
desc <- read.dcf("DESCRIPTION")[1, ]
root_direct <- unique(unlist(lapply(strong, function(f) parse_field(desc[f]))))
root_direct <- setdiff(root_direct, "R")

## recursive strong deps of each of lme4's direct dependencies, taken from CRAN
rec <- unique(unlist(c(
  root_direct,
  tools::package_dependencies(root_direct, db = db, which = strong,
                               recursive = TRUE)
)))
allpkgs <- c(root, rec)

## direct (non-recursive) edges among allpkgs, restricted to this package set
direct <- tools::package_dependencies(setdiff(allpkgs, root), db = db,
                                       which = strong, recursive = FALSE)
direct[[root]] <- root_direct
edges <- do.call(rbind, lapply(names(direct), function(p) {
  deps <- direct[[p]]
  deps <- deps[deps %in% allpkgs]
  if (length(deps) == 0) return(NULL)
  data.frame(from = p, to = deps, stringsAsFactors = FALSE)
}))

## package priority (base / recommended / NA) from the CRAN package db
priority <- setNames(db$Priority[match(allpkgs, db$Package)], allpkgs)
## base packages aren't listed on CRAN; fill those in from the local library
base_pkgs <- rownames(installed.packages(priority = "base"))
priority[allpkgs %in% base_pkgs] <- "base"

is_base <- !is.na(priority) & priority == "base"
is_recommended <- !is.na(priority) & priority == "recommended"
is_root <- allpkgs == root

node_style <- function(p) {
  if (is_root[allpkgs == p]) return('fillcolor="gold", penwidth=2')
  if (is_base[allpkgs == p]) return('fillcolor="gray92"')
  if (is_recommended[allpkgs == p]) return('fillcolor="palegreen3"')
  return('fillcolor="lightsteelblue"')
}

dot_lines <- c(
  "digraph deps {",
  "  rankdir=TB;",
  "  graph [nodesep=0.25, ranksep=0.6, splines=true];",
  '  node [style=filled, shape=ellipse, fontname="Helvetica", fontsize=11];',
  "  edge [color=gray60, arrowsize=0.6];",
  sprintf('  "%s" [%s];', allpkgs, vapply(allpkgs, node_style, character(1))),
  sprintf('  "%s" -> "%s";', edges$from, edges$to),
  "}"
)

dot_file <- "misc/lme4_deps.dot"
png_file <- "misc/lme4_deps.png"
writeLines(dot_lines, dot_file)

if (nzchar(Sys.which("dot"))) {
  system2("dot", c("-Tpng", shQuote(dot_file), "-o", shQuote(png_file)))
  cat("wrote", png_file, "via Graphviz\n")
} else {
  ## fall back to igraph's own layout if Graphviz isn't available
  library(igraph)
  g <- graph_from_data_frame(edges, directed = TRUE,
                              vertices = data.frame(name = allpkgs))
  V(g)$color <- vapply(allpkgs, function(p) {
    if (is_root[allpkgs == p]) "gold"
    else if (is_base[allpkgs == p]) "gray92"
    else if (is_recommended[allpkgs == p]) "palegreen3"
    else "lightsteelblue"
  }, character(1))
  root_id <- which(V(g)$name == root)
  lay <- layout_as_tree(g, root = root_id, mode = "out")
  lay[, 2] <- -lay[, 2]
  png(png_file, width = 1400, height = 1000, res = 130)
  par(mar = c(0.5, 0.5, 2, 0.5))
  plot(g, layout = lay, edge.arrow.size = 0.25, vertex.label.color = "black",
       main = "lme4: recursive strong dependencies (Depends/Imports/LinkingTo)")
  dev.off()
  cat("wrote", png_file, "via igraph (install Graphviz for a nicer layout)\n")
}

cat("\nlegend: gold = lme4 (root), gray = base R package,",
    "green = Recommended package, blue = ordinary CRAN dependency\n")
