source("scripts/runtime.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: cohort_flow.R POOL OUTPUT SITE_RETURN [...]")
run <- args[[1L]]
out <- args[[2L]]
returns <- args[-c(1L, 2L)]

flow <- do.call(rbind, lapply(returns, function(path) {
  utils::read.csv(file.path(path, "review/cohort_flow.csv"), stringsAsFactors = FALSE)
}))
assignments <- utils::read.csv(file.path(run, "combined_assignment.csv"), stringsAsFactors = FALSE)
counts <- tapply(flow$n, flow$category, sum)
phenotypes <- tapply(assignments$n, assignments$phenotype, sum)
eligible <- sum(counts)
landmark <- counts[["primary_landmark"]]
excluded <- sum(counts[c("early_death", "early_hospice", "early_live_discharge", "unclassifiable_outcome")])
assigned <- sum(phenotypes[c("P1", "P2", "P3", "P4")])
insufficient_n <- if ("I" %in% names(phenotypes)) phenotypes[["I"]] else 0
stopifnot(
  length(unique(flow$site)) == length(returns),
  eligible == landmark + excluded,
  landmark == assigned + insufficient_n
)

n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
label <- function(...) paste0("<", paste(c(...), collapse = "<BR/>"), ">")
nodes <- list(
  upstream = label("All ICU admissions at each participating site",
                   "January 2018 to December 2024"),
  eligible = label("Eligible lung cancer ICU hospitalizations",
                   "Adults with primary lung cancer", "Coded as present on hospital admission",
                   paste0("n = ", n(eligible))),
  landmark = label("72-hour landmark cohort",
                   "Alive and still hospitalized", "Classifiable outcome",
                   paste0("n = ", n(landmark))),
  assigned = label("Assigned to a shared phenotype",
                   paste0("n = ", n(assigned)),
                   paste0("P1, ", n(phenotypes[["P1"]]), "   |   P2, ", n(phenotypes[["P2"]])),
                   paste0("P3, ", n(phenotypes[["P3"]]), "   |   P4, ", n(phenotypes[["P4"]]))),
  excluded = label("<B>Excluded from landmark cohort</B>",
                   paste0("n = ", n(excluded)),
                   paste0("Early death, ", n(counts[["early_death"]])),
                   paste0("Early hospice discharge, ", n(counts[["early_hospice"]])),
                   paste0("Early live discharge without hospice, ", n(counts[["early_live_discharge"]])),
                   paste0("Unclassifiable outcome, ", n(counts[["unclassifiable_outcome"]]))),
  insufficient = label("<B>Retained without assignment</B>",
                       "Fewer than 12 hours of known support",
                       paste0("n = ", n(insufficient_n)))
)

dot <- sprintf('digraph cohort {
  graph [layout=dot, rankdir=TB, bgcolor="white", splines=ortho, nodesep=0.32, ranksep=0.17, pad=0.09]
  node [shape=box, style=filled, fillcolor="#F3F5F6", color="#657987", penwidth=0.9,
        fontcolor="#181D20", fontname="Helvetica", fontsize=12, margin="0.14,0.10"]
  edge [color="#3B454D", penwidth=0.9, arrowsize=0.65]
  upstream [label=%s, group=main, width=4.50, height=0.62]
  eligible [label=%s, group=main, width=4.50, height=0.84]
  landmark [label=%s, group=main, width=4.50, height=0.84]
  assigned [label=%s, group=main, width=4.50, height=0.87]
  split_ex [shape=point, label="", width=0, height=0, style=invis, group=main]
  split_info [shape=point, label="", width=0, height=0, style=invis, group=main]
  excluded [label=%s, shape=box, style=solid, color=white, margin="0.11,0.04", fontsize=10.5]
  insufficient [label=%s, shape=box, style=solid, color=white, margin="0.11,0.04", fontsize=10.5]
  upstream:s -> eligible:n [weight=20]
  eligible:s -> split_ex [arrowhead=none, weight=20]
  split_ex -> landmark:n [weight=20]
  landmark:s -> split_info [arrowhead=none, weight=20]
  split_info -> assigned:n [weight=20]
  { rank=same; split_ex; excluded; }
  { rank=same; split_info; insufficient; }
  split_ex -> excluded:w [constraint=false]
  split_info -> insufficient:w [constraint=false]
}', nodes$upstream, nodes$eligible, nodes$landmark, nodes$assigned, nodes$excluded, nodes$insufficient)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
writeLines(svg, file.path(out, "Figure S1 - Cohort Flow.svg"), useBytes = TRUE)
rsvg::rsvg_png(charToRaw(svg), file = file.path(out, "Figure S1 - Cohort Flow.png"), width = 2350)
