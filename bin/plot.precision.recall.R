suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(dplyr))

option_list = list(
  make_option(c("--performances"), action = "store", type = "character", default = NA,
              help = "List of all deletions result files for all platforms"),
  make_option(c("--population"), action = "store", type = "character", default = NA,
              help = "Population panel file"),
  make_option(c("--output"), action = "store", type = "character", default = NA,
              help = "Prefix of output plots")
)

# define set colors for all callers
light_rhg_cols <- c("PennCNV"="#8c9cd5", "QuantiSNP"="#e68d8d", "iPattern"="#8bcf8c", "EnsembleCNV"="#f5f4a0", "EnsembleCNV.no.gq"="#f5f4a0", "EnsembleCNV.gq.15"="#f5f4a0", "EnsembleCNV.gq.20"="#f5f4a0", "R-GADA"="#c7b497", "Gold Standard"="#e9f3fc")
dark_rhg_cols <- c("PennCNV"="#0e31a6", "QuantiSNP"="#cd1b1b", "iPattern"="#089908", "EnsembleCNV"="#eae83b", "EnsembleCNV.no.gq"="#eae83b", "EnsembleCNV.gq.15"="#eae83b", "EnsembleCNV.gq.20"="#eae83b", "R-GADA"="#8b6228", "Gold Standard"="#2fa9e2")
shapes <- c(21:25)

# create textbox backgrounds for title
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_classic()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_classic())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

plotPerformances <- function(performances, color.palette){
  performances$Caller <- factor(performances$Caller, levels=c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA"))
  performances[is.na(performances)] <- 0
  performances <- performances[order(sapply(performances$Caller, function(x) which(x == c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA")))), ]
  performances$x <- seq(1, nrow(performances))
  print(performances)
  p1 <- ggplot(performances, aes(x=Overlap.t, y=Average.Sensitivity, colour=Caller, fill=Caller, shape=Category)) + 
    geom_line() +
    geom_point(size=3) +
    labs(x="Percentage of reciprocal overlap", y="Sensitivity", shape="Category") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols) 
  
  p2 <- ggplot(performances, aes(x=Overlap.t, y=Average.Precision, colour=Caller, fill=Caller, shape=Category)) +
    geom_point(size=3) +
    geom_line() +
    labs(x="Percentage of reciprocal overlap", y="Precision") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols)
  
  p3 <- ggplot(performances, aes(x=Overlap.t, y=Average.F1.score, colour=Caller, fill=Caller, shape=Category)) +
    geom_point(size=3) +
    geom_line() +
    labs(x="Percentage of reciprocal overlap", y="F1 Score") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols)
  
  tiff(paste0(opt$output, "_average_precision_recall.tiff"), units="cm", width=40, height=20, res=300)
  g <- ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = T, legend="right", labels=c("A", "B", "C"), label.y=0.04)
  grid.draw(g)
  dev.off()
  
  p4 <- ggplot(performances, aes(x=Overlap.t, y=Overall.Sensitivity, colour=Caller, fill=Caller, shape=Category)) + 
    geom_line() +
    geom_point(size=3) +
    labs(x="Percentage of reciprocal overlap", y="Sensitivity", shape="Category") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols) 
  
  p5 <- ggplot(performances, aes(x=Overlap.t, y=Overall.Precision, colour=Caller, fill=Caller, shape=Category)) +
    geom_point(size=3) +
    geom_line() +
    labs(x="Percentage of reciprocal overlap", y="Precision") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols)
  
  p6 <- ggplot(performances, aes(x=Overlap.t, y=Overall.F1.score, colour=Caller, fill=Caller, shape=Category)) +
    geom_point(size=3) +
    geom_line() +
    labs(x="Percentage of reciprocal overlap", y="F1 Score") +
    theme_classic() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), legend.position = "right") +
    scale_shape_manual(values = shapes) + 
    scale_colour_manual(values = color.palette) +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = light_rhg_cols)
  
  tiff(paste0(opt$output, "_overall_precision_recall.tiff"), units="cm", width=40, height=20, res=300)
  g <- ggarrange(p4, p5, p6, ncol=3, nrow=1, common.legend = T, legend="right", labels=c("A", "B", "C"), label.y=0.04)
  grid.draw(g)
  dev.off()
}

opt <- parse_args(OptionParser(option_list = option_list))

# Execute analysis for all samples
# Collect all performance files for each unique caller in one dataframe
performance.files <- unlist(strsplit(opt$performances, "\\s+"))
caller.names <- lapply(performance.files, function(x) gsub(paste0(".", "performances", ".txt"), "", x))
performance.dfs <- lapply(performance.files, function(x) read.table(x, sep="\t", header=T))
performance.dfs <- mapply(cbind, performance.dfs, "Caller" = caller.names, SIMPLIFY=F)
all.performances <- do.call(rbind, performance.dfs)

# Determine all unique categories and overlap thresholds
categories <- unique(all.performances$Category)
thresholds <- unique(all.performances$Overlap.t)

# Create df with all overall and sample average scores for each caller, category and overlap threshold
overall.performances <- data.frame(Caller = character(), Category = character(), Overlap.t = numeric(),
                                   Overall.Sensitivity = numeric(), Overall.Precision = numeric(), Overall.F1.score = numeric(),
                                   Average.Sensitivity = numeric(), Average.Precision = numeric(), Average.F1.score = numeric())
# For each caller, for each category and for each overlap threshold, determine overall and sample average scores
for(caller in caller.names){
  for(category in categories){
    for(threshold in thresholds){
      # subset dataframe
      current.performances <- all.performances[all.performances$Caller == caller &
                                                 all.performances$Category == category &
                                                 all.performances$Overlap.t == threshold,]
      # determine overall scores
      overall.predictions <- colSums(current.performances[,c("TP", "FP", "FN")])
      Overall.Sensitivity = overall.predictions[['TP']] / (overall.predictions[['TP']] + overall.predictions[['FN']])
      Overall.Precision = overall.predictions[['TP']] / (overall.predictions[['TP']] + overall.predictions[['FP']])
      Overall.F1.score = 2 * ((Overall.Precision * Overall.Sensitivity) / (Overall.Precision + Overall.Sensitivity))
      # determine average scores
      average.scores <- colMeans(current.performances[,c("Sensitivity", "Precision", "F1.score")])
      # combine all current performance scores
      current.performances <- c(caller, category, threshold, Overall.Sensitivity, Overall.Precision, Overall.F1.score, average.scores)
      names(current.performances) <- c("Caller", "Category", "Overlap.t",
                                          "Overall.Sensitivity", "Overall.Precision", "Overall.F1.score",
                                          "Average.Sensitivity", "Average.Precision", "Average.F1.score")
      overall.performances[nrow(overall.performances)+1,] <- current.performances
    }
  }
}
overall.performances <- overall.performances %>% mutate_at(c("Overlap.t", "Overall.Sensitivity", "Overall.Precision",
                                                             "Overall.F1.score", "Average.Sensitivity", "Average.Precision",
                                                             "Average.F1.score"), as.numeric)
plotPerformances(overall.performances, dark_rhg_cols)

# # Execute analysis per super population
# if(!is.na(opt$population)){
#   population.info <- read.table(opt$population, sep="\t", header=T)
#   for(population in unique(population.info$super_pop)){
#     population.samples <- population.info[population.info$super_pop == population, 'sample']
# 
#     pop.deletion.dfs <- lapply(deletion.dfs, function(x) x[rownames(x) %in% population.samples, ])
#     pop.duplication.dfs <- lapply(duplication.dfs, function(x) x[rownames(x) %in% population.samples, ])
#     pop.longcnv.dfs <- lapply(longcnv.dfs, function(x) x[rownames(x) %in% population.samples, ])
#     pop.shortcnv.dfs <- lapply(shortcnv.dfs, function(x) x[rownames(x) %in% population.samples, ])
#     pop.allcnv.dfs <- lapply(allcnv.dfs, function(x) x[rownames(x) %in% population.samples, ])
#     
#     pop.deletion.means <- lapply(names(pop.deletion.dfs), function(x) colMeans(pop.deletion.dfs[[x]]))
#     pop.duplication.means <- lapply(names(pop.duplication.dfs), function(x) colMeans(pop.duplication.dfs[[x]]))
#     pop.longcnv.means <- lapply(names(pop.longcnv.dfs), function(x) colMeans(pop.longcnv.dfs[[x]]))
#     pop.shortcnv.means <- lapply(names(pop.shortcnv.dfs), function(x) colMeans(pop.shortcnv.dfs[[x]]))
#     pop.allcnv.means <- lapply(names(pop.allcnv.dfs), function(x) colMeans(pop.allcnv.dfs[[x]]))
#     names(pop.deletion.means) <- unlist(platformnames)
#     names(pop.duplication.means) <- unlist(platformnames)
#     names(pop.longcnv.means) <- unlist(platformnames)
#     names(pop.shortcnv.means) <- unlist(platformnames)
#     names(pop.allcnv.means) <- unlist(platformnames)
#     
#     pop.deletion.means.df <- as.data.frame(do.call(rbind, pop.deletion.means))
#     pop.duplication.means.df <- as.data.frame(do.call(rbind, pop.duplication.means))
#     pop.longcnv.means.df <- as.data.frame(do.call(rbind, pop.longcnv.means))
#     pop.shortcnv.means.df <- as.data.frame(do.call(rbind, pop.shortcnv.means))
#     pop.allcnv.means.df <- as.data.frame(do.call(rbind, pop.allcnv.means))
#     pop.deletion.means.df$Type <- "Deletions"
#     pop.duplication.means.df$Type <- "Duplications"
#     pop.longcnv.means.df$Type <- "Long CNVs"
#     pop.shortcnv.means.df$Type <- "Short CNVs"
#     pop.allcnv.means.df$Type <- "All CNVs"
#     pop.deletion.means.df$Caller <- rownames(pop.deletion.means.df)
#     pop.duplication.means.df$Caller <- rownames(pop.duplication.means.df)
#     pop.longcnv.means.df$Caller <- rownames(pop.longcnv.means.df)
#     pop.shortcnv.means.df$Caller <- rownames(pop.shortcnv.means.df)
#     pop.allcnv.means.df$Caller <- rownames(pop.allcnv.means.df)
#     
#     plotPerformances(pop.deletion.means.df, pop.duplication.means.df, pop.longcnv.means.df, pop.shortcnv.means.df, pop.allcnv.means.df, dark_rhg_cols, filename=paste0(opt$output, "_precision_recall_", population, ".tiff"))
#   }
# }
