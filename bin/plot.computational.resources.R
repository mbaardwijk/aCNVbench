library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

# define set colors for all callers
light_rhg_cols <- c("PennCNV"="#8c9cd5", "QuantiSNP"="#e68d8d", "iPattern"="#8bcf8c", "EnsembleCNV"="#f5f4a0", "EnsembleCNV.no.gq"="#f5f4a0", "EnsembleCNV.gq.15"="#f5f4a0", "EnsembleCNV.gq.20"="#f5f4a0", "R-GADA"="#c7b497", "Gold Standard"="#e9f3fc")
dark_rhg_cols <- c("PennCNV"="#0e31a6", "QuantiSNP"="#cd1b1b", "iPattern"="#089908", "EnsembleCNV"="#eae83b", "EnsembleCNV.no.gq"="#eae83b", "EnsembleCNV.gq.15"="#eae83b", "EnsembleCNV.gq.20"="#eae83b", "R-GADA"="#8b6228", "Gold Standard"="#2fa9e2")

computational.resources <- read.csv("C:/Users/myrth/Erasmus MC/Clinical Bioinformatics - Myrthe/Projects/CNV_Benchmark/Results/Resources/Computational_time_and_resources.csv", header=T)
computational.resources$X <- factor(computational.resources$X, levels=c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Preprocessing"))

p1 <- ggplot(computational.resources, aes(x=X, y=Processing.time.per.sample.s., colour=X, fill=X)) + 
  geom_bar(stat='identity', width=0.6) +
  labs(x="", y="Processing time per sample (s)", fill="Workflow", colour="Workflow") +
  theme_classic() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.text.x=element_text(angle=90)) +
  scale_colour_manual(values = dark_rhg_cols) +
  scale_fill_manual(values = light_rhg_cols)

p2 <- ggplot(computational.resources, aes(x=X, y=Peak.of.real.memory.GB., colour=X, fill=X)) + 
  geom_bar(stat='identity', width=0.6) +
  labs(x="", y="Peak of real memory (GB)", fill="Workflow", colour="Workflow") +
  theme_classic() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.text.x=element_text(angle=90)) +
  scale_colour_manual(values = dark_rhg_cols) +
  scale_fill_manual(values = light_rhg_cols)

ggarrange(p1, p2, ncol=2, nrow=1, common.legend = T, legend="top", labels=c("A", "B"), label.y=0.05)


tiff("C:/Users/myrth/Erasmus MC/Clinical Bioinformatics - Myrthe/Projects/CNV_Benchmark/Results/Figures/Computational_resources.tiff", units="cm", width=24, height=16, res=300)
g <- ggarrange(p1, p2, ncol=2, nrow=1, common.legend = T, legend="top", labels=c("A", "B"), label.y=0.05)
grid.draw(g)
dev.off()