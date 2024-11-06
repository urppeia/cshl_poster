library(clusterProfiler)
library(org.At.tair.db)
library(topGO)
library(GOstats)
library(ggplot2)
library(plyr)

cpg_data <- read.table(gzfile("data/cpg.txt.gz"), header = TRUE, sep = "\t")
cpg_data <- split(x = cpg_data, f = cpg_data$class)

chg_data <- read.table(gzfile("data/chg.txt.gz"), header = TRUE, sep = "\t")
chg_data <- split(x = chg_data, f = chg_data$class)

chh_data <- read.table(gzfile("data/chh.txt.gz"), header = TRUE, sep = "\t")
chh_data <- split(x = chh_data, f = chh_data$class)

run_enrichment <- function(data) {
  c10 <- unique(data[!(is.na(data$JSD_bit_10C)), "nearestTSS.gene_id"])
  c16 <- unique(data[!(is.na(data$JSD_bit_16C)), "nearestTSS.gene_id"])
  c22 <- unique(data[!(is.na(data$JSD_bit_22C)), "nearestTSS.gene_id"])
  
  gene_list <- unique(data$nearestTSS.gene_id)
  gene_universe <- bitr(gene_list, fromType = "TAIR", toType = "GO", OrgDb = org.At.tair.db)
  try(
    go_10 <- enrichGO(
      gene          = c10, # nearestTSS.gene_id
      OrgDb         = org.At.tair.db, # tair
      keyType       = "TAIR",
      ont           = "BP", # Ontology --> BP Biological Process, MF Molecular Function, CC Cellular Component
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  )

  try(
    go_10 <- simplify(go_10)
  )

  try(
    go_16 <- enrichGO(
      gene          = c16, # nearestTSS.gene_id
      OrgDb         = org.At.tair.db, # tair
      keyType       = "TAIR",
      ont           = "BP", # Ontology --> BP Biological Process, MF Molecular Function, CC Cellular Component
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  )
  
  try(
    go_16 <- simplify(go_16)
  )
  
  try(
    go_22 <- enrichGO(
      gene          = c22, # nearestTSS.gene_id
      OrgDb         = org.At.tair.db, # tair
      keyType       = "TAIR",
      ont           = "BP", # Ontology --> BP Biological Process, MF Molecular Function, CC Cellular Component
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  )
  
  try(
    go_22 <- simplify(go_22)
  )

  res <- list(go_10@result, go_16@result, go_22@result)
  res <- lapply(res, function(x){
    a <- x[x$p.adjust <= 0.05,][1:20,]
    a
  })
  names(res) <- c("10C", "16C", "22C")
  
  res %>% ldply(.id = "temp") -> res
  
  return(res)
}


cpg_bp <- lapply(cpg_data, run_enrichment)
chg_bp <- lapply(chg_data, run_enrichment)
chh_bp <- lapply(chh_data, run_enrichment)

cpg_bp %>% ldply(.id = "class") -> cpg_bp
chg_bp %>% ldply(.id = "class") -> chg_bp
chh_bp %>% ldply(.id = "class") -> chh_bp


cpg_bp <- cpg_bp[cpg_bp$p.adjust <= 0.05,]
chg_bp <- chg_bp[chg_bp$p.adjust <= 0.05,]

ggplot(data = cpg_bp, aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) +
  geom_point()+
  facet_grid(class~temp, scales = "free")

ggplot(data = chg_bp, aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) +
  geom_point()+
  facet_grid(class~temp, scales = "free")


c <- ggplot(data = arrange(go_df, Hyper_Fold_Enrichment), mapping = aes(
  x = name,
  y = Hyper_Fold_Enrichment,
  color = -log10(Hyper_Adjp_BH),
  size = Total_Genes_Annotated
)) +
  geom_point(show.legend = T, alpha = 0.7, stroke = 1.5) +
  coord_flip(expand = T) +
  # scale_fill_distiller(palette = "Spectral") +
  scale_colour_gradientn(
    colors = myPalette(n = nrow(go_df)),
    # colors = rev(cividis(n = nrow(go_df))),
    # guide = guide_colorbar(reverse = TRUE),
    breaks = c(ceiling(min(-log10(go_df$Hyper_Adjp_BH))), floor(max(-log10(go_df$Hyper_Adjp_BH))))
  ) +
  scale_size_continuous(breaks = c(200, 800)) +
  xlab("") +
  ylab("Fold enrichment") +
  facet_wrap(~Group, scales = "free_y", labeller = label_wrap_gen(width = 22), nrow = 1) +
  theme_bw(base_family = "Helvetica") +
  th +
  guides(
    size = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      nrow = 1, size = 10,
      label.hjust = 0.5, label.vjust = 0.5,
      order = 2
    ),
    colour = guide_colorbar(
      nrow = 1, order = 1,
      title.position = "left", title.hjust = 0.5, title.vjust = 1,
      barwidth = unit(2, "cm"), barheight = unit(0.3, "cm"),
      label.hjust = 0.5, label.vjust = 0.5, size = 10
    )
  ) +
  labs(
    color = expression(-Log[10] ~ "adjusted" ~ italic(P) ~ " "),
    size = "Annotated genes"
  ) +
  theme(axis.text.y = element_text(lineheight = 0.8), legend.spacing = unit(1, "cm"))

cg <- gridExtra::grid.arrange(egg::set_panel_size(p = c, width = unit(3.5, "cm"), height = unit(10, "cm")))

ggsave(plot = cg, filename = "output/c.pdf", width = unit(8.5, "in"), height = unit(11, "in"), dpi = 320)
ggsave(plot = cg, filename = "output/c.png", width = unit(8.5, "in"), height = unit(11, "in"), dpi = 320)
