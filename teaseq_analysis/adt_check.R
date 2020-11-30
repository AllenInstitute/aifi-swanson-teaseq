library(data.table)
library(dplyr)
library(ggplot2)

all_adt <- fread("/mnt/x060-multiome/barcounter/X061_full/X061-AP0C1W1-C_Tag_Counts.csv")
cellranger_barcodes <- fread("/mnt/x060-multiome/X061-Well1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE)

all_adt$arc_cell <- all_adt$cell_barcode %in% sub("-1","",cellranger_barcodes[[1]])

hist(log10(all_adt$total))

sum(all_adt$total > 500)
# ~8.7k barcodes have > 500

sum(all_adt$total > 5000)
# only 95 barcodes are high outliers

gt100_adt <- all_adt[total > 500]
gt100_adt <- gt100_adt[total < 5000]

hist(gt100_adt$total, breaks = 200,
     xlim = c(0, 5000))
# looks like there's another break above 500

hist(log10(gt100_adt$total), breaks = 50)

hist(log10(gt100_adt$total[gt100_adt$arc_cell]), breaks = 50)
hist(log10(gt100_adt$total[!gt100_adt$arc_cell]), breaks = 50)

ggplot(data = gt100_adt) +
  geom_histogram(aes(x = log10(total),
                     fill = arc_cell),
                 binwidth = 0.02)

# Total input reads
total_reads <- 2166370728 / 4
total_barcode_umis <- sum(all_adt$total)

reads_per_umi <- total_reads / total_barcode_umis

# ADT-based UMAP

arc_adt <- all_adt %>%
  filter(arc_cell)

arc_mat <- arc_adt %>%
  select(-cell_barcode, -total, -arc_cell) %>%
  as.matrix()

rownames(arc_mat) <- paste0(arc_adt$cell_barcode, "-1")

saveRDS(arc_mat, "adt_count_mat.rds")

norm_arc_mat <- arc_mat / rowSums(arc_mat) * 1e3

adt_umap <- uwot::umap(log10(norm_arc_mat + 1))
plot(adt_umap)

adt_umap_df <- as.data.frame(adt_umap)
names(adt_umap_df) <- c("adt_umap_1", "adt_umap_2")
adt_umap_df$original_barcode = arc_adt$cell_barcode

tile_umap <- read.csv("tile_umap_coords.csv")

adt_umap_df <- adt_umap_df %>%
  mutate(original_barcode = paste0(original_barcode, "-1")) %>%
  left_join(tile_umap)

ggplot(adt_umap_df) +
  geom_point(aes(x = adt_umap_1,
                 y = adt_umap_2,
                 color = archr_cluster))

all_adt_melt <- arc_adt %>%
  select(-total, -arc_cell) %>%
  reshape2::melt(id.vars = "cell_barcode",
                 variable.name = "target",
                 value.name = "count") %>%
  filter(!is.na(count)) %>%
  mutate(logcount = log10(count + 1))

ggplot() +
  geom_histogram(data = all_adt_melt,
                 aes(x = logcount),
                 bins = 50) +
  scale_y_continuous(limits = c(0, 3000)) +
  facet_wrap(~target,
             ncol = 8)

# Not so good: CD10, CD127, CD197, CD25, CD278, CD279, CD304, CD319, CD80, TCR-Va24-Ja18, TCR-Va7.2