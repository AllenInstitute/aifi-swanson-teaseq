library(data.table)
library(purrr)
library(GenomicRanges)
library(Matrix)
library(UpSetR)
options(stringsAsFactors = FALSE)

type_pos <- data.frame(type=c("B.Naive",
                              "B.Active",
                              "Plasmablast",
                              "T.CD4.Naive",
                              "T.CD4.Memory",
                              "T.CD4.CD45RO.CD45RA",
                              "T.CD4.Exhausted",
                              "T.CD8.Naive",
                              "T.CD8.Effector",
                              "T.CD8.Exhausted",
                              "T.DoublePositive",
                              "T.DoubleNegative",
                              "NK",
                              "DC.Myeloid",
                              "Mono.CD14",
                              "Mono.CD16"))
type_pos$xpos <- 1:nrow(type_pos)
type_pos$ypos <- nrow(type_pos):1

peak_bed_files <- list.files("X044_type_peaks/",
                             pattern = ".narrowPeak",
                             full.names = TRUE)
peak_beds <- map(peak_bed_files,
                 fread)
names(peak_beds) <- sub(".+/","",sub("_peaks.narrowPeak","",peak_bed_files))

n_peaks <- 2500

peak_beds <- map(peak_beds,
                 function(peak_bed) {
                   peak_bed <- peak_bed[order(peak_bed[[9]], decreasing = TRUE),]
                   keep_peaks <- min(nrow(peak_bed), n_peaks)
                   peak_bed[1:keep_peaks,]
                 })

peak_beds <- peak_beds[type_pos$type]

bed_to_gr <- function(bed) {
  GRanges(seqnames = bed[[1]],
          IRanges(bed[[2]],
                  bed[[3]]),
          strand = "+")
}

peak_grs <- map(peak_beds,
                bed_to_gr)

peak_grl <- GRangesList(peak_grs)

all_peaks <- GenomicRanges::reduce(unlist(peak_grl))

gr_overlaps <- function(query_gr,target_gr) {
  ol <- findOverlaps(query_gr, target_gr,
                     ignore.strand = TRUE)
  unique(subjectHits(ol))
}

peak_overlaps <- map(peak_grs,
                     gr_overlaps,
                     all_peaks)

peak_names <- paste(seqnames(all_peaks), start(all_peaks), end(all_peaks), sep = "_")

peak_mat <- sparseMatrix(i = unlist(peak_overlaps),
                         p = c(0,cumsum(map_int(peak_overlaps, length))),
                         x = rep(1, sum(map_int(peak_overlaps, length))),
                         index1 = TRUE,
                         dims = c(length(all_peaks),
                                  length(peak_overlaps)),
                         dimnames = list(peak_names,
                                         type_pos$type))

unique_peaks <- peak_mat[rowSums(peak_mat) == 1,]
unique_counts <- colSums(unique_peaks)

peak_fullmat <- as(peak_mat, "matrix")

# Output for Supp Data
peak_fulldf <- as.data.frame(peak_fullmat) %>%
  mutate(hg38_region = rownames(.)) %>%
  select(hg38_region, everything())
fwrite(peak_fulldf,
       "Figure-3_Supp-Data-3_Peak_Intersections.csv.gz")

peak_sets <- apply(peak_fullmat, 1, paste0, collapse = "")

peak_df <- as.data.frame(table(peak_sets))
peak_df <- peak_df[order(peak_df$Freq, decreasing = TRUE),]

peak_set_mat <- matrix(
  as.integer(
    unlist(
      sapply(
        as.character(peak_df[["peak_sets"]]), 
        strsplit, split = "")
    )
  ),
  ncol = ncol(peak_fullmat),
  byrow = TRUE)

peak_df$n_types <- rowSums(peak_set_mat)

type_count_df <- type_pos %>%
  mutate(total_peaks = map_int(peak_beds, nrow),
         unique_peaks = unique_counts)

inter_max <- 200
inter_ymax <- 10
n_intersections <- 20

intersection_sets <- peak_df[peak_df$n_types != 1,] %>%
  head(n_intersections) %>%
  mutate(xpos = 1:n(),
         ymin = 0,
         ymax = Freq / inter_max * inter_ymax)


inter_mat <- peak_set_mat[peak_df$n_types != 1,][1:n_intersections,]
inter_points <- reshape2::melt(inter_mat) %>%
  rename(xpos = Var1,
         ypos = Var2)

inter_lines <- inter_points %>%
  filter(value == 1) %>%
  group_by(xpos) %>%
  summarise(ymin = min(ypos),
         ymax = max(ypos))

inter_axis <- data.frame(vals = seq(50, 200, 50)) %>%
  mutate(ypos = vals / inter_max * inter_ymax)
inter_baseline <- data.frame(ypos = 0)

ggplot() +
  # y-axis for bars
  geom_segment(data = inter_axis,
               aes(x = 0, xend = 20.5,
                   y = ypos, yend = ypos),
               linetype = "dashed") +
  # y-axis labels
  geom_text(data = inter_axis,
            aes(x = 0.4, y = ypos + 0.3,
                label = vals),
            hjust = 1, vjust = 0) +
  # intersection bars
  geom_rect(data = intersection_sets,
            aes(xmin = xpos - 0.4, xmax = xpos + 0.4,
                ymin = 0, ymax = ymax)) +
  # y-axis baseline
  geom_segment(data = inter_baseline,
               aes(x = 0, xend = 20.5,
                   y = ypos, yend = ypos)) +
  geom_text(data = inter_baseline,
            aes(x = 0.4, y = ypos + 0.3,
                label = 0),
            hjust = 1, vjust = 0) +
  geom_point(data = inter_points %>%
               filter(value == 1),
             aes(x = xpos, y = -1 * ypos),
             size = 2) +
  geom_segment(data = inter_lines,
            aes(x = xpos, xend = xpos,
                y = -1 * ymin, yend = -1 * ymax)) +
  geom_point(data = inter_points %>%
               filter(value == 0),
             aes(x = xpos, y = -1 * ypos),
             alpha = 0.2) +
  geom_text(data = type_count_df,
            aes(x = -0.5,
                y = -1 * xpos,
                label = unique_peaks),
            hjust = 1,
            vjust = 0.3) + 
  geom_text(data = type_count_df,
            aes(x = -2.5,
                y = -1 * xpos,
                label = type),
            hjust = 1,
            vjust = 0.3) +
  scale_x_continuous("",
                     limits = c(-10, n_intersections + 1)) +
  theme_void()

ggsave("Figure_3_upset.pdf",
       width = 3,
       height = 3)

# Total is 2500 except plasmablasts, which only have 592 peaks.