

box::use(
  ./dependencies[...],
  data = ./datasets[se_long, se_short],
  purrr[...],
  tidyr[pivot_wider],
  ggplot2[...]
)



se_short <- filter(data$se_short,
                   rows(transcript_type == "protein_coding"))
.common_trans <- intersect(
  rowData(data$se_long)$transcript_id,
  rowData(se_short)$transcript_id
)

se_long <- data$se_long[match(.common_trans,
                              rowData(data$se_long)$transcript_id),] |>
  mutate(cols(seq_source = "long_read"))
se_short <- data$se_short[match(.common_trans,
                                rowData(data$se_short)$transcript_id),] |>
  mutate(cols(seq_source = "short_read"),
         counts = as.integer(round(rsem)))


## transcripts -----

#' @export
plot_data_transcripts <- list(short = se_short,
                              long = se_long) |>
  lapply(
    mutate,
    rows(trans_width = lapply(width(exons), sum, 1) |>
           map_dbl(identity),
         #kilobases
         trans_kb = trans_width/1000),
    RPK = counts/.rows$trans_kb,
    cols(
      scaleFactor = colSums(.assays_asis$RPK)/1e6
    ),
    TPM = RPK/.cols$scaleFactor
  ) |>
  lapply(
    select,
    TPM,
    cols(sample_name, tissue, seq_source)
  ) |> {
    \(plot_data) {
      com <- lapply(plot_data,
             mutate,
             cols(i = interaction(sample_name, tissue, drop = TRUE))) |>
        lapply(pull, cols(i))
      common <- intersect(com$short, com$long)
      rownames(se(plot_data[[2]])) <- rownames(se(plot_data[[1]]))
      map2(
        plot_data,
        com,
        ~ .x[,match(common, .y)])
    }
  }() |>
  lapply(as_tibble) |>
  bind_rows() |>
  pivot_wider(
    id_cols = c(.features, sample_name, tissue),
    names_from = seq_source,
    values_from = TPM
  ) |>
  filter(if_all(ends_with("read"),
                ~ .x > 1)) |>
  mutate(across(ends_with("read"),
                ~ log10(.x), 
                .names = "{.col}_log10"))

## genes -----


se_short <- filter(data$se_short,
                   rows(transcript_type == "protein_coding")) |>
  plyxp_on(
    mutate, .on = rowRanges,
    combo = paste(gene_id, seqnames, strand)
  )
se_long <- plyxp_on(
  data$se_long, mutate, .on = rowRanges,
  combo = paste(gene_id, seqnames, strand)
)

.common_genes <- intersect(
  rowData(se_long)$combo,
  rowData(se_short)$combo
)

se_long <- filter(se_long, rows(combo %in% .common_genes)) |>
  mutate(cols(seq_source = "long_read"))
se_short <- se_short |>
  filter(rows(combo %in% .common_genes)) |>
  mutate(cols(seq_source = "short_read"),
         counts = as.integer(round(rsem)))

# these should be aligned in terms of the transcripts!
#' @export
plot_data_genes <- list(short = se_short,
                        long = se_long) |> 
  lapply(
    arrange,
    rows(gene_id)
  ) |>
  lapply(
    mutate,
    rows(
      .start = lapply(start(exons), min) |> map_int(identity),
      .end = lapply(end(exons), max) |> map_int(identity),
      .seqnames = { seqs <- exons@unlistData@seqnames |>
        as.character()
      map2(start(exons@partitioning),
           end(exons@partitioning),
           function(i,j,x) x[c(i,j)], x = seqs) |>
        map_chr(unique)
      },
      .strand = {
        strnd <- exons@unlistData@strand |> as.character()
        map2(
          start(exons@partitioning),
          end(exons@partitioning),
          function(i,j,x) x[c(i,j)], x = strnd) |>
          map_chr(unique)
      }
    )
  ) |>
  lapply(
    group_by,
    rows(gene_id)
  ) |>
  lapply(
    summarise,
    counts = colSums(counts),
    rows(
      across(matches("\\.[se]"), list),
      .features = unique(gene_id)
    )
  ) |>
  lapply(ungroup) |>
  lapply(
    mutate,
    rows(
      .start = map_int(.start, min),
      .end = map_int(.end, max),
      across(c(.seqnames, .strand), ~map_chr(.x, unique)),
      gene_range = GRanges(.seqnames, IRanges::IRanges(start = .start,
                                                       end = .end),
                           strand = .strand),
      across(matches("\\.[se]"), ~ NULL)
    )
  ) |>
  lapply(
    mutate,
    rows(
      #kilobases
      trans_kb = width(gene_range)/1000),
    RPK = counts/.rows$trans_kb,
    cols(
      scaleFactor = colSums(.assays_asis$RPK)/1e6
    ),
    TPM = RPK/.cols$scaleFactor
  ) |> {
    \(x) {
      .common <- lapply(x, rowData) |>
        lapply(`[[`, "gene_id")
      .common <- intersect(.common[[1]], .common[[2]])
      lapply(
        x,
        function(x, g) x[match(g, rownames(se(x))),],
        g = .common
      )
    } 
  }() |>
  lapply(
    select,
    TPM,
    cols(sample_name, tissue, seq_source)
  ) |> {
    \(plot_data) {
      com <- lapply(plot_data,
                    mutate,
                    cols(i = interaction(sample_name, tissue, drop = TRUE))) |>
        lapply(pull, cols(i))
      common <- intersect(com$short, com$long)
      rownames(se(plot_data[[2]])) <- rownames(se(plot_data[[1]]))
      map2(
        plot_data,
        com,
        ~ .x[,match(common, .y)])
    }
  }() |>
  lapply(as_tibble) |>
  bind_rows() |>
  pivot_wider(
    id_cols = c(.features, sample_name, tissue),
    names_from = seq_source,
    values_from = TPM
  ) |>
  filter(if_all(ends_with("read"),
                ~ .x > 1)) |>
  mutate(across(ends_with("read"),
                ~ log10(.x), 
                .names = "{.col}_log10"))

## plot function ----

#' @export
plot_figc <- function(plot_data) {
  ggplot(plot_data, aes(short_read_log10, long_read_log10)) +
    geom_hex(binwidth = .125) +
    scale_fill_gradient(low = "grey95", high = "black") +
    geom_abline(color = "limegreen", slope = 1, intercept = 0) +
    theme_classic(base_size = 22) +
    labs(y = "log10 ONT TPM", x = "log10 Illumina TPM") +
    guides(fill = guide_colorbar(
      legend.title = element_text(vjust = .7),
      position = "bottom", theme = theme(
        legend.title = element_text(vjust = .7),
        legend.text = element_text(size = 12, angle = 45, vjust = 1)))) +
    coord_equal(xlim = c(0,5), ylim = c(0,5))
}
  