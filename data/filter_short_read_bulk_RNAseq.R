#################################
### DO NOT RUN 
### ONLY HERE TO SHOW
### HOW TO GENERATE RDS FILE
#################################


.file <- box::file(
  "GTEx_Analysis_v10_RSEMv1.3.3_transcripts_expected_count.txt.gz"
)

x <- data.table::fread(.file) |>
  tibble::as_tibble()
col_meta <- readr::read_csv(box::file("tissue_meta_map.csv"))

sample_meta <- dplyr::select(col_meta, -sample_name) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(tissue_id))


# get highly expressed genes

sample_key <- colnames(x) |>
  stringr::str_match("GTEX-[^-]+-([^-]+-?.*)-SM") |>
  _[,2] 


x <- x[c(1:2, match(sample_key, sample_meta$tissue_id) |>
           {Negate(is.na)}() |>
           which())]
sample_key <- colnames(x) |>
  stringr::str_match("GTEX-[^-]+-([^-]+-?.*)-SM") |>
  _[,2] 
col_data <- sample_meta[match(sample_key[-c(1:2)], sample_meta$tissue_id),]
track <- rtracklayer::import(
  box::file("gencode.v26.chr_patch_hapl_scaff.annotation.gff3.gz"))
trans <- plyranges::filter(track, type == "transcript")

# GRanges does not support NA indexes, So I will
# instead subset x by the indices that contain a
# match in the trans$ID column
x_sub <- x[x$transcript_id %in% trans$ID,]

trans <- trans[trans$ID %in% x_sub$transcript_id]
row_data <- trans[match(x_sub$transcript_id, trans$ID)]

## read in minimal ranges from long-read data

lr <- readr::read_csv(
  box::file("ranges_meta_map.csv")
) |> {
  \(x) {
    GenomicRanges::GRanges(
      seqnames = x$seqnames,
      ranges = IRanges::IRanges(
        start = x$start, end = x$end
      ),
      strand = x$strand,
      gene_id = x$gene_id
    )
  }
}()

# retain set that is "old"
lr_old <- filter(lr, grepl("^ENSG", gene_id))
row_data_old <- row_data |>
  filter(gene_id %in% .env$lr_old$gene_id)
# find possible genes arround start site
lr_new <- filter(lr, !grepl("^ENSG", gene_id)) |>
  GenomicRanges::resize(width = 1L) |>
  GenomicRanges::resize(width = 2000L, fix = "center")
g_id <- lr_new$gene_id
lr_new$gene_id <- NULL
row_data_new <- plyranges::find_overlaps(row_data, lr_new)
# merge these two together
row_data <- unique(c(row_data_old, row_data_new))

# re-match rows of x to new GRanges
x_sub <- x_sub[match(row_data$ID, x_sub$transcript_id),]

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    rsem = {\(x) {
      rn <- x$transcript_id
      x <- x[-(1:2)]
      x <- as.matrix(x)
      rownames(x) <- rn
      x
    }}(x_sub)
  ),
  rowData = row_data,
  colData = col_data
) 
se$sample_name <- stringr::str_extract(colnames(se), "GTEX-[^-]+")
saveRDS(se, box::file("bulk_RNAseq.rds"))
se_small <- se[,se$sample_name %in% col_meta$sample_name]
# google download link ~ 150MB:
# https://drive.google.com/file/d/1GMOuWT3ndUYDLRW5rML8TGn2Zg7vCpgq/view?usp=sharing
saveRDS(se_small, box::file("bulk_RNAseq_small.rds"))
