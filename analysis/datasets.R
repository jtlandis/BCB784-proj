
box::use(
  ./dependencies[...]
)

box::use(
  rtracklayer[import],
  readr[read_table],
  readxl[read_xlsx],
  stringr[str_extract]
)

is_hash <- function(x) {
  grepl(
    "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
    x)
}

#' @description
#' longread GTEX entire gtf file asis
#' 
#' @export
flair_gtf <- import(
  box::file("..", "data",
            "flair_filter_transcripts.gtf.gz")
)

#' @description
#' GTEX data as a data.frame
#' @export
gtex_data <- read_table(
  box::file("..", "data",
            "quantification_flair_filter.counts.txt.gz"),
  col_names = TRUE, show_col_types = FALSE
) |>
  mutate(
    key = str_extract(transcript, "^[^_]+"),
    .before = transcript
  )

#' @description
#' GRanges object of the `flair_gtf` object that only contains
#' the transcripts. The exons for each transcript is within the
#' `mcols(trans_data)$exons` and represents another GRanges as
#' a `CompressedGRanges` object
#' 
#' This becomes the rowRanges of `se_long` object
#' 
#' @export
trans_data <- local({
  .rle <- flair_gtf$type |>
    as.character() |>
    rle()
  transcripts <- filter(flair_gtf, type == "transcript")
  .exons <- filter(flair_gtf, type == "exon")
  # make grouping numbers
  .exons$transcript_group <- .rle$lengths[.rle$values=="exon"] |> {
    \(x) rep(seq_along(x), x)
  }()
  # split gene identifier in this subset
  exons <- split(.exons$gene_id, .exons$transcript_group)
  # quickly check that the order is the same
  stopifnot(
    all(transcripts$gene_id == purrr::map_chr(exons, unique))
  )
  transcripts$exons <- split(.exons, .exons$transcript_group)
  transcripts
}) |>
  mutate(
    key = str_extract(transcript_id, "^[^_]+")
  )

dup_key <- anyDuplicated(gtex_data$key)
dup_key <- gtex_data$key[dup_key]

ind <- which(trans_data$key %in% dup_key) # get the index
ord <- match(gtex_data$key, trans_data$key)
ind1 <- which(ord == ind[1]) #first index
trans <- trans_data[ord]
trans[ind1[2]] <- trans_data[ind[2]]

stopifnot(
  identical(trans$key, gtex_data$key)
)


#' @description
#' sample meta data. This is the colData of the `se_long` object
#' @export
col_data <- read_xlsx(
  box::file("..", "data", "metadata.xlsx"),
  sheet = "Suppl.Table1",
  range = "A1:T97"
)

col_slice <- match(select(gtex_data, -key, -transcript) |> colnames(),
                   col_data$sample_id)

stopifnot(
  identical(colnames(gtex_data)[-c(1:2)], col_data$sample_id[col_slice])
)


#' @name se_long
#' @description
#' The main dataset of long reads represented in a `PlySummarizedExperiment`,
#' which is a wrapper class of the `SummarizedExperiment` object that enables
#' `dplyr` syntax on the class. To access the `SummarizedExperiment` object,
#' use the `plyxp::se()` function.
#' 
#' @export
se_long <- SummarizedExperiment(
  assays = list(
    counts = gtex_data |> select(-key) %>% {
      rnms <- .$transcript
      x <- select(., -transcript)
      x <- as.matrix(x)
      rownames(x) <- rnms
      x
    }
  ),
  rowRanges = trans,
  colData = as(col_data[col_slice,], "DFrame")
) |>
  new_plyxp() |>
  group_by(cols(sample_name)) |>
  mutate(
    cols(RBP_kd = case_when(
      !grepl("Cells - Cultured fibroblasts", tissue) ~ NA, #only fibroblasts
      grepl("exp", sample_id) ~ TRUE, # exp for experimental purtebation
      any(grepl("exp", sample_id)) ~ FALSE, # all other items in group are false
      TRUE ~ NA # everything else should be NA
    ))
  ) |>
  ungroup()


#' @name se_long_subsets
#' @description
#' Subset of RBP knockdowns
#' @export
se_RBP_kd <- filter(se, cols( !is.na(RBP_kd) ))
# the above is a subset of the below.

#' @describeIn se_long_subsets only fibroblasts cells
#' @export
se_fibroblasts <- filter(se, cols(tissue == "Cells - Cultured fibroblasts"))

#' @describeIn se_long_subsets The 4 k562 cells used for quality control
#' @export
se_k562 <- filter(se, cols(sample_name == "K562"))

#' @describeIn se_long_subsets Donor tissues
#' @export
se_donor <- filter(se,
                cols(
                  !tissue %in% c("Cells - Cultured fibroblasts", "K562")
                )
)



#' @description
#' A subset summary (by rows) of unique TSS. 
#' TSS were grouped by transcription start site position and seqname
#' 
#' new transcripts were identified if they contained a 
#' hash identifier
#' 
#' The `encode_type` column of the rowData expresses the following
#' mannor that is similar to file permissions:
#' 1 -> novel: TSS group contains transcript with no associated gene parent
#'      and has hash
#' 2 -> ensbl novel: TSS group contains transcript with associated gene parent
#'      and has hash
#' 4 -> not novel: contains transcript that has no hash (assumed 
#'      to have a gene parent)
#'
#' Thus values between 1-7 are possible for the `encode_type`
#' 
#' @export
se_long_tss_groups <- local({
  .se <- plyxp_on(
    se_long, .on = rowRanges,
    mutate,
    .start = start,
    .seqnames = as.character(seqnames),
    new_transcript = is_hash(transcript_id),
    has_ensbl_parent = grepl("^ENSG", gene_id)
  ) |>
    group_by(rows(.start, .seqnames))
  rr <- .se |> {
    function(se) {
      row_groups <- group_data(se)$row_groups
      groups <- vctrs::list_unchop(
        as.list(row_groups$.indices_group_id),
        indices = row_groups$.indices)
      rr <- rowRanges(se(se))
      mcols(rr) <- NULL
      mcols(rr)$.exons <- rowData(se)$exons
      split(rr, groups)
    }
  }()
  
  .se <- .se |> 
    summarise(
      rows(
        n = n(),
        n_novel = sum(!has_ensbl_parent & new_transcript),
        n_novel_ensbl = sum(has_ensbl_parent & new_transcript),
        n_not_novel = sum(!new_transcript),
        across(c(gene_id, transcript_id), list)
      ),
      counts = list(counts)
    )
  rowData(.se)$.ranges <- rr
  ungrouop(.se) |>
    mutate(
      rows(
        encode_type = (n_novel>0) + (2 * (n_novel_ensbl>0)) + (4 * (n_not_novel>0))
      )
    )
})

