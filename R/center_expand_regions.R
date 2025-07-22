#' Modifies genomic regions by centering and then expanding them
#'
#' @description
#' [peakCombiner::center_expand_regions] is an optional step that re-defines the
#' genomic regions by expanding them from their center. The center information
#' has to be stored in the input data column `center`, while the information for
#' the expansion can either be user provided or input data derived. The accepted
#' input is a data frame created from [peakCombiner::prepare_input_regions].
#' Please see [peakCombiner::prepare_input_regions] for more details.
#'
#' @details
#' This is an optional function that resizes the genomic regions based on the
#' input peakCombiner standard data frame and the options you select. An
#' expected input data foam contains the following columns with the names:
#' `chrom`, `start`, `end`, `name`, `score`, `strand`, `center`, `sample_name`.
#' Such a data frame is created by the script
#' [peakCombiner::prepare_input_regions]. This step is useful if you want all of
#' your peaks to be the same size for your downstream analyses. In addition, if
#' you want to use the "summit" information, normally obtained by some peak
#' callers (e.g., Macs2), this function allows you to automatically center your
#' regions of interest on these summits. This enables you to capture
#' information about the most important region within a genomic region (e.g.,
#' TF-binding site or highest peak) and put that region in the center of your
#' downstream analyses (e.g., applicable to motif-finding or "heatmaps"
#' summarizing multiple genomic regions).
#'
#' There are two concepts that are relevant for
#' [peakCombiner::center_expand_regions]: how to define the center, and how much
#' to expand from the center.
#'
#' ## How to define the center?
#'
#' When you prepared your input regions, it is recommended to use the function
#' [peakCombiner::prepare_input_regions] provided by this package. This pre-
#' populated the `center` column with the absolute genomic coordinate of the
#' center of the peak region. You can either choose to define the center by
#' using pre-defined summit information (e.g., obtained from a peak caller like
#' MACS2) or re-compute the arithmetic mean and save that value in the column
#' `center`. (For details see the help for
#' [peakCombiner::prepare_input_regions]).
#'
#' ## How much to expand from the center
#' You can choose to expand the genomic region from the center either
#' symmetrically or asymmetrically (different lengths before and after the
#' center position).
#'
#' In the symmetrical case, if you want to choose the size of your genomic
#' region based on the input data, this function can also calculate the median
#' peak size across all of your genomic regions and use that value (`expand_by`
#' = NULL). Alternatively, the user is free to provide a numeric vector to
#' define the expansion. A numeric vector with one value is used to
#' symmetrically expand, while a vector with two values allows to expand
#' asymmetrically.
#'
#'
#' @param data        PeakCombiner data frame structure with required columns
#'                      named `chrom`, `start`, `end`, `name`,
#'                      `score`, `strand`, `center`, `sample_name`. Additional
#'                      columns will be maintained.
#' @param center_by   Allowed values are 'center_column' (default) or
#'                    'midpoint'.
#' * 'center_column' uses the value stored in the column `center` to center.
#' * 'midpoint' replaces the value stored in the column `center` based on the 
#' [GenomicRanges::resize] followed by the expansion from based on the user 
#' input.
#'
#' @param expand_by   Allowed values a numeric vector of length 1 or 2,
#'                      or 'NULL' (default).
#' * The value from the numeric vector of length 1
#'                          is expanded in both directions from center to define
#'                          the genomic region.
#'                          Thus, the size of the resulting genomic region is 2x
#'                          the provided value + 1 (for the center coordinate).
#' * The value of the numeric vector of length 2
#'                          subtracts the first value from the center and adds
#'                          the second value to the center to define the genomic
#'                          region. Thus, the size of the genomic regions is
#'                          the sum of the first value + the second value
#'                          + 1 (for the center coordinate).
#' * 'NULL' allows for data-driven definition of the
#'                          `expand_by` value. It calculates the median
#'                          genomic region size of the input data and uses this
#'                          value like a length 1 numeric vector for expansion.
#'            
#' @param output_format Character value to define format of output object. 
#'                      Accepted values are "GenomicRanges" (default), "tibble" 
#'                      or "data.frame".  
#'
#' @param show_messages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the [peakCombiner::prepare_input_regions] Details.
#' Use as input for functions [peakCombiner::filter_regions()] and
#' [peakCombiner::combine_regions()].
#'
#' @export
#'
#' @importFrom rlang .data
#' @import tidyr
#' @import here
#'
#' @examples
#' # Load in and prepare a an accepted tibble
#' utils::data(syn_data_bed)
#'
#' # Prepare input data
#' data_prepared <- prepare_input_regions(
#'   data = syn_data_bed,
#'   output_format = "tibble",
#'   show_messages = TRUE
#' )
#' # Run center and expand
#' data_center_expand <- center_expand_regions(
#'   data = data_prepared,
#'   center_by = "center_column",
#'   expand_by = NULL,
#'   output_format = "tibble",
#'   show_messages = TRUE
#' )
#'
#' data_center_expand
#'
#' # You can choose to use the midpoint and predefined values to expand
#'
#' data_center_expand <- center_expand_regions(
#'   data = data_prepared,
#'   center_by = "midpoint",
#'   expand_by = c(100, 600),
#'   output_format = "tibble",
#'   show_messages = FALSE
#' )
#'
#' data_center_expand
#'
center_expand_regions <- function(data,
                                  center_by = "center_column",
                                  expand_by = NULL,
                                  output_format = "GenomicRanges",
                                  show_messages = TRUE) {
  
  ### -----------------------------------------------------------------------###
  ### Check if output format is valid
  ### -----------------------------------------------------------------------###
  
  if (output_format %in% c("GenomicRanges", 
                           "GRanges", 
                           "tibble", 
                           "data.frame", 
                           "data.table")) {
    cli::cli_inform(c(
      "i" = "Argument {.arg output_format} is set to {.val {output_format}}."
    ))
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Argument {.arg output_format} has to be one of the following
      values: {.val GenomicRanges}, {.val tibble}, or {.val data.frame}.",
      "i" = "Provided value is {.val {output_format}}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Figure out what kind of input data was entered by the user and
  ### load the initial data for follow-up quality checks
  ### -----------------------------------------------------------------------###
  
  required_colnames <- c(
    "chrom", "start", "end", "sample_name"
  )
  
  if (inherits(data, "GRanges")) {
    cli::cli_inform(c(
      "!" = "Provided input {.arg data} is a class {.cls GRanges} and will be
      converted to class {.cls tibble}.",
      ">" = "Start converting and preparing data."
    ))
    
    data_filtered <-
      tibble::as_tibble(data) |>
      dplyr::rename(chrom = .data$seqnames) |>
      dplyr::mutate(
        start = as.numeric(.data$start),
        end = as.numeric(.data$end),
        strand = as.character(.data$strand)
      ) |>
      dplyr::mutate(strand = ifelse(.data$strand == "*", ".", .data$strand))
  } else if (all(required_colnames %in% colnames(data))) {
    cli::cli_inform(c(
      "i" = "Provide input {.arg data} is a {.cls data.frame} with three or four
      columns and paths to existing files.",
      ">" = "Start loading and preparing data."
    ))
    
    data_filtered <- data
    
  } else if (all(required_colnames %in% colnames(data))) {
    data_filtered <- data
    
    cli::cli_inform(c(
      "i" = "Provide input {.arg data} is a pre-loaded {.cls data.frame}  with
      the required column names.",
      ">" = "Start preparing data."
    ))
  } else {
    # show error independend of show_messages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Provide input {.arg data} does not have the required format.",
      "!" = "Please check your column names in {.arg data}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Show or hide messages
  ### -----------------------------------------------------------------------###

  if (!is.logical(show_messages)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg show_messages} has to be {.cls logical}."
    ))
  } else if (isTRUE(show_messages)) {
    options("rlib_message_verbosity" = "default")
  } else if (isFALSE(show_messages)) {
    options("rlib_message_verbosity" = "quiet")
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg show_messages} is a non-accepted {.cls logical}
      value.",
      "i" = "Argument {.arg show_messages} is {.val {show_messages}}."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Prepare parameters
  ### -----------------------------------------------------------------------###
  center_values <- c("center_column", "midpoint")

  ## Check parameter value correctness and calculate if needed
  expansion_value <- define_expansion(
    data = data,
    expand_by = expand_by
  )

  ## Calculate the values to expand the regions
  length_expansion_value <- length(expansion_value)
  expand_1 <- expansion_value[1]
  expand_2 <- expansion_value[length_expansion_value]

  ### -----------------------------------------------------------------------###
  ### Check input parameters
  ### -----------------------------------------------------------------------###

  if (is.null(center_by)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg center_by} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg center_by} is {.val NULL}."
    ))
  } else if (length(center_by) != 1) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg center_by} has a length of {length(center_by)}.",
      "i" = "{.arg center_by} allowed length is 1."
    ))
  } else if (!tolower(center_by) %in% center_values) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg center_by} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg center_by} is {.val {center_by}}."
    ))
  } else if (tolower(center_by) %in% center_values) {
    ## good values!
    center_by <- tolower(center_by)
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg center_by} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg center_by} is {.val NULL}."
    ))
  }

  ## Check the validity of the peakCombiner input data format

  data <- check_data_structure(
    data = data
  )

  ### -----------------------------------------------------------------------###
  ### Center and expand
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    ">" = "Genomic regions will be centered and expanded.",
    " " = " "
  ))

  if (center_by == "center_column") {
    cli::cli_inform(c(
      ">" = "Starting with expanding genomic regions from the column {.field
      center}."
    ))

    data_center_expand <-
      data |>
      dplyr::mutate(
        start = .data$center - !!expand_1,
        end = .data$center + !!expand_2
      ) |>
      dplyr::ungroup()
  } else if (center_by == "midpoint") {
    cli::cli_inform(c(
      ">" = "Starting with defining the {.field center} of the regions from the
      {.field start} and {.field end} coordinates."
    ))
  
    data 
    
    # Convert to GRanges
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df = data,
      keep.extra.columns = TRUE
    )
    
    # Resize (e.g., width = 500, fix = "center")
    gr_resized <- GenomicRanges::resize(gr, width = 2, fix = "center")
    
    # Convert back to tibble
    data_center_expand <- dplyr::as_tibble(as.data.frame(gr_resized)) |>
      dplyr::select(chrom = .data$seqnames, 
                    .data$start, 
                    .data$end,
                    .data$name, 
                    .data$score, 
                    .data$strand, 
                    .data$center, 
                    .data$sample_name) |>
      dplyr::mutate(strand = ".",
                    start = .data$center - !!expand_1,
                    end = .data$center + !!expand_2) |>
      dplyr::ungroup()
      
    # View result
    print(data_center_expand)
    
  }

  cli::cli_inform(c(
    "v" = "Genomic regions were successfully centered and expanded.",
    " " = " "
  ))

  if (any(data_center_expand$start < 1)) {
    neg_starts <-
      data_center_expand |>
      dplyr::filter(.data$start < 1) |>
      dplyr::pull(.data$name)

    cli::cli_inform(c(
      "i" = "Some newly-defined genomic regions have a {.field start}
      coordinate below {.val 1}.",
      ">" = "Values of {.field name} for these  site{?s}: {.val {neg_starts}}."
    ))

    data_center_expand <-
      data_center_expand |>
      dplyr::mutate(start = ifelse((.data$start < 1), 1, .data$start))

    cli::cli_inform(c(
      "v" = "These genomic regions were truncated to get {.field start}
      coordinate {.val 1}."
    ))

    rm(neg_starts)
  }

  ### -----------------------------------------------------------------------###
  ### Return data frame
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    "v" = "Genomic regions were successfully centered and expanded.",
    " " = " "
  ))
  
  ### -----------------------------------------------------------------------###
  ### Adjust output format
  ### -----------------------------------------------------------------------###
  
  if (output_format %in% c("GenomicRanges", "GRanges")) {
    cli::cli_inform(c(
      "i" = "Output format is set to {.val {output_format}}.")
    )
    
    data_center_expand <- 
      data_center_expand |>
      GenomicRanges::makeGRangesFromDataFrame(
        keep.extra.columns = TRUE,
      )
    
  } else if (output_format %in% c("tibble", "data.frame", "data.table")) {
    cli::cli_inform(c(
      "i" = "Output format is set to {.val tibble}."
    ))
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Argument {.arg output_format} has to be one of the following
      values: {.val GenomicRanges}, {.val tibble}, or {.val data.frame}.",
      "i" = "Provided value is {.val {output_format}}."
    ))
  } 
  
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###

  if (isFALSE(show_messages)) {
    options("rlib_message_verbosity" = "default")
  }


  return(data_center_expand)
}
