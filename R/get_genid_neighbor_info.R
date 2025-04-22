# Load required utilities and packages
# install.packages("devtools")
# library(devtools)
# devtools::install_github("rdev-create/rutils@v0.0.0.9001")
rutils::library_install("jsonlite")
rutils::library_install("logger")
rutils::library_install("rentrez")

library(jsonlite)
library(logger)
library(rutils)
library(rentrez)

source("R/NCBI_utils.r")

get_geneid_neighbor_info <- function(shortname) {
  # Global Switches
  use_cache <- TRUE # Set to TRUE to enable caching, FALSE to disable caching

  # Dataset Configuration
  do_small_dataset <- FALSE # TRUE for small dataset, FALSE for large
  output_folder <- if (do_small_dataset) "data_small" else "data"
  num_allgeneids <- if (do_small_dataset) 10 else "all"

  # File paths for caching and outputs
  cached_neighbor_geneids_file <- file.path(output_folder, "neighbor_gene_ids_cache.json")
  neighbor_fasta_cache_file <- file.path(output_folder, "neighbor_gene_fasta_cache.json")
  geneid_sequence_list_file <- file.path(output_folder, "geneid_sequence_list.json")
  neighbor_fastas_filename <- file.path(output_folder, "neighbor_protein_fasta_sequences.fa")

  # Create output folder if it doesn't exist
  rutils:::make_dir_if_not_exist(dir_path = output_folder)

  # --- Retrieve or Cache Neighbor Gene IDs ---
  if (use_cache && file.exists(cached_neighbor_geneids_file)) {
    log_info("Loading cached neighbor gene ids from {cached_neighbor_geneids_file}")
    neighbor_gene_ids <- read_json(cached_neighbor_geneids_file, simplifyVector = TRUE)
  } else {
    all_geneids <- get_gene_ids_of_shortname(shortname = shortname)
    n <- if (identical(num_allgeneids, "all")) length(all_geneids) else min(num_allgeneids, length(all_geneids))
    gene_infos <- get_gene_list_info(all_geneids[seq_len(n)])

    neighbor_gene_ids <- list()
    for (i in seq_along(gene_infos)) {
      log_info("Processing gene info {i} of {length(gene_infos)}")
      gene_neighbor_infos <- get_gene_neighbors(gene_infos[[i]]$uid)
      for (info in gene_neighbor_infos) {
        neighbor_gene_ids <- c(neighbor_gene_ids, info$uid)
      }
    }

    write_json(neighbor_gene_ids, cached_neighbor_geneids_file)
  }

  # --- Load or Initialize FASTA Cache ---
  if (use_cache && file.exists(neighbor_fasta_cache_file)) {
    log_info("Loading cached FASTA sequences from {neighbor_fasta_cache_file}")
    fasta_cache <- read_json(neighbor_fasta_cache_file, simplifyVector = TRUE)
  } else {
    fasta_cache <- list()
  }

  # --- Retrieve and Cache FASTA Sequences for Neighbor Gene IDs ---
  for (i in seq_along(neighbor_gene_ids)) {
    gene_id <- neighbor_gene_ids[[i]]

    # If caching is enabled and the FASTA is already cached, skip fetching
    if (use_cache && !is.null(fasta_cache[[gene_id]])) {
      log_info("FASTA for gene {gene_id} already cached, skipping fetch")
      next
    }

    log_info("Fetching FASTA for gene {gene_id} ({i} of {length(neighbor_gene_ids)})")
    fasta <- get_protein_fasta_for_gene(gene_id)

    if (is.null(fasta)) {
      log_warn("No FASTA sequence found for gene ID: {gene_id}")
      next
    }
    if (!is.character(fasta)) {
      log_warn("Invalid FASTA sequence for gene ID: {gene_id}")
      next
    }

    # Cache the fetched FASTA sequence (update file immediately if caching is enabled)
    fasta_cache[[gene_id]] <- fasta
    write_json(fasta_cache, neighbor_fasta_cache_file)
  }

  # --- Write FASTA Sequences to File and Build Gene ID -> Sequence Mapping ---
  neighbor_fastas_file <- file(neighbor_fastas_filename, open = "w")
  geneid_sequence_list <- list()

  for (gene_id in neighbor_gene_ids) {
    fasta <- fasta_cache[[gene_id]]

    if (is.null(fasta)) {
      log_warn("No cached FASTA for gene ID: {gene_id}, skipping writing to file")
      next
    }

    # Write the FASTA sequence to the file
    writeLines(fasta, neighbor_fastas_file)

    # Extract the sequence (excluding the header) and add to mapping
    fasta_lines <- strsplit(fasta, "\n")[[1]]
    sequence <- paste(fasta_lines[-1], collapse = "")
    geneid_sequence_list[[sequence]] <- gene_id
  }
  close(neighbor_fastas_file)

  # Save the gene ID to sequence mapping as JSON
  write_json(geneid_sequence_list, geneid_sequence_list_file)

  log_info("Processing complete, run InterProScan then move to final logic pt 2")
}
