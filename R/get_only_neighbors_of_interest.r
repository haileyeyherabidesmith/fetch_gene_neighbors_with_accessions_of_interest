library(jsonlite)
source("utils.r")

full_details_genes_of_interest_file <- file.path(output_folder, "full_details_genes_of_interest.json")
final_genes_of_interest_file <- file.path(output_folder, "final_genes_of_interest.json")
geneid_sequence_list_file <- file.path(output_folder, "geneid_sequence_list.json")
accessions_of_interest_file <- file.path(output_folder, "interprolocal_accessions_nuclease.json")
neighbor_protein_ipr_info_file <- file.path(output_folder, "neighbor_protein_fasta_sequences.fa.json")

gene_id_sequences <- read_json(geneid_sequence_list_file, simplifyVector = TRUE)
accessions_of_interest <- read_json(accessions_of_interest_file, simplifyVector = TRUE)
arfb_neighbors_info <- read_json(neighbor_protein_ipr_info_file, simplifyVector = FALSE)

interesting_accessions_list <- list()

for (i in seq_along(arfb_neighbors_info$results)) {
  # Extract the sequence and matches for each result
  result <- arfb_neighbors_info$results[[i]]
  result_sequence <- result$sequence # Extract the sequence at the result level
  for (match_i in seq_along(result$matches)) {
    match <- result$matches[[match_i]]
    # Extract the signature accession
    if (!is.null(match$signature)) {
      if (!is.null(match$signature$accession)) {
        # If entry is not NULL, extract the accession from entry
        if (match$signature$accession %in% accessions_of_interest) {
          gene_id <- gene_id_sequences[[result_sequence]]
          # Build an object as a named list
          match_info <- list(geneid = gene_id,
                             accession = match$signature$accession,
                             sequence = result_sequence)
          # Append the object to the overall list
          interesting_accessions_list <- c(interesting_accessions_list, list(match_info))
        }
      }
      if (!is.null(match$signature$entry) && !is.null(match$signature$entry$accession)) {
        # If entry is not NULL, extract the accession from entry
        if (match$signature$entry$accession %in% accessions_of_interest) {
          gene_id <- gene_id_sequences[[result_sequence]]
          # Build an object as a named list
          match_info <- list(geneid = gene_id,
                             accession = match$signature$entry$accession,
                             sequence = result_sequence)
          # Append the object to the overall list
          interesting_accessions_list <- c(interesting_accessions_list, list(match_info))
        }
      }
    }
  }
}

write_json(interesting_accessions_list, full_details_genes_of_interest_file)
gene_ids_of_interest <- sapply(interesting_accessions_list, function(x) x$geneid)
write_json(gene_ids_of_interest, final_genes_of_interest_file)
