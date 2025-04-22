# Load required utilities and packages
# install.packages("devtools")
# library(devtools)
# devtools::install_github("rdev-create/rutils")
rutils::library_install("jsonlite")
rutils::library_install("logger")
rutils::library_install("rentrez")
rutils::library_install("httr")

library(jsonlite)
library(logger)
library(rutils)
library(rentrez)
library(httr)

# Function to get Pfam Accession Number from InterPro for a given exact short name
get_pfamid_shortname_pair <- function(short_name) {
  # Base URL for InterPro API search
  search_url <- "https://www.ebi.ac.uk/interpro/api/entry/pfam/"

  # Step 1: Search for potential matches
  retry_get <- function(search_url, short_name, retry_attempts = 3, delay = 1) {
    attempt <- 1
    response <- NULL

    while (attempt <= retry_attempts) {
      tryCatch(
        {
          response <- GET(search_url, query = list(search = short_name))
          if (status_code(response) == 200) {
            break
          } else {
          }
        },
        error = function(e) {
        }
      )

      attempt <- attempt + 1
      Sys.sleep(delay) # Wait before retrying
    }

    if (is.null(response)) {
      stop("All retry attempts for ", short_name, " failed.")
      next
    }

    return(response)
  }
  tryCatch({
    response <- retry_get(search_url, short_name, retry_attempts = 5, delay = 1)
    # Check for successful response
    if (status_code(response) == 200) {
      break
    } else {
      next
    }
  })
  # Parse the JSON response
  data <- fromJSON(content(response, "text", encoding = "UTF-8"))

  # Check if there are results
  if (is.null(data$results) || length(data$results) == 0) {
    return("No results found for the given short name.")
  }

  # Step 2: Loop through results and check `short_name` via additional API requests
  for (entry in data$results) {
    pfam_ids <- entry$accession # Extract accession(s), may be a list

    # Ensure pfam_ids is a character vector
    pfam_ids <- as.character(pfam_ids)

    for (pfam_id in pfam_ids) {
      # Construct API URL to fetch details for this specific Pfam entry
      details_url <- paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/", pfam_id)

      retry_get_details <- function(details_url, retry_attempts = 3, delay = 1) {
        attempt <- 1
        details_response <- NULL

        while (attempt <= retry_attempts) {
          tryCatch(
            {
              details_response <- GET(details_url)
              if (status_code(details_response) == 200) {
                break # Success, exit the loop
              } else {
              }
            },
            error = function(e) {
            }
          )

          attempt <- attempt + 1
          Sys.sleep(delay) # Wait before retrying
        }

        if (is.null(details_response)) {
          stop("All retry attempts failed.")
        }

        return(details_response)
      }
      details_response <- retry_get_details(details_url)

      # Check for valid response
      if (status_code(details_response) == 200) {
        details_data <- fromJSON(content(details_response, "text", encoding = "UTF-8"))

        # Extract short_name from correct nested structure
        extracted_short_name <- tryCatch(
          {
            details_data$metadata$name$short
          },
          error = function(e) NULL
        ) # Handle missing fields safely

        # Check if short_name matches exactly
        if (!is.null(extracted_short_name) && extracted_short_name == short_name) {
          pfamid_shortname_pair <- paste(pfam_id, ": ", short_name)
          return(pfamid_shortname_pair)
        }
      }
    }
  }
  return(NULL) # nolint: return_linter.
}

# uses the short names of a gene to get its pfam id
get_pfamid_shortname_pairs <- function(short_names) {
  pfamid_shortname_pairs <- list()
  for (short_name in short_names) {
    pfamid_shortname_pair <- get_pfamid_shortname_pair(short_name)
    pfamid_shortname_pairs <- c(pfamid_shortname_pairs, pfamid_shortname_pair)
  }
  return(pfamid_shortname_pairs) # nolint: return_linter.
}

get_gene_ids_of_shortname <- function(batch_size = 500, retmax = 10, shortname = shortname) {
  # Initial search to get the total number of results
  initial_search <- entrez_search(db = "gene", term = paste(shortname, "[Gene Name]"), retmax = retmax)
  total_results <- initial_search$count

  if (total_results == 0) {
    stop(c("No ", shortname, " genes found."))
  }


  # Initialize vector to store gene IDs
  all_gene_ids <- c()

  # Loop to fetch results in batches
  for (start in seq(0, total_results - 1, by = batch_size)) {

    batch_results <- entrez_search(
      db = "gene",
      term = shortname,
      retmax = batch_size,
      retstart = start
    )

    all_gene_ids <- c(all_gene_ids, batch_results$ids)

    # Sleep to avoid overloading NCBI servers
    Sys.sleep(0.1)
  }
  return(all_gene_ids)
}

get_fasta_for_gene <- function(gene_id) {
  # Get gene summary information
  entrez_summary <- entrez_summary(db = "gene", id = gene_id)

  # Extract gene's genomic information (start and end positions)
  genomic_info <- entrez_summary$genomicinfo
  chraccver <- genomic_info$chraccver # Chromosome accession version
  start <- genomic_info$chrstart # Start position of the gene
  end <- genomic_info$chrstop # End position of the gene
  description <- entrez_summary$description # Gene description


  # Fetch the sequence of the chromosome containing the gene
  fasta_sequence <- entrez_fetch(
    db = "nucleotide",
    id = chraccver,
    rettype = "fasta",
    retstart = 0, # Fetch the entire chromosome sequence
    retmax = 1000000 # Fetch a large enough portion to cover the gene
  )

  # Check if the sequence was fetched properly
  if (is.null(fasta_sequence)) {
    stop("Failed to fetch sequence for chromosome: ", chraccver)
  }

  # Remove the header line from the FASTA sequence and concatenate the sequence
  sequence <- strsplit(fasta_sequence, "\n")[]
  sequence <- paste(sequence[-1], collapse = "") # Remove header and collapse sequence

  # Ensure that the sequence length is large enough to include the gene
  if (end > nchar(sequence)) {
    stop("Gene end position exceeds sequence length for chromosome: ", chraccver)
  }

  # Slice the gene sequence from the chromosome
  gene_sequence <- substr(sequence, start, end)

  # Format the header as requested, including description
  fasta_header <- paste(">", chraccver, ":", start, "-", end, " ", description, sep = "")

  # Break the gene sequence into 60-character lines
  sequence_lines <- strsplit(gene_sequence, "")[] # Split sequence into individual characters
  sequence_chunks <- split(sequence_lines, ceiling(seq_along(sequence_lines) / 60)) # Group into chunks of 60

  # Combine the chunks into a properly formatted sequence string
  formatted_sequence <- paste(sapply(sequence_chunks, paste, collapse = ""), collapse = "\n")

  # Combine header and sequence to get the full FASTA result
  formatted_fasta <- paste(fasta_header, "\n", formatted_sequence, sep = "")

  return(formatted_fasta) # nolint: return_linter.
}

get_protein_fasta_for_gene <- function(gene_id) {
  # Step 1: Retrieve the protein accession linked to the gene using entrez_link
  link_result <- tryCatch(
    {
      entrez_link(dbfrom = "gene", id = gene_id, db = "protein")
    },
    error = function(e) {
      return(NULL) # nolint: return_linter.
    }
  )

  # Check if there are any linked protein accessions
  if (is.null(link_result) || length(link_result$links) == 0) {
    return(NULL)
  }

  # Extract the protein accession (first one in the list)
  protein_accession <- link_result$links$gene_protein[]

  # Step 2: Fetch the protein sequence using the protein accession
  protein_fasta <- tryCatch(
    {
      entrez_fetch(
        db = "protein",
        id = protein_accession,
        rettype = "fasta"
      )
    },
    error = function(e) {
      return(NULL) # nolint: return_linter.
    }
  )
  # Return the protein sequence in FASTA format
  return(protein_fasta) # nolint: return_linter.
}

get_gene_list_info <- function(gene_ids, batch_size = 300) {
  total_results <- length(gene_ids)
  # Initialize a list to store all gene information
  gene_infos <- list()

  # Loop through each batch
  for (start in seq(1, total_results, by = batch_size)) {
    # Get the current batch of gene IDs
    current_batch <- gene_ids[start:min(start + batch_size - 1, total_results)]

    # Fetch summaries for the current batch
    retry_attempts <- 3
    attempt <- 1
    summaries <- NULL

    while (attempt <= retry_attempts) {
      tryCatch(
        {
          summaries <- entrez_summary(db = "gene", id = current_batch)
          if (!is.null(summaries)) {
            break
          }
        },
        error = function(e) {
        }
      )

      attempt <- attempt + 1
    }

    if (is.null(summaries)) {
      stop("Failed to fetch gene summaries.")
    }

    # Extract desired fields with error handling
    batch_gene_infos <- lapply(summaries, function(summary) {
      if (!is.list(summary) || summary$status == "discontinued") {
        return(NULL)
      }
      # Check that all required fields exist
      if (is.null(summary$uid) || is.null(summary$name) || is.null(summary$organism)) {
        return(NULL)
      }

      # Return the extracted fields
      list(
        uid = summary$uid,
        name = summary$name,
        organism = summary$organism # Organism info (contains taxid, scientificname, commonname)
      )
    })

    gene_infos <- c(gene_infos, Filter(Negate(is.null), batch_gene_infos))
  }

  return(gene_infos)
}

get_gene_neighbors <- function(gene_id) {
  retry_attempts <- 3
  attempt <- 1
  neighbors <- NULL

  while (attempt <= retry_attempts) {
    tryCatch(
      {
        neighbors <- entrez_link(dbfrom = "gene", db = "gene", id = gene_id, linkname = "gene_gene_neighbors")
        if (!is.null(neighbors$links)) {
          break
        }
      },
      error = function(e) {
      }
    )

    attempt <- attempt + 1
  }

  if (is.null(neighbors) || length(neighbors$links) == 0) {
    return(list())
  }

  # Extract neighbor gene IDs
  neighbor_ids <- neighbors$links$gene_gene_neighbors
  neighbor_ids <- Reduce(c, neighbor_ids)
  neighbor_infos <- get_gene_list_info(neighbor_ids)

  return(neighbor_infos) # nolint: return_linter.
}
