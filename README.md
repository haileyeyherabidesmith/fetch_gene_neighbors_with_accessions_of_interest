# R tool to fetch gene neighbors with accessions of interest

## interproscan testing

### resources:

- <https://interproscan-docs.readthedocs.io/en/v5/HowToUseViaContainer.html>

### notes:

- install docker desktop
- to execute interproscan:
  - docker run --rm -v ~/interproscan:/data interpro/interproscan
- download interproscan data
  - curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/alt/interproscan-data-5.73-104.0.tar.gz

### use InterProScan to get all accessions from that file

#### for small data:

docker run --rm \
-v $PWD/interproscan-5.73-104.0/data:/opt/interproscan/data \
-v $PWD/data_small:/input \
-v $PWD/temp:/temp \
-v $PWD/data_small:/output \
interpro/interproscan:5.73-104.0 \
--input /input/neighbor_protein_fasta_sequences.fa \
--output-dir /output \
--tempdir /temp \
--cpu 6 \
-exclappl Coils,MobiDBLite

#### for large data:

docker run --rm \
-v $PWD/interproscan-5.73-104.0/data:/opt/interproscan/data \
-v $PWD/data:/input \
-v $PWD/temp:/temp \
-v $PWD/data:/output \
interpro/interproscan:5.73-104.0 \
--input /input/neighbor_protein_fasta_sequences.fa \
--output-dir /output \
--tempdir /temp \
--cpu 6 \
-exclappl Coils,MobiDBLite
