## DADA2 16S rRNA-gene meta-taxonomy
# F1000 tutorial (minor adaptation by Roberto Siani), 060721
# Original at https://doi.org/10.12688/f1000research.8986.2

## setup env

p_load(dada2)

## trim adaptors with cutadapt

## parallel -j 4 --link cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed/{1} -p trimmed/{2} {1} {2} ::: *R1_001.fastq.gz ::: *R2_001.fastq.gz


## set path of fastq files

raw_path =
  "~/02_david/01_inData/david_16S/trimmed"

fns =
  sort(
    list.files(raw_path,
               full.names = T)
  )

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs = fns[grepl("R1", fns)]

fnRs = fns[grepl("R2", fns)]

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names =
  sapply(
    strsplit(
      basename(fnFs), "_"), `[`, 1)

# inspect quality profiles of the reads

quality_profiles =
  map(
  c(fnFs, fnRs),
  ~plotQualityProfile(.x)
)

# trim and filter low quality reads

dir.create("~/02_david/01_inData/david_16S/trimmed/filtered")

filtFs =
  file.path(
    raw_path,
    "filtered",
    paste0(
      sample.names,
      "_F_filt.fastq.gz"))

filtRs =
  file.path(
    raw_path,
    "filtered",
    paste0(
      sample.names,
      "_R_filt.fastq.gz"))

for (i in seq_along(fnFs)) {
  fastqPairedFilter(
    c(fnFs[[i]], fnRs[[i]]),
    c(filtFs[[i]], filtRs[[i]]),
    trimLeft = 20,
    trimRight = 20,
    truncLen = c(280,
                 220),
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    rm.phix = T,
    compress = T,
    multithread = T
  )
}

# dereplicate reads to reduce redundancy

derepFs = derepFastq(filtFs)

derepRs = derepFastq(filtRs)

# learn error rates

errF = learnErrors(derepFs)

errR = learnErrors(derepRs)

# check estimated error rates

plotErrors(errF, nominalQ = T)

plotErrors(errR, nominalQ = T)

# sample inference

dadaFs =
  dada(derepFs,
       err = errF,
       multithread = T,
       pool = T)

dadaRs =
  dada(derepRs,
       err = errR,
       multithread = T,
       pool = T)

# merge reads

mergedPs =
  mergePairs(
    dadaFs,
    derepFs,
    dadaRs,
    derepRs
  )

# construct sequence table and remove chimeras

seqTab =
  mergedPs %>%
  makeSequenceTable() %>%
  removeBimeraDenovo()

# assign taxonomy (download reference database here)

ref_fasta =
  "~/00_lib/silva_nr99_v138.1_wSpecies_train_set.fa.gz"

ref_species =
  "~/00_lib/silva_species_assignment_v138.1.fa.gz"

set.seed(3)
taxTab =
  assignTaxonomy(
    seqTab,
    ref_fasta,
    multithread = T,
    minBoot = 80)

speciesTab =
  assignSpecies(
    seqTab,
    ref_species,
    tryRC = T)

taxTab_plus =
  taxTab %>%
  as.data.frame() %>%
  mutate(Genus = coalesce(speciesTab[, "Genus"]),
         Species = coalesce(speciesTab[, "Species"]))

# get sequences

seqs =
  getSequences(seqTab)

names(seqs) = seqs

# save results for analyses

dir.create("~/02_david/01_inData/")

write_tsv(taxTab_plus %>%
            as.data.frame() %>%
            rownames_to_column("ASV"),
          "~/02_david/01_inData/taxonomy.tsv")

write_tsv(seqTab %>%
          as.data.frame() %>%
            rownames_to_column("Sample"),
          "~/02_david/01_inData/abundances.tsv")

Biostrings::writeXStringSet(
  seqs %>%
    Biostrings::DNAStringSet(),
  "~/02_david/01_inData/refSeq.fa")
