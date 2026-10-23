
library(GEOquery)

# Download GEO metadata
gse <- getGEO("GSE251911", GSEMatrix = TRUE, getGPL = FALSE)

# Convert metadata to data.frame
metadata_df <- as.data.frame(Biobase::pData(gse[[1]]))

# View metadata
View(metadata_df)
dim(metadata_df)
write.csv(metadata_df, file = "raw_metadata.csv", row.names = FALSE)


# Select only  experimental groups
metadata_clean <- metadata_df[
  metadata_df[["treatment:ch1"]] %in%
    c("24 Hr Untreated", "24 Hr Cytokines"),
]

# Keep relevant columns
metadata_clean <- data.frame(
  sample_id = metadata_clean$geo_accession,
  islet_id = sub(",.*", "", metadata_clean$title),
  donor_id = metadata_clean[["unos id:ch1"]],
  condition = ifelse(
    metadata_clean[["treatment:ch1"]] == "24 Hr Cytokines",
    "Cytokines",
    "Untreated"
  ),
  stringsAsFactors = FALSE
)

# Set reference group
metadata_clean$condition <- factor(
  metadata_clean$condition,
  levels = c("Untreated", "Cytokines")
)

# Sort by donor and condition
metadata_clean <- metadata_clean[
  order(metadata_clean$donor_id,
        metadata_clean$condition),
]

rownames(metadata_clean) <- metadata_clean$sample_id

# Verify matched pairs
stopifnot(
  nrow(metadata_clean) == 38,
  length(unique(metadata_clean$donor_id)) == 19,
  !anyDuplicated(metadata_clean$sample_id),
  all(table(metadata_clean$donor_id,
            metadata_clean$condition) == 1)
)

# Inspect cleaned metadata
dim(metadata_clean)
table(metadata_clean$condition)
head(metadata_clean)
View(metadata_clean)

# Save
write.csv(
  metadata_clean,
  "GSE251911_clean_metadata.csv",
  row.names = FALSE
)
