# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract arguments
variant_file <- args[1]
region_file <- args[2]
chrom_name <- args[3]
start_name <- args[4]
end_name <- args[5]
id_name <- args[6]
window_size <- as.numeric(args[7])
output_name <- args[8]

# Print arguments
cat("variant_file:", variant_file, "\n")
cat("region_file:", region_file, "\n")

library(data.table)
library(tidyverse)
library(stats)

# region file must have chrom column with 'chr' style seqnames, and region_id, region_start and region_end columns
regions_dt = fread(region_file)
setnames(regions_dt, chrom_name, 'chrom')
setnames(regions_dt, start_name, 'region_start')
setnames(regions_dt, end_name, 'region_end')
setnames(regions_dt, id_name, 'region_id')

# check chrom column
if (!grepl('chr', regions_dt$chrom[1])) {
  regions_dt[, chrom := paste0('chr', chrom)]
}

# reformat varint file from VCF format to be readable, and filter to SNVs
variants_dt = fread(variant_file)
variants_dt <- variants_dt %>% 
  mutate(INFO = gsub(';$', '', gsub('(?:.*?=(.*?))(?:;|$)', '\\1;', INFO))) %>%
  separate(INFO, into = c('AC', 'AN', 'NS_GT', 'NS_NOGT', 'NS_NODATA', 'FILE'), sep = ';') %>%
  mutate(across(c(AC, AN, NS_GT, NS_NOGT, NS_NODATA), as.numeric))
setDT(variants_dt)
variants_dt = variants_dt[nchar(REF) == 1 & nchar(ALT) == 1]

# allow some columns to allow foverlaps with region_file
variants_dt[, region_start := as.numeric(sapply(strsplit(ID, '-'), '[[', 2))]
variants_dt[, region_end := region_start]
setnames(variants_dt, '#CHROM', 'chrom')

# run foverlaps of SNVs with the gene/regions
setkey(regions_dt, chrom, region_start, region_end)
variants_regions=foverlaps(variants_dt, regions_dt, type = 'within')
variants_regions = variants_regions[!is.na(region_start)]

#1. calculate number of SNVs per position (per gene- all the rest of steps as well)
variants_snv_pos = variants_regions[, .(snv_count = .N), by = .(region_id, POS)]
variants_snv_pos[, total_snvs := sum(snv_count), by = region_id]

#2. create data frame of sliding windows spanning each gene
region_info = unique(variants_regions[, .(chrom, region_start, region_end, region_id)])
region_info[, region_length := region_end - region_start + 1]

windows = region_info[rep(seq_len(nrow(region_info)), each = max(region_info$region_length)), ]
windows[, relative_pos := rowid(region_id)-1]
windows[, start_pos := region_start + relative_pos]
windows[, end_pos := start_pos + window_size - 1]
windows = windows[end_pos <= region_end]
windows[, window_centre := start_pos + floor(window_size/2)]


#3. join with SNVs to get SNVs per window
windows <- variants_snv_pos[windows, on = .(POS >= start_pos, POS <= end_pos, region_id), nomatch = 0]
setnames(windows, old = c("POS", "POS.1"), new = c("start_pos", "end_pos"))


#4. Calculate proportion of observed SNVs, and median per window
windows <- windows[, .(snv_count = sum(snv_count)), 
                   by = .(region_id, region_length, total_snvs, relative_pos, window_centre)]
windows[, prop_snv := snv_count/(window_size*3)]
windows[, median_prop_snv := median(prop_snv, na.rm = TRUE), by = region_id]

windows[, dist_to_med := prop_snv-median_prop_snv]
windows[, min_dist_to_med := min(dist_to_med), by = region_id]

setorder(windows, min_dist_to_med)
fwrite(windows, output_name, sep = '\t')









