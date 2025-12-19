library(data.table)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-i", "--ibd"), action="store", default=NA, type='character',
              help="IBD file output by plink --genome function."),
  make_option(c("-m", "--metadata"), action="store", default=NA, type='character',
              help="Metadata from cohort."),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="Output file to save data.")
)
opt = parse_args(OptionParser(option_list=option_list))
### DEBUG
if(FALSE){
  opt = list(
    ibd = "2.Internal_files/RNA_genotype/IBD_DICE_cohort.genome",
    metadata = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/RNA_seq/mapping/BN/metadata_all.csv",
    outfile = "4.Output/QC_genotype/Genotype_predicted_DICE.csv"
  )
}

ifile <- data.table::fread(opt$ibd)
metadata <- data.table::fread(opt$metadata)

ifile <- ifile[IID1 %in% metadata$sample_ID & !(IID2 %in% metadata$sample_ID)]
write.csv(ifile, gsub(".csv", "_complete.csv", opt$outfile), row.names = FALSE)

donors_problm <- ifile[Z2 == 1, .N, by = IID1][N ==2]$IID1
problemdf <- ifile[Z2 == 1 & IID1 %in% donors_problm , .(donors_Z2_1 = paste(IID2, collapse=" ")), by = IID1]

ifile <- ifile[, .SD[which.max(DST)], by = IID1]
ifile[, c("FID1", "FID2", "RT", "EZ") := NULL]

if(nrow(problemdf) > 0){
   ifile <- merge(ifile, problemdf, by = "IID1", all = T)
}

setnames(ifile, "IID2", "Predicted_Donor")

ifile <- merge(metadata, ifile, by.x = "sample_ID", by.y = "IID1", all.x = T)

write.csv(ifile, opt$outfile, row.names = FALSE)
