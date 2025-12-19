library(ggplot2)
library(data.table)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-t", "--tpm"), action="store", default=NA, type='character',
              help="TPM counts file."),
  make_option(c("-m", "--metadata"), action="store", default=NA, type='character',
              help="Metadata with Gender in the columns."),
  make_option(c("-g", "--geneannot"), action="store", default=NA, type='character',
              help="Gene annotation."),
  make_option(c("-x", "--genes"), action="store", default="UTY,USP9Y,DDX3Y,XIST", type='character',
              help = "Genes to analyze. SPecifically sex-associated genes separated by comma. [default %default]"),
  make_option(c("-o", "--outpath"), action="store", default=NA, type='character',
              help="Output path to save data.")
)
opt = parse_args(OptionParser(option_list=option_list))

### DEBUG
if (FALSE){
  opt = list(
        tpm = "4.Output/counts/TPM_counts.csv",
        metadata = "metadata_all.csv",
        geneannot = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/reference/GRCh38-2020-A_build/GRCh38_annotation_filtered.csv",
        outpath = "4.Output/QC_sex/"
  )
}

genes <- strsplit(opt$genes, ",")[[1]]
## Loading TPM data and mean tab
tpm <- read.csv(opt$tpm, check.names = FALSE, row.names = 1)
annot <- fread(opt$metadata)
geneannotation <- fread(opt$geneannot)

tpm <- tpm[, annot$sample_ID]

sex.box <- function( gene.name, out.path = opt$outpath) {
    #gene.name <- "FOS"
    ens.id <- as.character( geneannotation[ gene_name == gene.name, "Geneid"])

    cat("Generating plot for gene: ", gene.name, "\n")
    plot.annot <- annot
    plot.annot$expr <- as.numeric( tpm[ ens.id,])
    d <- ggplot( plot.annot,aes(x=Gender,y=expr,color=Gender)) +
            geom_boxplot(outlier.colour = NA )+
            geom_point(position = position_jitter(width=0.2), size=1, color="black")+
            scale_color_manual(values=c("red", "blue") )+ #facet_wrap(~celltype, nrow=1)+
            theme_classic()+
            theme(axis.text.x=element_text(size=15),
                  axis.line=element_line(size=0.85, colour="black"),
                  axis.ticks=element_line(colour="black", size=0.85),
                  axis.ticks.length=unit(.25, "cm"),
                  axis.title.x=element_text(size=18),
                  legend.position="none",
                  plot.title=element_text(hjust=0.5, size=25),
                  axis.text.y=element_text(size=15),
                  axis.title.y=element_text(size=18) )+
        labs(title=gene.name, y="TPM", x="Sex")

    if (sum("celltype" == names(plot.annot)) > 0) {
      d <- d + facet_wrap(~ celltype, nrow = 1)
    }
    if (sum("Source" == names(plot.annot)) > 0) {
      d <- d + facet_wrap(~ Source, nrow = 1)
    }
    pdf( paste0(out.path, gene.name, "_Plot.pdf"), width=8, height= 8)
    print(d)
    dev.off()

    names(plot.annot)[ names(plot.annot) == "expr"] <- paste0( "tpm.", gene.name)
    write.csv( plot.annot, paste0( out.path, gene.name, "_Dat.csv"))
}


sapply( genes, sex.box)

for( gene.look in genes ){
    cat("analyzing:", gene.look, "\n")
    temp <- read.csv( paste0( opt$outpath, gene.look, "_Dat.csv"), row.names=1)
    length( which( rownames( temp) %in% rownames(annot)))
    temp <- temp[ rownames(annot),]
    annot[, paste0("tpm.", gene.look)] <- temp[,paste0("tpm.", gene.look)]
}

write.csv(annot, paste0( opt$outpath, "Sex-markers.csv"))
