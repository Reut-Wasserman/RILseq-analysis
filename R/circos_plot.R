install.packages("RCircos")
library(RCircos)

path <- "your_path/"
experiments <- c("Hfq_lambda_30", "Hfq_lambda_60")

genome1 <- "E. coli"
genome2 <- "lambda"
genome1_length <- 4641652
genome2_length <- 48502
multiply_factor <- round(genome1_length/genome2_length)
#multiply_factor <- 1

genes <- "" 
plot_type <- "half_circles"
#plot_type <- "real_proportions"

base_per_unit <- 45.137051
padding <- 700

for (experiment in experiments){
  my_cyto.info <- data.frame(Chromosome = c(genome1, genome2), chromStart = c(0,0), chromEnd=c(genome1_length,genome2_length*multiply_factor), Band=c(genome1, genome2), Stain=c("gpos25", "acen"))
  RCircos.Set.Core.Components(my_cyto.info, chr.exclude=NULL, tracks.inside=3, tracks.outside=0)
  params <- RCircos.Get.Plot.Parameters()
  params$base.per.unit <- base_per_unit
  params$chrom.paddings <- padding
  params$text.size <- 0.4
  RCircos.Reset.Plot.Parameters(params)
  
  file <- paste(path, experiment, "_", plot_type, ".pdf", sep="")
  pdf(file = file, height = 8, width = 8)
  RCircos.Set.Plot.Area()
  RCircos.Draw.Chromosome.Ideogram()
  RCircos.Label.Chromosome.Names()
  
  my_gene_label_data <- read.csv(paste(path, "two_genomes_scale_numbers_", plot_type, ".csv", sep=""))
  scale_data <- read.csv(paste(path, "two_genomes_scale_", plot_type, ".csv", sep=""))
  RCircos.Gene.Connector.Plot(scale_data, track.num = 1, side="in")
  RCircos.Gene.Name.Plot(my_gene_label_data, name.col = 4, track.num = 2, side="in")
  
  my_link_data = read.csv(paste(path, experiment, "_two_genomes_chimeras_", plot_type, genes, ".csv", sep=""))
  RCircos.Link.Plot(my_link_data, track.num = 3, lineWidth=rep(0.00001, nrow(my_link_data)))
  dev.off()
}


