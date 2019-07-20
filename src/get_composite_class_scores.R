#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list = list( 
  make_option(c("-d", "--data_transform"), type="character",default="quantile",
              help="how to transform the raw rnaposer scores [default: %default, others: max, zero-to-one]"),
  make_option(c("-o", "--output"), type="character",default="total_classification_scores.txt",
              help="name of output file [default %default]"),
  make_option(c("-q", "--quantile"), type="double",default=0.50,
              help="if data_transform=quantile, then raw rnaposer scores will be divide by this value [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)
parser = OptionParser(usage = "%prog [options] residue-wise-rnaposer-file", option_list=option_list)

arguments = parse_args(parser, positional_arguments = TRUE)
opt = arguments$options

if(length(arguments$args) != 1) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  if (opt$verbose){
    cat("Project: RNAPosers\n")
    cat("Purpose: Get Composite Pose Classifications\n")
    cat("Author: Aaron T. Frank\n")
    cat("Author: Sahil Chhabra\n")
    cat("Author: Jingru Xie\n")
    cat(sprintf("%s\n",date()))
  }
  
  # get arguments
  input_file = arguments$args[1]
  
  # get options
  data_transform = opt$data_transform
  output_file = opt$output
  quantile = opt$quantile
  verbose = opt$verbose
  
  # read in residue-wise class scores
  normalize = function(x) {((x - min(x)) / (max(x) - min(x)))}
  class_scores = read.table(input_file, col.names = c("pdb", "resid", "frame", "unk", "neg", "pos"))
  
  # 
  if (data_transform == "max"){
    scores = plyr::ddply(.data = class_scores, .variables = c("pdb", "resid"), .fun = function(x){data.frame(frame=x$frame, pos=x$pos/max(x$pos))})
  }
  if (data_transform == "quantile"){
    scores = plyr::ddply(.data = class_scores, .variables = c("pdb", "resid"), .fun = function(x){data.frame(frame=x$frame, pos=(x$pos)/quantile(x$pos, quantile))})
  }
  if (data_transform == "zero-to-one"){
    scores = plyr::ddply(.data = class_scores, .variables = c("pdb", "resid"), .fun = function(x){data.frame(frame=x$frame, pos=normalize(x$pos))})
  }
  
  scores = plyr::ddply(.data = scores, .variables = c("pdb", "frame"), .fun = function(x){-log(prod(x$pos, na.rm = TRUE))})
  scores = scores[order(scores$V1, decreasing = FALSE), ]

  # write out composite score
  scores = scores[1:10,]
  #scores$V1 = round(scores$V1/sum(scores$V1), 4)
  write.table(scores, file = output_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # print out
  print(head(scores, 10))
  # write pymol visualization string
  cat("\n")
  cat(sprintf("test_rna_poser_proteins_new_again(pdbid = \"%s\", offset = 0, poses = [%s, %s, %s, %s, %s, %s, %s, %s, %s, %s], render = False, scale = 1)\n", scores$pdb[1], scores$frame[1], scores$frame[2], scores$frame[3], scores$frame[4], scores$frame[5], scores$frame[6], scores$frame[7], scores$frame[8], scores$frame[9], scores$frame[10]))
  cat("\n")
}


