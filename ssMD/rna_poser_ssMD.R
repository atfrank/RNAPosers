# user functions

# test with R scripting front-end version 3.3.1 (2016-06-21)
# conda install -c r r=3.6.0 r-pls r-optparse -m -n my-r

library("optparse")
library("pls")

option_list = list(make_option(c("-o","--output"), type = "character", default = "ssMD_ML_profile.txt",
                help = "name of output file")
)

parser = OptionParser(usage = "%prog [options] ssMD_model feature_file", option_list = option_list)
arguments = parse_args(parser, positional_arguments = TRUE)
opt = arguments$options

if(length(arguments$args) != 2) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  cat("ML2ssMD: ML predicted ssMD unbinding profile\n")
  cat("Author: Yichen Liu\n")
  cat("Author: Aaron T. Frank\n")

  # get arguments
  model = arguments$args[1]
  data =  arguments$args[2]
  
  # get options
  output = opt$output
  
  # load model
  load(model)
    
  # load data
  features = read.table(data)
  features = features[, -1] # ignore the first column
  
  # predict
  ncomp = pls::selectNcomp(models$plsa, method = "onesigma", plot = FALSE) # figure out optimal number of components
  A = as.vector(predict(models$plsa, ncomp = ncomp, newdata = features)) # predict taus

  ncomp = pls::selectNcomp(models$plsb, method = "onesigma", plot = FALSE)
  tau = as.vector(predict(models$plsb, ncomp = ncomp, newdata = features)) # predict taus

  ncomp = pls::selectNcomp(models$plsc, method = "onesigma", plot = FALSE)
  B = as.vector(predict(models$plsc, ncomp = ncomp, newdata = features)) # predict taus
  result = data.frame(pose_number = 1:nrow(features), A = exp(A), tau = exp(tau), B = exp(B)) # store results
  
  # write output
  write.table(result, file = output, row.names = F, col.names = T, quote = F)  
  print(head(result))
}
