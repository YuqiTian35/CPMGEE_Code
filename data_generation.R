# data generation
rcont.clm <- function(clsize = clsize,  # cluster size
                      betas = betas, # nolint
                      xformula = formula(xdata),
                      xdata = parent.frame(),
                      cor.matrix = cor.matrix){
  # beta*x
  xmat <- model.matrix(xformula, data = xdata)[,-1]
  linear_predictor <- matrix(xmat %*% betas,
                             ncol = clsize, byrow=TRUE) # changed from cluster_size to clsize
  
  
  sample_size <- nrow(linear_predictor)
  
  # e
  distr <- 'qlogis' # logit link
  simulated_latent_variables <- rnorta(
    sample_size, cor.matrix, rep(distr, nrow(cor.matrix))
  )
  
  # y
  y <- c(t(simulated_latent_variables + linear_predictor)) 
  
  # data
  id <- rep(1:sample_size, each = clsize)
  time <- rep(1:clsize, sample_size)
  sim_model_frame <- model.frame(formula = xformula, 
                                 data = xdata)
  
  simdata <- data.frame(y, sim_model_frame, id, time)
  
  return(list(simdata=simdata, latent=simulated_latent_variables))
}
