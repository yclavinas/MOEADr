rm(list = ls(all = TRUE))

suppressPackageStartupMessages(library(smoof))
suppressPackageStartupMessages(library(MOEADps))
suppressPackageStartupMessages(library(CAISEr))


source('~/MOEADr/R/linPF.R')
source("~/MOEADr/R/utils.R")
source('~/MOEADr/R/moead.R')
source('~/MOEADr/R/linPF_ParetoFronts.R')
source('~/MOEADr/R/perform_variation.R')
source('~/MOEADr/R/resource_allocation_select_random.R')
source('~/MOEADr/R/resource_allocation_update.R')

source('~/MOEADr/R/variation_diffmut.R')
### Build function names (instances: UF1 - UF7, dimensions 10 - 40)
fname   <- paste0("UF_", 1:7)
dims    <- c(100)
allfuns <- expand.grid(fname, dims, stringsAsFactors = FALSE)

# Assemble instances list
instances <- vector(nrow(allfuns), mode = "list")
for (i in 1:length(instances)){
  instances[[i]]$FUN <- paste0(allfuns[i,1], "_", allfuns[i,2])
}

### Build the functions listed above (so that they can be properly used)
for (i in 1:nrow(allfuns)){
  assign(x = instances[[i]]$FUN,
         value = MOEADr::make_vectorized_smoof(prob.name  = "UF",
                                               dimensions = allfuns[i, 2],
                                               id = as.numeric(strsplit(allfuns[i, 1], "_")[[1]][2])))
}

# Prepare algorithm function to be used in run_experiment():
myalgo <- function(type, instance){
  # Input parameters:
  #     - type (variant to use: "original", "original2", "moead.de" or "moead.de2")
  #     - instance (instance to be solved, e.g., instance = instances[[i]])
  # All other parameters are set internally
  
  ## Extract instance information to build the MOEADr problem format
  fdef  <- unlist(strsplit(instance$FUN, split = "_"))
  uffun <- smoof::makeUFFunction(dimensions = as.numeric(fdef[3]),
                                 id         = as.numeric(fdef[2]))
  fattr    <- attr(uffun, "par.set")
  prob.dim <- fattr$pars$x$len
  
  ## Build MOEADr problem list
  problem <- list(name = instance$FUN,
                  xmin = fattr$pars$x$lower,
                  xmax = fattr$pars$x$upper,
                  m    = attr(uffun, "n.objectives"))
  
  ## Load presets for the algorithm provided in input 'type' and 
  ## modify whatever is needed for this particular experiment
  algo.preset <- MOEADr::preset_moead("moead.de")
  algo.preset$decomp$H <- 349 # <-- set population size
  algo.preset$stopcrit[[1]]$name <-
    "maxeval" # <-- type of stop criterion
  algo.preset$stopcrit[[1]]$maxeval <- 30000 # stop crit.
  algo.preset$update$UseArchive = TRUE
  poly.ind <- which(sapply(algo.preset$variation,
                           function(x) {
                             x$name == "polymut"
                           }))
  algo.preset$variation[[poly.ind]]$pm <-
    1 / prob.dim # <--- pm = 1/d
  algo.preset$scaling <- list(name = "simple")
  
  ## Run algorithm on "instance"
  if (type == "moead.ps") {
    
    resource.allocation <-
      list(
        name = "random",
        dt = 0,
        selection = "n",
        n = 35
      )
    out <-
      moeadps(
        preset = algo.preset,
        problem = problem,
        resource.allocation = resource.allocation,
        showpars = list(show.iters = "none")
      )
    
  }
  else if(type == "moead.ps2"){
    resource.allocation <-
      list(
        name = "random",
        dt = 20,
        selection = "n",
        n = 35
      )
    out <-
      moeadps(
        preset = algo.preset,
        problem = problem,
        resource.allocation = resource.allocation,
        showpars = list(show.iters = "none")
      )
  }
  else{
    out <-
      moeadps(
        preset = algo.preset,
        problem = problem,
        showpars = list(show.iters = "none")
      )
  }
  
  ## Read reference data to calculate the IGD
  Yref  <- as.matrix(read.table(paste0("./inst/extdata/pf_data/",
                                       fdef[1], fdef[2], ".dat")))
  IGD = MOEADr::calcIGD(Y = out$Y, Yref = Yref)
  
  ## Return IGD as field "value" in the output list
  return(list(value = IGD))
}

algorithms <- list(
  list(FUN   = "myalgo",
       alias = "MOEAD.PS",
       type  = "moead.ps"),
  list(FUN   = "myalgo",
       alias = "PS-dt-20",
       type  = "moead.ps2"),
  list(FUN   = "myalgo",
       alias = "MOEAD.DE",
       type  = "moead.de")
)

my.results <- run_experiment(instances  = instances,
                             algorithms = algorithms,
                             power = 0.8,      # Desired power: 80%
                             power.target = "mean", # on average,
                             d = 0.5,          # to detect differences greater
                             # than 0.5 standard deviations
                             sig.level = 0.05, # at a 95% confidence level. 
                             se.max = 0.05,    # Measurement error: 5% 
                             dif = "perc",     # on the paired percent 
                             # differences of means,
                             method = "param", # calculated using parametric
                             # formula.
                             comparisons = "all.vs.all", # Compare all algorithms 
                             # vs all others,
                             nstart = 15,      # start with 15 runs/algo/inst
                             nmax   = 400,     # and do no more than 400 runs/inst.
                             # NOTICE: Using all but 1 cores. Change if needed
                             ncpus  = 1)




plot(my.results)




algopairs <- paste(my.results$data.summary$Alg1,
                   my.results$data.summary$Alg2,
                   sep = " - ")


par(mfrow = c(1, 1))
df <- cbind(Comparison = algopairs, my.results$data.summary)

mp <- ggplot(df, aes(x = Comparison, y = Phi, fill = Comparison))
mp +
  geom_violin(alpha = 0.6,
              show.legend = FALSE,
              scale = "width") +
  geom_boxplot(
    alpha = 0,
    show.legend = FALSE,
    outlier.shape = NA,
    width = .15
  ) +
  geom_point(
    shape = 16,
    col = "black",
    fill = "black",
    alpha = 0.6,
    position = position_jitter(width = .15)
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    col = "red",
    lty = 2
  ) +
  ylab("Percent difference in IGD") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary(my.results, test = "wilcoxon")

print(ggplot(df, aes(x = Instance, y = Phi, colour = Comparison,
                     ymin = Phi - SE, ymax = Phi + SE)) + 
        geom_pointrange(show.legend = FALSE) + 
        geom_abline(slope = 0, intercept = 0, col = 1, lty = 2) + 
        facet_grid(Comparison ~ .) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
        xlab("")
)

print(aggregate(my.results$data.raw, by = list(my.results$data.raw$Algorithm,my.results$data.raw$Instance), mean))

# linPF_CEC_results.rds <- my.results
# 
# save(linPF_CEC_results.rds, "UF_results.rds")

