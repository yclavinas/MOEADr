rm(list = ls(all = TRUE))

library(smoof)
library(emoa)
library(mco)
library(feather)
library(parallel)
library(MOEADps)

source("~/MOEADr/R/utils.R")

repetitions <-  10
dimensions <- 40
lambda <- 50
pop.size <- 250
cores <-  16
cl <- makeCluster(cores)
maxevals <- 100000

decomp.small    <-
  list(name       = "sld", H = pop.size/10-1)


decomp.big    <-
  list(name       = "sld", H = pop.size - 1 )

variation = preset_moead("moead.de")$variation
variation[[2]]$pm = 1 / dimensions
scaling <- list()
scaling$name <- "simple"
neighbors.big <- preset_moead("moead.de")$neighbors
neighbors.big$T <- pop.size * 0.2
neighbors.small <- preset_moead("moead.de")$neighbors
neighbors.small$T <- (pop.size /10) * 0.2
update <- preset_moead("moead.de")$update
update$UseArchive = TRUE

stopcrit  <- list(list(name    = "maxeval",
                       maxeval = maxevals))



n.obj <- 2

resource.allocation.1 <-
  list(
    name = "random",
    dt = 0,
    selection = "n",
    n = 1 
  )

resource.allocation.small <-
  list(
    name = "random",
    dt = 0,
    selection = "n",
    n = (pop.size /10) - n.obj
  )



problem.to.solve <- c("DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7")
for (fun_name in problem.to.solve) {
  problem <-
    load.DTLZ.function(fun_name, dimension = dimensions, n.obj = n.obj)
  problem.smoof.DTLZ <- problem$fn
  problem.DTLZ <- function(X) {
    t(apply(X, MARGIN = 1,
             FUN = problem.smoof.DTLZ))
  }
  par.set = ParamHelpers::getParamSet(problem.smoof.DTLZ)
  fun <- list(
    name       = "problem.DTLZ",
    xmin       = as.numeric(getLower(par.set)),
    xmax       = as.numeric(getUpper(par.set)),
    m          = n.obj
  )
  
  
  dir.name <- paste0("~/tec/",fun_name, "/")
  if (!dir.exists(dir.name))
    dir.create(dir.name)
  
  
  for (j in 1:repetitions) {
    cat("rep", j, "\n")
    
    seed <- sample(1:1000)[1]
    
    
    dir.name <- paste0("~/tec/",fun_name,"/moead_small_iter_",j,"/")
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    moead.small <- moeadps(
      problem  = fun,
      preset   = preset_moead("moead.de"),
      decomp = decomp.small,
      variation = variation,
      stopcrit = stopcrit,
      scaling = scaling,
      neighbors = neighbors.small,
      showpars = list(show.iters = "none", showevery = 1000),
      seed = seed,
      update = update,
      saving.dir = dir.name
    )
    
    rm(moead.small)
    
    dir.name <- paste0("~/tec/",fun_name,"/moead_big_iter_",j,"/")
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    moead.big <- moeadps(
      problem  = fun,
      preset   = preset_moead("moead.de"),
      decomp = decomp.big,
      variation = variation,
      stopcrit = stopcrit,
      scaling = scaling,
      neighbors = neighbors.big,
      showpars = list(show.iters = "none", showevery = 1000),
      seed = seed,
      update = update,
      saving.dir = dir.name
    )
    
    rm(moead.big)
    
    dir.name <- paste0("~/tec/",fun_name,"/moead_steady_iter_",j,"/")
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    moead.ps.1 <- moeadps(
      problem  = fun,
      preset   = preset_moead("moead.de"),
      decomp = decomp.big,
      variation = variation,
      stopcrit = stopcrit,
      scaling = scaling,
      neighbors = neighbors.big,
      showpars = list(show.iters = "none", showevery = 1000),
      update = update,
      resource.allocation = resource.allocation.1,
      seed = seed,
      saving.dir = dir.name
    )
    rm(moead.ps.1)
    
    dir.name <- paste0("~/tec/",fun_name,"/moead_ps_25_iter_",j,"/")
    if (!dir.exists(dir.name))
      dir.create(dir.name)
    
    moead.ps.small <- moeadps(
      problem  = fun,
      preset   = preset_moead("moead.de"),
      decomp = decomp.big,
      variation = variation,
      stopcrit = stopcrit,
      scaling = scaling,
      neighbors = neighbors.big,
      showpars = list(show.iters = "none", showevery = 1000),
      update = update,
      resource.allocation = resource.allocation.small,
      seed = seed,
      saving.dir = dir.name
    )
    
    rm(moead.ps.small)
    
    
    
  }
  
  
}
