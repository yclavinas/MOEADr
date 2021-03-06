rm(list = ls(all = TRUE))
library(smoof)
library(emoa)
library(feather)
library(compiler)
library(ggplot2)
library(ggthemes)
library(MOEADps)
library(eaf)


number.fun <- 1
repetitions <- 9

checkpoints <- (0:40) * 2500
checkpoints[1] <- 1
collab_results_uf9 <-
  read_feather("~/france_data/uf9_collab_results")
collab_results_dtlz7 <-
  read_feather("~/france_data/dtlz7_collab_results")

collab_results <- rbind(collab_results_uf9, collab_results_dtlz7)


fun.names1 <- list()

for (i in 9:9) {
  fun.names1[[length(fun.names1) + 1]] = paste0("UF", i)
}

for (i in 7:7) {
  fun.names1[[length(fun.names1) + 1]] = paste0("DTLZ", i)
}



source("~/MOEADr/R/summary_moead.R")
source("~/MOEADr/R/utils.R")
source("~/MOEADr/R/loadPlotData.R")
#
strategy <-
  c("DS",
    "RI",
    "MOEA/D-PS")

names <- c("moead.norm",
           "moead.RI",
           "moead.random")



fun_hv <- data.frame()
fun_igd <- data.frame()
results <- data.frame()
for (fun in fun.names1) {
  ref1 <- data.frame()
  print(fun)
  
  benchmark <- strsplit(fun, "[0-9]")[[1]][1]
  number <- strsplit(fun, "[A-Z]")[[1]][3]
  if (benchmark == "DTLZ") {
    Yref <-
      as.matrix(read.table(paste0(
        "../inst/extdata/pf_data/", fun, ".2D.pf"
      )))
    colnames(Yref) <- c("f1", "f2")
    ref.point <- c(1, 1)
    number_subproblems <-
      c(3, 4, 6, 8, 10, 30, 50, 100, 150, 250)
  }
  else {
    Yref <-
      as.matrix(read.table(paste0(
        "../inst/extdata/pf_data/", fun, ".dat"
      )))
    if (as.numeric(number) == 8 ||
        as.numeric(number) == 9 || as.numeric(number) == 10) {
      colnames(Yref) <- c("f1", "f2", "f3")
      ref.point <- c(1, 1, 1)
      number_subproblems <-
        c(4, 6, 8, 10, 30, 50, 100, 150, 250)
    }
    else{
      colnames(Yref) <- c("f1", "f2")
      ref.point <- c(1, 1)
      number_subproblems <-
        c(3, 4, 6, 8, 10, 30, 50, 100, 150, 250)
    }
  }
  
  
  
  
  for (j in 1:repetitions) {
    for (lambda in number_subproblems) {
      moead.random <-
        loadPlotData(
          name = paste0(fun, "_moead.random_", lambda, "_"),
          j = j,
          wd = "~/france_data/"
        )
      
      
      moead.ds <-
        loadPlotData(
          name = paste0(fun, "_moead.norm_", lambda, "_"),
          j = j,
          wd = "~/france_data/"
        )
      
      
      moead.RI <-
        loadPlotData(
          name = paste0(fun, "_moead.RI_", lambda, "_"),
          j = j,
          wd = "~/france_data/"
        )
      ref1 <-
        rbind(ref1,
              moead.random$Y,
              moead.ds$Y,
              moead.RI$Y)
    }
  }
  
  total_hv <- data.frame()
  total_igd <- data.frame()
  
  my_hv <- data.frame()
  my_igd <- data.frame()
  stg_idx <- 1
  for (name in names) {
    print(name)
    for (my_rep in 1:repetitions) {
      for (lambda in number_subproblems) {
        moea <-
          loadPlotData(
            name = paste0(fun, "_", name, "_", lambda, "_"),
            j = my_rep,
            wd = "~/france_data/"
          )
        moea$n.iter <- as.integer(moea$n.iter)
        
        
        ck_idx <- 1
        
        for (iter in 1:moea$n.iter) {
          if (iter == 1) {
            nfe <- dim(moea$X)[1]
          }
          else{
            nfe <- iter * lambda + 500
          }
          if (nfe >= checkpoints[ck_idx] || iter == moea$n.iter) {
            if (number == 9) {
              PF <-
                data.frame(
                  cbind(
                    moea$plot.paretofront[moea$plot.paretofront$stage == iter, ]$V1,
                    moea$plot.paretofront[moea$plot.paretofront$stage == iter, ]$V2,
                    moea$plot.paretofront[moea$plot.paretofront$stage == iter, ]$V3
                  )
                )
            }
            else{
              PF <-
                data.frame(cbind(
                  moea$plot.paretofront[moea$plot.paretofront$stage == iter, ]$V1,
                  moea$plot.paretofront[moea$plot.paretofront$stage == iter, ]$V2
                ))
            }
            
            if (ck_idx == length(checkpoints)) {
              my.iter <- 100000
            }
            else{
              my.iter <- checkpoints[ck_idx]
            }
            
            
            # calc HV
            colnames(PF) <- colnames(ref1)
            
            my_igd <-
              rbind(my_igd,
                    cbind(
                      igd = igd(PF, Yref),
                      iter = my.iter,
                      fun = fun,
                      Strategy = paste0(strategy[stg_idx], "_", lambda)
                    ))
            
            PF <- scaling_Y(PF, ref1)
            hv <- dominated_hypervolume(t(PF), ref = ref.point)
            
            
            my_hv <-
              rbind(my_hv,
                    cbind(
                      hv = hv,
                      iter = my.iter,
                      fun = fun,
                      Strategy = paste0(strategy[stg_idx], "_", lambda)
                    ))
            ck_idx <- ck_idx + 1
          }
        }
      }
      nfe <- 0
      total_hv <- rbind(total_hv, my_hv)
      total_igd <- rbind(total_igd, my_igd)
    }
    
    
    fun_hv <-
      rbind(fun_hv, total_hv)
    
    fun_igd <-
      rbind(fun_igd, total_igd)
    
    stg_idx <- stg_idx + 1
    
  }
  
  fun_hv$hv <- as.numeric(as.character(fun_hv$hv))
  fun_igd$igd <- as.numeric(as.character(fun_igd$igd))
  
}

write_feather(fun_hv, "fun_hv")
write_feather(fun_igd, "fun_igd")