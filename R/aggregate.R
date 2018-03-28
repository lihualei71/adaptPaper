## Aggregate estrogen_random
NumRej <- matrix(0, 13, 30)
for (seed in 0:49){
    filename <- paste0("../data/estrogen_random_", seed, ".RData")
    load(filename)
    NumRej <- NumRej + NumRej_list[[1]] + NumRej_list[[2]]
}
NumRej <- NumRej / 100
filename <- "../data/estrogen_random.RData"
save(file = filename, NumRej)

## Aggregate results of Simulation 1
simul1.sum <- list()
for (k in 1:3){
    simul1.sum[[k]] <- list(FDP = list(), power = list())
}
for (k in 1:3){
    for (j in 1:7){
        simul1.sum[[k]]$FDP[[j]] <- rep(0, 30)
        simul1.sum[[k]]$power[[j]] <- rep(0, 30)
    }
}

for (seed in 0:49){
    filename <- paste0("../data/simul1_seed_", seed, ".RData")
    load(filename)
    for (k in 1:3){
        for (j in 1:7){
            simul1.sum[[k]]$FDP[[j]] <-
                simul1.sum[[k]]$FDP[[j]] +
                result[[k]]$FDP[[j]]
            simul1.sum[[k]]$power[[j]] <-
                simul1.sum[[k]]$power[[j]] +
                result[[k]]$power[[j]]
        }
    }
}
for (k in 1:3){
    for (j in 1:7){
        simul1.sum[[k]]$FDP[[j]] <-
            simul1.sum[[k]]$FDP[[j]] / 50
        simul1.sum[[k]]$power[[j]] <-
            simul1.sum[[k]]$power[[j]] / 50
    }
}

result <- list()
for (k in 1:3){
    tmp <- list()
    tmp$FDP <- Reduce(rbind, simul1.sum[[k]]$FDP)
    rownames(tmp$FDP) <- NULL
    tmp$power <- Reduce(rbind, simul1.sum[[k]]$power)
    rownames(tmp$power) <- NULL
    result[[k]] <- tmp
}

filename <- "../data/simul1.RData"
save(file = filename, result)

## Aggregate results of Simulation 2
simul2.sum <- list(FDP = list(), power = list())
for (j in 1:5){
    simul2.sum$FDP[[j]] <- rep(0, 30)
    simul2.sum$power[[j]] <- rep(0, 30)
}

for (seed in 0:49){
    filename <- paste0("../data/simul2_seed_", seed, ".RData")
    load(filename)
    for (j in 1:5){
        simul2.sum$FDP[[j]] <-
            simul2.sum$FDP[[j]] +
            result$FDP[[j]]
        simul2.sum$power[[j]] <-
            simul2.sum$power[[j]] +
            result$power[[j]]
    }
}

for (j in 1:5){
    simul2.sum$FDP[[j]] <-
        simul2.sum$FDP[[j]] / 50
    simul2.sum$power[[j]] <-
        simul2.sum$power[[j]] / 50
}

result <- list()
result$FDP <- Reduce(rbind, simul2.sum$FDP)
rownames(result$FDP) <- NULL
result$power <- Reduce(rbind, simul2.sum$power)
rownames(result$power) <- NULL

filename <- "../data/simul2.RData"
save(file = filename, result)
