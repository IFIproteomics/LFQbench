species <- vector(mode="list", length=3)
names(species) <- c("HUMAN", "YEAST", "ECOLI")
species[[1]] <- "_HUMAN"
species[[2]] <- "_YEAS"
species[[3]] <- "_ECOLI"

experiments <- vector(mode="list", length=4)
names(experiments) <- c("5600-32w", "5600-64w", "6600-32w", "6600-64w")
# for each experiment:
#   three first ones are sample A, other three ones are sample B
experiments[[1]] <- c("lgillet_L150206_001", "lgillet_L150206_003", "lgillet_L150206_005",   # A
                      "lgillet_L150206_002", "lgillet_L150206_013", "lgillet_L150206_014")   # B 

experiments[[2]] <- c("lgillet_L150206_007", "lgillet_L150206_009", "lgillet_L150206_011",   # A
                      "lgillet_L150206_008", "lgillet_L150206_010", "lgillet_L150206_012")   # B

experiments[[3]] <- c("lgillet_I150211_002", "lgillet_I150211_004", "lgillet_I150211_006",   # A
                      "lgillet_I150211_003", "lgillet_I150211_005", "lgillet_I150211_007")   # B

experiments[[4]] <- c("lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012",   # A
                      "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013")   # B



