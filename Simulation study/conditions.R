# conditions #
# this code makes the conditions for the simulation study #

conditions <- list(setup = c("both", "weights", "loadings"),
                   dimensions = c("low", "high"),
                   ones = c(0),
                   vafx = c(1, 0.9, 0.50), 
                   sparsity = c(0.9, 0.7, 0.5),
                   overlap = c(0, 0.4),
                   reps = c(1:50))

condition_df <- data.frame(setup = NA, 
                           dimensions = NA,
                           ones = NA,
                           vafx = NA,
                           sparsity = NA,
                           overlap = NA,
                           reps = NA)

counts <- 0

for (setupz in 1:3){
  for (dimensionsz in 1:2){
    for (onesz in 1:1){
      for (vafxz in 1:3){
        for (sparsityz in 1:3){
          for (overlapz in 1:2){
            for (repsz in 1:50){
              
              counts <- counts + 1
              
              setup_now <- conditions$setup[setupz]
              dimension_now <- conditions$dimensions[dimensionsz]
              ones_now <- conditions$ones[onesz]
              sparsity_now <- conditions$sparsity[sparsityz]
              vafx_now <- conditions$vafx[vafxz]
              overlap_now <- conditions$overlap[overlapz]
              reps_now <- conditions$reps[repsz]
              
              condition_df[counts,] <- c(setup_now, 
                                         dimension_now,
                                         ones_now,
                                         vafx_now,
                                         sparsity_now,
                                         overlap_now,
                                         reps_now)
              
              print(counts)
            }
          }
        }
      }
    }
  }
}


# save("condition_df", file = "C://Users//park//Desktop//project_spca_models//industry//sim_29_july_2020/condition_df.Rdata")

# getwd()

