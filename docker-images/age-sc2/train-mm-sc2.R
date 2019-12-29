#!/usr/bin/env Rscript

age.cutoff <- 50
model.state <- list("age" = age.cutoff)

model.state.metadata.file <- "./model-state-metadata.Rd"
save(model.state, file=model.state.metadata.file)

cat(paste0("Successfully saved model metadata to: ", model.state.metadata.file, "\n"))

