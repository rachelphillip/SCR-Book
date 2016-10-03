## Code to fit a simple example SCR model, taken from the secr.fit() help file.
library(secr)

## Always load in data with a relative (rather than full) file specification.
## That is, use "../data/data-file" instead of "/home/ben/GitHub/analysis/data/data-file"
load("../data/example-data.RData")

## Fitting model.
example.fit <- secr.fit (detections, buffer = 100, trace = FALSE)

## Saving results. Again, use a relative file specification.
save(example.fit, file = "../objects/example.RData")
