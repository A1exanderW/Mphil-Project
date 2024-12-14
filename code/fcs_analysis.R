if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowCore")

library(flowCore)

file.name <- "C:/Users/lxwg3/OneDrive - The University of Queensland/2024/Mphil/Experiments/Fluorescence activated sorting/FACS data/QBI data/Session 1/fcs files/S4_dualstain.fcs"
x <- read.FCS(file.name, transformation=FALSE)
summary(x)

metadata <- keyword(x)

# Check for flow rate information
flow_rate <- metadata[["$TIMESTEP"]]  # This is a common keyword for flow rate,unit is in seconds
print(flow_rate)

# Extract timestamps and number of events
timestamps <- exprs(x)[, "Time"]
num_events <- nrow(exprs(x))

# Calculate flow rate (events per second)
total_time <- max(timestamps) - min(timestamps)
flow_rate <- num_events / total_time
print(flow_rate)