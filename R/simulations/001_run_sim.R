## Use this script to run a set of simulatios in a saver mode when using a server
## via ssh. The script will run the simulations in parallel and save the results in a file.

## Source the setup scripts

source("R/000_setup.R")
source("R/simulations/000_setup_sim.R")

## Create a vector of files to be run. For exsample to run the
## simulations in Scenario 1B for L1 statistic and montecarlo method.

files_to_run <- character(0)

# files_to_run <- list.files(
#   "R/simulations/Scenario1B",
#   pattern = "^l1std*.montecarlo.*\\.R$",
#   full.names = TRUE
# )

# files_to_run <- list.files(
#   "R/simulations/Scenario1A",
#   pattern = "l2std_permutations_200_100_1A.R",
#   full.names = TRUE
# )

## Uncomment and adapt the pattern to run the simulations 
## to your choice.

## Run the scripts

for (f in files_to_run) {
  source(f)
}

## Close the cluster and stop the script when all simulations are done

stopCluster(cl)

## Go to the repository folder in a terminal window and run the script with the
## following command:

# nohup Rscript R/simulations/001_run_sim.R > simulation.log 2>&1 &
# echo $!

## Note the assigned PID,e.g., 12345, which can be monitored with the command:

# ps -p 12345

## Check the log file for the progress of the simulations with the command:

# tail -f simulation.log