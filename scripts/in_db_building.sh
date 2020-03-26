#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set max wallclock time
#SBATCH --time=50:00:00

# set name of job
#SBATCH --job-name=test_db_querying

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=david.hemprich-bennett@zoo.ox.ac.uk

# run the application

module purge
module load R/3.5.2

Rscript db_building.R