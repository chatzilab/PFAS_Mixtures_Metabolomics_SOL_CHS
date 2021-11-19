# Log in to HPC 
# Performed analysis on the epyc-64	partition on Discovery
#ssh jagoodri@endeavour.usc.edu
ssh jagoodri@discovery.usc.edu

# Change Working Directory
cd /project/dconti_624/Users/jagoodri/sol


# Run SOLAR analysis
dos2unix 1_1_solar_mixtures_mwas_hpc.sh
sbatch 1_1_solar_mixtures_mwas_hpc.sh


# Change Working Directory
cd /project/dconti_624/Users/jagoodri/chs

# Run CHS analysis 
dos2unix 1_2_chs_mixtures_mwas_hpc.sh
sbatch 1_2_chs_mixtures_mwas_hpc.sh

exit

#############################################################
# Commands for troubleshooting and testing ------------------
ssh jagoodri@endeavour.usc.edu
cd /project/dconti_624/Users/jagoodri/sol
dos2unix test_pbdmpi.sh
sbatch test_pbdmpi.sh


# To rerun code test

# Different Partitions:
# To run on epic:
#SBATCH --partition=epyc-64
# To run on endeavour:
#SBATCH --partition=conti


# TO SOLVE ERROR WHERE JAGS IS NOT RECOGNIZED ON WORKER NODES:
# (this only needs to be run once)
module load gcc/11.2.0 patchelf
patchelf --set-rpath /spack/apps2/linux-centos7-x86_64/gcc-11.2.0/jags-4.3.0-d373ytkqeamusbww7n2qjxyfwgikw2w5/lib /home1/jagoodri/R/x86_64-pc-linux-gnu-library/4.1/rjags/libs/rjags.so



# To install R2jags
module spider jags

module load usc r
module load jags
module load gcc/11.2.0
module load openblas/0.3.18


export PKG_CONFIG_PATH=/packages/jags/4.3.0/lib/pkgconfig
pkg-config ––modversion jags
R
install.packages("rjags", configure.args="––enable-rpath")
install.packages("R2jags", configure.args="––enable-rpath")