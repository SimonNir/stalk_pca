#!/bin/bash


#SBATCH -N 1

#SBATCH -n 32

##SBATCH --mem-per-cpu=1G     # for most cases, default memory expectations should be fine

#SBATCH -t 4:00:00
#SBATCH -o hessian.out
##SBATCH --account=brubenst-condo
#SBATCH --job-name=stalk_pca
#SBATCH --mail-user=simon_nirenberg@brown.edu
#SBATCH --mail-type=END

#SBATCH --partition=batch

module load openblas
module load hpcx-mpi

# Unload old hdf5 only if it is loaded
if module list 2>&1 | grep -q 'hdf5/1.12.2-s6aacp3'; then
    module unload hdf5/1.12.2-s6aacp3
fi

module load hdf5/1.14.1-2-rdd6y6v

module load boost-mpi/1.82.0-atgr6tz
module load cmake/3.26.3-xi6h36u
module load libarchive/3.6.2-mnc5shn
module load intel-oneapi-mkl/2023.1.0-xbcd2g

export MPI_BIND_PROCESSOR=0
export OMPI_MCA_rmaps_base_oversubscribe=1

export PATH=/users/snirenbe/gamess:$PATH
export PATH=/users/snirenbe/qmcpack_v4.1.0/qmcpack/build/bin:$PATH
export PATH=/users/snirenbe/qmcpack_v4.1.0/qmcpack/nexus/bin:$PATH

module load python/3.11.0s-ixrhc3q

export PYTHONPATH=/users/snirenbe/stalk_venv/lib64/python3.11/site-packages:$PYTHONPATH

export PYTHONPATH=/users/snirenbe/qmcpack_v4.1.0/qmcpack/nexus/lib:$PYTHONPATH
export PYTHONPATH=/users/snirenbe/qmcpack_v4.1.0/qmcpack/src/QMCTools:$PYTHONPATH
export PYTHONPATH=/users/snirenbe/stalk_pca:$PYTHONPATH
export PYTHONPATH=/users/snirenbe/surrogate_hessian_relax:$PYTHONPATH
export PYTHONPATH=/users/snirenbe/surrogate_hessian_relax/legacy:$PYTHONPATH

source /users/snirenbe/stalk_venv/bin/activate
export PATH=/users/snirenbe/stalk_venv/bin:$PATH

export LD_LIBRARY_PATH=/users/snirenbe/OpenBLAS:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/hdf5-1.12.2-olqx5gzjtts4hcun542ahnq5zyczugzm/lib:$LD_LIBRARY_PATH

python -u run2_tpw_hessian.py > output_hessian.txt

