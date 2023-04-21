#!/bin/sh
#SBATCH --account=nn9497k
#SBATCH --job-name=C2H5_HtoC2H6
#SBATCH --time=2-00:00:00
#SBATCH --nodes=5
#SBATCH --output=slurm-%j.out
NAME="C2H5_HtoC2H6"
# Load Vasp modules
module --force purge
module load StdEnv
module load VASPModules
module load VaspExtra
module load VASP/6.3.2-intel-2021b-std-wannier90-libxc-hdf5-beef-vtst-097e6cb5a78f237dc588ba9c7877f23b # for the VTST parts
# Loads personalized enviorment on fram
module load Anaconda3/2022.05
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/home/chrisafa/.conda/envs/master
WORKDIR=$USERWORK/temp/$SLURM_JOB_ID
mkdir -p $WORKDIR
cp -r job/* $WORKDIR
cd ${WORKDIR}
python *.py
cat *.py
STOREDIR=$USERWORK/temp/out/$NAME
mkdir -p $STOREDIR
# Copy only relevant folders into dat folder
for folder in $(ls -d $WORKDIR*);do
    f="$(basename -- $folder)"
    mkdir -p $USERWORK/temp/out/"$NAME"/"$f"
    for file in OUTCAR OSZICAR vasprun.xml; do
        cp "$folder"/"$file" $USERWORK/temp/out/"$NAME"/"$f"/"$file"
    done;
done
cp $SLURM_SUBMIT_DIR/slurm-$SLURM_ARRAY_JOB_ID.out $STOREDIR/"out.out