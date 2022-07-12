#$ -V -cwd
#$ -l h_rt=48:00:00
#$ -m be
#$ -l np=1024
#$ -l h_vmem=4.8G
#$ -j y
#$ -N DNS

module list
module load anaconda
source activate base
conda activate dedalus

cd /nobackup/scbd/PhD/Year1/Dedalus/3D/RotatingRBC/Perturbation/DNS
#export OMP_NUM_THREADS=1
mpiexec python3 3d-rrbc.py --ra=2.8e8 --ek=1e-5 --N=256 --max_dt=5e-4 --init_dt=1e-8 --mesh=32,32
