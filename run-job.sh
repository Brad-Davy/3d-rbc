#$ -V -cwd
#$ -l h_rt=48:00:00
#$ -m be
#$ -l np=64
#$ -l h_vmem=3g
#$ -j y
#$ -N DNS

module list
module load anaconda
source activate base
conda activate dedalus

cd /nobackup/scbd/PhD/Year1/Dedalus/3D/RotatingRBC/Perturbation/DNS
#export OMP_NUM_THREADS=1
mpiexec python3 3d-rrbc.py --ra=1e4 --ek=1 --N=64 --max_dt=5e-5 --init_dt=1e-8 --mesh=8,8
