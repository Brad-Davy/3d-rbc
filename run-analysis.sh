#$ -V -cwd
#$ -l h_rt=1:00:00
#$ -m be
#$ -l np=1
#$ -l h_vmem=4.6G
#$ -j y
#$ -N analysis
#$ -P tartarus

module list
module load anaconda
source activate base
conda activate dedalus

directory='results/Ra_8-00e+08_Ek_1-00e-06_Pr_7-0_N_256_Asp_2-0/'
vtrFileName=$directory'/fields.vtr'
echo 'Running analysis on directory:'
echo $directory
echo 'The vtr file name is:'
echo $vtrFileName

cd /nobackup/scbd/PhD/Year1/Dedalus/3D/rotatingRBC/Perturbation/DNS
#export OMP_NUM_THREADS=1

mpiexec python3 analysis.py --dir=$directory --snap_t=5
mpiexec python3 forces-spectrums.py --dir=$directory --snap_t=5 --mask=4
python3 convert-to-vtr.py --dir=$directory --fields=T,u,v,w --output_file=$vtrFileName
