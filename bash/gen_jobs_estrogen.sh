#!/bin/bash

params_filename='corr_params_estrogen.txt'

while read times seed
do
    filename='../results/estrogen_"'$seed'".sh'
    echo $filename
    Rout_filename='../results/estrogen_"'$seed'".Rout'
    echo $Rout_filename
    touch $filename
    echo '#!/bin/bash' > $filename
    echo '#SBATCH --cpus-per-task 1' >> $filename
    echo 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK' >> $filename
    echo '' >> $filename
    echo 'export seed="'$seed'"' >> $filename
    echo 'export times="'$times'"' >> $filename    
    echo "R CMD BATCH --no-save estrogen_random_cluster.R "$Rout_filename >> $filename
    chmod 755 $filename
    sbatch $filename
done < $params_filename
