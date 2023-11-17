snakemake -s batch_s.smk --unlock
#snakemake --jobs 8  --printshellcmds \
snakemake --jobs 40 \
 --cluster-config slurm.json \
 --cluster "sbatch  \
-A "zps5164_sc" --qos "burst4x" -p "burst" \
-t {cluster.time} --output {cluster.output} --error {cluster.error} \
 --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --mem {cluster.mem}" \
-s batch_s.smk --latency-wait 300

#--mail-type=all --mail-user=aur1111@psu.edu \

