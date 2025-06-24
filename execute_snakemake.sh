#!/bin/bash

snakemake --executor cluster-generic \
  --jobs 100 \
  --default-resources qos=normal \
  --default-resources mem_mb=5000 \
  --default-resources time="00:10:00" \
  --default-resources requeue=0 \
  --default-resources trigger=0 \
  --latency-wait 60 \
  --keep-going \
  --printshellcmds \
  --scheduler greedy \
  --use-conda \
  --local-cores 1 \
  --max-jobs-per-second 10 \
  --max-status-checks-per-second 1 \
  --cluster-generic-submit-cmd 'mkdir -p logs/slurm/{rule} && sbatch \
    --account=applicata \
    --nodes=1 \
    --qos={resources.qos} \
    --cpus-per-task={threads} \
    --mem={resources.mem_mb} \
    --job-name=smk-{rule}-{jobid} \
    --output=logs/slurm/{rule}/{rule}-{jobid}-%j.out \
    --error=logs/slurm/{rule}/{rule}-{jobid}-%j.err \
    $( [ "{resources.requeue}" -eq "1" ] && echo "--requeue" ) \
    $( [ "{resources.trigger}" -eq "1" ] && echo "--reservation=prj-trigger --nodelist=mtx30" ) \
    --time={resources.time} \
'
  "$@"
