#!/bin/bash
ml R
R CMD INSTALL -l ~/R/x86_64-pc-linux-gnu-library/4.0 ~/GenITR/GenITR_0.1.0.tar.gz
for m in 3 4 5
do
  for k in 350 500 800 
    do
    for j in 3200
      do
      for p in `seq 1 500`
      do
      sbatch -t 01-00:00:00 R CMD BATCH --no-save "--args n=${k}; N=${j}; case=as.character(${m});index=${p};delta=0" sim5.R res_${k}_${j}_${m}_${p}.out
      done
    done
  done  
done
squeue -A zhao_y | awk '{print $4}' | sort | uniq -c

for n in 0 1 2 3 4 5
do
for m in 3 4 5
do
  for k in 350 
    do
    for j in 3200
      do
      for p in `seq 1 500`
      do
      sbatch -t 01-00:00:00 R CMD BATCH --no-save "--args n=${k}; N=${j}; case=as.character(${m});index=${p}; delta=${n}" sim5.R res_${k}_${j}_${m}_${p}_${n}.out
      done
    done
  done  
done
done
squeue -A zhao_y | awk '{print $4}' | sort | uniq -c

for m in 4
do
  for k in 350 500 800 
    do
      for p in `seq 1 500`
      do
      sbatch -t 01-00:00:00 R CMD BATCH --no-save "--args n=${k}; case=as.character(${m});index=${p};delta=0" sim4.R res_${k}_${j}_${m}_${p}.out
      done
  done  
done
squeue -A zhao_y | awk '{print $4}' | sort | uniq -c

for p in `seq 1 500`
do
sbatch -t 01-00:00:00 R CMD BATCH --no-save "--args index=${p}" real_data.R res_${p}.out
done