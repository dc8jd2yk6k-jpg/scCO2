#!/bin/bash
# Minimal sequential job queue: each line of runs/queue.txt is
#   <script.py> [args...]
# and is executed with 4 MPI ranks once the previous job has finished.
# Append lines to queue.txt at any time; touch runs/queue.stop to exit when idle.
cd /home/user/scCO2/zrncl_dft
touch runs/queue.txt
i=$(cat runs/queue.pos 2>/dev/null || echo 0)
while true; do
  line=$(sed -n "$((i+1))p" runs/queue.txt)
  if [ -z "$line" ]; then
    [ -f runs/queue.stop ] && break
    sleep 5; continue
  fi
  if [[ "$line" == \#* ]]; then i=$((i+1)); echo $i > runs/queue.pos; continue; fi
  name=$(basename ${line%% *} .py)
  echo "$(date '+%F %T') START [$i] $line" >> runs/queue.log
  mpiexec --allow-run-as-root -n 4 /opt/dftenv/bin/python $line > runs/logs/${i}_${name}.out 2>&1
  rc=$?
  echo "$(date '+%F %T') END   [$i] rc=$rc $line" >> runs/queue.log
  i=$((i+1)); echo $i > runs/queue.pos
done
