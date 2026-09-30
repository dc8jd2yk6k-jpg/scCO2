#!/bin/bash
cd /home/user/scCO2/zrncl_dft
until grep -q "^ecut 800" runs/convergence.log; do sleep 5; done
pid=$(pgrep -x mpiexec)
[ -n "$pid" ] && kill $pid
sleep 5
nohup runs/queue.sh > runs/queue.out 2>&1 &
echo "queue started"
