#!/bin/bash
# All stages of compute_certificate.jl, ending with hexagon_certificate.mrdi. Setups and cover in 3 processes, greedy in 4
# (as for the certificate of 2026-11-07). Each process is restarted after a crash or a kill by its memory cap (at most 30
# times); every stage continues where it stopped, so the script can also be stopped and started again.
# From the julia folder (or any folder):   bash studies/hexagon/certificate/run_compute_overnight.sh
# Progress: tail -n 2 studies/hexagon/certificate/work/*.log

cd "$(dirname "$0")/../../.." || exit 1
WORK=studies/hexagon/certificate/work
mkdir -p $WORK

run() {   # $1: name of the log, then the arguments of compute_certificate.jl
    LOG=$WORK/$1.log
    shift
    for attempt in $(seq 1 30); do
        echo "[$(date -Iseconds)] attempt $attempt" >> "$LOG"
        # the heap size hint makes Julia's GC collect before the memory cap kills the process
        systemd-run --user --scope -p MemoryMax=4500M -p MemorySwapMax=0 \
            julia --heap-size-hint=3G --project=. studies/hexagon/certificate/compute_certificate.jl "$@" >> "$LOG" 2>&1 && return 0
        sleep 10
    done
    echo "[$(date -Iseconds)] gave up" >> "$LOG"
    return 1
}

# the dimensions of the cones of Trop, once, before the parallel stages
run cones --stage=cones || exit 1

run setups_1 --stage=setups --part=1 --parts=3 &
run setups_2 --stage=setups --part=2 --parts=3 &
run setups_3 --stage=setups --part=3 --parts=3 &
wait

run greedy_1 --stage=greedy --part=1 --parts=4 &
run greedy_2 --stage=greedy --part=2 --parts=4 &
run greedy_3 --stage=greedy --part=3 --parts=4 &
run greedy_4 --stage=greedy --part=4 --parts=4 &
wait

[ -f $WORK/greedy/R.mrdi ] || run merge --stage=merge || exit 1

# the coordinates in the order of 2026-11-07: first one of each D12-orbit, so these are done early
COORDINATES=15,24,20,17,54,38,39,41,14,11,7,1,16,18,19,21,22,23,25,26,27,28,29,30,31,32,33,34,35,36,37,40,42,43,44,45,46,47,48,49,50,51,52,53,55,56,8,9,10,12,13,2,3,4,5,6
run cover_1 --stage=cover --coordinates=$COORDINATES --part=1 --parts=3 &
run cover_2 --stage=cover --coordinates=$COORDINATES --part=2 --parts=3 &
run cover_3 --stage=cover --coordinates=$COORDINATES --part=3 --parts=3 &
wait

run assemble --stage=assemble
