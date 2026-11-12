#!/bin/bash
# The verification of hexagon_certificate.mrdi (verify_certificate.jl) in 3 processes, then the summary.
# Each process is restarted after a crash or a kill by its memory cap (at most 30 times); checked cones are kept, so the
# script can also be stopped and started again.
# From the julia folder (or any folder):   bash studies/hexagon/certificate/run_verification_overnight.sh
# Progress: tail -n 2 studies/hexagon/certificate/verification/run_*.log

cd "$(dirname "$0")/../../.." || exit 1
mkdir -p studies/hexagon/certificate/verification

run_part() {   # $1: number of the part
    LOG=studies/hexagon/certificate/verification/run_$1.log
    for attempt in $(seq 1 30); do
        echo "[$(date -Iseconds)] attempt $attempt" >> "$LOG"
        # the heap size hint makes Julia's GC collect before the memory cap kills the process
        systemd-run --user --scope -p MemoryMax=4500M -p MemorySwapMax=0 \
            julia --heap-size-hint=3G --project=. studies/hexagon/certificate/verify_certificate.jl --part=$1 --parts=3 >> "$LOG" 2>&1 && break
        sleep 10
    done
    echo "[$(date -Iseconds)] part $1 ends" >> "$LOG"
}

run_part 1 &
run_part 2 &
run_part 3 &
wait

julia --project=. studies/hexagon/certificate/verify_certificate.jl --summary >> studies/hexagon/certificate/verification/summary.log 2>&1
