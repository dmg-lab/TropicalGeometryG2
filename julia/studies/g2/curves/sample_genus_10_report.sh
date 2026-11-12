#!/bin/bash
# Summary of a run of sample_genus_10.jl: experiments by status, curves by
# genus (balanced or not, core trivalent or not), candidates and verifications.
#   bash studies/g2/curves/sample_genus_10_report.sh [DIR]
DIR=${1:-results/curves/genus_10_sampling}
cat "$DIR"/samples/task_*.csv 2>/dev/null | grep -v '^experiment,' | awk -F, '
    { n++; status[$2 ~ /^failed/ ? "failed" : $2]++; seconds += $13 }
    $2 == "curve" {
        trivalent = ($8 ~ /^3\^[0-9]+$/) ? "trivalent core" : "other core"
        curves[$5 " " ($10 == "true" ? "balanced" : "unbalanced") ", " trivalent]++
    }
    END {
        printf "%d experiments, %.0f s each on average\n", n, n ? seconds / n : 0
        for (s in status) printf "  %s: %d\n", s, status[s]
        print "curves by genus:"
        for (c in curves) printf "  genus %s: %d\n", c, curves[c] | "sort -n -k2"
    }'
echo "candidates (genus 10, trivalent core):"
cat "$DIR/candidates.csv" 2>/dev/null || echo "  none"
echo "verified (certified census rows):"
cat "$DIR/verified.csv" 2>/dev/null || echo "  none"
