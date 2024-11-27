#!/bin/zsh
MOD=SignalAnalysisBenchmarks

DT=$(date '+%Y-%m-%d_%H:%M:%S')
TMP=/tmp/SignalAnalysis
[ -d $TMP ] || mkdir -p $TMP
FNAME="$TMP/regressions_$DT.json"
SCRIPT=$(cat <<EOF
using $MOD 
using BenchmarkTools
$MOD.load!("dsp";tune=false)
@info "running benchmark suite..."
results = run($MOD.SUITE["dsp"])
BenchmarkTools.save("$FNAME", results)
EOF
)
julia --project=. -e $SCRIPT --color=no
