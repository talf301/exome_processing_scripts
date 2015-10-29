#!/usr/bin/env bash
set -eu
set -o pipefail

function usage {
    cat << EOF
Usage: $0 vcf outdir

Splits vcf into chunks, submits jobs for running cadd on each chunk, 
and leaves a script for combining in the outdirectory. Can be run 
again to resubmit failed jobs
EOF
    exit 1
}

splitcount=100000
vcf=$1
v=`basename $vcf`
out=$2
logs=/hpf/largeprojects/agoldenb/tal/cadd/logs/"$v"_logs

mkdir -pv "$out"/scripts
mkdir -pv "$logs"

# Split
split -l $splitcount "$vcf" "$out"/"$v".

# Do submissions
for file in "$out"/"$v".*; do
    f=`basename $file`
    if [[ -s "$out"/"$f".cadd ]]; then
        echo "VCF part already processed: $f" >&2
        continue
    fi
    script="$out"/scripts/dispatch_"$f".sh
    cat > "$script" << EOF
#!/usr/bin/env bash
#PBS -S /bin/bash
#PBS -N $f
#PBS -l vmem=5g
#PBS -joe $logs

cat $out/$f | /hpf/tools/centos6/python/2.7.6/bin/python /hpf/largeprojects/agoldenb/tal/scripts/extractCADDscores.py -p /hpf/largeprojects/agoldenb/tal/cadd/whole_genome_SNVs.tsv.gz > /localhd/\$PBS_JOBID/$f.temp

# Big move
mv /localhd/\$PBS_JOBID/"$f".temp "$out"/
# Small move
mv "$out"/"$f".temp "$out"/"$f".cadd
EOF

    qsub "$script"
done
