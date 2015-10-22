#!/usr/bin/env bash
set -eu
set -o pipefail
function usage {
	cat <<EOF
Usage: $0 dir

Dispatch the given VCFs for running CADD on
EOF
	exit 1
}

if [ $# -ne 1 ]; then
	usage
fi

dir=$1
dir_name=`basename $dir`
memory=5g
processors=1

mkdir -pv "$dir/scripts"

for file in $dir/*.vcf; do
	f=`basename $file .vcf`
	script="$dir/scripts/dispatch_$f.sh"
	cat > "$script" <<EOF
#!/usr/bin/env bash
#PBS -N "cadd_$f"
#PBS -l vmem="$memory"

set -eu

echo "Current directory: \$(pwd)" >&2
echo "Input file: $file" >&2


cat $file | /hpf/tools/centos6/python/2.7.6/bin/python /hpf/largeprojects/agoldenb/tal/scripts/extractCADDscores.py -p /hpf/largeprojects/agoldenb/tal/cadd/whole_genome_SNVs.tsv.gz | gzip -c > $dir/$f.scores.tsv.gz
EOF

	sleep 1
	qsub -S /bin/sh "$script"
done
