#!/usr/bin/env bash
set -eu
set -o pipefail
function usage {
	cat <<EOF
Usage: $0 dir

Dispatch the given annovar input files as separate jobs to run.
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


logdir="/hpf/largeprojects/agoldenb/tal/logs/annovar/$dir_name"
mkdir -pv "$logdir"
mkdir -pv "$dir/scripts"

for file in $dir/*.avinput; do
	f=`basename $file .avinput`
	script="$dir/scripts/dispatch_$f.sh"
	cat > "$script" <<EOF
#!/usr/bin/env bash
#PBS -N "anno_$f"
#PBS -l vmem="$memory"
#PBS -joe "$logdir"

set -eu
#temp=\$TMPDIR/$f
echo "Log directory: $logdir" >&2
echo "Current directory: \$(pwd)" >&2
#echo "Temp file: \$temp" >&2
echo "Input file: $file" >&2

/hpf/tools/centos6/annovar/2013.08.23/table_annovar.pl --remove -protocol refGene,ljb23_all,snp137NonFlagged,esp6500si_all -operation g,f,f,f -buildver hg19  "$file" /hpf/tools/centos6/annovar/2013.08.23/humandb
EOF

	sleep 1
	qsub -S /bin/sh "$script"
done
