#!/bin/bash

##
## VELOUR SINGLE-END ASSEMBLY
##
##   apply user provided node k-mer coverage cutoff
##

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

set -o errexit
set -o nounset

set +o nounset
if [[ -n "$1" ]] && [[ -n "$2" ]] && [[ -n "$3" ]] && [[ -n "$4" ]] && [[ -n "$5" ]] ; then
  true
else
    echo "Usage: covcutoff.sh outputDir oddKmerLength numParts gigabytesMemory covCutoff" >&2
    exit 1
fi
set -o nounset

WORK=$1
FULLK=$2
VELOUR_PARTS=$3
VELOUR_MEMORY=$4
COV_CUTOFF=$5

set +o nounset
if [ "$VELOUR_PARTS" = 4 ] ; then
    VELOUR_PARTITIONS=$VELOUR_PARTS
elif [ "$VELOUR_PARTS" = 16 ] ; then
    VELOUR_PARTITIONS=$VELOUR_PARTS
elif [ "$VELOUR_PARTS" = 64 ] ; then
    VELOUR_PARTITIONS=$VELOUR_PARTS
elif [ "$VELOUR_PARTS" = 256 ] ; then
    VELOUR_PARTITIONS=$VELOUR_PARTS
elif [ "$VELOUR_PARTS" = 1024 ] ; then
    VELOUR_PARTITIONS=$VELOUR_PARTS
else
    echo "User Error: VELOUR_PARTS environment variable must be equal to 4, 16, 64, 256, or 1024." >&2
    exit 1
fi
set -o nounset

set +o nounset
if [ -n "$VELOUR_MEMORY" ] ; then
    OPTS="-mem $(expr $VELOUR_MEMORY \* 1024 \* 80 \/ 100)"
else
    echo "User Error: VELOUR_MEMORY variable not defined." >&2
    exit 1
fi
set -o nounset

VELOUR_ROOT=$(dirname $(readlink -f $0))
VELOUR=$VELOUR_ROOT/velour

if [ ! -d "$VELOUR_ROOT/minikmer_ptables" ] ; then
    echo "User Error: please install the mini-kmer tables in '$VELOUR_ROOT/minikmer_ptables'" >&2
    exit 1
fi

if [ -f "$WORK/2SUCCESS" ] ; then
  rm -f "$WORK/2SUCCESS"
fi

# get number of actual partitions created
VELOUR_PARTITIONS=`cat $WORK/work/common.partitions`

echo "VELOUR: Finishing single-end assembly with coverage cutoff of $COV_CUTOFF..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -bubble_removal -cov_cutoff $COV_CUTOFF -quilt -bucket $FINALS >& "$WORK/covcutoff.log"
#$VELOUR "$WORK/work" $FULLK $OPTS -cov_cutoff $COV_CUTOFF -quilt -bucket $FINALS >& "$WORK/covcutoff.log"
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour single-end assembly finishing failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

mv "$WORK/work/PreGraph" "$WORK/work/PreGraph.cov"
cat "$WORK/work/PreGraph.cov" | sed '1d' | sed 's/NODE/>NOD/g' | sed 's/[\t]/_/g' > $WORK/contigs.fa

# strip the coverage information to produce a Velvet-compatible format
echo "VELOUR: Producing Velvet-compatible PreGraph file..."
cat "$WORK/work/PreGraph.cov" | sed 's/\(^NODE\t[0-9]\+\t[0-9]\+\).*/\1/g' > "$WORK/PreGraph"

echo "VELOUR: Computing single-end assembly statistics, min contig length of 100..."
$VELOUR_ROOT/contig_stats.pl -k $FULLK -m 100 "$WORK/contigs.fa"

echo "VELOUR: Done."

#
# DONE!
#

touch "$WORK/2SUCCESS"

