#!/bin/bash

##
## VELOUR SINGLE-END ASSEMBLY
##
##   partitioned single end assembly
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
    echo "Usage: assemble.sh outputDir oddKmerLength numParts gigabytesMemory fastaFiles" >&2
    exit 1
fi
set -o nounset

WORK=$1
FULLK=$2
VELOUR_PARTS=$3
VELOUR_MEMORY=$4
shift
shift
shift
shift
INPUT="$@"

#echo "VELOUR: k-mer length is $FULLK"

MINIK=13

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
    echo "User Error: VELOUR_PARTS variable must be equal to 4, 16, 64, 256, or 1024." >&2
    exit 1
fi
set -o nounset

set +o nounset
if [ -n "$VELOUR_MEMORY" ] ; then
    OPTS="-mem $(expr $VELOUR_MEMORY \* 1024)"
else
    echo "User Error: VELOUR_MEMORY variable not defined." >&2
    exit 1
fi
set -o nounset

set +o nounset
if [ -z "$RETRY" ] ; then
    RETRY=""
fi
set -o nounset

VELOUR_ROOT=$(dirname $(readlink -f $0))
VELOUR=$VELOUR_ROOT/velour

if [ ! -d "$VELOUR_ROOT/minikmer_ptables" ] ; then
    echo "User Error: please install the mini-kmer tables in '$VELOUR_ROOT/minikmer_ptables'" >&2
    exit 1
fi

if [ -d "$WORK" ] ; then
    rm -rf "$WORK/work"
    rm -f "$WORK/*.log"
    rm -f "$WORK/*.txt"
    rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK/work"

#echo "VELOUR: desired maximum physical memory use is $VELOUR_MEMORY gigabytes"

echo "VELOUR: Partitioning input $VELOUR_PARTITIONS ways..."
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -part $VELOUR_PARTITIONS $MINIK $INPUT >& "$WORK/partitioning.log"
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour partitioner failed.  Exit code $RETVAL." >&2
  tail -n 4 "$WORK/partitioning.log" >&2
  exit $RETVAL
fi

# get number of actual partitions created
VELOUR_PARTITIONS=`cat "$WORK/work/common.partitions"`

echo "VELOUR: Flowing each partition..."
for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do

  while true ; do

  echo "VELOUR:  partition $p of $VELOUR_PARTITIONS"
  INBOXES=""
  if [ $p -ne 1 ] ; then
    for ((i=1; i < $p ; i++)) ; do
        NEXT_INBOX="$WORK/work/inbox-for-$p/InboxBucket-from-$i.bucket"
        INBOXES="${INBOXES} ${NEXT_INBOX}"
    done
    set +o errexit
    $VELOUR "$WORK/work" $FULLK $OPTS -flow $VELOUR_PARTITIONS "$WORK/work/Subsequences-$p.loom" -bucket $INBOXES >& "$WORK/flowing-$p.log"
    RETVAL=$?
    set -o errexit
  else
    set +o errexit
    $VELOUR "$WORK/work" $FULLK $OPTS -flow $VELOUR_PARTITIONS "$WORK/work/Subsequences-$p.loom" >& "$WORK/flowing-$p.log"
    RETVAL=$?
    set -o errexit
  fi
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour flowing of partition $p failed.  Exit code $RETVAL." >&2
    if [ -z "$RETRY" ] ; then
        exit $RETVAL
    fi
  else
    break # success, don't retry
  fi

  done # end while loop

  # flowing success for partition.  delete inputs.
  rm -f "$WORK/work/Subsequences-$p.loom"  $INBOXES
done

echo "VELOUR: Finishing single-end assembly, with no coverage cutoff..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
set +o errexit
#$VELOUR "$WORK/work" $FULLK $OPTS -bubble_removal -quilt -bucket $FINALS >& $WORK/nocovcutoff.log
$VELOUR "$WORK/work" $FULLK $OPTS -quilt -bucket $FINALS >& "$WORK/nocovcutoff.log"
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour single-end assembly finishing failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

mv "$WORK/work/PreGraph" "$WORK/work/PreGraph.wcov"

echo "VELOUR: Generating node k-mer coverage database..."
grep NODE "$WORK/work/PreGraph.wcov" | cut -f 3,5 | sort -n > "$WORK/nodekmercoverage.txt"

echo "VELOUR: Computing length weighted median contig coverage..."
HALF_LENGTH=`awk '{sum+=$1}END{printf "%d", sum/2}' "$WORK/nodekmercoverage.txt"`
awk '{sum+=$1; if(sum >= $HALF_LENGTH) { printf "lwmcc           = %10.2f\n", $2; exit }}' "$WORK/nodekmercoverage.txt"

# strip the coverage information to produce a Velvet-compatible format
echo "VELOUR: Producing Velvet-compatible PreGraph file..."
cat "$WORK/work/PreGraph.wcov" | sed 's/\(^NODE\t[0-9]\+\t[0-9]\+\).*/\1/g' > "$WORK/PreGraph"

echo "VELOUR: Computing single-end assembly statistics, min contig length of 100..."
$VELOUR_ROOT/contig_stats.pl -k $FULLK -m 100 "$WORK/PreGraph"

# finishing success.  delete final buckets.
#rm -f "$WORK"/work/FinalBucket-from-*.bucket
rmdir "$WORK"/work/inbox-for-*

echo "VELOUR: Done."

#
# DONE!
#

touch "$WORK/SUCCESS"

