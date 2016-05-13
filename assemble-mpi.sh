#!/bin/bash

##
## VELOUR SINGLE-END ASSEMBLY
##
##   MPI-enabled partitioned single end assembly
##

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

# NOTE: Only uncomment (define) the NTHREADS variable if Velour was compiled with TBB threading support
# Optional value: specify maximum number of threads; default behavior is for TBB to choose
#NTHREADS=

#KEEP_ALL_INTERMEDIATE_FILES=true
KEEP_PARTITIONED_SUBSEQUENCES=true

#RETRY=true
#RESTART=1

if [[ ${NTHREADS+_} ]] ; then
    RETRY_NTHREADS=1
    #RETRY_NTHREADS=$NTHREADS
fi

set -o errexit
set -o nounset

set +o nounset
if [[ -n "$1" ]] && [[ -n "$2" ]] && [[ -n "$3" ]] && [[ -n "$4" ]] && [[ -n "$5" ]] ; then
  true
else
    echo "Usage: assemble-mpi.sh outputDir oddKmerLength numParts gigabytesMemory fastaFiles" >&2
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

set +o nounset
if [ -z "$RESTART" -a -d "$WORK" ] ; then
    echo "Usage Error:  Output directory exists:  '$WORK'"
    echo "              Please remove it or choose a different location."
    exit 1
fi
set -o nounset

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
    echo "Usage Error: VELOUR_PARTS variable must be equal to 4, 16, 64, 256, or 1024." >&2
    exit 1
fi
set -o nounset

set +o nounset
if [ -z "$KEEP_ALL_INTERMEDIATE_FILES" ] ; then
    KEEP_ALL_INTERMEDIATE_FILES=""
fi
if [ -z "$KEEP_PARTITIONED_SUBSEQUENCES" ] ; then
    KEEP_PARTITIONED_SUBSEQUENCES=""
fi
if [ -z "$RETRY" ] ; then
    RETRY=""
fi
if [ -z "$RESTART" ] ; then
    RESTART=""
fi
set -o nounset

VELOUR_ROOT=$(dirname $(readlink -f $0))
VELOUR=$VELOUR_ROOT/velour
OPTS="-mem $(expr $VELOUR_MEMORY \* 1024)"

set +o nounset
if [ -z "$NTHREADS" ] ; then
    THR=""
else
    THR="-thr $NTHREADS"
fi
set -o nounset

if [ ! -d "$VELOUR_ROOT/minikmer_ptables" ] ; then
    echo "Usage Error: please install the mini-kmer tables in '$VELOUR_ROOT/minikmer_ptables'" >&2
    exit 1
fi

if [ -n "$RETRY" ] ; then
    echo "VELOUR: RETRY is enabled. Warning: retry can infinite loop if failure is not transient."
fi
if [ -n "$RESTART" ] ; then
    echo "VELOUR: RESTART=$RESTART  Partitioning skipped and flowing restarted at partition $RESTART."
fi

if [ -n "$RESTART" ] ; then
    rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK"

#echo "VELOUR: desired maximum physical memory use is $VELOUR_MEMORY gigabytes"

if [ -z "$RESTART" ] ; then
    echo "VELOUR: Partitioning input $VELOUR_PARTITIONS ways..."
    set +o errexit
    $VELOUR "$WORK/work" $FULLK $OPTS $THR -part $VELOUR_PARTITIONS $MINIK $INPUT >& "$WORK/partitioning.log"
    RETVAL=$?
    set -o errexit
    if [ $RETVAL -ne 0 ] ; then
        echo "Velour partitioner failed.  Exit code $RETVAL." >&2
        tail -n 4 "$WORK/partitioning.log" >&2
        exit $RETVAL
    fi
fi

# get number of actual partitions created
VELOUR_PARTITIONS=`cat "$WORK/work/common.partitions"`

#if [ -n "$RESTART" -a "$RESTART" -gt $VELOUR_PARTITIONS ] ; then
#    echo "VELOUR: ERROR, RESTART=$RESTART is larger than the maximum partition index $VELOUR_PARTITIONS." >&2
#    exit 1
#fi

# if restarting, delete stale inbox buckets
#if [ -n "$RESTART" ] ; then
#    echo "VELOUR: Restarting, deleting stale inbox buckets..."
#    for ((p=$RESTART; p <=$VELOUR_PARTITIONS; p++)) ; do
#        rm -f "$WORK/work/quilt/FinalBucket-from-$p.bucket"
#        for ((i=($p+1) ; i <=$VELOUR_PARTITIONS; i++)) ; do
#            rm -f "$WORK/work/inbox/$i/InboxBucket-from-$p.bucket"
#        done
#    done
#fi

#if [ -n "$RESTART" ] ; then
#    START_INDEX=$RESTART
#else
#    START_INDEX=1
#fi

echo "VELOUR: Flowing partitions MPI parallel..."
#for ((p=$START_INDEX; p <=$VELOUR_PARTITIONS; p++)) ; do
#
#  while true ; do
#
#  echo "VELOUR:  partition $p of $VELOUR_PARTITIONS"
#  INBOXES=""
#  if [ $p -ne 1 ] ; then
#    for ((i=1; i < $p ; i++)) ; do
#        NEXT_INBOX="$WORK/work/inbox/$p/InboxBucket-from-$i.bucket"
#        INBOXES="${INBOXES} ${NEXT_INBOX}"
#    done
#  fi
  set +o errexit
  mpirun -l -n $VELOUR_PARTITIONS $VELOUR "$WORK/work" $FULLK $OPTS $THR -flow $VELOUR_PARTITIONS >& "$WORK/flowing-mpi.log"
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour MPI flowing failed.  Exit code $RETVAL." >&2
    if [ -z "$RETRY" ] ; then
        exit $RETVAL
    fi
  fi
#    if [[ ${NTHREADS+_} ]] ; then
#        THR="-thr $RETRY_NTHREADS"
#    fi
#  else
#    if [[ ${NTHREADS+_} ]] ; then
#        if [ -n "$NTHREADS" ] ; then
#            THR="-thr $NTHREADS"
#        else
#            THR=""
#        fi
#    fi
#    break # success, don't retry
#  fi
#
#  done # end while loop

#  # flowing success for partition.  delete inputs.
#  if [ -z "$KEEP_ALL_INTERMEDIATE_FILES" -a -z "$KEEP_PARTITIONED_SUBSEQUENCES" ] ; then
#    rm -f "$WORK/work/loom/Subsequences-$p.loom"
#  fi
#  if [ -z "$KEEP_ALL_INTERMEDIATE_FILES" ] ; then
#    rm -f $INBOXES
#  fi
#done

# TODO: use the coverage cutoff script like assemble.sh does

echo "VELOUR: Finishing single-end assembly, with no coverage cutoff..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/quilt/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
set +o errexit
#echo "        NOTE: Bubble removal ENABLED.  Single-end assembly does NOT require further finishing."
#$VELOUR "$WORK/work" $FULLK $OPTS $THR -bubble_removal -quilt $VELOUR_PARTITIONS -bucket $FINALS >& $WORK/nocovcutoff.log
echo "        WARNING: Bubble removal DISABLED.  Single-end assembly should be consumed by Velvet for further finishing."
$VELOUR "$WORK/work" $FULLK $OPTS $THR -quilt $VELOUR_PARTITIONS -bucket $FINALS >& "$WORK/nocovcutoff.log"
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
#if [ -z "$KEEP_ALL_INTERMEDIATE_FILES" ] ; then
#    rm -f "$WORK"/work/quilt/FinalBucket-from-*.bucket
#    rmdir "$WORK"/work/inbox/* "$WORK"/work/inbox
#fi

echo "VELOUR: Done."

#
# DONE!
#

touch "$WORK/SUCCESS"

