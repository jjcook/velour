#!/bin/bash
#
# partitioning velour run script

# required parameter environment variables: $PARTS $MINIK

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

set -o errexit
set -o nounset

WORK=$1
FULLK=$2
shift
shift
INPUT="$@"

if [ $PARTS -ge 1 ] ; then
    VELOUR_PARTITIONS=$PARTS
else
    echo "Error: PARTS variable not valid." >&2
    exit 1
fi

if [ $MINIK -gt $((2 * $FULLK)) ] ; then
    echo "Error: MINIK variable not valid." >&2
    exit 1
fi

set +o nounset
PREGRAPH_PARTITIONS=0
if [ -n "$PGPARTS" ] ; then
    if [ $PGPARTS -gt 1 ] ; then
        PREGRAPH_PARTITIONS=$PGPARTS
        if [ ! -x $HMETIS ] ; then
            echo "Error: HMETIS variable not valid executable." >&2
            exit 1
        fi
    else
        echo "Error: PGPARTS variable not valid." >&2
        exit 1
    fi
fi
set -o nounset

set +o nounset
if [ -z "$OPTS" ] ; then
    OPTS=""
fi
if [ -z "$RETRY" ] ; then
    RETRY=""
fi
set -o nounset

VELOUR_HOME=${VELOUR_HOME-.}
VELOUR=${VELOUR-$VELOUR_HOME/velour}

if [ -d "$WORK" ] ; then
    rm -rf "$WORK/work"
    rm -f "$WORK/*.log"
    rm -f "$WORK/*.txt"
    rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK"

echo "PART: Partitioning input..."
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -part $VELOUR_PARTITIONS $MINIK $INPUT >& $WORK/partitioning.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour partitioner failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

#du -ms $WORK > "$WORK/du-partitioning.txt"

# get number of actual partitions created
VELOUR_PARTITIONS=`cat $WORK/work/common.partitions`

echo "PART: Flowing partitions MPI parallel..."

  while true ; do

  set +o errexit
  mpirun -l -n $VELOUR_PARTITIONS $VELOUR "$WORK/work" $FULLK $OPTS -flow $VELOUR_PARTITIONS >& $WORK/flowing-mpi.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour MPI flowing failed.  Exit code $RETVAL." >&2
    if [ -z "$RETRY" ] ; then
        exit $RETVAL
    fi
  else
    break # success, don't retry
  fi

  done # end while loop

for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do
  INBOXES=""
  if [ $p -ne 1 ] ; then
    for ((i=1; i < $p ; i++)) ; do
        NEXT_INBOX="$WORK/work/inbox/$p/InboxBucket-from-$i.bucket"
        INBOXES="${INBOXES} ${NEXT_INBOX}"
    done
  fi
  # flowing success for partition.  delete inputs.
  #rm -f "$WORK/work/loom/Subsequences-$p.loom"  $INBOXES
done

echo "PART: Quilting final bucket..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/quilt/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
#du -ms $WORK/work > "$WORK/du-final.txt"
set +o errexit
if [ $PREGRAPH_PARTITIONS -ne 0 ] ; then
  $VELOUR "$WORK/work" $FULLK $OPTS -pgpart -quilt -bucket $FINALS >& $WORK/quilting.log
else
  $VELOUR "$WORK/work" $FULLK $OPTS -quilt $VELOUR_PARTITIONS -bucket $FINALS >& $WORK/quilting.log
fi
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour quilting of final bucket failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

# quilting success.  delete final buckets.
#rm -f "$WORK"/work/FinalBucket-from-*.bucket
#rmdir "$WORK"/work/inbox/* "$WORK"/work/inbox

# optionally, partition pregraph
if [ $PREGRAPH_PARTITIONS -ne 0 ] ; then
  echo "VELOUR: Partitioning PreGraph using hMetis..."
  set +o errexit
  $HMETIS "$WORK/work/PreGraph.metis" $PREGRAPH_PARTITIONS >& $WORK/hmetis.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour hMetis partitioning of Pregraph failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi

  echo "VELOUR: Distributing nodes to partitioned PreGraphs..."
  set +o errexit
  #TODO $VELOUR "$WORK/work" $FULLK $OPTS -pgdist $PREGRAPH_PARTITIONS -pgfilter 27 $WORK/work/PreGraph.quilt >& $WORK/distribution.log
  $VELOUR "$WORK/work" $FULLK $OPTS -pgdist $PREGRAPH_PARTITIONS $WORK/work/PreGraph.quilt >& $WORK/distribution.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour distribution of partitioned PreGraphs failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi
fi

#
# DONE!
#

touch $WORK/SUCCESS

