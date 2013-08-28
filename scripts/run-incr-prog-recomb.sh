#!/bin/bash
#
# incremental recombination
#

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

if [ $PARTS -ge 2 ] ; then
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
if [ -z "$OPTS" ] ; then
    OPTS=""
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

mkdir -p "$WORK/work"

echo "PART: Partitioning input..."
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -part $VELOUR_PARTITIONS $MINIK $INPUT >& $WORK/partitioning.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour partitioner failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

# get number of actual partitions created
VELOUR_PARTITIONS=`cat $WORK/work/common.partitions`

echo "PART: Looming each partition..."
for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do
  echo "PART:  partition $p of $VELOUR_PARTITIONS"
  set +o errexit
  $VELOUR "$WORK/work" $FULLK $OPTS -loom "$WORK/work/Subsequences-$p.loom" >& $WORK/looming-$p.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour looming of partition $p failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi
done

echo "PART: Bucket'ing each partition..."
for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do
  echo "PART:  partition $p of $VELOUR_PARTITIONS"
  set +o errexit
  $VELOUR "$WORK/work" $FULLK $OPTS -makebucket "$WORK/work/Partition-$p.quilt" >& $WORK/makebucket-$p.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour makebucket of partition $p failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi
done


echo "PART: Incremental progressive recombination..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/Partition-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -noslice -incr $FINALS >& $WORK/incr.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour incremental progressive recombination failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

touch $WORK/SUCCESS

