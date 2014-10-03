#!/bin/bash
#
# partitioning velour run script

# required parameter environment variables: $PARTS $MINIK

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

# TODO not updated or deleted

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
$VELOUR "$WORK/work" $FULLK -part $VELOUR_PARTITIONS $MINIK $INPUT >& $WORK/partitioning.log
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
  $VELOUR "$WORK/work" $FULLK -loom "$WORK/work/Subsequences-$p.loom" >& $WORK/looming-$p.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour looming of partition $p failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi
done

echo "PART: Quilting partitions..."
QUILTS=
for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do
    QUILTS="${QUILTS} $WORK/work/Partition-$p.quilt"
done
set +o errexit
$VELOUR "$WORK/work" $FULLK -quilt $QUILTS >& $WORK/quilting.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour quilting of partitions failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

touch $WORK/SUCCESS

