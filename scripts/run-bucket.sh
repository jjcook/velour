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

echo "PART: Splitting each partition..."
for ((p=1; p <=$VELOUR_PARTITIONS; p++)) ; do
  echo "PART:  split: partition $p of $VELOUR_PARTITIONS"
  set +o errexit
  $VELOUR "$WORK/work" $FULLK -split $VELOUR_PARTITIONS "$WORK/work/Partition-$p.quilt" >& $WORK/split-$p.log
  RETVAL=$?
  set -o errexit
  if [ $RETVAL -ne 0 ] ; then
    echo "Velour splitting of partition $p failed.  Exit code $RETVAL." >&2
    exit $RETVAL
  fi
  if [ $p -ne 1 ] ; then
    echo "PART:  split-bucket: partition $p of $VELOUR_PARTITIONS"
    INBOXES=""
    for ((i=1; i < $p ; i++)) ; do
        NEXT_INBOX="$WORK/work/inbox-for-$p/InboxBucket-from-$i.bucket"
        INBOXES="${INBOXES} ${NEXT_INBOX}"
    done
    set +o errexit
    $VELOUR "$WORK/work" $FULLK -split $VELOUR_PARTITIONS -bucket "$WORK/work/inbox-for-$p/SelfBucket-$p.bucket" $INBOXES >& $WORK/splitbucket-$p.log
    RETVAL=$?
    set -o errexit
    if [ $RETVAL -ne 0 ] ; then
        echo "Velour split-bucket of partition $p failed.  Exit code $RETVAL." >&2
        exit $RETVAL
    fi
  fi
done

echo "PART: Quilting final bucket..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
set +o errexit
$VELOUR "$WORK/work" $FULLK -quilt -bucket $FINALS >& $WORK/quilting.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour quilting of final bucket failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

touch $WORK/SUCCESS

