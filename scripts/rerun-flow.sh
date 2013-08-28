#!/bin/bash
#
# partitioning velour run script

# required parameter environment variables: $PARTS $MINIK $START

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
#    rm -f "$WORK/*.log"
#    rm -f "$WORK/*.txt"
    rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK/work"

# already partitioned into loom files

# get number of actual partitions created
VELOUR_PARTITIONS=`cat $WORK/work/common.partitions`

set +o nounset
if [ -z $START ] ; then
    echo "Error: START variable not set." >&2
    exit 1
fi
set -o nounset

if [ $START -gt $VELOUR_PARTITIONS ] ; then
    echo "Error: START variable not valid." >&2
    exit 1
fi

# delete stale inbox buckets
echo "PART: Deleting stale inbox buckets..."
for ((p=$START; p <=$VELOUR_PARTITIONS; p++)) ; do
    rm -f "$WORK/FinalBucket-from-$p.bucket"
    for ((i=($p+1) ; i <=$VELOUR_PARTITIONS; i++)) ; do
        rm -f "$WORK/work/inbox-for-$i/InboxBucket-from-$p.bucket"
    done
done

echo "PART: Flowing each partition..."
for ((p=$START; p <=$VELOUR_PARTITIONS; p++)) ; do

  while true ; do

  echo "PART:  partition $p of $VELOUR_PARTITIONS"
  INBOXES=""
  if [ $p -ne 1 ] ; then
    for ((i=1; i < $p ; i++)) ; do
        NEXT_INBOX="$WORK/work/inbox-for-$p/InboxBucket-from-$i.bucket"
        INBOXES="${INBOXES} ${NEXT_INBOX}"
    done
    #du -ms $WORK/work/inbox-for-$p > "$WORK/du-inbox-for-$p.txt"
    set +o errexit
    $VELOUR "$WORK/work" $FULLK $OPTS -flow $VELOUR_PARTITIONS "$WORK/work/Subsequences-$p.loom" -bucket $INBOXES >& $WORK/flowing-$p.log
    RETVAL=$?
    set -o errexit
  else
    set +o errexit
    $VELOUR "$WORK/work" $FULLK $OPTS -flow $VELOUR_PARTITIONS "$WORK/work/Subsequences-$p.loom" >& $WORK/flowing-$p.log
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

echo "PART: Quilting final bucket..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
#du -ms $WORK/work > "$WORK/du-final.txt"
set +o errexit
$VELOUR "$WORK/work" $FULLK $OPTS -quilt -bucket $FINALS >& $WORK/quilting.log
RETVAL=$?
set -o errexit
if [ $RETVAL -ne 0 ] ; then
  echo "Velour quilting of final bucket failed.  Exit code $RETVAL." >&2
  exit $RETVAL
fi

# quilting success.  delete final buckets.
#rm -f "$WORK/work/FinalBucket-from-*.bucket"

touch $WORK/SUCCESS

