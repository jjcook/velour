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

set +o nounset
if [ -z "$OPTS" ] ; then
    OPTS=""
fi
set -o nounset

VELOUR_HOME=${VELOUR_HOME-.}
VELOUR=${VELOUR-$VELOUR_HOME/velour}

if [ ! -d "$WORK/work" ] ; then
  echo "Missing directory." >&2
  exit 1
fi

rm -f "$WORK/SUCCESS"

# get number of actual partitions created
VELOUR_PARTITIONS=`cat $WORK/work/common.partitions`

echo "PART: Quilting final bucket..."
FINALS=""
for ((p=1; p <= $VELOUR_PARTITIONS; p++)) ; do
  NEXT_FINAL="$WORK/work/FinalBucket-from-$p.bucket"
  FINALS="${FINALS} ${NEXT_FINAL}"
done
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

