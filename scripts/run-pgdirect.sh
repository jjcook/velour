#!/bin/bash
#
# direct (non-partitioned) velour run script

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

if [ $PGPART -ge 2 ] ; then
    VELVET_PARTITIONS=$PGPART
else
    echo "Error: PGPART variable not valid." >&2
    exit 1
fi

VELOUR=${VELOUR-./velour}

if [ -d "$WORK" ] ; then
	rm -rf "$WORK/work"
	rm -f "$WORK/*.log"
	rm -f "$WORK/*.txt"
	rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK/work"

set +o errexit
$VELOUR $WORK/work $FULLK -pgpart $INPUT >& $WORK/direct.log
RETVAL=$?
set -o errexit

if [ $RETVAL -ne 0 ] ; then
  echo "Velour direct construction failed.  Exit code $RETVAL." >& 2
  exit $RETVAL
fi

touch $WORK/SUCCESS

