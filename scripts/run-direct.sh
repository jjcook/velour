#!/bin/bash
#
# direct (non-partitioned) velour run script

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

set +o nounset
if [ -z "$OPTS" ] ; then
    OPTS=""
fi
set -o nounset

VELOUR=${VELOUR-./velour}

if [ -d "$WORK" ] ; then
	rm -rf "$WORK/work"
	rm -f "$WORK/*.log"
	rm -f "$WORK/*.txt"
	rm -f "$WORK/SUCCESS"
fi

mkdir -p "$WORK/work"

set +o errexit
$VELOUR $WORK/work $FULLK $OPTS $INPUT >& $WORK/direct.log
RETVAL=$?
set -o errexit

if [ $RETVAL -ne 0 ] ; then
  echo "Velour direct construction failed.  Exit code $RETVAL." >& 2
  exit $RETVAL
fi

touch $WORK/SUCCESS

