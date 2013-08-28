#!/bin/bash

# regression.sh
#
#   initiates regression suites and reports results
#
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

if [ -z "$TESTBASEDIR" ] ; then
    echo "ERROR: Do not run this script $BASH_SOURCE directly." >& 2
    exit 1
fi

set -o errexit
set -o nounset

##
## status counters
##

PASS_COUNT=0
FAIL_COUNT=0

EXPECT_PASS_COUNT=0
EXPECT_FAIL_COUNT=0
EXPECT_MUSTFAIL_COUNT=0

FIXED_COUNT=0
REGRESSION_COUNT=0

DIRECTORIES=pregraph

function finish ()
{
    local RETVAL=0
    
    echo ""
    echo "# # # # # # # # # # # # # # # # # # # #"
    echo ""
    echo "Passed           : $PASS_COUNT"
    echo "Failed           : $FAIL_COUNT"
    echo "Expected pass    : $EXPECT_PASS_COUNT"
    echo "Expected fail    : $EXPECT_FAIL_COUNT"
    echo "Expected mustfail: $EXPECT_MUSTFAIL_COUNT"
    echo ""

    if [ $FIXED_COUNT -gt 0 ] ; then
        echo -ne "\\033[0;32m" ;
    fi
    echo -e "***        FIXED : $FIXED_COUNT \\033[0m"

    if [ $REGRESSION_COUNT -gt 0 ] ; then
        echo -e "\\033[0;31m***  REGRESSIONS : $REGRESSION_COUNT \\033[0m"
        RETVAL=2
    else
        echo -e "\\033[0;32m***  REGRESSIONS : $REGRESSION_COUNT \\033[0m"
    fi

    return $RETVAL
}

COMPARE_PREGRAPH=$TESTSCRIPTS/compare-pregraph-simple.py
COMPARE_PASS=0
COMPARE_ERROR=1
COMPARE_FAIL=2

function check_for_regression ()
{
    local RESULT="$1"
    local PREFIX="$2"

    if [ $RESULT -eq $COMPARE_PASS ] ; then
        PASS_COUNT=$(( $PASS_COUNT + 1 ))
        if [ -e "$PREFIX.pass" ] ; then
            EXPECT_PASS_COUNT=$(( $EXPECT_PASS_COUNT + 1 ))
            echo "EXPECT: pass $PREFIX"
            return 0
        elif [ -e "$PREFIX.mustfail" ] ; then
            EXPECT_MUSTFAIL_COUNT=$(( $EXPECT_MUSTFAIL_COUNT + 1 ))
            echo "BAD: pass but mustfail $PREFIX"
            return 2
        elif [ -e "$PREFIX.FAIL" ] ; then
            EXPECT_FAIL_COUNT=$(( $EXPECT_FAIL_COUNT + 1 ))
            echo "GOOD: expected fail but didn't! $PREFIX"
            FIXED_COUNT=$(( $FIXED_COUNT + 1 ))
            return 0
        else
            echo "NEW TEST PASSED: $PREFIX"
            #touch $PREFIX.pass
            return 0
        fi
    fi

    if [ $RESULT -eq $COMPARE_ERROR ] ; then
        echo "ERROR: testing $PREFIX" >& 2
        exit 1
    fi

    if [ $RESULT -eq $COMPARE_FAIL ] ; then
        FAIL_COUNT=$(( $FAIL_COUNT + 1 ))
        if [ -e "$PREFIX.pass" ] ; then
            EXPECT_PASS_COUNT=$(( $EXPECT_PASS_COUNT + 1 ))
            echo "BAD: regression on $PREFIX"
            REGRESSION_COUNT=$(( $REGRESSION_COUNT + 1 ))
            return 2
        elif [ -e "$PREFIX.mustfail" ] ; then
            EXPECT_MUSTFAIL_COUNT=$(( $EXPECT_MUSTFAIL_COUNT + 1 ))
            echo "OK: mustfail $PREFIX"
            return 0
        elif [ -e "$PREFIX.FAIL" ] ; then
            EXPECT_FAIL_COUNT=$(( $EXPECT_FAIL_COUNT + 1 ))
            echo "EXPECT: fail $PREFIX"
            return 0
        else
            echo "NEW TEST FAILED: $PREFIX"
            #touch $PREFIX.FAIL
            return 2 # TODO: return 2
        fi
    fi

    echo "ERROR: script error" >& 2
    exit 1
}

# runtest(): executes single regression test
function runtest()
{
    local PREGRAPH=$1
    local WORK=$2
    local FULLK=$3
    local INPUT=$4

    if [ ! -r "$PREGRAPH" ] ; then
        echo "ERROR: Invalid reference PreGraph file: $PREGRAPH" >& 2
        exit 1
    fi

    echo "TEST: $WORK"
    set +o errexit
    $RUN $WORK $FULLK $INPUT
    local RETVAL=$?
    set -o errexit

    if [ $RETVAL -ne 0 ] ; then
        echo "ERROR: Run script $RUN had an error." >& 2
        exit 1
    fi

    if [ ! -r "$WORK/work/PreGraph" ] ; then
        echo "ERROR: Missing resultant PreGraph file: $WORK/work/Pregraph" >& 2
        exit 1
    fi

    set +o errexit
    export BESILENT=silent
    $COMPARE_PREGRAPH "$WORK/work/PreGraph" "$PREGRAPH"
    RETVAL=$?
    unset BESILENT

    check_for_regression $RETVAL "$PREGRAPH"
    set -o errexit
}

# runsuite(): executes tests in current directory
function runsuite ()
{
    for fna in `ls *.fna` ; do
        local ID=`basename $fna | tr '.' '\t' | cut -f1`
        local KMER=`basename $fna | tr '.' '\t' | cut -f2`
        for pg in `ls $ID.*.pregraph` ; do # TODO: handle case where no pregraph solutions exist
            if [ -d "test.$ID" ] ; then
                rm -rf test.$ID
            fi
            set +o errexit
            runtest $pg test.$ID $KMER $fna
            local RETVAL=$?
            set -o errexit
            #if [ $RETVAL -eq 0 -a -d "test.$ID" ] ; then
            #    rm -rf test.$ID
            #fi
        done
    done
    return 0
}

# lastly, actually initiate the regression tests
if [ -f "./tests.sh" ] ; then
    RUN=${RUN-$BASEDIR/scripts/run-direct.sh}
    source ./tests.sh
    finish
else
    echo "ERROR: Missing './tests.sh' script!" >& 2
    exit 1
fi

