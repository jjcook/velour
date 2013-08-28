if [ -z "$TESTBASEDIR" ] ; then
    echo "ERROR: Do not run this script $BASH_SOURCE directly." >& 2
    exit 1
fi

#DIRECTORIES="stage1-unit"
#DIRECTORIES="stage2-synthetic"
DIRECTORIES="stage1-unit stage2-synthetic"
#DIRECTORIES="stage1-unit stage2-synthetic stage3-genome"

for dir in $DIRECTORIES ; do
    echo "   [[ $dir ]]"
    pushd $dir > /dev/null
    runsuite
    popd > /dev/null
done
