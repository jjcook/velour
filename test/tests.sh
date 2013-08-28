if [ -z "$TESTBASEDIR" ] ; then
    echo "ERROR: Do not run this script $BASH_SOURCE directly." >& 2
    exit 1
fi

DIRECTORIES="pregraph"

for dir in $DIRECTORIES ; do
    echo "[[[ $dir ]]]"
    pushd $dir > /dev/null
    source ./tests.sh
    popd > /dev/null
done
