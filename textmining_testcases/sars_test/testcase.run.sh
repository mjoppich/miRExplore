

SCRIPTPATH=$(dirname "$0")

python3 $SCRIPTPATH/../../python/textmining/textmineDocument.py --input $SCRIPTPATH/testcase.sent --synonyms $SCRIPTPATH/testcase.syn --output - --process-level single