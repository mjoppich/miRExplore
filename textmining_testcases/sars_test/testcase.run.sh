

SCRIPTPATH=$(dirname "$0")

# report test1 test3
python3 $SCRIPTPATH/../../python/textmining/textmineDocument.py --input $SCRIPTPATH/testcase.sent --synonyms $SCRIPTPATH/testcase.syn --output - --process-level single

# report test1 test2 test3
python3 $SCRIPTPATH/../../python/textmining/textmineDocument.py --input $SCRIPTPATH/testcase.sent --synonyms $SCRIPTPATH/testcase.syn --test-is-word --output - --process-level single

# report test3
python3 $SCRIPTPATH/../../python/textmining/textmineDocument.py --input $SCRIPTPATH/testcase.sent --synonyms $SCRIPTPATH/testcase.syn --submatch-exclude $SCRIPTPATH/testcase.exclude --output - --process-level single