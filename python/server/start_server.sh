#python3 mirexplore_server.py --textmine /home/mjoppich/owncloud/data/miRExplore/textmine/ --obodir /home/mjoppich/owncloud/data/miRExplore/obodir/ --sentdir /home/mjoppich/dev/data/pmid/ --feedback /home/mjoppich/owncloud/data/miRExplore/obodir/feedback --port 65500

CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

BASE=/mnt/w/miRExplore_pmid_pmc/
PMIDSENTS=/mnt/w/PubMed/

python3 $CURDIR/mirexplore_server_fast.py --textmine $BASE --obodir $BASE/obodir/ --sentdir $PMIDSENTS --feedback $BASE/obodir/feedback --port 65500