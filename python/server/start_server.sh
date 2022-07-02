#python3 mirexplore_server.py --textmine /home/mjoppich/owncloud/data/miRExplore/textmine/ --obodir /home/mjoppich/owncloud/data/miRExplore/obodir/ --sentdir /home/mjoppich/dev/data/pmid/ --feedback /home/mjoppich/owncloud/data/miRExplore/obodir/feedback --port 65500

CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

BASE=${1:-/mnt/w/miRExplore_pmid_pmc/}
PMIDSENTS=${2:-/mnt/w/PubMed/}
PMCSENTS=${3:-/mnt/w/PubMedCentral/}
SERVICEPORT=${4:-65500}

echo $BASE
echo $PMIDSENTS
echo $PMCSENTS
echo $SERVICEPORT

sudo service mongodb start

python3 $CURDIR/mirexplore_server_fast.py --textmine $BASE --obodir $BASE/obodir/ --sentdir $PMIDSENTS --feedback $BASE/obodir/feedback --port $SERVICEPORT --sentdir-pmc $PMCSENTS --load-pmc