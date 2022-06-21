#python3 mirexplore_server.py --textmine /home/mjoppich/owncloud/data/miRExplore/textmine/ --obodir /home/mjoppich/owncloud/data/miRExplore/obodir/ --sentdir /home/mjoppich/dev/data/pmid/ --feedback /home/mjoppich/owncloud/data/miRExplore/obodir/feedback --port 65500


#ln -s aggregated_litcovid/ aggregated_pmid

CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

BASE=/mnt/biocluster/projekte/Corona2020/Texts/
SENTDIR=$BASE/output/LitCovid/singledocs/

python3 $CURDIR/mirexplore_server.py --textmine $BASE --obodir $BASE/obodir/ --sentdir $SENTDIR --feedback $BASE/obodir/feedback --port 65500