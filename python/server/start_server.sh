#python3 mirexplore_server.py --textmine /home/mjoppich/owncloud/data/miRExplore/textmine/ --obodir /home/mjoppich/owncloud/data/miRExplore/obodir/ --sentdir /home/mjoppich/dev/data/pmid/ --feedback /home/mjoppich/owncloud/data/miRExplore/obodir/feedback --port 65500


BASE=/mnt/d/dev/data/pmid_jun2020/

python3 mirexplore_server.py --textmine $BASE --obodir $BASE/obodir/ --sentdir $BASE/pmid/ --feedback $BASE/obodir/feedback --port 65500