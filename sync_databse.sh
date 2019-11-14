#!/bin/bash
# Just to syncronize the database
echo "Syncronizing database"
rsync -rhaP -e ssh database/* verdugm9@login.univie.ac.at:/home/verdugm9/html/database/


