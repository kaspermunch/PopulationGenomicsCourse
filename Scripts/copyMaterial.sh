#!/bin/bash

fld_size=(`du -hs /usr/$COURSE`)
job_id=(`echo $HOSTNAME | tr "-" "\n"`)

while true; do
    read -p "Do you wish to copy ${fld_size[0]} of data? " yn
    case $yn in
        [Yy]* ) unlink $COURSE 2>/dev/null; \
                cp -fa /usr/$COURSE/ /work; \
                echo "----> COPY done"; \
                echo "----> FIND your files in uCloud at USER_FOLDER/Jobs/Genomics Sandbox/JOB_NAME (${job_id[1]})"/; \
                break;;
        [Nn]* ) echo "----> ABORT copy!!!"; exit;;
        * ) echo "----> PLEASE answer Y/y or N/n.";;
    esac
done