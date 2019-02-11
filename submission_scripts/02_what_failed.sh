#!/bin/bash

# A script that checks all output folders ("out_*") and checks to see which
# jobs were successful. If they were successful, then the script echos the
# successful run folder. If there are multiple successful folders, it echos
# a crap error message

IRWIN_SHARE='/nbi/group-data/ifs/JIC/Judith-Irwin'
READ_DIR_1=$IRWIN_SHARE'/2017_11_21_eleri_sequencing/C101HW17071262/raw_data'
READ_DIR_2=$IRWIN_SHARE'/Judith/C101HW17071262/raw_data/'

UNIQUE_IDENTIFIERS=($(find $READ_DIR_1 $READ_DIR_2 -name "*.fq.gz" | \
                      awk '{split($0, a, "/");
                          gsub(/_[1-2]\.fq\.gz/, "", a[length(a)]);
                          print(a[length(a)])}' | \
                      sort -u))

OUTPUT_FOLDERS=($(find ./ -type d -name "out_*"))

for (( ident_idx=0; ident_idx<${#UNIQUE_IDENTIFIERS[@]}; ident_idx++ ));
do
    IDENTIFIER=${UNIQUE_IDENTIFIERS[$ident_idx]}
    SUCCESS_FOLDER=''
    for (( out_idx=0; out_idx<${#OUTPUT_FOLDERS[@]}; out_idx++ ));
    do
        if [[ $(sacct --jobs $(find ${OUTPUT_FOLDERS[$out_idx]}/run_logs/ \
            -name *$IDENTIFIER*.err | \
            awk '{split($0, a, "."); print(a[4])}') | grep -c COMPLETED) == "9" ]];
        then
            if [[ "$SUCCESS_FOLDER" == "" ]];
            then
                SUCCESS_FOLDER=${OUTPUT_FOLDERS[$out_idx]}
            else
                echo "ERROR $IDENTIFIER"
            fi
        fi
    done
    if [[ "$SUCCESS_FOLDER" == "" ]];
    then
        SUCCESS_FOLDER='failed'
    fi
    echo "$IDENTIFIER $SUCCESS_FOLDER"
done
