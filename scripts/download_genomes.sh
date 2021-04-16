#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset


function download_genome() {
    dest_file=data/seqs/$1
    url=$2
    check_sum=$3
    if [ ! -f ${dest_file} ];
    then
        wget --no-check-certificate -O ${dest_file} ${url} 
    fi
    echo ${check_sum} ${dest_file} | md5sum -c        
}

while read file location md5
do
     download_genome ${file} ${location} ${md5}
done < data/seq_info.tsv






