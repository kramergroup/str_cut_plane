#!/bin/bash
# This script is rather ugly. Should do a Python version written better.

# Use Stdin instead of filename
if ! [ $1 ] 
then
    echo "Using stdin" >&2
    awk 'NR==1{s=$0; print "#from stdin"} NR>1&&NR<6 {print}NR==6{print s; print} NR>6' 
    exit 0
fi

# Change file in place
if [ "$2" == "-i" ] 
then
   tmpfile=$(mktemp /tmp/add_elem.XXXXXX)
   echo "Changing file $1 in place using tempfile $tmpfile">&2
   comment=$(echo $1 | sed -E "s:.+r_:: ; s:.vasp::")
   comment="r_$comment"
   awk -v com="$comment" 'NR==1{s=$0; print "#",com} NR>1&&NR<6 {print}NR==6{print s; print} NR>6' $1 > $tmpfile 
   mv $tmpfile $1
   exit 0
fi

# Standard file-name case
comment=$(echo $1 | sed -E "s:.+r_:: ; s:.vasp::")
comment="r_$comment"
awk -v com="$comment" 'NR==1{s=$0; print "#",com} NR>1&&NR<6 {print}NR==6{print s; print} NR>6' $1 
exit 0

