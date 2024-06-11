#!/bin/bash
# This script joins several files with names $1XX.gff into one
# It keeps the header of 13 lines from the first file and removes them from the others
# $1 is the first argument that is passed when running the script, which must be the path with the prefix of the files to join
# $2 is the second argument that is passed when running the script, which must be the name of the organism

# We create the output file in the desired folder with the name {organism}.gff
touch $2

# We get the name of the first file
first_file=$(ls $1*.gff | head -n 1)

# We copy the header and content of the first file to the output file
cat $first_file > $2

# We loop through all the files with a for loop
for file in $1*.gff
do
    # We skip the first file, since we have already copied its header
    if [ $file != $first_file ]
    then
        # We remove the first 13 lines of the current file and add them to the output file
        cat $file >> $2
    fi
done

# Confirmation message
echo "The file $2 has been created successfully"
