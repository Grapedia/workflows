#!/bin/bash

pattern="$2"

if [[ $pattern == "default" ]]
then
ls -l ${1}/*transcriptome.gtf | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'
elif [[ $pattern == "alt" ]]
then
ls -l ${1}/*transcriptome.AltCommands.gtf | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'
fi