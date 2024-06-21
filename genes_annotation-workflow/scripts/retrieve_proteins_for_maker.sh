#!/bin/bash

cat $1 | awk -F, '$3 ~ "yes" {print $2}' | sed "s/^/\/protein_path\//" | tr '\n' ',' | sed 's/,$//'
