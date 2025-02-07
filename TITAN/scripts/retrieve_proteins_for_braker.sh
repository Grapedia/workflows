#!/bin/bash

cat $1 | sed 1d | awk -F, '{print $2}' | sed "s/^/\/protein_path\//" | tr -d '\r' | tr '\n' ',' | sed 's/,$//'