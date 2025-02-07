#!/bin/bash

if [ -d "$1" ]; then
    ls -l ${1}/* | grep "work" | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'
else
    echo ""
fi