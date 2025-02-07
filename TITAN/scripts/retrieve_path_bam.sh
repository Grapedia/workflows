#!/bin/bash

ls -l ${1}/*bam | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'
