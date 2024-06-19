#!/bin/bash

ls -l ${1}/glimmerhmm_prediction_*gff | sed 's/.*-> .*work\//\/work\//' | tr '\n' ' ' | sed 's/ $//'
