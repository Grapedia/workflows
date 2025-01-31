#!/bin/bash

ls -l ${1}/*.gtf | sed 's/.*-> .*work\//\/work\//' | tr '\n' ' ' | sed 's/ $//'