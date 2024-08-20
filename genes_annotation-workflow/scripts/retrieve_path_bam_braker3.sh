#!/bin/bash

ls -l ${1}/* | grep "work" | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'