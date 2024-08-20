#!/bin/bash

find ${1}/ -maxdepth 2 -type l -printf '%f\n' | sed 's/.*-> .*work\//\/work\//' | tr '\n' ',' | sed 's/,$//'