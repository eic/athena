#!/bin/bash
npsim --runType run --compactFile compact/subsystem_views/drich_only.xml --macro test.mac --outputFile test.root --enableG4GPS
root -b -q drawHits.C
