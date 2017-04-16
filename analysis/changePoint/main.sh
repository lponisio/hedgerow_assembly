#!/usr/bin/env bash
RScript dataPrep.R
python hedegrows.py
RScript prepCptOutput.R
python postChangePoint.py paradigms.txt
