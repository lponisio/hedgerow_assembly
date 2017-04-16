#!/usr/bin/env bash
python hedegrows.py
RScript prepCptOutput.R
python postChangePoint.py paradigms.txt
