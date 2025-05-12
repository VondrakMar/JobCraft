#!/bin/bash
script_path=/home/mvondrak/development/JobCraft/jobcraft/restart_script.py 
python ${script_path} raven \
    --usedN 5 \
    -N 5 \
    --per_submit 2 \
    --wall_time="4:00:00"
