#!/bin/bash
banjo="/home/samanthagonzales/software/banjo_2.2.0/banjo.jar"

analysis_name="nov17.mouse.liver"
duration="32"

settings_file="${analysis_name}.banjo.${duration}h.settings.txt"

echo "Running banjo on data in: ${PWD}"
nohup java -jar $banjo settingsFile=$settings_file > ./logs/banjo_${analysis_name}_${duration}h_$(date +%F).log 2>&1 &
