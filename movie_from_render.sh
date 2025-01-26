#!/bin/bash

# ffmpeg -framerate 24 -i rendered_results/render.%04d.png -c:v libx264 -pix_fmt yuv420p -crf 23 output.mp4

ffmpeg -framerate 24 -i "rendered_results/Final Result/40-60-fast-final.%04d.png" -c:v libx264 -pix_fmt yuv420p -crf 23 output.mp4

