#!/bin/bash

ffmpeg -framerate 24 -i rendered_results/render.%04d.png -c:v libx264 -pix_fmt yuv420p -crf 23 output.mp4

