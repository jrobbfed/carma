#!/bin/bash

rm carmaorion12co.mpeg
~/Desktop/ffmpeg -f image2 -pattern_type glob -i 'channel*.png' -framerate 1 carmaorion12co.mpeg
