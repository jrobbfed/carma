#!/bin/bash

rm carmaorion13co.mpeg
~/Desktop/ffmpeg -f image2 -pattern_type glob -i 'channel*.png' -framerate 1 carmaorion13co.mpeg
