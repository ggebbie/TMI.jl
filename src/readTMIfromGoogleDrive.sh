#!/bin/sh
# Extract data from Google Drive using your favorite method.

# 1) Download directly in Julia using `downloadTMI()` in `src/TMI.jl`

# 2) Download TMI_4deg_data.nc from Google Drive using command line
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7" -O TMI_4deg_data.nc && rm -rf /tmp/cookies.txt

# 3) download manually at:
#   https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7
#   https://drive.google.com/file/d/1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7/view?usp=sharing
