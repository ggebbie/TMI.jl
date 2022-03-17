#!/bin/sh
# Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

# Or, download manually at:
#   https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7
#   https://drive.google.com/file/d/1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7/view?usp=sharing

# Download TMI_4deg_data.nc from Google Drive
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=11zD1nOfT6V7G0qIHdjK2pDGHFk-ExXwU' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=11zD1nOfT6V7G0qIHdjK2pDGHFk-ExXwU" -O TMI_modern_180x90x33_GH11_GH12.mat.gz && rm -rf /tmp/cookies.txt
