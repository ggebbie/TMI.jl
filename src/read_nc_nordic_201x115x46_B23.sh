#!/bin/sh
# Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

# Old Google ID = 1mAmJKS3xeNFjKLLTMBpgSS-aWZ1WHNm-
# Google ID = 1aLCoNkAMSujC-ImX6Xw-mMK5qQeWwZHd

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1aLCoNkAMSujC-ImX6Xw-mMK5qQeWwZHd' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1aLCoNkAMSujC-ImX6Xw-mMK5qQeWwZHd" -O TMI_nordic_201x115x46_B23.nc && rm -rf /tmp/cookies.txt
