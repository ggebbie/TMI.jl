#!/bin/sh
# Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

# Old Google ID = 1VUkucxO7aKVsVD8XpHUT9GCN-_CHz9aU
# Google ID = 1LddFEsGVCor1IOr4sYZ1S8SAVR9nDQ4R
#wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1YTVjvHNQ1fAGjGyBH2-sLHk6gOrg0igK' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1YTVjvHNQ1fAGjGyBH2-sLHk6gOrg0igK" -O TMI_nordic_271x205x46_B23.mat.gz && rm -rf /tmp/cookies.txt
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1YTVjvHNQ1fAGjGyBH2-sLHk6gOrg0igK' -O TMI_nordic_271x205x46_B23.mat.gz