#!/usr/bin/env bash
# Install Miniconda3 for Travis CI.

MINCONDA_DIR="$HOME/miniconda3"
MINICONDA_BIN="$MINICONDA_DIR/bin"

if [[ -d $MINICONDA_BIN ]]; then
	echo "INFO: Using Miniconda from cache"
else
	echo "INFO: Installing Miniconda"
	wget https://repo.anaconda.com/miniconda/Miniconda3-4.5.12-Linux-x86_64.sh
	bash Miniconda3-4.5.12-Linux-x86_64.sh -b -p -u $MINICONDA_DIR
	rm -f Miniconda3-4.5.12-Linux-x86_64.sh
fi

