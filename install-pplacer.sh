#!/usr/bin/env bash

if [[ "$OSTYPE" == "darwin"* ]]; then
  DOWNLOAD_URL="https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Darwin-v1.1.alpha17.zip"
elif [[ "$OSTYPE" == "linux"* ]]; then
  DOWNLOAD_URL="https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip"
else
  echo "Detected OS version (${OSTYPE}) is not supported. Aborting."
  exit 1
fi

echo "Fetching pplacer binaries from ${DOWNLOAD_URL}..."
curl -L "${DOWNLOAD_URL}" > pplacer.zip

echo "Extracting..."
unzip pplacer.zip
rm pplacer.zip

echo "Installing pplacer in $CONDA_PREFIX..."
if [[ ! -d "$CONDA_PREFIX/bin/" ]]; then
  mkdir $CONDA_PREFIX/bin/
fi
mv pplacer*/guppy $CONDA_PREFIX/bin/
mv pplacer*/pplacer $CONDA_PREFIX/bin/
mv pplacer*/rppr $CONDA_PREFIX/bin/

mkdir $CONDA_PREFIX/bin/scripts/
mv pplacer*/scripts/* $CONDA_PREFIX/bin/scripts/
rm -r pplacer*

echo "Testing installation..."
if [[ $(which pplacer) == "$CONDA_PREFIX/bin"* ]]; then
  echo "Success!"
# TODO: make sure later that this really is not necessary (try to install with conda on macOS and Ubuntu)
#  pplacer --version
#  retVal=$?
#  if [[ $retVal -ne 0 ]]; then
#    echo "pplacer was installed in the correct location but there was a problem with the installation."
#    echo "See below for the potential error message:"
#    pplacer --version
#    exit 1
#  else
#    echo "Success!"
#  fi
else
  echo "Installation failed."
  exit 1
fi
