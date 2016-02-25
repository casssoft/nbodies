#!/bin/bash

set -e

BASE_DIR=$(dirname $0)
EXE=${1:-lab09-serial}
CURR_DIR=$pwd

function returnDir {
  cd $CURR_DIR
}

trap returnDir EXIT

cd "$BASE_DIR/.."
echo -e "\033[0;36mMaking $EXE\033[0m"
make $EXE

echo
echo
echo -e "\033[0;36mRunning $EXE\033[0m"
./$EXE 100 resources/patrick.txt | scripts/diffPatrick.pl

echo -e "\033[0;32mSuccess!\033[0m"
