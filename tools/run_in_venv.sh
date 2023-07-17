#!/usr/bin/env bash
# setup venv
CURRENT_DIR="$(cd "$(dirname "$0")" && pwd)"

ROOT_DIR=$CURRENT_DIR/../
VENV_DIR=$ROOT_DIR/venv/
SRC_DIR=$ROOT_DIR/src
source ${VENV_DIR}/bin/activate
export PYTHONPATH=${SRC_DIR}:$PYTHONPATH
python $1 "${@:2}"
