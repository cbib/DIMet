# setup file for local venv
# setup venv
CURRENT_DIR="$(cd "$(dirname "$0")" && pwd)"

ROOT_DIR=$CURRENT_DIR/../
VENV_DIR=$CURRENT_DIR/../venv/
export LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
export CPPFLAGS="-I/opt/homebrew/opt/openblas/include"
export PKG_CONFIG_PATH="/opt/homebrew/opt/openblas/lib/pkgconfig"

# some pb with the default virtualenv install not working, had to force upgrade
# create the venv
cd ${ROOT_DIR} && virtualenv -p python3.9 venv --prompt="(DIMet)"

# setup in the venv

source ${VENV_DIR}/bin/activate && pip3.9 install -r ${ROOT_DIR}/requirements.txt
