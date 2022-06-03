#!/usr/bin/env bash
# Forked from: point-process-simple-example/run_script.bash: https://github.com/sohale/point-process-simple-example/blob/82a62d013d909f365a391aa254dc598d62a0c2d4/run_script.bash
# Forked from https://github.com/sosi-org/scientific-code/blob/main/timescales-state/run-script.bash
# which was in turn forked from https://github.com/sosi-org/primsoup/blob/master/actn/run-actn.bash

set -xu

cd /Users/sohail/cs/implisolid/sandbox/sympy-experiment

mkdir -p temp
source ./temp/my-bash-utils.sh || curl -k \
    https://raw.githubusercontent.com/sohale/implisolid/revival-sohale/scripts/bash-utils.sh \
    >./temp/my-bash-utils.sh

source ./temp/my-bash-utils.sh

set -e

# if behind a firewall
export PIPFLAGS="\
    --trusted-host pypi.python.org \
    --trusted-host files.pythonhosted.org \
    --trusted-host pypi.org"

# not behind a firewall
export PIPFLAGS=""

echo "PIPFLAGS>>>> $PIPFLAGS"


function chk_venv(){
    # a solution based on `venv` as opposed to `virutalenv`

    #set -ex
    if [[  -d ./p3-for-me ]]
    then
    echo "venv exists"
    return 0
    fi

    echo "INSTALLING THEM"
    rm -rf p3-for-me || :

    # venv is shipped with python3
    #python3 -m venv -v --python=python3 p3-for-me
    python3 -m venv p3-for-me
    source ./p3-for-me/bin/activate

    python --version
    # Python 3.9.12

    # For trusted sources: see  https://stackoverflow.com/questions/49324802/pip-always-fails-ssl-verification

    python -m \
        pip install \
            $PIPFLAGS \
            --upgrade pip

}

# to refresh: `rm -rf ./p3-for-me`
MAKE_HAPPEN "./p3-for-me/bin/activate" || {
    chk_venv
}


source ./p3-for-me/bin/activate
# export PYTHON_VER="python3.9"
export PYTHON_VER="$(ls -1t ./p3-for-me/lib/|grep -i "python"|head -n 1)"

export VENV_PACKAGES="./p3-for-me/lib/$PYTHON_VER/site-packages"

MAKE_HAPPEN "$VENV_PACKAGES/numpy/LICENSE.txt" || {
    pip install \
        $PIPFLAGS \
        numpy

    pip install \
        $PIPFLAGS \
        matplotlib
}

MAKE_HAPPEN "$VENV_PACKAGES/scipy/LICENSE.txt" || {
  pip install $PIPFLAGS scipy
}
MAKE_HAPPEN "$VENV_PACKAGES/sympy/__init__.py" || {
  pip install $PIPFLAGS sympy
}
MAKE_HAPPEN "$VENV_PACKAGES/yaml/__init__.py" || {
  pip install $PIPFLAGS PyYAML
}
#MAKE_HAPPEN "$VENV_PACKAGES/pdb/__init__.py" || {
#  pip install $PIPFLAGS pdb
#}

# python -m pip install -U autopep8

MAKE_HAPPEN "$VENV_PACKAGES/graphviz/__init__.py" || {
  pip install $PIPFLAGS graphviz
}

######################################################
# Attempts to run mayavi-based script
######################################################

echo > /dev/null '
pip install numexpr
brew install vtk  # installs vtk@9.1
pip install vtk
'

pip install wheel
pip install ipdb
pip install numexpr
pip install vtk
pip install mayavi
pip install pyqt6 # ?
# pip install traitsui  # no need, already installed


#problem:
# -->
#ImportError: Could not import backend for traitsui.  Make sure you
#        have a suitable UI toolkit like PyQt/PySide or wxPython
#        installed.
#

# with:
#export ETS_TOOLKIT=pyqt6
#export QT_API=pyqt6
#---->
# problem:
# RuntimeError: No pyface.toolkits plugin found for toolkit qt6

# PySide2 does not work

# with:
#   export ETS_TOOLKIT=
#   export QT_API=pyqt6
# --->
#     raise ve_exc
# ValueError: 65534 is not a valid QEvent.Type

# suggested: online
#   export ETS_TOOLKIT=qt4
#   export QT_API=pyqt5
#
# suggested to yuse pyqt4: https://stackoverflow.com/questions/44501987/using-mayavi-on-macos-with-pyqt5


######################################################
# End of attempts to run mayavi-based script
######################################################



###############################################
# Attempts to install a different backend for matplotlib (not for implisolid)
###############################################
#brew install pkg-config
#brew link pkg-config
#brew install pygtk
#brew install freetype
#brew install libpng

true || MAKE_HAPPEN "$VENV_PACKAGES/mplcairo/__init__.py" || {
    # brew install llvm
    brew info llvm # keg-only

    export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
    export PATH="/opt/homebrew/opt/llvm/bin:$PATH"

    # https://github.com/matplotlib/mplcairo#macos
    # pip install mplcairo  # from PyPI
    # The wheel package is not available.

    # pip install git+https://github.com/matplotlib/mplcairo

echo '
==> llvm
To use the bundled libc++ please add the following LDFLAGS:
  LDFLAGS="-L/opt/homebrew/opt/llvm/lib -Wl,-rpath,/opt/homebrew/opt/llvm/lib"

llvm is keg-only, which means it was not symlinked into /opt/homebrew,
because macOS already provides this software and installing another version in
parallel can cause all kinds of trouble.

If you need to have llvm first in your PATH, run:
  echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc

For compilers to find llvm you may need to set:
  export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
  export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
'

}

#MAKE_HAPPEN "$VENV_PACKAGES/mpl_interactions/__init__.py" || {
#  pip install mpl_interactions
#}

# for llvm (failed attempt)
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"

######################################################
# End of attempts for new backends for matplotlib
######################################################


echo "Main script"

source ./p3-for-me/bin/activate

python --version

echo '
source ./p3-for-me/bin/activate
export ETS_TOOLKIT=qt5
# export ETS_TOOLKIT=pyqt5   ?

cd ../../python_implicit
python demoImplicitObject.py
'

cd ../../python_implicit
python demoImplicitObject.py


# failed

