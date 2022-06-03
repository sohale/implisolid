
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR
REPO_ROOT=$(git rev-parse --show-toplevel)
cd $REPO_ROOT

ERROR_MESSAGE_SKIP_ERROR="\n ERROR:  NEEDS ATTENTION \n\nMoving on to the next test.\n"

source sandbox/sympy-experiment/p3-for-me/bin/activate

# first run the `run_script.bash`
# python -c 'import vtk'
# python -c 'import numexpr'
python -c 'import mayavi' || echo 'First run `sandbox/sympy-experiment/run_script.bash` and then source ..../p3-for-me'
python -c 'import mayavi'  # to make sure run_script is executed already, stop if eg mayavi is not there


cd $REPO_ROOT/python_implicit
# cd experimentation
# python demoImplicitObject.py
# python twisted_cube_object.py

# needs to get fixed
python implicit_test.py   || printf "$ERROR_MESSAGE_SKIP_ERROR"
python mesh_examples_for_tests.py

# pip install numpy-stl       # not this! #pip install stl
# pip install scikit-image
python stl_tests.py  || printf "$ERROR_MESSAGE_SKIP_ERROR"

python test_meshutils.py || printf "$ERROR_MESSAGE_SKIP_ERROR"
python test_optimize_dual_mesh.py || printf "$ERROR_MESSAGE_SKIP_ERROR"
python test_screw.py || printf "$ERROR_MESSAGE_SKIP_ERROR"
