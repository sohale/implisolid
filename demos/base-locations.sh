# use with `source` only

# export USER=a9858770
export USER_HOME=/Users/$USER

# Parameters
#export BASELOC1=$USER_HOME/cs/mp5
#export BASELOC2=$USER_HOME/cs/mp5
export BASELOC1=$USER_HOME/cs
export BASELOC2=$USER_HOME/cs
export BASELOC3=$USER_HOME/cs

# $BASELOC1 -> implisolid
# $BASELOC2 -> mp5-private
# $BASELOC3 -> sohale.github.io

# /Users/$USER/cs/mp5/implisolid
# mycomputer-specific

[[ $OSTYPE == 'darwin'* ]] || "Error: MacOS-specific code"

# Base locations for implisolic, mp5-private, sohale.github.io
# Other scrips use this as the seed parameters
