#############################################
#    Sensagri Environment Variables         #
#    To make it active, type the command:   #
#    source sensagri.profile                #
#############################################

echo "*** Sensagri Chain Environment Variables Loaded ***"

export PATH=~/local/bin:$PATH
export LD_LIBRARY_PATH=~/local/lib:$LD_LIBRARY_PATH

export SENSAGRICHAIN_HOME=~/SensagriChain
export SSOTB_HOME=${SENSAGRICHAIN_HOME}/bin/OTB-6.2.0-Linux64
export PATH=${SSOTB_HOME}/bin:$PATH
export PATH=${SENSAGRICHAIN_HOME}/chdb:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/cpp/Executables:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/scripts:$PATH
export LD_LIBRARY_PATH=${SSOTB_HOME}/lib:$LD_LIBRARY_PATH
export ITK_AUTOLOAD_PATH=${SSOTB_HOME}/otb/applications/
 
PATH=$PATH:~/SenSAgriChain/bin/scripts
PATH=$PATH:~/SenSAgriChain/bin/OTB/Executables

alias log="cat $1 | grep 'LOG'"
alias dbg="cat $1 | grep 'DBG'"
