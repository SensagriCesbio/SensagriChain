#############################################
#    Sensagri Environment Variables         #
#    To make it active, type the command:   #
#    source sensagri.profile                #
#############################################

echo "*** Sensagri Chain Environment Variables Loaded ***"

export SENSAGRICHAIN_HOME=/datalocal1/home/arnaudl/SensagriChain
export PATH=/datalocal1/home/arnaudl/local/bin:$PATH
export LD_LIBRARY_PATH=/datalocal1/home/arnaudl/local/lib:$LD_LIBRARY_PATH

export SSOTB_HOME=${SENSAGRICHAIN_HOME}/bin/OTB-6.2.0-Linux64
export PATH=${SSOTB_HOME}/bin:$PATH
export PATH=${SENSAGRICHAIN_HOME}/chdb:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/cpp/Executables:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/scripts:$PATH
export LD_LIBRARY_PATH=${SSOTB_HOME}/lib:$LD_LIBRARY_PATH
export ITK_AUTOLOAD_PATH=${SSOTB_HOME}/otb/applications/
 
PATH=$PATH:${SENSAGRICHAIN_HOME}/bin/scripts
PATH=$PATH:${SENSAGRICHAIN_HOME}/bin/OTB/Executables

alias log="cat $1 | grep 'LOG'"
alias dbg="cat $1 | grep 'DBG'"
