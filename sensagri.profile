#############################################
#    Sensagri Environment Variables         #
#    To make it active, type the command:   #
#    source sensagri.profile                #
#############################################

# TODO: Fill those field automatically.


# Add path of the user local bin directory
export PATH=/datalocal1/home/arnaudl/local/bin:$PATH

# Add path of the SensagriChain installation directory
export SENSAGRICHAIN_HOME=/datalocal1/home/arnaudl/SensagriChain

echo "*** Sensagri Chain Environment Variables Loaded ***"
export LD_LIBRARY_PATH=/datalocal1/home/arnaudl/local/lib:$LD_LIBRARY_PATH

export SSOTB_HOME=${SENSAGRICHAIN_HOME}/bin/OTB-6.2.0-Linux64
export PATH=${SSOTB_HOME}/bin:$PATH
export PATH=${SENSAGRICHAIN_HOME}/chdb:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/cpp/Executables:$PATH
export PATH=${SENSAGRICHAIN_HOME}/bin/scripts:$PATH
export LD_LIBRARY_PATH=${SSOTB_HOME}/lib:$LD_LIBRARY_PATH
export ITK_AUTOLOAD_PATH=${SSOTB_HOME}/otb/applications/

export GDAL_DATA=${SENSAGRICHAIN_HOME}/bin/OTB-6.2.0-Linux64/share/gdal/
export PROJ_LIB=${SENSAGRICHAIN_HOME}/bin/OTB-6.2.0-Linux64/share/proj/

PATH=$PATH:${SENSAGRICHAIN_HOME}/bin/scripts
PATH=$PATH:${SENSAGRICHAIN_HOME}/bin/OTB/Executables

alias log="cat $1 | grep 'LOG'"
alias dbg="cat $1 | grep 'DBG'"
