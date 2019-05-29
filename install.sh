#!/bin/bash
############################################################
### SenSAgri Classification Chain Prototype Installation ###
###                CESBIO-TOULOUSE (2018)                ###
###                Update 29/05.2017                     ### 
############################################################

echo "************************************************************"
echo "*** SenSAgri Classification Chain Prototype Installation ***"
echo "***          v1.0 CESBIO-TOULOUSE (2018-2019)            ***"
echo "************************************************************"

echo ""
echo "*** INSTALLATION CHDB ***"
tar zxvf chdb.tar.gz
cd chdb
make
ln -s chdb.exe chdb
cd ..

#echo ""
#echo "*** INSTALLATION OTB 6.2.0 ***"
#echo "*** This installation will not be in conflict with existing OTB installation ***"
#
#chmod 744 OTB-6.2.0-Linux64.run
#./OTB-6.2.0-Linux64.run >> install.log
#
#echo ""
#echo "*** OTB PATCHES ***"
#cp patches/otbRequiresOpenCVCheck.h OTB-6.2.0-Linux64/include/OTB-6.2/otbRequiresOpenCVCheck.h
#cp patches/otbOpenCVUtils.h OTB-6.2.0-Linux64/include/OTB-6.2/otbOpenCVUtils.h
#cp patches/otbMachineLearningModelFactory.txx OTB-6.2.0-Linux64/include/OTB-6.2/otbMachineLearningModelFactory.txx
#
#echo ""
#echo "*** CHAIN COMPILATION ***"
#cd cpp/build
#rm -r *
#cmake ..
#make
#cd ../CreateProbabilityMap/build
#rm -r *
#cmake ../src/
#make
#echo "# SENSAGRI CHAIN ENVIRONMENT VARIABLES #" >> ~/.bashrc
#echo "export SENSAGRICHAIN_HOME=~/SenSAgriChain" >> ~/.bashrc
#echo "export SSOTB_HOME=\$SENSAGRICHAIN_HOME/bin/OTB-6.2.0-Linux64" >> ~/.bashrc
#echo "export PATH=\${SSOTB_HOME}/bin:\$PATH" >> ~/.bashrc
#echo "export PATH=\${SENSAGRICHAIN_HOME}/bin/chdb-master:\$PATH" >> ~/.bashrc
#echo "export PATH=\${SENSAGRICHAIN_HOME}/bin/cpp/Executables:\$PATH" >> ~/.bashrc
#echo "export PATH=\${SENSAGRICHAIN_HOME}/bin/scripts:\$PATH" >> ~/.bashrc
#echo "export LD_LIBRARY_PATH=\${SSOTB_HOME}/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
#
#echo ""
#echo "*** Installation Almost Done ***"
#echo "To terminate the installation process,"
#echo "please source the user .bashrc file "
#echo "with the command:"
#echo "  source ~/.bashrc"
#echo " " >> ~/.bashrc 
#
