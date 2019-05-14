//*********************************************//
// Join Radar and Optical samples in streaming //
//      Ludo 06/11/2017                        //
//*********************************************//
#include <iostream>
#include <stdexcept>
#include "itkListSample.h"
#include "itkVariableLengthVector.h"
#include <fstream>
#include <sstream>



void StreamJoin(std::string RadarFile, std::string OpticalFile, std::string RadarOpticalFile, std::string DateOrderFile, int RadarChannelNbr, int OpticalChannelNbr);
//unsigned int ComputeXDim (std::string OpticalFile, std::string Radarfile, RadarOpticalFile DateOrderFile);
unsigned int ComputeXDim (std::string filename);

int main(int argc, char * argv[])
{
    std::string RadarFile = argv[1];
    std::string OpticalFile = argv[2];
    std::string RadarOpticalFile = argv[3];
    std::string DateOrderFile = argv[4];
    int RadarChannelNbr = static_cast<int>(atoll(argv[5]));
    int OpticalChannelNbr = static_cast<int>(atoll(argv[6]));

    std::cout << RadarFile << std::endl; 
    std::cout << OpticalFile << std::endl; 
    std::cout << RadarOpticalFile << std::endl; 
    std::cout << DateOrderFile << std::endl; 
    
    // Streaming //
    StreamJoin(RadarFile, OpticalFile, RadarOpticalFile, DateOrderFile, RadarChannelNbr, OpticalChannelNbr);
    return 0;
}

void StreamJoin(std::string RadarFile, std::string OpticalFile, std::string RadarOpticalFile, std::string DateOrderFile, int RadarChannelNbr, int OpticalChannelNbr)
{
    std::ifstream RadarStream;
    std::ifstream OpticalStream;
    std::ifstream DateOrderStream;
    std::ofstream RadarOpticalStream;

    RadarStream.open(RadarFile.c_str());
    OpticalStream.open(OpticalFile.c_str());
    DateOrderStream.open(DateOrderFile.c_str());
    RadarOpticalStream.open(RadarOpticalFile.c_str());

    if(!RadarStream.is_open())
    {
        perror(RadarFile.c_str());
        exit(0);
    }
    if(!OpticalStream.is_open())
    {
        perror(OpticalFile.c_str());
        exit(0);
    }
    if(!DateOrderStream.is_open())
    {
        perror(DateOrderFile.c_str());
        exit(0);
    }
    if(!RadarOpticalStream.is_open())
    {
        perror(RadarOpticalFile.c_str());
        exit(0);
    }  
    
    auto RadarDim = ComputeXDim(RadarFile);
    auto OpticalDim = ComputeXDim(OpticalFile);
    auto counter = 0;

    // Vector where the sample is store 
    itk::VariableLengthVector<float> RadarSample;
    itk::VariableLengthVector<float> OpticalSample;
    itk::VariableLengthVector<int> DateOrder;
    RadarSample.SetSize(RadarDim);
    OpticalSample.SetSize(OpticalDim);
    DateOrder.SetSize(RadarDim + OpticalDim);
    // Get Date Order data consistant with number of channel//
    std::string line;
    int dateidx;
    int ChannelNbr;
    while (getline(DateOrderStream,line))
    {
      dateidx = std::stoi(line);
      if(dateidx == 1) ChannelNbr = RadarChannelNbr; 
      else if(dateidx == 2) ChannelNbr = OpticalChannelNbr; 
      for(int i = 0; i < ChannelNbr; i++)
      {
        DateOrder[counter] = dateidx;
        counter++;
      }
    }
    DateOrderStream.clear();
    DateOrderStream.close();

    //std::cout << "dateData size " << DateOrder.Size() << std::endl;
    //std::cout << "dateData = " << DateOrder << std::endl;

    int rc = 0;
    int oc = 0;
    std::string field;
    std::string opticalline;
    for (std::string radarline; std::getline(RadarStream, radarline); ) {
        // Sparsing
        std::getline(OpticalStream , opticalline);
        std::stringstream radarss(radarline);
        std::stringstream opticalss(opticalline);
        counter = 0;
        // Get Radar Sample
        while (getline( radarss, field, ',' ))
        {
            std::stringstream fs(field);
            fs >> RadarSample[counter];
            counter = counter + 1;
        }
        // Get Optical Sample
        counter = 0;
        while (getline( opticalss, field, ',' ))
        {
            std::stringstream fs(field);
            fs >> OpticalSample[counter];
            counter = counter + 1;
        }
        rc = 0;
        oc = 0;
        // Creation Joined Radar/Optical Sample then write in file //
        for(int i = 0; i < (RadarDim + OpticalDim); i++)
        {
          if(DateOrder[i] == 1)
          {
            RadarOpticalStream << RadarSample[rc];
            if((rc+oc) != (RadarDim + OpticalDim) - 1) RadarOpticalStream << ",";
            rc++;
          } 
          else if(DateOrder[i] == 2)
          {
            RadarOpticalStream << OpticalSample[oc];
            if((rc+oc) != (RadarDim + OpticalDim) - 1) RadarOpticalStream << ",";
            oc++;
          }
        }
        RadarOpticalStream << std::endl;

        //std::cout << "Radar Sample size:\t" << RadarDim << std::endl;
        //std::cout << RadarSample  << std::endl;
        //std::cout << "Optical Sample size:\t" << OpticalDim << std::endl;
        //std::cout << OpticalSample  << std::endl;

        // Debug //
        /* 
        std::string str; 
        std::cin >> str;
        if (str.compare("q") == 0)
        {
           RadarStream.close();
           OpticalStream.close();
           RadarOpticalStream.close();
           break;
        }
        */

    }
    RadarStream.clear();
    OpticalStream.clear();
    RadarOpticalStream.clear();
    RadarStream.close();
    OpticalStream.close();
    RadarOpticalStream.close();
}

unsigned int ComputeXDim (std::string filename)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.is_open())
    {
        perror(filename.c_str());
        exit(10);
    }        
    // computing the number of columns of the txt file
    auto XDim = 0;
    std::string lineCompteur;
    std::getline(myfile , lineCompteur);
    std::stringstream ssCompteur(lineCompteur);
    std::string fieldCompteur;
    while (getline( ssCompteur, fieldCompteur, ',' ))
    {
        XDim = XDim + 1;
    }
    // std::cout << XDim <<std::endl;
    myfile.clear();
    myfile.close();
    return XDim;    
}
