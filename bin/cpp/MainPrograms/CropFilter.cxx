//*********************************************//
//      Isolate crop profile in streaming      //
//      Last Version: Ludo 21/03/2018          //
//                    Ludo 06/11/2017          //
//*********************************************//
#include <iostream>
#include <stdexcept>
#include "itkListSample.h"
#include "itkVariableLengthVector.h"
#include <fstream>
#include <sstream>
#include <set>


void StreamFilter(std::string InProfilesFile, std::string InCoordinatesFile, std::string InCropTypeLabelsFile, std::string InBinaryLabelsFile,std::string OutProfilesFile, std::string OutCoordinatesFile, std::string OutCropTypeLabelsFile);

int main(int argc, char * argv[])
{
    std::string InProfilesFile = argv[1];
    std::string InCoordinatesFile = argv[2];
    std::string InCropTypeLabelsFile = argv[3];
    std::string InBinaryLabelsFile = argv[4];
    std::string OutProfilesFile = argv[5];
    std::string OutCoordinatesFile = argv[6];
    std::string OutCropTypeLabelsFile = argv[7];
    
    /*
    std::cout << InProfilesFile << std::endl; 
    std::cout << InCoordinatesFile << std::endl; 
    std::cout << InCropTypeLabelsFile << std::endl; 
    std::cout << OutProfilesFile << std::endl; 
    std::cout << OutCoordinatesFile << std::endl; 
    std::cout << OutCropTypeLabelsFile << std::endl; 
    */

    // Streaming //
    StreamFilter(InProfilesFile, InCoordinatesFile, InCropTypeLabelsFile, InBinaryLabelsFile,OutProfilesFile, OutCoordinatesFile, OutCropTypeLabelsFile);
    return 0;
}

void StreamFilter(std::string InProfilesFile, std::string InCoordinatesFile, std::string InCropTypeLabelsFile, std::string InBinaryLabelsFile,std::string OutProfilesFile, std::string OutCoordinatesFile, std::string OutCropTypeLabelsFile)
{
    // In Streams
    std::ifstream InProfilesStream;
    std::ifstream InCoordinatesStream;
    std::ifstream InCropTypeLabelsStream;
    std::ifstream InBinaryLabelsStream;
   
    InProfilesStream.open(InProfilesFile.c_str());
    InCoordinatesStream.open(InCoordinatesFile.c_str());
    InCropTypeLabelsStream.open(InCropTypeLabelsFile.c_str());
    InBinaryLabelsStream.open(InBinaryLabelsFile.c_str());
  
    if(!InProfilesStream.is_open())
    {
        perror(InProfilesFile.c_str());
        exit(10);
    }
    if(!InCoordinatesStream.is_open())
    {
        perror(InCoordinatesFile.c_str());
        exit(10);
    }
    if(!InCropTypeLabelsStream.is_open())
    {
        perror(InCropTypeLabelsFile.c_str());
        exit(10);
    }
    if(!InBinaryLabelsStream.is_open())
    {
        perror(InBinaryLabelsFile.c_str());
        exit(10);
    }
  
    // Out Streams
    std::ofstream OutProfilesStream;
    std::ofstream OutCoordinatesStream;
    std::ofstream OutCropTypeLabelsStream;

    OutProfilesStream.open(OutProfilesFile.c_str());
    OutCoordinatesStream.open(OutCoordinatesFile.c_str());
    OutCropTypeLabelsStream.open(OutCropTypeLabelsFile.c_str());

    if(!OutProfilesStream.is_open())
    {
        perror(OutProfilesFile.c_str());
        exit(10);
    }
    if(!OutCoordinatesStream.is_open())
    {
        perror(OutCoordinatesFile.c_str());
        exit(10);
    }
    if(!OutCropTypeLabelsStream.is_open())
    {
        perror(OutCropTypeLabelsFile.c_str());
        exit(10);
    }
   
    std::string InCoordinatesLine;
    std::string InProfilesLine;
    std::string InBinaryLabelsLine;
    for (std::string InCropTypeLabelsLine; std::getline(InCropTypeLabelsStream, InCropTypeLabelsLine);)
    {
        // Sparsing
        std::getline(InCoordinatesStream, InCoordinatesLine);
        std::getline(InProfilesStream, InProfilesLine);
        std::getline(InBinaryLabelsStream, InBinaryLabelsLine);
        
        /*
        // Debug //
        std::cout << InCropLabelsLine << std::endl;
        std::cout << InCropTypeLabelsLine << std::endl;
        std::cout << InCoordinatesLine << std::endl;
        std::cout << InProfilesLine << std::endl;
        std::cout << "*************************" << std::endl;
        */

        // Check if a crop is pointed at // 
        if(InBinaryLabelsLine=="1")
        {   
            //std::cout << "*** CROP ***" << std::endl;
            OutProfilesStream << InProfilesLine << std::endl;
            OutCoordinatesStream << InCoordinatesLine << std::endl;
            OutCropTypeLabelsStream << InCropTypeLabelsLine << std::endl; 
            //std::string str; 
            //std::cin >> str;
            //if (str.compare("q") == 0)
            //{
            //   break;
            //}
        } 
    }

    // Clean streams and close files //
    InProfilesStream.clear();InProfilesStream.close();
    InCoordinatesStream.clear();InCoordinatesStream.close();
    InCropTypeLabelsStream.clear();InCropTypeLabelsStream.close();
    InBinaryLabelsStream.clear();InBinaryLabelsStream.close();

    OutProfilesStream.clear();OutProfilesStream.close();
    OutCoordinatesStream.clear();OutCoordinatesStream.close();
    OutCropTypeLabelsStream.clear();OutCropTypeLabelsStream.close();
}
