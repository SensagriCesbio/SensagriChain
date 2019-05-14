//***********************************************//
//             Samples extraction                //
//                                               //
//  Extract a given number of samples randomly   //
//                                               //
//             Ludo 01/02/2017                   //
//***********************************************//
#include <iostream>
#include <stdexcept>
#include "itkListSample.h"
#include "itkVariableLengthVector.h"
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>

std::set<int> GetSamplesPositions(std::string InCropTypeLabelsFile, int SamplesNumber);
void ExtractSamples(std::vector<std::string> InName, std::vector<std::string> OutName, std::set<int> positions);
unsigned int ComputeXDim (std::string filename);

int main(int argc, char * argv[])
{

    if(argc < 3)
    {
      // Input Parameters
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<
      " NumberOfSamples input1.txt input2.txt ... inputN.txt" << std::endl;
      std::cerr << "NumberOfSamples \t Number of samples to extract randomly" << std::endl;
      std::cerr << "input2 ... inputN.txt \t List of files that should be sampling randomly. The first file should be a crop type label file." << std::endl;
    return EXIT_FAILURE;
    }
    //std::cout << "test1" << std::endl;
    std::vector<std::string> InName;
    std::vector<std::string> OutName;
    int SamplesNumber = static_cast<unsigned int>(atoll(argv[1]));
    // Handle files name, especially output file name
    for(int i = 2; i < argc; i++)
    {
      std::string name = argv[i];
      InName.push_back(name);
      OutName.push_back(name.insert(name.size()-4,"_"+std::to_string(SamplesNumber)));
    }
 
    // Obtain reshuffled samples positions according
    std::set<int> SamplesPos; 
    std::set<int>::iterator it;
    SamplesPos = GetSamplesPositions(InName[0], SamplesNumber);

    // Extract lines from input file following samples positions
    ExtractSamples(InName, OutName, SamplesPos); 
    return 0;
}

std::set<int> GetSamplesPositions(std::string InCropTypeLabelsFile, int SamplesNumber)
{
    // Set of Crop classes //
    //std::set<std::string> cropSet({"10", "12", "14", "41", "438", "4351"});
    //std::set<std::string> croplabelSet({"0", "1"});

    // Test if output file is not filled with class index.
    //std::ifstream InCropLabelsStream;
    //InCropLabelsStream.open(OutCropLabelsFile.c_str());
    //std::string InCropLabelsLine;
    //std::getline(InCropLabelsStream, InCropLabelsLine);
    
    // In Streams
    std::ifstream InCropTypeLabelsStream;
    InCropTypeLabelsStream.open(InCropTypeLabelsFile.c_str());
    if(!InCropTypeLabelsStream.is_open())
    {
        perror(InCropTypeLabelsFile.c_str());
        exit(10);
    }

    int idx = 0;
    int label = 0;
    std::map<int, std::vector<int> > IdxList;  
    std::set<int> CropSet;
    std::set<int> NewIdx;
    std::set<int>::iterator it;
    for (std::string InCropTypeLabelsLine; std::getline(InCropTypeLabelsStream, InCropTypeLabelsLine);)
    {
      idx++;
      label = std::stoi(InCropTypeLabelsLine);
      //std::cout << idx << "," << label << std::endl;
      //label = static_cast<unsigned int>(atoll(InCropTypeLabelsLine));
      if(CropSet.find(label) == CropSet.end() )
      {  
         CropSet.insert(label); 
      } 
      IdxList[label].push_back(idx);
    }

    // shuffle the selected postion
    for (it=CropSet.begin(); it!=CropSet.end(); ++it)
    {
      std::random_shuffle(IdxList[*it].begin(), IdxList[*it].end() );    
    }

    // Gather the shuffled positions in set
    for (it=CropSet.begin(); it!=CropSet.end(); ++it)
    {
      for(int i=0;i<SamplesNumber;i++)
      {
        NewIdx.insert(IdxList[*it][i]); 
      }
    }
  
    std::ofstream OutStream;
    OutStream.open("Sampling-Statistics.txt");
    if(!OutStream.is_open())
    {
      perror("Sampling-Statistics.txt");
      exit(10);
    }

    OutStream << "#Class\t\t\tCount" << std::endl;
    OutStream << "------------------------------\n";
    for (it=CropSet.begin(); it!=CropSet.end(); ++it)
    {
      OutStream << *it << "\t\t\t" << IdxList[*it].size()  << std::endl;
    }

    OutStream << std::endl;
    OutStream << "#Number of samples\t" << std::to_string(SamplesNumber) << std::endl;	
    // Clean streams and close files //
    InCropTypeLabelsStream.clear();InCropTypeLabelsStream.close();
    OutStream.clear();OutStream.close();
    return NewIdx;
}


void ExtractSamples(std::vector<std::string> InName, std::vector<std::string> OutName, std::set<int> positions)
{

    // Create all input and output streams
    
    std::ifstream InStream[InName.size()];
    for(int i = 0; i < InName.size(); i++)
    {
      InStream[i].open(InName[i].c_str());
      if(!InStream[i].is_open())
      {
        perror(InName[i].c_str());
        exit(10);
      }
    }

    std::ofstream OutStream[OutName.size()];
    for(int i = 0; i < OutName.size(); i++)
    {
      OutStream[i].open(OutName[i].c_str());
      if(!OutStream[i].is_open())
      {
        perror(OutName[i].c_str());
        exit(10);
      }
    }

    // Write selected samples in corresponding files
    int idx = 0;
    std::string InLine[InName.size()];
    for(std::string line; std::getline(InStream[0], line);)
    {
      idx++;
      InLine[0] = line;
      for(int i = 1; i < InName.size(); i++)
      {
        std::getline(InStream[i], InLine[i]);
      }

      if(positions.find(idx) != positions.end() )
      {
        for(int i = 0; i < InName.size(); i++)
        {
          OutStream[i] << InLine[i] << std::endl;
        }
      } 
   }
}

