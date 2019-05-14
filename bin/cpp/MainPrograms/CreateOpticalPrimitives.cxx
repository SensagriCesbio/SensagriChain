//***********************************************//
//             Create Optical Primitive          //
//             in streaming                      //
//             Update 22/0é20169                 //
//             Ludo 14/03/2018                   //
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
#include "ReadWriteSamples.hpp"

typedef itk::VariableLengthVector<float> SampleType;
typedef itk::Statistics::ListSample<SampleType> ListSampleType;

SampleType CreatePrimitive(SampleType Profile, int NbDates, int NbPrimitives, int type);

ListSampleType::Pointer PCAMatrix = ListSampleType::New();
ListSampleType::Pointer PCAMoyList = ListSampleType::New();
SampleType PCAMoy;

//unsigned int ComputeXDim (std::string filename);

// Global

int main(int argc, char * argv[])
{
    if(argc < 5)
    {
      // Input Parameters
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<
      " ProfilesIn.txt ProfilesOut.txt NbDates type" << std::endl;
      std::cerr << "ProfilesIn.txt\tFile that contains the optical bands for all dates." << std::endl;
      std::cerr << "ProfilesOut.txt\tOutput file with the primitives added." << std::endl;
      std::cerr << "NbDates\t\tNumber of optical dates." << std::endl;
      std::cerr << "type\t\tType of primitives calculated (0: OSO, 1: Red Edge, 2: OSO and Red Edge)" << std::endl;
      return EXIT_FAILURE;
    }

    std::string matrixFile;
    std::string moyFile;
    std::string FileIn = argv[1];
    std::string FileOut = argv[2];
    int NbDates = std::stoi(argv[3]);
    int type = std::stoi(argv[4]);
    if(argc == 7)
    {
     matrixFile = argv[5];
     moyFile = argv[6];
    }
    //std::cout << FileIn << std::endl;
    //std::cout << FileOut << std::endl;
    //std::cout << NbDates<< std::endl;
    //std::cout << type<< std::endl;

    ListSampleType::Pointer ListSampleIn = ListSampleType::New();
    ListSampleType::Pointer ListSampleOut = ListSampleType::New();

    //ReadProfiles//
    SampleType Profile;
    SampleType newProfile;
    int NbBands = ComputeXDim(FileIn);
    int NbPrimitives;

    if(type == 0) NbPrimitives = 13;
    if(type == 1) NbPrimitives = 12;
    if(type == 2) NbPrimitives = 15;
    if(type == 3)
    {
      std::cout << "PCA with "<< matrixFile << std::endl;
      std::cout << "         "<< moyFile << std::endl;
      NbPrimitives = 3;
  
      // Get PCA rotation matrix
      PCAMatrix->SetMeasurementVectorSize(10*NbDates);
      ReadListSample(PCAMatrix, matrixFile);

      // Get mean PCA vector
      PCAMoyList->SetMeasurementVectorSize(10*NbDates);
      ReadListSample(PCAMoyList, moyFile);
      PCAMoy = PCAMoyList->GetMeasurementVector(0);
      //std::cout << PCAMoy << std::endl;
    }

    // Construct Sample Primitive
    // Stream loop //
    
    std::ifstream instream;
    std::ofstream outstream;
    instream.open(FileIn.c_str());
    if(!instream.is_open())
    {
        perror(FileIn.c_str());
        exit(10);
    }

    outstream.open(FileOut.c_str());
    if(!outstream.is_open())
    {
        perror(FileOut.c_str());
        exit(10);
    }

    auto XDim = ComputeXDim(FileIn);
    Profile.SetSize(XDim);
    for (std::string line; std::getline(instream , line); )
    {
      // now we'll use a stringstream to separate the fields out of the line
      std::stringstream ss( line );
      std::string field;
      auto i = 0;
      while (getline( ss, field, ',' ))
      {
        // for each field we wish to convert it to a float
        std::stringstream fs( field );
        // float singleValue = 0.0;  // (default value is 0.0)
        fs >> Profile[i];
        i = i  + 1;
      }   

      newProfile = CreatePrimitive(Profile,NbDates,NbPrimitives,type);
   
      for(auto j = 0; j < newProfile.Size() - 1; j++)
      {
        outstream << newProfile[j] << ",";
      }      
      outstream << newProfile[newProfile.Size()-1];
      outstream << std::endl;
      
    }
    instream.close();
    outstream.close();

    return 0;
}

SampleType CreatePrimitive(SampleType Profile, int NbDates, int NbPrimitives, int type)
{
  SampleType newProfile;
  SampleType uvec;
  float A;
  float B;
  float C;
  if(type==0)
  {
    // OSO case //
    int NbBands = 13;
    newProfile.SetSize(NbDates*NbPrimitives);
    int j = 0;  
    for (int i = 0; i < NbDates*10; i++)
    {   
      newProfile[j] =  static_cast<float>(Profile[i]);
      j++;
      if( (i+1)%10 == 0)
      {
        A = static_cast<float>(Profile[i - 11 + 8]); //B8
        B = static_cast<float>(Profile[i - 11 + 4]); //B4
        if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDVI
        else newProfile[j] = (A-B)/0.00000001;
        j++;

        //A = static_cast<float>(Profile[i - 11 + 3]); //B3
        //B = static_cast<float>(Profile[i - 11 + 8]); //B8
        //if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDWI
       
        A = static_cast<float>(Profile[i - 11 + 9]); //B8A
        B = static_cast<float>(Profile[i - 11 + 11]); //B11
        if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDWI SWIR
        else newProfile[j] = (A-B)/0.00000001;
        j++; 
      
        A = 0.0;
        for (int k = 0; k<10; k++) A = A + (static_cast<float>(Profile[i - k])*static_cast<float>(Profile[i - k]));
        newProfile[j] = std::sqrt(A);                 //Brightness
        j++; 
      }
    }
  }
  if(type==1)
  {
    // RE case //
    int NbBands = 12;
    newProfile.SetSize(NbDates*NbPrimitives);
    int j = 0;  
    for (int i = 0; i < NbDates*10; i++)
    {   
      newProfile[j] = static_cast<float>(Profile[i]);
      j++;
      if( (i+1)%10 == 0)
      {
        A = static_cast<float>(Profile[i - 11 + 8]); //B8
        B = static_cast<float>(Profile[i - 11 + 4]); //B4
        C = static_cast<float>(Profile[i - 11 + 3]); //B3
        if ( B != 0) newProfile[j] = (B-C)/A; //PSRI_NIR
        else newProfile[j] = (B-C)/0.00000001;
        j++;
       
        A = static_cast<float>(Profile[i - 11 + 5]); //B5 
        B = static_cast<float>(Profile[i - 11 + 8]); //B8
        if ( B != 0) newProfile[j] = A/B; //CHL_RE
        else newProfile[j] = A/0.00000001;
        j++; 
      }
    }
  }

  if(type==2)
  {
    // OSO + RE case //
    int NbBands = 15;
    newProfile.SetSize(NbDates*NbPrimitives);
    int j = 0;  
    for (int i = 0; i < NbDates*10; i++)
    {   
      newProfile[j] = static_cast<float>(Profile[i]);
      j++;
      if( (i+1)%10 == 0)
      {
        A = static_cast<float>(Profile[i - 11 + 8]); //B8
        B = static_cast<float>(Profile[i - 11 + 4]); //B4
        if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDVI
        else newProfile[j] = (A-B)/0.00000001;
        j++;
       
        //A = static_cast<float>(Profile[i - 11 + 3]); //B3
        //B = static_cast<float>(Profile[i - 11 + 8]); //B8
        A = static_cast<float>(Profile[i - 11 + 9]); //B8A
        B = static_cast<float>(Profile[i - 11 + 11]); //B11

        if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDWI SWIR
        else newProfile[j] = (A-B)/0.00000001;
        j++; 
      
        A = 0.0;
        for (int k = 0; k<10; k++) A = A + (static_cast<float>(Profile[i - k])*static_cast<float>(Profile[i - k]));
        newProfile[j] = std::sqrt(A);                 //Brightness
        j++; 

        A = static_cast<float>(Profile[i - 11 + 8]); //B8
        B = static_cast<float>(Profile[i - 11 + 4]); //B4
        C = static_cast<float>(Profile[i - 11 + 3]); //B3
        if ( B != 0) newProfile[j] = (B-C)/A; //PSRI_NIR
        else newProfile[j] = (B-C)/0.00000001;
        j++;
       
        A = static_cast<float>(Profile[i - 11 + 5]); //B5 
        B = static_cast<float>(Profile[i - 11 + 8]); //B8
        if ( B != 0) newProfile[j] = A/B; //CHL_RE
        else newProfile[j] = A/0.00000001;
        j++; 
      }
    } 
  }
  if(type==3)
  {
    // CAUTION: Experimental using fixed pca bases   
    newProfile.SetSize(NbDates*NbPrimitives);
    for(int i = 0; i < NbDates*NbPrimitives; i++)
    {
      uvec = PCAMatrix->GetMeasurementVector(i);
      newProfile[i] = 0.0;
      for(int j = 0; j < uvec.GetSize() ;j++)
      {
        newProfile[i] = newProfile[i] + (Profile[j]-PCAMoy[j])*uvec[j];
      }
    }
  }
  return newProfile;
}



