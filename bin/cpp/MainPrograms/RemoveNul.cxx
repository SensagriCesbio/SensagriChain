//****************************//
// Remove zero sample         //
//      Ludo 11/10/2017       //
//****************************//
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <stdio.h>

double sumString(std::string str);

int main(int argc, char * argv[])
{

    std::string filename_CropTypeLabelsA = argv[1];
    std::string filename_CropTypeLabelsB= argv[2];

    // Sparsing input files //
    std::string StringToReplace = "CropTypeLabels";
    unsigned positionA = filename_CropTypeLabelsA.find(StringToReplace);
    unsigned positionB = filename_CropTypeLabelsB.find(StringToReplace);

    std::string filename_CropLabelsA = filename_CropTypeLabelsA;
    std::string filename_CropLabelsB = filename_CropTypeLabelsB;

    std::string filename_CoordinatesA = filename_CropTypeLabelsA;
    std::string filename_CoordinatesB = filename_CropTypeLabelsB;

    std::string filename_ProfilesA = filename_CropTypeLabelsA;
    std::string filename_ProfilesB = filename_CropTypeLabelsB;
    
    filename_CropLabelsA = filename_CropLabelsA.replace(positionA, 14, "BinaryCropLabels");
    filename_CropLabelsB = filename_CropLabelsB.replace(positionB, 14, "BinaryCropLabels");
    
    filename_CoordinatesA = filename_CoordinatesA.replace(positionA, 14, "Coordinates");
    filename_CoordinatesB = filename_CoordinatesB.replace(positionB, 14, "Coordinates");

    filename_ProfilesA = filename_ProfilesA.replace(positionA, 14, "Profiles");
    filename_ProfilesB = filename_ProfilesB.replace(positionB, 14, "Profiles");

    std::string out_CropTypeLabelsA = filename_CropTypeLabelsA.substr(0,filename_CropTypeLabelsA.find_last_of('.'))+".dat";
    std::string out_CropTypeLabelsB = filename_CropTypeLabelsB.substr(0,filename_CropTypeLabelsB.find_last_of('.'))+".dat";
    
    std::string out_CropLabelsA = filename_CropLabelsA.substr(0,filename_CropLabelsA.find_last_of('.'))+".dat";
    std::string out_CropLabelsB = filename_CropLabelsB.substr(0,filename_CropLabelsB.find_last_of('.'))+".dat";

    std::string out_CoordinatesA = filename_CoordinatesA.substr(0,filename_CoordinatesA.find_last_of('.'))+".dat";
    std::string out_CoordinatesB = filename_CoordinatesB.substr(0,filename_CoordinatesB.find_last_of('.'))+".dat";
    
    std::string out_ProfilesA =  filename_ProfilesA.substr(0,filename_ProfilesA.find_last_of('.'))+".dat";
    std::string out_ProfilesB =  filename_ProfilesB.substr(0,filename_ProfilesB.find_last_of('.'))+".dat";
 
    std::cout << "************** Input files Optical *************" << std::endl; 
    std::cout << filename_CropTypeLabelsA << std::endl; 
    std::cout << filename_CropLabelsA << std::endl; 
    std::cout << filename_CoordinatesA << std::endl; 
    std::cout << filename_ProfilesA << std::endl; 
    std::cout << "************** Output files Optical************" << std::endl; 
    std::cout << out_CropTypeLabelsA << std::endl; 
    std::cout << out_CropLabelsA << std::endl; 
    std::cout << out_CoordinatesA << std::endl; 
    std::cout << out_ProfilesA << std::endl; 
    std::cout << "************************************************" << std::endl; 
    std::cout << "************** Input files Radar ***************" << std::endl; 
    std::cout << filename_CropTypeLabelsB << std::endl; 
    std::cout << filename_CropLabelsB << std::endl; 
    std::cout << filename_CoordinatesB << std::endl; 
    std::cout << filename_ProfilesB << std::endl; 
    std::cout << "************** Output files Optical************" << std::endl; 
    std::cout << out_CropTypeLabelsB << std::endl; 
    std::cout << out_CropLabelsB << std::endl; 
    std::cout << out_CoordinatesB << std::endl; 
    std::cout << out_ProfilesB << std::endl; 
    std::cout << "***********************************************" << std::endl; 
 

    // Streaming //
 
    std::ifstream CropTypeLabelsInA; std::ofstream CropTypeLabelsOutA;
    std::ifstream CropTypeLabelsInB; std::ofstream CropTypeLabelsOutB;
    
    std::ifstream CropLabelsInA; std::ofstream CropLabelsOutA;
    std::ifstream CropLabelsInB; std::ofstream CropLabelsOutB;
   
    std::ifstream CoordinatesInA; std::ofstream CoordinatesOutA;
    std::ifstream CoordinatesInB; std::ofstream CoordinatesOutB;
  
    std::ifstream ProfilesInA; std::ofstream ProfilesOutA;
    std::ifstream ProfilesInB; std::ofstream ProfilesOutB;
   
    CropTypeLabelsInA.open(filename_CropTypeLabelsA); CropTypeLabelsOutA.open(out_CropTypeLabelsA);
    CropTypeLabelsInB.open(filename_CropTypeLabelsB); CropTypeLabelsOutB.open(out_CropTypeLabelsB);

    CropLabelsInA.open(filename_CropLabelsA); CropLabelsOutA.open(out_CropLabelsA);
    CropLabelsInB.open(filename_CropLabelsB); CropLabelsOutB.open(out_CropLabelsB);

    CoordinatesInA.open(filename_CoordinatesA); CoordinatesOutA.open(out_CoordinatesA);
    CoordinatesInB.open(filename_CoordinatesB); CoordinatesOutB.open(out_CoordinatesB);
 
    ProfilesInA.open(filename_ProfilesA); ProfilesOutA.open(out_ProfilesA);
    ProfilesInB.open(filename_ProfilesB); ProfilesOutB.open(out_ProfilesB);

    std::string crop_labelsA;
    std::string crop_labelsB;

    std::string coordinatesA;
    std::string coordinatesB;
       
    std::string profilesA;
    std::string profilesB;

    std::string crop_type_labelsB; 

    int in = 0;
    int out = 0; 
    double sum = 0;
    std::cout << "Streaming inputs files" << std::endl;
    for (std::string crop_type_labelsA; std::getline(CropTypeLabelsInA, crop_type_labelsA);)
    {

       std::getline(CropTypeLabelsInB, crop_type_labelsB);
      
       std::getline(CropLabelsInA, crop_labelsA);
       std::getline(CropLabelsInB, crop_labelsB);
  
       std::getline(CoordinatesInA, coordinatesA);
       std::getline(CoordinatesInB, coordinatesB);
     
       std::getline(ProfilesInA, profilesA);
       std::getline(ProfilesInB, profilesB);
   
       sum = sumString(profilesA);
        
       if(sum > 0)
       {
           in++;
           CropTypeLabelsOutA << crop_type_labelsA << std::endl;
           CropTypeLabelsOutB << crop_type_labelsB << std::endl;
         
           CropLabelsOutA << crop_labelsA << std::endl; 
           CropLabelsOutB << crop_labelsB << std::endl; 
         
           CoordinatesOutA << coordinatesA << std::endl;
           CoordinatesOutB << coordinatesB << std::endl;
       
           ProfilesOutA << profilesA << std::endl;
           ProfilesOutB << profilesB << std::endl;
       }
       else
       {
           out++;
       }
    }
    
    std::cout << "Number of deleted profile: " << out << std::endl; 
    std::cout << "Number of kept profile: " << in << std::endl; 
   
    CropTypeLabelsInA.close(); CropTypeLabelsOutA.close();
    CropTypeLabelsInB.close(); CropTypeLabelsOutB.close();

    CropLabelsInA.close(); CropLabelsOutA.close();
    CropLabelsInB.close(); CropLabelsOutB.close();

    CoordinatesInA.close(); CoordinatesOutA.close();
    CoordinatesInB.close(); CoordinatesOutB.close();

    ProfilesInA.close(); ProfilesOutA.close();
    ProfilesInB.close(); ProfilesOutB.close();

    std::cout << "Rename files: " << std::endl;
    int result = 0;
    result = result + rename(out_CropTypeLabelsA.c_str(),filename_CropTypeLabelsA.c_str());
    result = result + rename(out_CropTypeLabelsB.c_str(),filename_CropTypeLabelsB.c_str());

    result = result + rename(out_CropLabelsA.c_str(),filename_CropLabelsA.c_str());
    result = result + rename(out_CropLabelsB.c_str(),filename_CropLabelsB.c_str());

    result = result + rename(out_CoordinatesA.c_str(),filename_CoordinatesA.c_str());
    result = result + rename(out_CoordinatesB.c_str(),filename_CoordinatesB.c_str());

    result = result + rename(out_ProfilesA.c_str(),filename_ProfilesA.c_str());
    result = result + rename(out_ProfilesB.c_str(),filename_ProfilesB.c_str());
    std::cout << "Rename files code: " << result << std::endl;

    return 0;
}

double sumString(std::string chaine)
{
    double sum = 0;
    std::string::size_type stTemp = chaine.find(",");
    while(stTemp != std::string::npos)
    {
        sum = sum + std::stod(chaine.substr(0, stTemp));
	chaine = chaine.substr(stTemp + 1);
	stTemp = chaine.find(",");
    }
    return sum;
}
