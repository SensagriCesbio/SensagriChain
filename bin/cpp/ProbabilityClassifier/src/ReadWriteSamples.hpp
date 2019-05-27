/*
 * Outils.hpp
 *
 *  Created on: 18/3/2015
 *      Author: silvia
 */

#ifndef OUTILS_READ_WRITE_SAMPLES_HPP_
#define OUTILS_READ_WRITE_SAMPLES_HPP_

#include "itkListSample.h"
#include "itkVariableLengthVector.h"
#include <fstream>
#include <sstream>


void writeListSample (itk::Statistics::ListSample<itk::VariableLengthVector<double>>::Pointer ListSample, std::string filename)
{
   std::ofstream myfile;
    myfile.open (filename);
    myfile.precision(6);
   
    
    for ( auto i = 0; i < ListSample->Size(); ++i )
    {
         itk::VariableLengthVector<double> Sample = ListSample->GetMeasurementVector(i);
        
        for(auto j = 0; j < Sample.Size() -1 ; j++)
        {
            myfile<<Sample[j] << ",";
        }
        
        myfile<<Sample[Sample.Size()-1];
        
        myfile<< std::endl;
    }
    myfile.close();
    
    
}



void readListLabel (itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer ListLabel, std::string filename)
{
    std::ifstream myfile;
    myfile.open (filename);
    
    
    
    for (std::string line; std::getline(myfile , line); )
    {
        
       
        itk::FixedArray<int, 1> Sample;
        
        // now we'll use a stringstream to separate the fields out of the line
        std::stringstream ss( line );
        std::string field;
        
        auto i = 0;
        while (getline( ss, field, ',' ))
        {
            // for each field we wish to convert it to a float
            std::stringstream fs( field );
            // float singleValue = 0.0;  // (default value is 0.0)
            fs >> Sample[i];      
            i = i  + 1;
        }
	     
        ListLabel->PushBack(Sample);
    }
    

    
    
   }


void writeListLabel (itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer ListLabel, std::string filename)
{
    std::ofstream myfile;
    myfile.open (filename);
    myfile.precision(3);
    
    for ( auto i = 0; i < ListLabel->Size(); ++i )
    {
        itk::FixedArray<int, 1> Sample = ListLabel->GetMeasurementVector(i);
        myfile<<Sample[0];
        myfile<< std::endl;
    }
    
    myfile.close();
    
    
}


void writeListCoord (itk::Statistics::ListSample<itk::FixedArray<unsigned int, 2>>::Pointer ListLabel, std::string filename)
{
    std::ofstream myfile;
    myfile.open (filename);
    myfile.precision(3);
    
    for ( auto i = 0; i < ListLabel->Size(); ++i )
    {
        itk::FixedArray<int, 2> Sample = ListLabel->GetMeasurementVector(i);
        myfile<<Sample[0]<< ",";
 	myfile<<Sample[1];
        myfile<< std::endl;
    }
    
    myfile.close();
    
    
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
        XDim = XDim  +1;
    }
    // std::cout << XDim <<std::endl;
    myfile.clear();
    myfile.close();
    
    return XDim;
    
}


void ReadListSample (itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer ListSample, std::string filename)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.is_open())
    {
        perror(filename.c_str());
        exit(10);
    }
    
    auto XDim = ComputeXDim (filename);
    
    for (std::string line; std::getline(myfile , line); ) {
        
        //std::cout << line <<std::endl;
        
        itk::VariableLengthVector<float> Sample;
        Sample.SetSize(XDim);
        
        // now we'll use a stringstream to separate the fields out of the line
        std::stringstream ss( line );
        std::string field;
        
        auto i = 0;
        while (getline( ss, field, ',' ))
        {
            // for each field we wish to convert it to a float
            std::stringstream fs( field );
            // float singleValue = 0.0;  // (default value is 0.0)
            fs >> Sample[i];
            
            i = i  + 1;
        }

        ListSample->PushBack(Sample);
    }
    
    myfile.close();
}


void ReadListSampleInstantT(itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer ListSample, std::string filename, int InstantT)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.is_open())
    {
        perror(filename.c_str());
        exit(10);
    }
    
    auto XDim = ComputeXDim (filename);
    
    for (std::string line; std::getline(myfile , line); ) {
        
        //std::cout << line <<std::endl;
        
        itk::VariableLengthVector<float> Sample;
        Sample.SetSize(InstantT);
        
        // now we'll use a stringstream to separate the fields out of the line
        std::stringstream ss( line );
        std::string field;
        
        auto i = 0;
        while (getline( ss, field, ',' ) && i < InstantT )
        {
            // for each field we wish to convert it to a float
            std::stringstream fs( field );
            // float singleValue = 0.0;  // (default value is 0.0)
            fs >> Sample[i];
            
            i = i  + 1;
        }

        ListSample->PushBack(Sample);
    }
    
    myfile.close();
}

void ReadListSampleUnsignedInt (itk::Statistics::ListSample<itk::VariableLengthVector<unsigned int>>::Pointer ListSample, std::string filename)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.is_open())
    {
        perror(filename.c_str());
        exit(10);
    }
    
    auto XDim = ComputeXDim (filename);
    
    for (std::string line; std::getline(myfile , line); ) {
        
        //std::cout << line <<std::endl;
        
        itk::VariableLengthVector<unsigned int> Sample;
        Sample.SetSize(XDim);
        
        // now we'll use a stringstream to separate the fields out of the line
        std::stringstream ss( line );
        std::string field;
        
        auto i = 0;
        while (getline( ss, field, ',' ))
        {
            // for each field we wish to convert it to a float
            std::stringstream fs( field );
            // float singleValue = 0.0;  // (default value is 0.0)
            fs >> Sample[i];
            
            i = i  + 1;
        }
        ListSample->PushBack(Sample);
    }
    
    myfile.close();
}

void ReadMaskListSample (itk::Statistics::ListSample<itk::VariableLengthVector<int>>::Pointer ListSample, std::string filename)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.is_open())
    {
        perror(filename.c_str());
        exit(10);
    }
    
    auto XDim = ComputeXDim (filename);
    
    for (std::string line; std::getline(myfile , line); ) {
        
        //std::cout << line <<std::endl;
        
        itk::VariableLengthVector<int> Sample;
        Sample.SetSize(XDim);
        
        // now we'll use a stringstream to separate the fields out of the line
        std::stringstream ss( line );
        std::string field;
        
        auto i = 0;
        while (getline( ss, field, ',' ))
        {
            // for each field we wish to convert it to a float
            std::stringstream fs( field );
            fs >> Sample[i];
            
            i = i  + 1;
        }
        ListSample->PushBack(Sample);
    }
    
    myfile.close();
}



#endif /* READ_WRITE_SAMPLES_HPP_ */
