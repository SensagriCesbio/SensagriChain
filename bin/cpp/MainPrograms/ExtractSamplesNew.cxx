/*
 * main_Hello.cxx
 *
 *  Created on: 02/02/2018
 *  Modification: Ludo
 *
 */


#include "otbVectorImage.h"
#include "otbVectorData.h"

#include "otbVectorDataFileReader.h"
#include "otbImageFileReader.h"

#include <iostream>
#include <stdlib.h>

#include "ExtractSampleFromVectorFile.hpp"
#include "ReadWriteSamples.hpp"

using PixelType = float;
using ImageType = otb::VectorImage<PixelType, 2>;
using VectorDataType = otb::VectorData<unsigned int,2>;


typedef  unsigned int ClassLabelType;
typedef typename ImageType::PixelType SampleType;
typedef itk::Statistics::ListSample<SampleType> ListSampleType;


typedef itk::FixedArray<ClassLabelType, 1> CropTypeLabel; 
typedef itk::Statistics::ListSample<CropTypeLabel> ListCropLabelType;

typedef itk::FixedArray<ClassLabelType, 1> BinaryCropLabel; 
typedef itk::Statistics::ListSample<BinaryCropLabel> ListBinaryCropLabelType;

typedef itk::FixedArray<unsigned int, 2> CoordinateType; 
typedef itk::Statistics::ListSample<CoordinateType> CoordinateList;


int main(int argc, char * argv[])
{
    
    std::string ShpFileName  = argv[1];
    std::string image_name= argv[2];
    std::string path = argv[3];
    std::string type;   
 
    if(ShpFileName.find("learn") != std::string::npos)
    {
      type="_learn";
    }
    else if(ShpFileName.find("val") != std::string::npos)
    {
      type="_val";
    }
    else
    {
      std::cout << "No learn or val in shapefile found." << '\n';
      return EXIT_FAILURE;
    } 

    if(argc == 5)
    {
      type = type + "_" + argv[4];
    }


    std::string Outfilename_Samples = path + "Profiles" + type  + ".txt";
    std::string Outfilename_CropTypeLabels = path + "CropTypeLabels"+type+".txt";
    std::string Outfilename_BinaryCropLabels = path + "BinaryCropLabels"+type+".txt";
    std::string Outfilename_Coordinates = path + "Coordinates"+type+".txt";
 
    //std::cout << Outfilename_Samples << std::endl;
    //std::cout << Outfilename_CropTypeLabels << std::endl;
    //std::cout << Outfilename_BinaryCropLabels << std::endl;
    //std::cout << Outfilename_Coordinates << std::endl;

    auto readerImage = otb::ImageFileReader<ImageType>::New();
    readerImage->SetFileName(image_name);
    readerImage->GetOutput()->UpdateOutputInformation();
    std::cout << "Image defined"  << std::endl;
    
    auto vectorReader =  otb::VectorDataFileReader<VectorDataType>::New();
    vectorReader->SetFileName( ShpFileName);
    vectorReader->Update( );
    std::cout << "Shapefile defined"  << std::endl;
    
    ExtractSampleFromVectorFile<ImageType, VectorDataType>::Pointer  SampleGenerator;
    SampleGenerator = ExtractSampleFromVectorFile<ImageType, VectorDataType>::New();
    
    SampleGenerator->SetInputImage(readerImage->GetOutput());
    SampleGenerator->SetInputVectorData(vectorReader->GetOutput());
    std::cout << "Sample generator defined"  << std::endl;
    
    //std::cout<<" Generate List  "  << std::endl;
    SampleGenerator->GenerateLists( );

    auto NbClasses = SampleGenerator->GetNumberOfClasses();
    
    std::cout<<" The NbClasses is  " << NbClasses << std::endl;
    
    SampleGenerator->PrintInformation( );
    
    
    ListSampleType::Pointer ListSample   = SampleGenerator->GetListSample();
    writeListSample (ListSample,  Outfilename_Samples);

    ListCropLabelType::Pointer CropTypeListSample   = SampleGenerator->GetCropTypeListLabel();
    writeListLabel (CropTypeListSample, Outfilename_CropTypeLabels);

    ListBinaryCropLabelType::Pointer BinaryListSample   = SampleGenerator->GetBinaryCropTypeListLabel(); 
    writeListLabel (BinaryListSample,  Outfilename_BinaryCropLabels);

    CoordinateList::Pointer CoordListSample = SampleGenerator->GetCoordinatesList();
    writeListCoord(CoordListSample,Outfilename_Coordinates);
    
    
    return 0;
}
