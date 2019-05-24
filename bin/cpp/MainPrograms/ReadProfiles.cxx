/////////////////////////////////////////////
//      Read Profile                       //
//      Ludo 07/2017                    //
/////////////////////////////////////////////
#include "otbImage.h"
#include "otbVectorImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbStandardFilterWatcher.h"
#include "otbRAMDrivenAdaptativeStreamingManager.h"

#include "ReadWriteSamples.hpp"
#include "otbImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionSplitter.h"
#include <iostream>
#include <stdexcept>
#include <math.h>

//#include <TemporalAxisLinearInterpolator.hpp>
#include <tuple>
#include <itkVariableLengthVector.h>
#include <vector>
#include <fstream>

typedef float PixelType;
//typedef unsigned int LabelType;

using ImageType = otb::VectorImage<PixelType, 2>;
using RegionType = ImageType::RegionType;

typedef otb::RAMDrivenAdaptativeStreamingManager<ImageType> StreamingManagerType;
typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
typedef itk::ImageRegionIterator<ImageType> OutputIteratorType;
typedef itk::ImageRegionSplitter<2> ImageSplitter;
typedef itk::VariableLengthVector<PixelType> SampleType;

int main(int argc, char * argv[])
{
     
    if (argc != 3 ) 
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile pixelNumber" << std::endl;
      return EXIT_FAILURE;
    }
    // Input Parameters //
    std::string filename_InputImage = argv[1];
    int PixelNumber = std::stoi(argv[2]); 
    
    // Load input Images //
    auto reader = otb::ImageFileReader<ImageType>::New();
    const int availableRAM = 128;
    RegionType   streamingRegion;
    reader->SetFileName(filename_InputImage);
    reader->GetOutput()->UpdateOutputInformation();
    reader->GenerateOutputInformation();
   
    RegionType largestRegion = reader->GetOutput()->GetLargestPossibleRegion();
    int NbComponents = reader->GetOutput()->GetNumberOfComponentsPerPixel();
    StreamingManagerType::Pointer
    streamingManager = StreamingManagerType::New();
    streamingManager->SetAvailableRAMInMB(availableRAM);
    streamingManager->PrepareStreaming(reader->GetOutput(), largestRegion);
    const unsigned long numberOfStreamDivisions = streamingManager->GetNumberOfSplits(); 
    std::cout << "largest region " << largestRegion << std::endl;
    std::cout << "Number of features per pixel for input image = " << NbComponents << std::endl;
    std::cout << "Number of splits = " << numberOfStreamDivisions << std::endl;

    itk::VariableLengthVector<PixelType> pixel;
    //std::string str;
    std::cout << "*** Image Reading ***" << std::endl;
    int counter = 0;
    for (auto piece = 0;piece < numberOfStreamDivisions; piece++)
    {
        std::cout << "Piece = " << piece << std::endl;
        streamingRegion = streamingManager->GetSplit(piece);
        reader->GetOutput()->SetRequestedRegion(streamingRegion);
        reader->GetOutput()->PropagateRequestedRegion();
        reader->GetOutput()->UpdateOutputData();
        IteratorType it(reader->GetOutput(), streamingRegion);

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {   
            pixel = it.Get();
            if(counter == PixelNumber) // Set test to true (with a condition if needed) to see Debug information //
            {
              // Print if non zero
              std::cout << "Vector Size: " << pixel.Size() << std::endl;
              std::cout << "Vector: " << std::endl;
              std::cout << it.Get() << std::endl;
              //std::cin >> str;
              //if (str.compare("q") == 0) return 0;
              return 0;
            }
            counter ++;
        }
     }      

   return 0;
}
