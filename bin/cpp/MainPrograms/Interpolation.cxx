/////////////////////////////////////////////
//         Ludo 07/2017                    //
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

#include <TemporalAxisLinearInterpolator.hpp>
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

// Prototypes //
std::tuple<std::vector<unsigned int>,std::vector<unsigned int>,int,int> LoadDates(std::string infilename, std::string outfilename);
int dateconversion(int date, int referenceYear);

int main(int argc, char * argv[])
{
     
    if (argc != 6 ) 
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile inputDatesFile outputImageFile targetDatesFile NbPrimitives" << std::endl;
      return EXIT_FAILURE;
    }
    // Input Paramters //
    std::string filename_InputImage = argv[1];
    std::string filename_InputDates = argv[2]; 
    std::string filename_OutputImage = argv[3];
    std::string filename_OutputDates = argv[4];
    int NbPrim = static_cast<int>(atoll(argv[5]));
    
   // std::cout << " Test dates." << std::endl;
   // std::cout << "20151104 -> " << dateconversion(20151104, 16)<< std::endl;
   // std::cout << "20151231 -> " << dateconversion(20151231, 16)<< std::endl;
   // std::cout << "20160101 -> " << dateconversion(20160101, 16)<< std::endl;
   // std::cout << "20160201 -> " << dateconversion(20160201, 16)<< std::endl;
    

    std::cout << "*** Dates conversions ***" << std::endl;
    std::tuple<std::vector<unsigned int>,std::vector<unsigned int>,int,int> dates;
    dates = LoadDates(filename_InputDates,filename_OutputDates);
    std::vector<unsigned int> indates; indates = std::get<0>(dates);
    std::vector<unsigned int> outdates; outdates = std::get<1>(dates);
    int NbInDates = std::get<2>(dates);
    int NbOutDates = std::get<3>(dates);
    std::cout << "There is " << NbInDates << " input dates." << std::endl;
    std::cout << "There is " << NbOutDates << " output dates." << std::endl;
    
    //std::cout << indates[0] << std::endl;
    //std::cout << indates[1] << std::endl; 
    //std::cout << indates[2] << std::endl; 
    //std::cout << outdates[0] << std::endl;
    //std::cout << outdates[1] << std::endl; 
    //std::cout << outdates[3] << std::endl;
    
 
    // mask creation
    itk::VariableLengthVector<int> mask;
    mask.SetSize(indates.size());
    mask.Fill(0);
    std::cout << "Mask =  " << mask << std::endl;

    // Load input Images //
    auto reader = otb::ImageFileReader<ImageType>::New();
    const int availableRAM = 60;
    //const int availableRAM = 32;
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
    std::cout << "Largest region " << largestRegion << std::endl;
    std::cout << "Number of features per pixel for input image = " << NbComponents << std::endl;
    std::cout << "Number of splits = " << numberOfStreamDivisions << std::endl;

    // Set output image to input image size //
    ImageType::Pointer OutputImage = ImageType::New();
    OutputImage->SetRegions(largestRegion);
  
    OutputImage->CopyInformation(reader->GetOutput());
    //OutputImage->SetNumberOfComponentsPerPixel(10*outdates.size()); // CORRECT ONLY FOR S2
    OutputImage->SetNumberOfComponentsPerPixel(NbPrim*outdates.size()); // CORRECT ONLY FOR S2
    OutputImage->Allocate();
    //std::cout << "DEBUG"<< std::endl;

    
    itk::VariableLengthVector<PixelType> pixel;
    itk::VariableLengthVector<float> pixelIn_b;
    itk::VariableLengthVector<float> pixelOut_b;
    itk::VariableLengthVector<PixelType> pixelOut;
    
    pixelIn_b.SetSize(NbInDates);
    pixelOut_b.SetSize(NbOutDates);
    pixelOut.SetSize(NbPrim*NbOutDates);

    std::cout << "*** Image Construction ***" << std::endl;
    for (auto piece = 0;piece < numberOfStreamDivisions; piece++)
    {  
        if(piece%100 == 0) std::cout << "Piece = " << piece << std::endl;
        streamingRegion = streamingManager->GetSplit(piece);
        //std::cout << "\nProcessing region Image: " << streamingRegion << std::endl;
        reader->GetOutput()->SetRequestedRegion(streamingRegion);
        reader->GetOutput()->PropagateRequestedRegion();
        reader->GetOutput()->UpdateOutputData();
        IteratorType it(reader->GetOutput(), streamingRegion);

        OutputIteratorType itOut(OutputImage, streamingRegion);

        itOut.GoToBegin();
        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {   
            std::string str; 
            //itOut.Set(interpolation(it,indates,outdates,mask));
            pixel = it.Get();
            //std::cout << "pixel = " << pixel << std::endl;
            for(int b = 0; b<NbPrim; b++)
            {
              // Construct input pixel per band//
              for(int t=0;t<NbInDates;t++)
              {
                pixelIn_b[t] = float(pixel[(NbPrim)*t+b]);
              }
              //std::cout << "pixelIn_b = " << pixelIn_b << std::endl;
              // Interpolate pixel per band  
              pixelOut_b = DoYLinearInterpolate(pixelIn_b, mask, indates, outdates); // Inter S2
              //std::cout << "pixelOut_b = " << pixelOut_b << std::endl;

              // Construct output pixel from all interpolated bands
              for(int t=0;t<NbOutDates;t++)
              {
                pixelOut[(NbPrim)*t+b] = static_cast<PixelType>(pixelOut_b[t]);
              }
              //std::cout << "pixelOut = " << pixelOut << std::endl;
              //std::cout << std::endl;
            }

            itOut.Set(pixelOut); // Inter. Orbite 0
            if(false) // Set test to true (with a condition if needed) to see Debug information //
            {
              // Print if non zero
              std::cout << it.Get() << std::endl;
              std::cout << std::endl;
              std::cout << itOut.Get() << std::endl;  
              std::cin >> str;
              if (str.compare("q") == 0) return 0;
            }
            ++itOut;
        }
     }      

   // Load Output Image //
   //otb::StandardFilterWatcher watcher(writer, "Writing Classification Map...");
   std::cout << "*** Image Creation ***"<< std::endl;
   auto writer = otb::ImageFileWriter<ImageType>::New();
   writer->SetFileName(filename_OutputImage);
   writer->SetInput(OutputImage);
   try 
   { 
       writer->Update(); 
   } 
   catch (itk::ExceptionObject& err) 
   { 
       std::cerr << "ExceptionObject caught !" << std::endl; 
       std::cerr << err << std::endl; 
       return EXIT_FAILURE; 
   }

   return 0;
}

std::tuple<std::vector<unsigned int>,std::vector<unsigned int>,int,int> LoadDates(std::string infilename, std::string outfilename)
{
    std::ifstream file;
    file.open(infilename);    
    std::string line;
    int sizex = 0;
    int first;
    int date;
    // Get infile size and first date
    while (getline(file,line))
    {
      if (sizex == 0)
      {
          //date = std::stoi(line.substr(17,25));
          date = std::stoi(line.substr(0,8));
          first = dateconversion(date, 16);
      }
      sizex = sizex + 1;
    }
    file.close();
    
    // Get infile size and first date if in the second date file
    int sizey = 0;
    int firsty;
    file.open(outfilename);    
    while (getline(file,line))
    {
      if (sizey == 0)
      {
          //date = std::stoi(line.substr(17,25));
          date = std::stoi(line.substr(0,8));
          firsty = dateconversion(date, 16);
          if (firsty < first) first = firsty;
      }
      sizey = sizey + 1;
    }
    file.close();

    // Load infil in vector x with shift according to first date
    std::vector<unsigned int> x(sizex);
    file.open(infilename);    
    int i = 0;
    while (getline(file,line))
    {
      //date = std::stoi(line.substr(17,25));
      date = std::stoi(line.substr(0,8));
      x[i] = dateconversion(date, 16) - first;
      i++;
    }
    file.close();

    // Get outfile size
    // Load outfile in y vector with shift according to first date
    std::vector<unsigned int> y(sizey);
    file.open(outfilename);    
    i = 0;
    while (getline(file,line))
    { 
      //date = std::stoi(line.substr(17,25));
      date = std::stoi(line.substr(0,8));
      y[i] = dateconversion(date, 16) - first;
      i++;
    }
    file.close();
    return std::tie(x,y,sizex,sizey);
}

int dateconversion(int date, int referenceYear){
    int d,m,mt,y,stamp;
    int monthlength[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    d = date % 100;
    m = (date/100) % 100;
    y = (date/10000) % 100;
    
    mt = 0;
    for(int i = 0; i < m-1; i++){
      mt = mt + monthlength[i];
    } 
    stamp = (y - referenceYear)*365 + mt  + (d-1);
    return stamp;
}


