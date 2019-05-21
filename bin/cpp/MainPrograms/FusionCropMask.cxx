/////////////////////////////////////////////
//          Create Fusion Crop Mask        //
//          Ludo 17/05/2018                //
/////////////////////////////////////////////
#include "otbImage.h"
#include "otbVectorImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbStandardFilterWatcher.h"
#include "otbRAMDrivenAdaptativeStreamingManager.h"
#include "otbImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionSplitter.h"
#include <ostream>
#include <iostream>
#include <stdexcept>
#include <math.h>

typedef unsigned int PixelType;
using InputImageType = otb::VectorImage<PixelType, 2>;
//using OutputImageType = otb::Image<unsigned int, 2 >;
using OutputImageType = otb::Image<PixelType, 2>;
using RegionType = InputImageType::RegionType;
typedef otb::RAMDrivenAdaptativeStreamingManager<InputImageType> StreamingManagerType;
typedef itk::ImageRegionConstIterator<InputImageType> InIteratorType;
typedef itk::ImageRegionIterator<OutputImageType> OutIteratorType;
typedef itk::ImageRegionSplitter<2> ImageSplitter;

using Vector = itk::VariableLengthVector<unsigned int>;
using FusionVariables = std::tuple<unsigned int, unsigned int>;

// Prototypes

FusionVariables CalculateFusionCropMask(Vector px1, Vector px2, std::vector<unsigned int> labels);
FusionVariables CalculateFusionCropType(Vector px1, Vector px2, std::vector<unsigned int> labelsCT);
void ReadLabels(std::string filename, std::vector<unsigned int> * labels);
void printProgress (double percentage);
int main(int argc, char * argv[])
{
    // Input Paramters //
    std::string filename_ProbaMap1;
    std::string filename_ProbaMap2;
    std::string filename_Classes;
    std::string filename_BinaryClasses;
    std::string filename_FusionMap;
    std::string filename_ConfidenceMap;
    std::string filename_FusionMapCT;
    std::string filename_ConfidenceMapCT;
 

    Vector pix1;
    Vector pix2;
    unsigned int Fus;
    unsigned int Con;
    unsigned int FusCT;
    unsigned int ConCT;
 

    FusionVariables FusVar;  
    FusionVariables FusVarCT;  

    if (argc == 10)
    {
      filename_ProbaMap1 = argv[1];
      filename_ProbaMap2 = argv[2];
      filename_Classes = argv[3];
      filename_BinaryClasses = argv[4];
      filename_FusionMap = argv[5];
      filename_ConfidenceMap = argv[6];
      filename_FusionMapCT = argv[7];
      filename_ConfidenceMapCT = argv[8];
 
    }
    else
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<
      " ProbabilityMap1 ProbabilityMap2 BinaryClassesList FusionMap ConfidenceMap" << std::endl;
      return EXIT_FAILURE;
    }

   std::vector<unsigned int> labelsCT;
   ReadLabels(filename_Classes, &labelsCT);

   std::vector<unsigned int> labels;
   ReadLabels(filename_BinaryClasses, &labels);


   //std::cout << filename_Map1 << std::endl;
   //std::cout << filename_Map2 << std::endl;
   //std::cout << filename_ValidationMap << std::endl;
   //std::cout << filename_OutputMap << std::endl;
   
   // Load input Images //
   auto reader1 = otb::ImageFileReader<InputImageType>::New();
   auto reader2 = otb::ImageFileReader<InputImageType>::New();
   RegionType streamingRegion;
 
   // Map 1
   reader1->SetFileName(filename_ProbaMap1);
   reader1->GetOutput()->UpdateOutputInformation();
   reader1->GenerateOutputInformation();

   // Map 2
   reader2->SetFileName(filename_ProbaMap2);
   reader2->GetOutput()->UpdateOutputInformation();
   reader2->GenerateOutputInformation();

   // Assigne largest region
   RegionType largestRegion1 = reader1->GetOutput()->GetLargestPossibleRegion();
   RegionType largestRegion2 = reader2->GetOutput()->GetLargestPossibleRegion();
   std::cout << "Largest region Map 1" << largestRegion1 << std::endl;
   std::cout << "Largest region Map 2" << largestRegion2 << std::endl;
   std::cout << "Number of features per pixel for Map 1 = " << reader1->GetOutput()->GetNumberOfComponentsPerPixel() << std::endl;
   std::cout << "Number of features per pixel for Map 2 = " << reader2->GetOutput()->GetNumberOfComponentsPerPixel() << std::endl;
 
   const int availableRAM = 128;
   // Streaming Manager fo Map1
   StreamingManagerType::Pointer streamingManager1 = StreamingManagerType::New();
   streamingManager1->SetAvailableRAMInMB(availableRAM);
   streamingManager1->PrepareStreaming(reader1->GetOutput(), largestRegion1);
   const unsigned long numberOfStreamDivisions1 = streamingManager1->GetNumberOfSplits(); 
   
   // Streaming Manager fo Map2
   StreamingManagerType::Pointer streamingManager2 = StreamingManagerType::New();
   streamingManager2->SetAvailableRAMInMB(availableRAM);
   streamingManager2->PrepareStreaming(reader2->GetOutput(), largestRegion2);
   const unsigned long numberOfStreamDivisions2 = streamingManager2->GetNumberOfSplits(); 

   std::cout << std::endl;
   std::cout << "Number of splits Map 1 = " << numberOfStreamDivisions1 << std::endl;
   std::cout << "Number of splits Map 2 = " << numberOfStreamDivisions1 << std::endl;
   
   if (numberOfStreamDivisions1 == numberOfStreamDivisions2)
   {
     std::cout << "Maps files are compatible" << std::endl;
   }
   else {
     std::cout << "Problem in maps files. Check if the number of channels are the same." << std::endl;
     return EXIT_FAILURE;
   }  

   // Set output image for fusion //
   OutputImageType::Pointer FusionImage = OutputImageType::New();
   FusionImage->SetRegions(largestRegion1);
   FusionImage->CopyInformation(reader1->GetOutput());   
   FusionImage->Allocate();

   // Set output image for confidenc //
   OutputImageType::Pointer ConfidenceImage = OutputImageType::New();
   ConfidenceImage->SetRegions(largestRegion1);
   ConfidenceImage->CopyInformation(reader1->GetOutput());   
   ConfidenceImage->Allocate();

   // Set output image for fusion CT//
   OutputImageType::Pointer FusionImageCT = OutputImageType::New();
   FusionImageCT->SetRegions(largestRegion1);
   FusionImageCT->CopyInformation(reader1->GetOutput());   
   FusionImageCT->Allocate();

   // Set output image for confidenc CT //
   OutputImageType::Pointer ConfidenceImageCT = OutputImageType::New();
   ConfidenceImageCT->SetRegions(largestRegion1);
   ConfidenceImageCT->CopyInformation(reader1->GetOutput());   
   ConfidenceImageCT->Allocate();
 

   std::cout << "*** Map construction *** " << std::endl;
   for (auto piece = 0;piece < numberOfStreamDivisions1; piece++)
   {
       std::cout << "#" << std::flush;
       if((piece+1) % 50 == 0) std::cout << std::endl;
       streamingRegion = streamingManager1->GetSplit(piece);

       // Stream 1
       reader1->GetOutput()->SetRequestedRegion(streamingRegion);
       reader1->GetOutput()->PropagateRequestedRegion();
       reader1->GetOutput()->UpdateOutputData();
       // Stream 2
       reader2->GetOutput()->SetRequestedRegion(streamingRegion);
       reader2->GetOutput()->PropagateRequestedRegion();
       reader2->GetOutput()->UpdateOutputData();
   
       InIteratorType it1(reader1->GetOutput(), streamingRegion);
       InIteratorType it2(reader2->GetOutput(), streamingRegion);

       OutIteratorType itFus(FusionImage, streamingRegion);
       OutIteratorType itCon(ConfidenceImage, streamingRegion);

       OutIteratorType itFusCT(FusionImageCT, streamingRegion);
       OutIteratorType itConCT(ConfidenceImageCT, streamingRegion);

       //itOut.GoToBegin();
       for(it1.GoToBegin(); !it1.IsAtEnd(); ++it1)
       {   
            pix1 = it1.Get();
            pix2 = it2.Get();
            //std::cout << "p1" << pix1 << std::endl;
            //std::cout << "p2" << pix2 << std::endl;
            FusVar = CalculateFusionCropMask(pix1,pix2,labels);
            Fus = std::get<0>(FusVar);
            Con = std::get<1>(FusVar);

            FusVarCT = CalculateFusionCropType(pix1,pix2,labelsCT);
            FusCT = std::get<0>(FusVar);
            ConCT = std::get<1>(FusVar);

 
            ++it2;

            // Fill fusion and confidence pixels
            //std::cout << "Fus = "<< Fus << std::endl; 
            //std::cout << "Con = "<< Con << std::endl; 
            itFus.Set(Fus);
            itCon.Set(Con);
	    ++itFus;
            ++itCon;

            ++itFusCT;
            ++itConCT;
        }
     }      
   
   // Write Fusion Map //
   ////otb::StandardFilterWatcher watcher(writer, "Writing Classification Map...");
   std::cout << "*** Creation Fusion Map ***"<< std::endl;
   auto writer = otb::ImageFileWriter<OutputImageType>::New();
   writer->SetFileName(filename_FusionMap);
   writer->SetInput(FusionImage);
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

   // Write Fusion Map //
   std::cout << "*** Creation Confidence Map ***"<< std::endl;
   writer->SetFileName(filename_ConfidenceMap);
   writer->SetInput(ConfidenceImage);
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

   std::cout << "*** Creation Fusion Map CT***"<< std::endl;
   writer = otb::ImageFileWriter<OutputImageType>::New();
   writer->SetFileName(filename_FusionMapCT);
   writer->SetInput(FusionImageCT);
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

   // Write Fusion Map //
   std::cout << "*** Creation Confidence Map CM***"<< std::endl;
   writer->SetFileName(filename_ConfidenceMapCT);
   writer->SetInput(ConfidenceImage);
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

FusionVariables CalculateFusionCropMask(Vector pix1, Vector pix2, std::vector<unsigned int> labels)
{
  FusionVariables val;
  unsigned int Fus;
  unsigned int Con;
  unsigned int argmax = 0;
  float max = 0;
  float norm = 0.0; 
  float p0 = 0.0; 
  float p1 = 0.0; 
  itk::VariableLengthVector<float> Pc;   
  Pc.SetSize(pix1.GetSize());

  // normalize PIX1 //
  norm = 0.0;
  for(int i=0; i<pix1.GetSize(); i++)
  {
    norm = norm + pix1[i];
  }
  for(int i=0; i<pix1.GetSize(); i++)
  {
    pix1[i] = (int)1000*pix1[i]/norm;
  }
 
  // normalize PIX2 //
  norm = 0.0;
  for(int i=0; i<pix2.GetSize(); i++)
  {
    norm = norm + pix2[i];
  }
  for(int i=0; i<pix2.GetSize(); i++)
  {
    pix2[i] = (int)1000*pix2[i]/norm;
  }
 
  // Calculate new probability vector thanks to the fusion method // 
  for(int i=0; i<pix1.GetSize(); i++)
  {
    Pc[i] = pix1[i]*pix2[i];
    norm = norm + Pc[i];
  }

  // Treat pathological case // 
  if(norm == 0.0)
  {
    for(int i=0; i<pix1.GetSize(); i++)
    {
      Pc[i] = (pix1[i] + pix2[i])*0.5;
      norm = norm + Pc[i];
    }
  }
 
  // Calculate total probability of Crop and NoCrop //
  for(int i=0; i<Pc.GetSize(); i++)
  {
    Pc[i] = Pc[i]/norm;
    p0 = p0 + Pc[i]*(int)(labels[i]==0);
    p1 = p1 + Pc[i]*(int)(labels[i]==1);
  }
  
  // Compare total probability of Crop and NoCrop //
  Fus = (int)(p0<p1);
  Con = (1-Fus)*std::round(1000*p0/(p0+p1)) + Fus*std::round(1000*p1/(p0+p1));
  //if(Fus == 1 and Con>800)
  //{
  //std::cout << "bin = ";
  //for(int i = 0; i<18; i++) std::cout << labels[i] << ",";
  //std::cout << std::endl;
  //std::cout << "Pc = "<< Pc << std::endl;
  //std::cout << "(p0,p1) = " << p0 << "," << p1 << std::endl;
  //std::cout << "(Fus,Con) = " << Fus << "," << Con << std::endl;
  //std::cout << std::endl;
  //} 
  
  std::get<0>(val) = Fus;
  std::get<1>(val) = Con;

  return val;
}

FusionVariables CalculateFusionCropType(Vector pix1, Vector pix2, std::vector<unsigned int> labelsCT)
{
  
  FusionVariables val;
  return val;

}

void ReadLabels(std::string filename, std::vector<unsigned int> * labels)
{
  std::ifstream inputFile(filename);
  for (std::string line; std::getline(inputFile , line);)
  {
    labels->push_back(atoi(line.c_str()));
  }  
}
