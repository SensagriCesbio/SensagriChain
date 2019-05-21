/////////////////////////////////////////////
//          Create Fusion Products         //
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
using OutputImageType = otb::VectorImage<PixelType, 2>;
using RegionType = InputImageType::RegionType;
typedef otb::RAMDrivenAdaptativeStreamingManager<InputImageType> StreamingManagerType;
typedef itk::ImageRegionConstIterator<InputImageType> InIteratorType;
typedef itk::ImageRegionIterator<OutputImageType> OutIteratorType;
typedef itk::ImageRegionSplitter<2> ImageSplitter;

using Vector = itk::VariableLengthVector<unsigned int>;

// Prototypes

Vector CalculateFusionProducts(Vector px1, Vector px2, std::vector<unsigned int> classlabels, std::vector<unsigned int> binlabels);
void ReadLabels(std::string filename, std::vector<unsigned int> * labels);
void printProgress (double percentage);
int main(int argc, char * argv[])
{
    // Input Paramters //
    std::string filename_ProbaMap1;
    std::string filename_ProbaMap2;
    std::string filename_Classes;
    std::string filename_BinaryClasses;
    std::string filename_ProductsMap;

    Vector pix1;
    Vector pix2;
    Vector pixProd ;

    if (argc == 6)
    {
      filename_ProbaMap1 = argv[1];
      filename_ProbaMap2 = argv[2];
      filename_Classes = argv[3];
      filename_BinaryClasses = argv[4];
      filename_ProductsMap = argv[5];
    }
    else
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<
      " ProbabilityMap1 ProbabilityMap2 ClassesList BinaryClassesList ProductsMap" << std::endl;
      return EXIT_FAILURE;
    }

   std::vector<unsigned int> classlabels;
   std::vector<unsigned int> binlabels;
   ReadLabels(filename_Classes, &classlabels);
   ReadLabels(filename_BinaryClasses, &binlabels);

   // DEBUG: Print labels
   //for(int i=0; i<classlabels.size();i++)
   //{
   //  std::cout << classlabels[i] << ", " << binlabels[i]   << std::endl;  
   //}


   //std::cout << filename_Map1 << std::endl;
   //std::cout << filename_Map2 << std::endl;
   //std::cout << filename_ValidationMap << std::endl;
   //std::cout << filename_OutputMap << std::endl;
   
   // Load input Images //
   auto reader1 = otb::ImageFileReader<InputImageType>::New();
   auto reader2 = otb::ImageFileReader<InputImageType>::New();
   RegionType streamingRegion;
 
   // Map 1: A priori S1 probability map
   reader1->SetFileName(filename_ProbaMap1);
   reader1->GetOutput()->UpdateOutputInformation();
   reader1->GenerateOutputInformation();

   // Map 2: A priori S2 probability map
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
   OutputImageType::Pointer FusionProductsImage = OutputImageType::New();
   FusionProductsImage->SetRegions(largestRegion1);
   FusionProductsImage->CopyInformation(reader1->GetOutput());   
   FusionProductsImage->SetNumberOfComponentsPerPixel(4);
   FusionProductsImage->Allocate();

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

       OutIteratorType itProd(FusionProductsImage, streamingRegion);

       //itOut.GoToBegin();
       for(it1.GoToBegin(); !it1.IsAtEnd(); ++it1)
       {   
            pix1 = it1.Get();
            pix2 = it2.Get();
            //std::cout << "p1" << pix1 << std::endl;
            //std::cout << "p2" << pix2 << std::endl;

	    pixProd = CalculateFusionProducts(pix1,pix2,classlabels,binlabels);
            itProd.Set(pixProd);
            ++itProd;
            ++it2;
        }
     }      
   
   // Write Fusion Map //
   ////otb::StandardFilterWatcher watcher(writer, "Writing Classification Map...");
   std::cout << "*** Creation Products Map ***"<< std::endl;
   auto writer = otb::ImageFileWriter<OutputImageType>::New();
   writer->SetFileName(filename_ProductsMap);
   writer->SetInput(FusionProductsImage);
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

Vector CalculateFusionProducts(Vector pix1, Vector pix2, std::vector<unsigned int> classlabels, std::vector<unsigned int> binlabels)
{
  unsigned int CropMask;
  unsigned int CropMaskQuality;
  unsigned int CropType;
  unsigned int CropTypeQuality;

  unsigned int argmax = 0;
  float max = 0.0;
  float norm = 0.0; 
  float n0 = 0.0; 
  float n1 = 0.0; 
  float p0 = 0.0; 
  float p1 = 0.0; 

  unsigned int NbrBin[2]; 
  float ProbCM[2]; 
  Vector pixProd;
  pixProd.SetSize(4); // CropMask,CropMaskQuality,CropType,CropTypeQuality

  unsigned int Nbr = pix1.GetSize();
  itk::VariableLengthVector<float> floatpix1;   
  itk::VariableLengthVector<float> floatpix2;   
  itk::VariableLengthVector<float> Pc;
  floatpix1.SetSize(Nbr);
  floatpix2.SetSize(Nbr);
  Pc.SetSize(Nbr);

  // normalize PIX1 //
  norm = 0.0;
  for(int i=0; i<Nbr; i++)
  {
    norm = norm + pix1[i];
  }
  for(int i=0; i<Nbr; i++)
  {
    floatpix1[i] = (float)pix1[i]/(float)norm;
  }
 
  // normalize PIX2 //
  norm = 0.0;
  for(int i=0; i<Nbr; i++)
  {
    norm = norm + pix2[i];
  }
  for(int i=0; i<Nbr; i++)
  {
    floatpix2[i] = (float)pix2[i]/(float)norm;
  }
 
  // Calculate new probability vector thanks to the fusion method // 
  norm = 0.0;
  for(int i=0; i<Nbr; i++)
  {
    Pc[i] = floatpix1[i]*floatpix2[i];
    norm = norm + Pc[i];
  }

  // Treat pathological case // 
  if(norm == 0.0)
  {
    for(int i=0; i<Nbr; i++)
    {
      Pc[i] = (floatpix1[i] + floatpix2[i])*0.5;
      norm = norm + Pc[i];
    }
  }

  // Calculate number of Crop and NoCrop //
  // And found max index
  for(int i=0; i<Nbr; i++)
  {
    Pc[i] = Pc[i]/norm;
    if(Pc[i]>max)
    {
      max = Pc[i];
      argmax = i;
    }
    n0 = n0 + (int)(binlabels[i]==0);
    n1 = n1 + (int)(binlabels[i]==1);

    p0 = p0 + Pc[i]*(int)(binlabels[i]==0);
    p1 = p1 + Pc[i]*(int)(binlabels[i]==1);

  }
  // Vector number of NoCrop/Crop classes
  NbrBin[0] = n0;
  NbrBin[1] = n1;

  ProbCM[0] = p0;
  ProbCM[1] = p1;



//  // Calculate total probability of Crop and NoCrop //
//  // And found max index
//  for(int i=0; i<Pc.GetSize(); i++)
//  {
//    Pc[i] = Pc[i]/norm;
//    if(Pc[i]>max)
//    {
//      max = Pc[i];
//      argmax = i;
//    }
//    p0 = p0 + Pc[i]*(int)(labels[i]==0);
//    p1 = p1 + Pc[i]*(int)(labels[i]==1);
//  }
 
  CropMask = binlabels[argmax];
  CropType = classlabels[argmax]*CropMask;
  CropMaskQuality = (unsigned int)1000.0*(ProbCM[CropMask]*2.0-1.0);
  CropTypeQuality = (unsigned int)1000.0*((Pc[argmax]*(float)Nbr-1.0)/((float)Nbr-1.0))*CropMask;

  pixProd[0] = CropMask;
  pixProd[1] = CropMaskQuality;
  pixProd[2] = CropType;
  pixProd[3] = CropTypeQuality;

//  No use anymore. Has been show to be underperformant
//  // Compare total probability of Crop and NoCrop //
//  Fus = (int)(p0<p1);
//  Con = (1-Fus)*std::round(1000*p0/(p0+p1)) + Fus*std::round(1000*p1/(p0+p1));
//  if(Fus == 1 and Con>800)
//  {
//  std::cout << "bin = ";
//  for(int i = 0; i<18; i++) std::cout << labels[i] << ",";
//  std::cout << std::endl;
//  std::cout << "Pc = "<< Pc << std::endl;
//  std::cout << "(p0,p1) = " << p0 << "," << p1 << std::endl;
//  std::cout << "(Fus,Con) = " << Fus << "," << Con << std::endl;
//  std::cout << std::endl;
//  } 

  return pixProd;
}

void ReadLabels(std::string filename, std::vector<unsigned int> * labels)
{
  std::ifstream inputFile(filename);
  for (std::string line; std::getline(inputFile , line);)
  {
    labels->push_back(atoi(line.c_str()));
  }  
}
