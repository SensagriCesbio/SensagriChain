// Calculate statistics from probability vector //

#include "ReadWriteSamples.hpp"
#include <iostream>
#include <stdexcept>
#include "otbConfusionMatrixCalculator.h"
#include "WriteConfusionMatrix.hpp"

typedef unsigned int LabelType;
typedef itk::VariableLengthVector<float> Vector;
typedef itk::VariableLengthVector<unsigned int> uintVector;
typedef itk::Statistics::ListSample<Vector> ListProbType;
typedef itk::Statistics::ListSample<itk::FixedArray<LabelType, 1>> ListLabelSampleType;
typedef otb::ConfusionMatrixCalculator<ListLabelSampleType,ListLabelSampleType> ConfusionMatrixCalculatorType;
using FusionVariables = std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>;

// Prototypes
uintVector ReadLabels(std::string filename);
int argmax(Vector vector);
FusionVariables CalculateFusion(Vector p1, Vector p2, uintVector labels, uintVector binarylabels);   
unsigned int CalculateCropMask(Vector p, uintVector binarylabels);

int main(int argc, char * argv[])
{

      std::string filename_probabilitiesS1;
      std::string filename_probabilitiesS2;
      std::string filename_Validation_Labels;
      std::string filename_Classes;
      std::string filename_ConfusionMatrixFusion;
      std::string filename_Validation_BinaryLabels;
      std::string filename_Binary;
      std::string filename_S1CM;
      std::string filename_S2CM;
      std::string filename_FusCM;
      ListLabelSampleType::Pointer validationListBinaryLabel = ListLabelSampleType::New() ;

    if(argc == 6)
    {
      std::cout << "### Parameters Consistant with Crop Type ###" << std::endl;
      filename_probabilitiesS1 = argv[1];
      filename_probabilitiesS2 = argv[2];
      filename_Validation_Labels = argv[3];
      filename_Classes = argv[4];
      filename_ConfusionMatrixFusion = argv[5];
    }
    else if(argc == 11)
    {
      std::cout << "### Parameters Consistant with Crop Mask ###" << std::endl;
      filename_probabilitiesS1 = argv[1];
      filename_probabilitiesS2 = argv[2];
      filename_Validation_Labels = argv[3];
      filename_Classes = argv[4];
      filename_ConfusionMatrixFusion = argv[5];
      filename_Validation_BinaryLabels = argv[6];
      filename_Binary = argv[7];
      filename_S1CM = argv[8];
      filename_S2CM = argv[9];
      filename_FusCM = argv[10];
    }
    else
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] <<
      " prob 1 prob2 labelsVal Classes FusionConfMat (binarylabelsVal Binary S1CM S2CM FusionBinaryCM" << std::endl;
      return EXIT_FAILURE;
    }
    auto NbProbabilities =  ComputeXDim (filename_probabilitiesS1);
    auto NbLabels =  ComputeXDim (filename_Validation_Labels);
    
    ListProbType::Pointer ProbabilitiesS1 = ListProbType::New();
    ProbabilitiesS1->SetMeasurementVectorSize(NbProbabilities);

    ListProbType::Pointer ProbabilitiesS2 = ListProbType::New();
    ProbabilitiesS2->SetMeasurementVectorSize(NbProbabilities);
   
    ListLabelSampleType::Pointer validationListLabel = ListLabelSampleType::New() ;
    validationListLabel->SetMeasurementVectorSize(1);
    readListLabel (validationListLabel, filename_Validation_Labels);
    std::cout << "Validation labels loaded" << std::endl;
 
    uintVector labels;
    uintVector binarylabels;
    labels = ReadLabels(filename_Classes);
    if(argc == 11)
    {
      binarylabels = ReadLabels(filename_Binary);
      validationListBinaryLabel->SetMeasurementVectorSize(1);
      readListLabel (validationListBinaryLabel, filename_Validation_BinaryLabels);
      std::cout << "Validation binary labels loaded" << std::endl;

    }
    else
    {
      binarylabels = labels;
      binarylabels.Fill(1);
    }

    ReadListSample(ProbabilitiesS1, filename_probabilitiesS1);
  
    std::cout << "Proba S1 loaded" << std::endl;

    ReadListSample(ProbabilitiesS2, filename_probabilitiesS2);
    std::cout << "Proba S2 loaded" << std::endl;
   
    if(validationListLabel->Size()!= ProbabilitiesS1->Size())
        throw std::invalid_argument("Problem with the training data\n");
    
    if(validationListLabel->Size()!= ProbabilitiesS2->Size())
        throw std::invalid_argument("Problem with the validation data\n");

    int val;
    Vector p1;
    Vector p2;
    ListLabelSampleType::Pointer predictedFusion = ListLabelSampleType::New() ;
    ListLabelSampleType::Pointer predictedS1CM = ListLabelSampleType::New() ;
    ListLabelSampleType::Pointer predictedS2CM = ListLabelSampleType::New() ;
    ListLabelSampleType::Pointer predictedFusCM = ListLabelSampleType::New() ;

    predictedFusion->SetMeasurementVectorSize(1);
    predictedS1CM->SetMeasurementVectorSize(1);
    predictedS2CM->SetMeasurementVectorSize(1);
    predictedFusCM->SetMeasurementVectorSize(1);

    FusionVariables fusion;
    int in = 0;
    int out = 0;
    for (auto i = 0; i < ProbabilitiesS1->Size(); ++i){
      p1 = ProbabilitiesS1->GetMeasurementVector(i);
      p2 = ProbabilitiesS2->GetMeasurementVector(i);
      val = validationListLabel->GetMeasurementVector(i)[0];
  
      // Calculate Fusion
      fusion = CalculateFusion(p1, p2, labels, binarylabels);
      predictedFusion->PushBack(std::get<0>(fusion)); 
      predictedFusCM->PushBack(std::get<2>(fusion)); 
      in = in + std::get<3>(fusion);
      //std::cout << val << ", "  << std::get<0>(fusion) << std::endl;  
      //std::cout << std::endl;
       
      // Calculate Predicted Crop mask
      
      if(argc == 11)
      {
        predictedS1CM->PushBack(CalculateCropMask(p1,binarylabels));
        predictedS2CM->PushBack(CalculateCropMask(p2,binarylabels));
      }
    }

   // writeListLabel(predictedList,filename_PredictedLabels);
    ConfusionMatrixCalculatorType::Pointer confMatCalc = ConfusionMatrixCalculatorType::New();    
    confMatCalc->SetReferenceLabels(validationListLabel);
    confMatCalc->SetProducedLabels(predictedFusion);
    confMatCalc->Compute();
    writeConfusionMatrix(filename_ConfusionMatrixFusion,confMatCalc);
 
    if(argc == 11)
    {
      ConfusionMatrixCalculatorType::Pointer MatS1CM = ConfusionMatrixCalculatorType::New();    
      ConfusionMatrixCalculatorType::Pointer MatS2CM = ConfusionMatrixCalculatorType::New();    
      ConfusionMatrixCalculatorType::Pointer MatFusCM = ConfusionMatrixCalculatorType::New();    

      MatS1CM->SetReferenceLabels(validationListBinaryLabel);
      MatS2CM->SetReferenceLabels(validationListBinaryLabel);
      MatFusCM->SetReferenceLabels(validationListBinaryLabel);

      MatS1CM->SetProducedLabels(predictedS1CM);
      MatS2CM->SetProducedLabels(predictedS2CM);
      MatFusCM->SetProducedLabels(predictedFusCM);

      MatS1CM->Compute();
      MatS2CM->Compute();
      MatFusCM->Compute();

      writeConfusionMatrix(filename_S1CM,MatS1CM);
      writeConfusionMatrix(filename_S2CM,MatS2CM);
      writeConfusionMatrix(filename_FusCM,MatFusCM);
   }
   
   std::cout << in << "/" << ProbabilitiesS1->Size() << std::endl;

   // if(binarylabels.GetSize() > 0)
    //{
     // std::cout << "Calculate Crop Mask statistics" << std::endl;
        
   // }    

    return 0;
}


itk::VariableLengthVector<unsigned int> ReadLabels(std::string filename)
{
  int i = 1;
  int dim = 0;
  itk::VariableLengthVector<unsigned int> labels;
  std::ifstream inputFile;
  inputFile.open(filename);
  for (std::string line; std::getline(inputFile , line);)
  {        
     dim = dim + 1;
     i++;
  }
  inputFile.close();
  labels.SetSize(dim);

  i = 0;
  inputFile.open(filename);
  for (std::string line; std::getline(inputFile , line);)
  { 
    labels[i] = atoi(line.c_str());
    i++;
  }
  inputFile.close();
  return labels;
}

int argmax(itk::VariableLengthVector<float> vector)
{
  float max = 0;
  int arg = 0;
  for(int i = 0; i < vector.GetSize(); i++)
  {
    if(max<=vector[i])
    {
      max = vector[i];
      arg = i;
    }        
  }
 return arg;
}

FusionVariables CalculateFusion(Vector p1, Vector p2, uintVector labels, uintVector binarylabels)
{
  FusionVariables val;
  unsigned int Fus;
  unsigned int Con;
  unsigned int flag;
  unsigned int argmax = 0;
  float max = 0;
  float norm = 0.0; 
  itk::VariableLengthVector<int> Pc;   
  itk::VariableLengthVector<int> pix1;   
  itk::VariableLengthVector<int> pix2;   
  pix1.SetSize(p1.GetSize());
  pix2.SetSize(p2.GetSize());
  Pc.SetSize(p1.GetSize());

  // normalize P1 //
  norm = 0.0;
  for(int i=0; i<p1.GetSize(); i++)
  {
    norm = norm + p1[i];
  }
  for(int i=0; i<p1.GetSize(); i++)
  {
    pix1[i] = (int)1000*p1[i]/norm;
  }
 
  // normalize P2 //
  norm = 0.0;
  for(int i=0; i<p2.GetSize(); i++)
  {
    norm = norm + p2[i];
  }
  for(int i=0; i<p2.GetSize(); i++)
  {
    pix2[i] = (int)1000*p2[i]/norm;
  }

  //norm = 0.0;
  //flag = 0;
  //for(int i=0; i<pix1.GetSize(); i++)
  //{
  //  Pc[i] = (int)pix1[i]*pix2[i];
  //  if(pix1[i] == 0.0 or pix2[i] == 0.0 or pix1[i] > 900.0 or pix2[i] > 900.0)
  //  {
  //    flag = flag + 1;    
  //    Pc[i] = (int)(pix1[i] + pix2[i])*0.5;
  //  }
  //  norm = norm + Pc[i];
  //}

  norm = 0.0;
  for(int i=0; i<pix1.GetSize(); i++)
  {
    Pc[i] = (int)pix1[i]*pix2[i];
    norm = norm + Pc[i];
  }

  if(norm == 0.0)
  {
    flag = 1; 
    for(int i=0; i<pix1.GetSize(); i++)
    {
     Pc[i] = (int)(pix1[i] + pix2[i])*0.5;
     norm = norm + Pc[i];
    }
  }

  for(int i=0; i<pix1.GetSize(); i++)
  {
    Pc[i] = (int)std::round(1000*Pc[i]/norm);
    if(Pc[i]>max)
    {
      max = Pc[i];
      argmax = i;
    }
  }

  if (flag > 0)
  {
    flag = 1;
  }
  //std::cout << flag << std::endl;
  //std::cout << pix1 << std::endl;
  //std::cout << pix2 << std::endl;
  //std::cout << Pc << std::endl;

  // Calculate crop mask label
  int pcrop=0;
  for(int i=0; i<binarylabels.GetSize(); i++)
  {
    pcrop = pcrop + binarylabels[i]*Pc[i];
  }
  unsigned int CM = 0;
  if (pcrop>500) CM = 1; 

  Fus = labels[argmax];
  Con = max;
  std::get<0>(val) = Fus;
  std::get<1>(val) = Con;
  std::get<2>(val) = CM;
  std::get<3>(val) = flag;
  return val;

}

unsigned int CalculateCropMask(Vector p, uintVector binarylabels)
{
  itk::VariableLengthVector<int> Pnorm;   
  Pnorm.SetSize(p.GetSize());

  // normalize //
  int norm = 0.0;
  for(int i=0; i<p.GetSize(); i++)
  {
    norm = norm + p[i];
  }
  for(int i=0; i<p.GetSize(); i++)
  {
    Pnorm[i] = (int)1000*p[i]/norm;
  }

  int pcrop=0;
  for(int i=0; i<binarylabels.GetSize(); i++)
  {
    pcrop = pcrop + binarylabels[i]*Pnorm[i];
  }
  unsigned int CM = 0;
  if (pcrop>500) CM = 1; 
 
  return CM;
}
