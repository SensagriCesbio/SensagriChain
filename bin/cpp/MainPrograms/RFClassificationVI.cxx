

#include "ReadWriteSamples.hpp"

#include <iostream>
#include <stdexcept>

# include "otbRandomForestsMachineLearningModel.h"
#include "otbMachineLearningModelFactory.h"

#include "otbConfusionMatrixCalculator.h"

#include "WriteConfusionMatrix.hpp"

typedef float PixelType;
typedef unsigned int LabelType;

typedef itk::VariableLengthVector<PixelType> SampleType;
typedef itk::Statistics::ListSample<SampleType> ListSampleType;

typedef itk::Statistics::ListSample<itk::FixedArray<LabelType, 1>> ListLabelSampleType;

typedef itk::FixedArray<double,1>            ConfidenceSampleType;
typedef itk::Statistics::ListSample<ConfidenceSampleType> ConfidenceListSampleType;

typedef otb::RandomForestsMachineLearningModel<PixelType, LabelType> RandomForestType;

typedef otb::ConfusionMatrixCalculator<ListLabelSampleType,ListLabelSampleType> ConfusionMatrixCalculatorType;

int main(int argc, char * argv[])
{
      std::string filename_Training_Samples;
      std::string filename_Label_Training_Samples;
      std::string filename_Validation_Samples;
      std::string filename_Label_Validation_Samples;
    
      std::string filename_ClassifierModel;
      std::string filename_PredictedLabels;
      std::string filename_ConfusionMatrix;
      std::string filename_VariableImportance;
      unsigned int NbFeaturesForEachDate;
      unsigned int NbDates;

      ListSampleType::Pointer trainingListSample = ListSampleType::New() ;
      ListLabelSampleType::Pointer trainingListLabel = ListLabelSampleType::New() ;

    if (argc == 11 )
    {
      std::cout << "Classification: Learning + Validation" << std::endl;
      filename_Training_Samples = argv[1];
      filename_Label_Training_Samples = argv[2];
      filename_Validation_Samples = argv[3];
      filename_Label_Validation_Samples = argv[4];
    
      filename_ClassifierModel = argv[5];
      filename_PredictedLabels = argv[6];
      filename_ConfusionMatrix = argv[7];
      filename_VariableImportance = argv[8];
      NbFeaturesForEachDate = static_cast<unsigned int>(atoll(argv[9]));
      NbDates = static_cast<unsigned int>(atoll(argv[10]));

      auto NbTotalFeatures =  ComputeXDim (filename_Training_Samples );
      auto NbFeatures =  NbDates *NbFeaturesForEachDate ;

      trainingListSample->SetMeasurementVectorSize(NbFeatures);
    
      ReadListSampleInstantT(trainingListSample, filename_Training_Samples,NbFeatures );
      trainingListLabel ->SetMeasurementVectorSize(1);

      readListLabel (trainingListLabel  ,filename_Label_Training_Samples);

    }

    if (argc == 9 )
    {
      std::cout << "Classification: Validation only" << std::endl;
      filename_Validation_Samples = argv[1];
      filename_Label_Validation_Samples = argv[2];
    
      filename_ClassifierModel = argv[3];
      filename_PredictedLabels = argv[4];
      filename_ConfusionMatrix = argv[5];
      filename_VariableImportance = argv[6];
      NbFeaturesForEachDate = static_cast<unsigned int>(atoll(argv[7]));
      NbDates = static_cast<unsigned int>(atoll(argv[8]));

    }
    
    //std::cout << " NbFeatures for each sample :: " <<  NbTotalFeatures  << std::endl;
    //std::cout << " NbFeatures for each date :: " <<  NbFeaturesForEachDate  << std::endl;
    //std::cout << " Nb of total dates :: " <<  NbTotalFeatures/NbFeaturesForEachDate  << std::endl;

    auto NbFeatures =  NbDates *NbFeaturesForEachDate ;
    
    ListSampleType::Pointer validationListSample = ListSampleType::New() ;
    validationListSample->SetMeasurementVectorSize(NbFeatures);
    ReadListSampleInstantT(validationListSample,filename_Validation_Samples,NbFeatures );

    ListLabelSampleType::Pointer validationListLabel = ListLabelSampleType::New() ;
    validationListLabel->SetMeasurementVectorSize(1);
    
    readListLabel (validationListLabel ,filename_Label_Validation_Samples);
    
    RandomForestType::Pointer classifier = RandomForestType::New();
   
    if(argc == 11)
    {

      std::cout << "*** Training part ***"  << std::endl; 
      if(trainingListLabel->Size()!=  trainingListSample->Size())
          throw std::invalid_argument("Problem with the training data\n");
    
      if(validationListLabel->Size()!=  validationListSample->Size())
          throw std::invalid_argument("Problem with the validation data\n");
   

      classifier->SetInputListSample(trainingListSample);
      classifier->SetTargetListSample( trainingListLabel);
      classifier->SetMaxDepth(25);
      classifier->SetCalculateVariableImportance(true);
    
      // Number of Sample Count -> Choice between sqrt(N) or N/3
      unsigned int p = floor(std::sqrt(NbFeatures));
      classifier->SetMinSampleCount(p);
      classifier->SetMaxNumberOfTrees(100);
      std::cout << "STARTING LEARNING "<<std::endl;
    
      classifier->Train();
      std::cout << "LEARNING DONE"  <<std::endl;
      classifier->Save(filename_ClassifierModel);
    }
 
    std::cout << "*** Validation part ***"  << std::endl; 
    typedef otb::MachineLearningModelFactory<PixelType, LabelType> MachineLearningModelFactoryType;
    typedef MachineLearningModelFactoryType::MachineLearningModelTypePointer ModelPointerType;
    
    ModelPointerType model = MachineLearningModelFactoryType::CreateMachineLearningModel(filename_ClassifierModel ,MachineLearningModelFactoryType::ReadMode);
    model->Load(filename_ClassifierModel);
  
    ListLabelSampleType::Pointer predictedList =  ListLabelSampleType:: New();
    //ConfidenceListSampleType quality ;

    
    predictedList = model->PredictBatch(validationListSample);
    
    writeListLabel(predictedList,filename_PredictedLabels);
    ConfusionMatrixCalculatorType::Pointer confMatCalc = ConfusionMatrixCalculatorType::New();    
    confMatCalc->SetReferenceLabels(validationListLabel);
    confMatCalc->SetProducedLabels(predictedList);
    confMatCalc->Compute();
    
    std::cout<< "Precision of the different class: " << confMatCalc->GetPrecisions() << std::endl;
    std::cout<<"Recall of the different class: " << confMatCalc->GetRecalls() << std::endl;
    std::cout<< "F-score of the different class: " << confMatCalc->GetFScores() << std::endl;
    std::cout<< "Kappa index: " << confMatCalc->GetKappaIndex() << std::endl;
    
    writeConfusionMatrix(filename_ConfusionMatrix ,confMatCalc);

  //  auto VI = classifier->GetVariableImportance();
  //  std::ofstream outfile(filename_VariableImportance);                
  //  if (outfile.is_open())       
  //  {                            
  //    for(int i = 0; i < VI.Cols(); i++)
  //    {
  //      outfile << VI(0,i) << std::endl ;                          
  //    }                        
  //    outfile.close();         
  //  }                            
  //  else std::cout << "Unable to open variable importance file"; 

    return 0;
}
