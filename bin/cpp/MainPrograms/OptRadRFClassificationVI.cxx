

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

    std::string filename_Training_Samples = argv[1];
    std::string filename_Label_Training_Samples = argv[2];
    std::string filename_Validation_Samples = argv[3];
    std::string filename_Label_Validation_Samples = argv[4];
    
    std::string filename_ClassifierModel = argv[5];
    std::string filename_PredictedLabels = argv[6];
    std::string filename_ConfusionMatrix = argv[7];
    std::string filename_VariableImportance = argv[8];
    unsigned int NbFeatures = static_cast<unsigned int>(atoll(argv[9]));
    unsigned int NbDates = static_cast<unsigned int>(atoll(argv[10]));

    auto NbTotalFeatures =  ComputeXDim (filename_Training_Samples );
    
    //std::cout << " NbFeatures for each sample :: " <<  NbTotalFeatures  << std::endl;
    //std::cout << " NbFeatures for each date :: " <<  NbFeaturesForEachDate  << std::endl;
    //std::cout << " Nb of total dates :: " <<  NbTotalFeatures/NbFeaturesForEachDate  << std::endl;
    //std::cout << "Nb of NbFeatures at this step :: " <<  NbFeatures   << std::endl;
   
    ListSampleType::Pointer trainingListSample = ListSampleType::New() ;
    trainingListSample->SetMeasurementVectorSize(NbFeatures);
 
    ListSampleType::Pointer validationListSample = ListSampleType::New() ;
    validationListSample->SetMeasurementVectorSize(NbFeatures);
    
    ReadListSampleInstantT(trainingListSample, filename_Training_Samples,NbFeatures );
    ReadListSampleInstantT(validationListSample,filename_Validation_Samples,NbFeatures );
    
    ListLabelSampleType::Pointer trainingListLabel = ListLabelSampleType::New() ;
    trainingListLabel ->SetMeasurementVectorSize(1);
    
    ListLabelSampleType::Pointer validationListLabel = ListLabelSampleType::New() ;
    validationListLabel->SetMeasurementVectorSize(1);
    
    readListLabel (trainingListLabel  ,filename_Label_Training_Samples);
    readListLabel (validationListLabel ,filename_Label_Validation_Samples);
    
    if(trainingListLabel->Size()!=  trainingListSample->Size())
        throw std::invalid_argument("Problem with the training data\n");
    
    if(validationListLabel->Size()!=  validationListSample->Size())
        throw std::invalid_argument("Problem with the validation data\n");
  

    RandomForestType::Pointer classifier = RandomForestType::New();
    classifier->SetInputListSample(trainingListSample);
    classifier->SetTargetListSample( trainingListLabel);
    classifier->SetMaxDepth(25);
    classifier->SetCalculateVariableImportance(true);
    //unsigned int p = floor(std::sqrt(NbFeatures));
    unsigned int p = floor(float(NbFeatures)/3.0);
    classifier->SetMinSampleCount(p);
    classifier->SetMaxNumberOfTrees(100);
    std::cout << "STARTING LEARNING "<<std::endl;
    
    classifier->Train();
    std::cout << "LEARNING DONE"  <<std::endl;
    classifier->Save(filename_ClassifierModel);
  
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

    auto VI = classifier->GetVariableImportance();
    std::ofstream outfile(filename_VariableImportance);                
    if (outfile.is_open())       
    {                            
      for(int i = 0; i < VI.Cols(); i++)
      {
        outfile << VI(0,i) << std::endl ;                          
      }                        
      outfile.close();         
    }                            
    else std::cout << "Unable to open variable importance file"; 
 
    return 0;
}
