/*=========================================================================

Program:   ORFEO Toolbox
Language:  C++
Date:      $Date$
Version:   $Revision$


Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
See OTBCopyright.txt for details.


This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

#include "otbChangeLabelImageFilter.h"
#include "otbStandardWriterWatcher.h"
#include "otbStatisticsXMLFileReader.h"
#include "otbShiftScaleVectorImageFilter.h"
#include "ImageClassificationFilter.h"
#include "otbMultiToMonoChannelExtractROI.h"
#include "otbImageToVectorImageCastFilter.h"
#include "MachineLearningModelFactory.h"
#include <boost/algorithm/string.hpp>
#include "MachineLearningModel.h"
#include "RandomForestsMachineLearningModel.h"
#include "ReadWriteSamples.hpp"
#include "otbConfusionMatrixCalculator.h"
#include "WriteConfusionMatrix.hpp"


// #include "otbVectorImage.h"

namespace otb
{
namespace Wrapper
{

class RFClassification : public Application
{
public:
	/** Standard class typedefs. */
	typedef RFClassification		Self;
	typedef Application                   	Superclass;
	typedef itk::SmartPointer<Self>       	Pointer;
	typedef itk::SmartPointer<const Self> 	ConstPointer;

	/** Standard macro */
	itkNewMacro(Self);
	itkTypeMacro(RFClassification, otb::Application);

	/** Filters typedef */
	typedef UInt16ImageType                                                                     OutputImageType;
	typedef UInt8ImageType                                                                      MaskImageType;
	typedef itk::VariableLengthVector<FloatVectorImageType::InternalPixelType>                  MeasurementType;
	typedef otb::StatisticsXMLFileReader<MeasurementType>                                       StatisticsReader;
	typedef otb::ShiftScaleVectorImageFilter<FloatVectorImageType, FloatVectorImageType>        RescalerType;
	typedef myImageClassificationFilter<FloatVectorImageType, OutputImageType, MaskImageType>   ClassificationFilterType;
	
	typedef ClassificationFilterType::Pointer                                                   ClassificationFilterPointerType;
	typedef ClassificationFilterType::ModelType                                                 ModelType;
	typedef ModelType::Pointer                                                                  ModelPointerType;
	typedef ClassificationFilterType::ValueType                                                 ValueType;
	typedef ClassificationFilterType::LabelType                                                 LabelType;
	typedef myMachineLearningModelFactory<ValueType, LabelType>                                 MachineLearningModelFactoryType;
	typedef ClassificationFilterType::ConfidenceImageType                                       ConfidenceImageType;
	typedef ClassificationFilterType::ProbaImageType                                            ProbaImageType;
	typedef ClassificationFilterType::MapType                                                   MapType;

	typedef otb::myRandomForestsMachineLearningModel<float, unsigned int> 		 	    RandomForestType;
        typedef otb::ConfusionMatrixCalculator<ListLabelSampleType,ListLabelSampleType> ConfusionMatrixCalculatorType;

private:

	void DoInit()
	{
		SetName("ProbRFClassification");

		/* INPUTS */

		AddParameter(ParameterType_InputFilename, "trainprofiles", "Input samples for training");

		AddParameter(ParameterType_InputFilename, "valprofiles",  "Input samples for validation");

		AddParameter(ParameterType_Int, "bands",  "Features per date");
		AddParameter(ParameterType_StringList, "dateslist",  "List of dates indices use for the classification");	
		AddParameter(ParameterType_InputFilename, "model", "Model file to export training model");
                MandatoryOff("model");

		AddParameter(ParameterType_InputFilename, "usemodel", "Model file to import for validation");
                MandatoryOff("usemodel");

                AddParameter(ParameterType_InputFilename, "adapt", "Option to adapt legend (OBSOLETE)");
                MandatoryOff("adapt");

		AddParameter(ParameterType_InputFilename, "trainlabels", "Input labels for training");		

		AddParameter(ParameterType_InputFilename, "vallabels", "Input labels for validation");		
		
		/* OUTPUTS */
		// AddParameter(ParameterType_OutputFilename, "ooberror", "OOB error array filename");
		AddParameter(ParameterType_OutputFilename, "predlabels", "Predicted Labels for validation");
		AddParameter(ParameterType_OutputFilename, "prob", "Probabilities for each classes map");
		AddParameter(ParameterType_OutputFilename, "conf", "Confidence values");	
		AddParameter(ParameterType_OutputFilename, "cmatrix", "Confusion matrix");
		AddParameter(ParameterType_OutputFilename, "vi", "Variables importances");

                /* Add possibility of adding a date mask */
		AddParameter(ParameterType_InputFilename, "datesmask", "Date mask");
                MandatoryOff("datesmask");
		AddRAMParameter();

		AddParameter(ParameterType_InputFilename,"nomen", "CSV file containing the list of classes");

	}

	void DoUpdateParameters()
	{
		// Nothing to do here : all parameters are independent
	}


        itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer AllToDate(itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer AllSamples, unsigned int numFeatures)
        {
          int NbrSamples = AllSamples->Size();
          int AllSize = AllSamples->GetMeasurementVector(0).Size();
          itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer Samples
	       = itk::Statistics::ListSample<itk::VariableLengthVector<float>>::New();
          Samples->SetMeasurementVectorSize(numFeatures);
          itk::VariableLengthVector<float> sample_i;
          sample_i.SetSize(numFeatures);

          for(int i = 0; i<NbrSamples; i++)
          {
            for(int j = 0; j < numFeatures; j++)
            {
              sample_i[j] = AllSamples->GetMeasurementVector(i)[j];
            }
            Samples->PushBack(sample_i);     
          }
          return Samples;
        }

	void DoExecute()
	{
		// Get input parameter

                std::string listTrainLabelsFilename     	= GetParameterString("trainlabels");
		std::string listValLabelsFilename		= GetParameterString("vallabels");
                std::string listTrainSamplesFilename	        = GetParameterString("trainprofiles");
		std::string listValSamplesFilename 		= GetParameterString("valprofiles");
		// std::string listOOBeFilename	 		= GetParameterString("ooberror");

		std::string listlabelsTemplate		= GetParameterString("predlabels");
		std::string listconfidencesTemplate	= GetParameterString("conf");
		std::string listprobabilitiesTemplate 	= GetParameterString("prob");
		std::string confusionmatrixTemplate 	= GetParameterString("cmatrix");
		std::string variableimportanceTemplate 	= GetParameterString("vi");

                bool Qtrain = HasValue("model");
                bool Qusemodel = HasValue("usemodel");
                bool Qadapt = HasValue("adapt");

		if(Qtrain==false and Qusemodel==false)
		{
	          otbAppLogFATAL("The model or usemodel parameter should be provided");
                  exit(1);
		}

		std::string classifiermodelTemplate;
		std::string classifiermodelTemplateUse;
		classifiermodelTemplate = GetParameterString("model");
		if(Qusemodel==true) classifiermodelTemplateUse = GetParameterString("usemodel");

               	std::string ListLabelsFilename;
               	std::string ListConfidencesFilename;
		std::string ListProbabilitiesFilename;
		std::string ConfusionMatrixFilename;
		std::string VariableImportanceFilename;
		std::string ClassifierModelFilename;
		
		unsigned int numDates;
		unsigned int numBands 	= GetParameterInt("bands");

		// ATTENTION MAPPING //
	        // Conversion legende//
		// PARTICULIER  A NOTRE ANALYSE. A ENLEVER PLUS TARD.
		std::map<unsigned int, unsigned int> AdaptInToMod;
		std::map<unsigned int, unsigned int> AdaptModToOut;
		if(Qadapt==true) 
		{
		    otbAppLogINFO("*** Legend mapping activated ***");

		    AdaptInToMod[33] = 3;
		    AdaptInToMod[10] = 10;
		    AdaptInToMod[12] = 12;
		    AdaptInToMod[14] = 14;
		    AdaptInToMod[41] = 41;
		    AdaptInToMod[438] = 438;
		    AdaptInToMod[9115] = 791;
		    AdaptInToMod[2001] = 2001;
		    AdaptInToMod[3001] = 3001;
		    AdaptInToMod[4001] = 4001;
		    AdaptInToMod[435] = 4351;
		    AdaptInToMod[5002] = 5002;
		    AdaptInToMod[5003] = 5003; 
		    AdaptInToMod[7001] = 7001;
		    AdaptInToMod[8001] = 8001;
		    AdaptInToMod[3] = 10012;


		    AdaptModToOut[3] = 33;
		    AdaptModToOut[10] = 10;
		    AdaptModToOut[12] = 12;
		    AdaptModToOut[14] = 14;
		    AdaptModToOut[41] = 41;
		    AdaptModToOut[438] = 438;
		    AdaptModToOut[791] = 9115;
		    AdaptModToOut[2001] = 2001;
		    AdaptModToOut[3001] = 3001;
		    AdaptModToOut[4001] = 4001;
		    AdaptModToOut[4351] = 435;
		    AdaptModToOut[5002] = 5002;
		    AdaptModToOut[5003] = 5003; 
		    AdaptModToOut[7001] = 7001;
		    AdaptModToOut[8001] = 8001;
		    AdaptModToOut[10012] = 3;
		}
      	        else
		{
		    for(int cl = 0; cl < 10000; cl++)
		    {
   		        AdaptInToMod[cl] = cl;
		        AdaptModToOut[cl] = cl;
		    }
		    otbAppLogINFO("*** No legend mapping ***");
		} 
		otbAppLogINFO("*** RF classification Started***");
 
                // Construct dates list
                std::vector< std::string > DatesListString = GetParameterStringList("dateslist");
	        itk::VariableLengthVector<unsigned int> DatesList;
                DatesList.SetSize(DatesListString.size());
                for(int i=0; i<DatesListString.size(); i++)
                { 
                  DatesList[i] = static_cast<unsigned int>(std::stoi(DatesListString[i])); 
                }
	        otbAppLogINFO("List of dates: " << DatesList);

                // Prepare loading of input data
                itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer listTrainSamplesAll 
	          = itk::Statistics::ListSample<itk::VariableLengthVector<float>>::New();
                itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listTrainLabels 
	          = itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();
                itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer listValSamplesAll 
	          = itk::Statistics::ListSample<itk::VariableLengthVector<float>>::New();
                itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listValLabels 
	        = itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();


                unsigned int numFeaturesAll = ComputeXDim(listValSamplesFilename);
	        otbAppLogINFO("Total number of features: " << numFeaturesAll);

                if(Qtrain)
		{
	  	  otbAppLogINFO("Loading training labels");
                  listTrainLabels->SetMeasurementVectorSize(1);
	          readListLabel (listTrainLabels, listTrainLabelsFilename);
	          otbAppLogINFO("-> " << listTrainLabels->Size() << " samples loaded\n");
                }

	        otbAppLogINFO("Loading validation labels");
	        listValLabels->SetMeasurementVectorSize(1);
	        readListLabel (listValLabels, listValLabelsFilename);
	        otbAppLogINFO("-> " << listValLabels->Size() << " samples loaded\n");

		if(Qtrain)
	        {
      	          otbAppLogINFO("Loading training profiles");
                  listTrainSamplesAll->SetMeasurementVectorSize(numFeaturesAll);
	          //ReadListSampleInstantT(listTrainSamples, listTrainSamplesFilename, numFeatures);
	          ReadListSample(listTrainSamplesAll, listTrainSamplesFilename);
	          otbAppLogINFO("Loaded\n");
                }
		 
                otbAppLogINFO("Loading validation profiles");
                listValSamplesAll->SetMeasurementVectorSize(numFeaturesAll);
	        //ReadListSampleInstantT(VallistSamples, listValSamplesFilename, numFeatures);
	        ReadListSample(listValSamplesAll, listValSamplesFilename);
	        otbAppLogINFO("Loaded\n");

                for(int d = 0; d < DatesList.Size(); d++)
                {
                  unsigned int numDates = DatesList[d];
		  auto numFeatures = numDates * numBands;
                  otbAppLogINFO("Classification Started: date="<< numDates << ", Nbr of Features=" << numFeatures);
               
                  // Create Output File names for current Date 
                  ListLabelsFilename		= listlabelsTemplate         + "_Date" + DatesListString[d] + ".txt";
                  ListConfidencesFilename	= listconfidencesTemplate    + "_Date" + DatesListString[d] + ".txt";
		  ListProbabilitiesFilename 	= listprobabilitiesTemplate  + "_Date" + DatesListString[d] + ".txt";
		  ConfusionMatrixFilename 	= confusionmatrixTemplate    + "_Date" + DatesListString[d] + ".csv";
		  VariableImportanceFilename 	= variableimportanceTemplate + "_Date" + DatesListString[d] + ".txt";
      	          ClassifierModelFilename       = classifiermodelTemplate    + "_Date" + DatesListString[d] + ".txt";

                  //otbAppLogINFO("listconfidencesfilename: " << listconfidencesfilename);
                  //otbAppLogINFO("confusionmatrixfilename: " << confusionmatrixfilename);
	          //otbAppLogINFO("Skip training ...");
	          /* train the model to obtain the OOB error */

	          RandomForestType::Pointer classifier = RandomForestType::New();
                  std::vector<float> oob_array;
                  std::vector< std::vector<float> > oob_votes;

		  itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer listTrainSamples; 
		  if(Qtrain)
	          {
	            
                    otbAppLogINFO("Calculate date samples for training");
		    auto listTrainSamples = AllToDate(listTrainSamplesAll,numFeatures);
		    otbAppLogDEBUG("Done: First train sample size " << listTrainSamples->GetMeasurementVector(0).Size());

		    classifier->SetInputListSample(listTrainSamples);
	            classifier->SetTargetListSample(listTrainLabels);

	            classifier->SetMaxDepth(25);
	            classifier->SetCalculateVariableImportance(true);
	            unsigned int p = floor(std::sqrt(numFeatures));
	            classifier->SetMinSampleCount(p);
	            classifier->SetMaxNumberOfTrees(100);
	            	
	            otbAppLogINFO("Starting learning ...");
	            classifier->Train();
	            otbAppLogINFO("Learning done."); 
	            otbAppLogINFO("Saving model ...");
	            classifier->Save(ClassifierModelFilename);

                    oob_array = classifier->GetOOBErrors();
                    oob_votes = classifier->GetOOBVotes();

		  }

		  otbAppLogINFO("Calculate date samples for validation");
                  auto listValSamples = AllToDate(listValSamplesAll,numFeatures); 
                  otbAppLogINFO("Done");

       	          // Load model
		  if(Qusemodel==true) ClassifierModelFilename = classifiermodelTemplateUse + "_Date" + DatesListString[d] + ".txt";
	          otbAppLogINFO("Loading model from "<< ClassifierModelFilename);
	          m_Model = MachineLearningModelFactoryType::CreateMachineLearningModel(
	          	ClassifierModelFilename,
	          	MachineLearningModelFactoryType::ReadMode);

	          if (m_Model.IsNull())
	          {
	          	otbAppLogFATAL(<< "Error when loading model " << ClassifierModelFilename << " : unsupported model type");
	          }
	          m_Model->Load(ClassifierModelFilename);

	          otbAppLogINFO("Model loaded");

	          //Si parametre CSV actif initialiser map
	          if(IsParameterEnabled("nomen") && HasValue("nomen"))
	          {
	          	std::string in_nomen_file{""};
	          	in_nomen_file = GetParameterString("nomen");
	          	std::ifstream csvFile(in_nomen_file);
	          	
	          	if(!csvFile)
	          	{
	          		otbAppLogWARNING(<< "Could not open file" << in_nomen_file << "\n");
	          		this->DisableParameter("nomen");
	          	}
	          	
	          	std::string ligne;
	          	int nbLignes = 0;
	          	
	          	while(std::getline(csvFile, ligne))
	          		nbLignes++;

	          	csvFile.close();
	          	nbLignes = nbLignes-2; // supprimer la ligne des non decisions

	          	//int idx;
	          	std::ifstream csvFile2(in_nomen_file);
	          	MapType maplabels;
	          	char delim = '\n';
	          	
	          	for(int i = 0; i<=nbLignes; ++i)
	          	{
	          		std::getline(csvFile2, ligne, delim);
	          		std::vector<std::string> chaine;
	          		boost::split(chaine,ligne,boost::is_any_of(":"));
	          		//std::cout << stoi(chaine[1]) << " " << i << "\n";
	          		maplabels[stoi(chaine[1])] = i;
	          		//idx++;
	          	}

	          	csvFile2.close();
	          	
	          	if (maplabels.size() > 1)
	          	{
	          		otbAppLogINFO("Actif LabelMap.");
	          		m_Model->SetUseLabelMap(true);
	          		m_Model->SetLabelMap(maplabels);
	          	}
	          	else
	          	{
	          		otbAppLogWARNING(<<"Empty file : "<<in_nomen_file<<"\n");
	          		this->DisableParameter("nomen");
	          	}
	          }

	          // std::vector<unsigned int> classlabels;
	          // m_Model->GetLabels(&classlabels);

	          // for(unsigned int i = 0; i < classlabels.size(); ++i)
	          // {
	          // 	std::cout << "Label["<< i << "]: " << classlabels[i] << std::endl;
	          // }

	          otbAppLogINFO("Loading OOB Error array ...");
	          m_Model->SetOOBErrors(oob_array);
	          m_Model->SetOOBVotes(oob_votes);

	          // std::vector<float> test = m_Model->GetOOBErrors();

	          // for(unsigned int i = 0; i < 10; ++i)
	          // {
	          // 	std::cout << "oob_array["<< i << "]: " << test[i] << std::endl;
	          // }

	          otbAppLogINFO("Defining outputs ...");

	          // Outputs definitons 
	          itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listLabels 
	          = itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();


	          listLabels->SetMeasurementVectorSize(1);

	          itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listConfidences 
	          = itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();
	          
	          listConfidences->SetMeasurementVectorSize(1);

	          itk::Statistics::ListSample<itk::VariableLengthVector<double>>::Pointer listProbabilities 
	          = itk::Statistics::ListSample<itk::VariableLengthVector<double>>::New();

	          listProbabilities->SetMeasurementVectorSize(m_Model->GetNumberOfClasses());

	          /* CLASSIFICATION */
	          double confidenceIndex = 0.0;

	          itk::VariableLengthVector<double> probaVector{m_Model->GetNumberOfClasses()};
	          probaVector.Fill(0.0);

	          otbAppLogINFO("Classification starts");

	          for (unsigned int i = 0; i < listValSamples->Size(); ++i)
	          {
                        
	          	itk::FixedArray<unsigned int, 1> label
			= m_Model->Predict(listValSamples->GetMeasurementVector(i), &confidenceIndex, &probaVector)[0];
			
			// Old Version Without Mapping //
			//listLabels->PushBack(label);	
			itk::FixedArray<unsigned int, 1>  adapted_legend_label;
                 	adapted_legend_label[0]= AdaptModToOut[label[0]];
	          	listLabels->PushBack(adapted_legend_label);	
                        listConfidences->PushBack(confidenceIndex);
	          	listProbabilities->PushBack(probaVector);
			
                        //itk::FixedArray<unsigned int, 1>  adapted_legend_label = label;
                     	// otbAppLogINFO(<< "Predicted Label: " << label);
	          	// otbAppLogINFO(<< "Predicted Conf: "  << confidenceIndex);
	          	// otbAppLogINFO(<< "Predicted Proba: " << listProbabilities[i]);

	          }


                  // confusion matrix //
                  ConfusionMatrixCalculatorType::Pointer confMatCalc = ConfusionMatrixCalculatorType::New();
                  confMatCalc->SetReferenceLabels(listValLabels);
                  confMatCalc->SetProducedLabels(listLabels);
                  confMatCalc->Compute();
                 
                  //std::cout<< "Precision of the different class: " << confMatCalc->GetPrecisions() << std::endl;
                  //std::cout<<"Recall of the different class: " << confMatCalc->GetRecalls() << std::endl;
                  //std::cout<< "F-score of the different class: " << confMatCalc->GetFScores() << std::endl;
                  //std::cout<< "Kappa index: " << confMatCalc->GetKappaIndex() << std::endl;


	          otbAppLogINFO("Classification Finished.");

                  
                  auto VI = classifier->GetVariableImportance(); 
                  std::ofstream outfile(VariableImportanceFilename);                 
                  if (outfile.is_open())        
                  {                             
                    for(int i = 0; i < VI.Cols(); i++) 
                    { 
                      outfile << VI(0,i) << std::endl ;                           
                    }                         
                    outfile.close();
	            otbAppLogINFO("Variable Importance File created.");
                  }                             
                  else std::cout << "Unable to open variable importance file";  
                
                
                  writeListLabel(listLabels, ListLabelsFilename);
	          otbAppLogINFO("Predicted Labels File created.");

	          writeListLabel(listConfidences, ListConfidencesFilename);
	          otbAppLogINFO("Confidences File created.");

	          writeListSample(listProbabilities, ListProbabilitiesFilename);
	          otbAppLogINFO("Probabilities File created.");

                  writeConfusionMatrix(ConfusionMatrixFilename,confMatCalc);
                  otbAppLogINFO("Confusion matrix file created.");
               }

	}

	ClassificationFilterType::Pointer m_ClassificationFilter;
	ModelPointerType m_Model;
	RandomForestType::Pointer classifier;
	RescalerType::Pointer m_Rescaler;
};


}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::RFClassification)
