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
	typedef RFClassification			  	Self;
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
	typedef myImageClassificationFilter<FloatVectorImageType, OutputImageType, MaskImageType>	ClassificationFilterType;
	
	typedef ClassificationFilterType::Pointer                                                   ClassificationFilterPointerType;
	typedef ClassificationFilterType::ModelType                                                 ModelType;
	typedef ModelType::Pointer                                                                  ModelPointerType;
	typedef ClassificationFilterType::ValueType                                                 ValueType;
	typedef ClassificationFilterType::LabelType                                                 LabelType;
	typedef myMachineLearningModelFactory<ValueType, LabelType>                                 MachineLearningModelFactoryType;
	typedef ClassificationFilterType::ConfidenceImageType                                       ConfidenceImageType;
	typedef ClassificationFilterType::ProbaImageType                                            ProbaImageType;
	typedef ClassificationFilterType::MapType                                                   MapType;

	typedef otb::myRandomForestsMachineLearningModel<float, unsigned int> 						RandomForestType;
        typedef otb::ConfusionMatrixCalculator<ListLabelSampleType,ListLabelSampleType> ConfusionMatrixCalculatorType;

private:

	void DoInit()
	{
		SetName("ProbRFClassification");

		/* INPUTS */

		AddParameter(ParameterType_InputFilename, "trainprofiles", "Input samples for training");

		AddParameter(ParameterType_InputFilename, "valprofiles",  "Input samples for validation");

		AddParameter(ParameterType_Int, "bands",  "Features per date");
		AddParameter(ParameterType_Int, "dates",  "Number of dates");	
		AddParameter(ParameterType_InputFilename, "model", "Model file");

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

	void DoExecute()
	{
		// Load input Profiles
		std::string listSamplesFilename 		= GetParameterString("valprofiles");

		std::string listValLabelsFilename		= GetParameterString("vallabels");

		std::string listLabelsFilename	 		= GetParameterString("predlabels");

		// std::string listOOBeFilename	 		= GetParameterString("ooberror");

		std::string listConfidencesFilename		= GetParameterString("conf");
		std::string listProbabilitiesFilename 	= GetParameterString("prob");
		std::string ConfusionMatrixFilename 	= GetParameterString("cmatrix");
		std::string VariableImportanceFilename 	= GetParameterString("vi");

		std::string filename_ClassifierModel 	= GetParameterString("model");
		
		unsigned int numBands 	= GetParameterInt("bands");
		unsigned int numDates	= GetParameterInt("dates");

		auto numFeatures = numDates * numBands;
                
                itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer listTrainSamples 
		  = itk::Statistics::ListSample<itk::VariableLengthVector<float>>::New();
                itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listTrainLabels 
		  = itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();


                itk::Statistics::ListSample<itk::VariableLengthVector<float>>::Pointer listSamples 
		= itk::Statistics::ListSample<itk::VariableLengthVector<float>>::New();

        	listSamples->SetMeasurementVectorSize(numFeatures);

		ReadListSampleInstantT(listSamples, listSamplesFilename, numFeatures);
		otbAppLogINFO("Loaded Profiles samples for classification.");

 
                std::string listTrainSamplesFilename	= GetParameterString("trainprofiles");
        	std::string listTrainLabelsFilename		= GetParameterString("trainlabels");

	        otbAppLogINFO("Learning is going to be performed ...");
		              
	    	listTrainSamples->SetMeasurementVectorSize(numFeatures);

		ReadListSampleInstantT(listTrainSamples, listTrainSamplesFilename, numFeatures);
		otbAppLogINFO("Loaded Profiles samples for training.");

		  
            	listTrainLabels->SetMeasurementVectorSize(1);

	 	readListLabel (listTrainLabels, listTrainLabelsFilename);
	 	otbAppLogINFO("Loaded Reference labels for training.");
                

        	itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::Pointer listValLabels 
		= itk::Statistics::ListSample<itk::FixedArray<unsigned int, 1>>::New();

		listValLabels->SetMeasurementVectorSize(1);

		readListLabel (listValLabels, listValLabelsFilename);
		otbAppLogINFO("Loaded Reference labels for validation.");


		std::cout << "listSamples Size:"  << listSamples->Size() << std::endl;


		//otbAppLogINFO("Skip training ...");
		/* train the model to obtain the OOB error */
                 
	        RandomForestType::Pointer classifier = RandomForestType::New();
               
                // Modif start here so far
 

 
	        classifier->SetInputListSample( listTrainSamples );
	        classifier->SetTargetListSample( listTrainLabels );
	        classifier->SetMaxDepth(25);
	        classifier->SetCalculateVariableImportance(false);
	        unsigned int p = floor(std::sqrt( numFeatures ));
	        classifier->SetMinSampleCount( p );
	        classifier->SetMaxNumberOfTrees(100);
	        	
	        otbAppLogINFO("Starting learning ...");
	        
	        classifier->Train();

	        otbAppLogINFO("Learning done."); 
	        otbAppLogINFO("Saving model ...");

	        classifier->Save(filename_ClassifierModel);

                
         	std::vector<float> oob_array = classifier->GetOOBErrors();
	       	std::vector< std::vector<float> > oob_votes = classifier->GetOOBVotes();

		// Load model
		otbAppLogINFO("Loading model ...");

		m_Model = MachineLearningModelFactoryType::CreateMachineLearningModel(
			filename_ClassifierModel,
			MachineLearningModelFactoryType::ReadMode);

		if (m_Model.IsNull())
		{
			otbAppLogFATAL(<< "Error when loading model " << GetParameterString("model") << " : unsupported model type");
		}

		m_Model->Load(filename_ClassifierModel);

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

		otbAppLogINFO("Classification starts.");

		for (unsigned int i = 0; i < listSamples->Size(); ++i)
		{
                        
			itk::FixedArray<unsigned int, 1> label = 
			m_Model->Predict(listSamples->GetMeasurementVector(i), &confidenceIndex, &probaVector)[0];
			listLabels->PushBack(label);	
                        listConfidences->PushBack(confidenceIndex);
			listProbabilities->PushBack(probaVector);
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
                
               
              	writeListLabel(listLabels, listLabelsFilename);
		otbAppLogINFO("Predicted Labels File created.");

		writeListLabel(listConfidences, listConfidencesFilename);
		otbAppLogINFO("Confidences File created.");

		writeListSample(listProbabilities, listProbabilitiesFilename);
		otbAppLogINFO("Probabilities File created.");

                writeConfusionMatrix(ConfusionMatrixFilename,confMatCalc);
                otbAppLogINFO("Confusion matrix file created.");


	}

	ClassificationFilterType::Pointer m_ClassificationFilter;
	ModelPointerType m_Model;
	RandomForestType::Pointer classifier;
	RescalerType::Pointer m_Rescaler;
};


}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::RFClassification)
