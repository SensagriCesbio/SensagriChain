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

private:
	
	void DoInit()
	{
		SetName("RFClassification");
		SetDescription("Performs a classification of the input image according to a model file.");

		// Documentation
		SetDocName("RFClassification");
		SetDocLongDescription("This application performs an image classification based on a model file produced by the TrainImagesClassifier application. Pixels of the output image will contain the class labels decided by the classifier (maximal class label = 65535). The input pixels can be optionally centered and reduced according to the statistics file produced by the ComputeImagesStatistics application. An optional input mask can be provided, in which case only input image pixels whose corresponding mask value is greater than 0 will be classified. The remaining of pixels will be given the label 0 in the output image.");

		SetDocLimitations("The input image must have the same type, order and number of bands than the images used to produce the statistics file and the SVM model file. If a statistics file was used during training by the TrainImagesClassifier, it is mandatory to use the same statistics file for classification. If an input mask is used, its size must match the input image size.");
		SetDocAuthors("OTB-Team");
		SetDocSeeAlso("TrainImagesClassifier, ValidateImagesClassifier, ComputeImagesStatistics");

		AddDocTag(Tags::Learning);

		AddParameter(ParameterType_InputImage, "in",  "Input Image");
		SetParameterDescription( "in", "The input image to classify.");

		AddParameter(ParameterType_InputImage,  "mask",   "Input Mask");
		SetParameterDescription( "mask", "The mask allows to restrict classification of the input image to the area where mask pixel values are greater than 0.");
		MandatoryOff("mask");

		AddParameter(ParameterType_InputFilename, "model", "Model file");
		SetParameterDescription("model", "A model file (produced by TrainImagesClassifier application, maximal class label = 65535).");

		AddParameter(ParameterType_InputFilename, "imstat", "Statistics file");
		SetParameterDescription("imstat", "A XML file containing mean and standard deviation to center and reduce samples before classification (produced by ComputeImagesStatistics application).");
		MandatoryOff("imstat");

		AddParameter(ParameterType_OutputImage, "out",  "Output Image");
		SetParameterDescription( "out", "Output image containing class labels");
		SetDefaultOutputPixelType( "out", ImagePixelType_uint32);

		AddParameter(ParameterType_OutputImage, "confmap",  "Confidence map");
		SetParameterDescription( "confmap", "Confidence map of the produced classification. The confidence index depends on the model : \n"
		"  - LibSVM : difference between the two highest probabilities (needs a model with probability estimates, so that classes probabilities can be computed for each sample)\n"
		"  - OpenCV\n"
		"    * Boost : sum of votes\n"
		"    * DecisionTree : (not supported)\n"
		"    * GradientBoostedTree : (not supported)\n"
		"    * KNearestNeighbors : number of neighbors with the same label\n"
		"    * NeuralNetwork : difference between the two highest responses\n"
		"    * NormalBayes : (not supported)\n"
		"    * RandomForest : Confidence (proportion of votes for the majority class). Margin (normalized difference of the votes of the 2 majority classes) is not available for now.\n"
		"    * SVM : distance to margin (only works for 2-class models)\n");
		SetDefaultOutputPixelType( "confmap", ImagePixelType_double);
		MandatoryOff("confmap");
		// Proba output
		AddParameter(ParameterType_OutputImage, "proba", "Probabilities for each classes map");
		SetParameterDescription("proba", "Return the confidence for each classes\n RF only");
		SetDefaultOutputPixelType( "proba", ImagePixelType_double);
		MandatoryOff("proba");
		// List of classes
		// AddParameter(ParameterType_OutputFilename,"labels","Labels output file");
		// SetParameterDescription("labels","Filename to store the list of classes in the same order of the probabilities image");
		// MandatoryOff("labels");
		AddRAMParameter();

		AddParameter(ParameterType_InputFilename,"nomen","CSV file containing the list of classes");
		SetParameterDescription("nomen","CSV file containing the list of classes");
		MandatoryOff("nomen");

               	AddParameter(ParameterType_Int, "dates",  "Number of dates");
		SetParameterDescription( "dates", "Number of dates in the images.");

                AddParameter(ParameterType_Int, "prim",  "Type of primitive");
		SetParameterDescription( "prim", "Type of primitive\n"
                                                "0:                    10 Optical bands\n"  
                                                "1:              10 Optical bands + OSO\n"  
                                                "2:         10 Optical bands + Red Edge\n"  
                                                "3:   10 Optical bands + OSO + Red Edge\n"
                                                "4:       VH and VV Radar Polarizations\n"
                                                "5:            26 primitive of Sen2agri\n"
                                                "6: Crop Type sen2agri(DO NOT WORK YET)\n"
                                      );



		// Doc example parameter settings
		SetDocExampleParameterValue("in", "QB_1_ortho.tif");
		SetDocExampleParameterValue("imstat", "EstimateImageStatisticsQB1.xml");
		SetDocExampleParameterValue("model", "clsvmModelQB1.svm");
		SetDocExampleParameterValue("out", "clLabeledImageQB1.tif");
	}

	void DoUpdateParameters()
	{
		// Nothing to do here : all parameters are independent
	}

	void DoExecute()
	{
		// Load input image
		FloatVectorImageType::Pointer inImage = GetParameterImage("in");
		inImage->UpdateOutputInformation();

		// Load svm model
		otbAppLogINFO("Loading model");
		m_Model = MachineLearningModelFactoryType::CreateMachineLearningModel(
			GetParameterString("model"),
			MachineLearningModelFactoryType::ReadMode);

		if (m_Model.IsNull())
		{
			otbAppLogFATAL(<< "Error when loading model " << GetParameterString("model") << " : unsupported model type");
		}

		m_Model->Load(GetParameterString("model"));
		otbAppLogINFO("Model loaded");

		// Normalize input image (optional)
		StatisticsReader::Pointer statisticsReader = StatisticsReader::New();
		MeasurementType meanMeasurementVector;
		MeasurementType stddevMeasurementVector;
		m_Rescaler = RescalerType::New();
		
		// Classify
		m_ClassificationFilter = ClassificationFilterType::New();
		m_ClassificationFilter->SetModel(m_Model);


                // Primitives
		m_ClassificationFilter->SetNbDates(GetParameterInt("dates"));
		m_ClassificationFilter->SetPrimitives(GetParameterInt("prim"));

		// Normalize input image if asked
		if(IsParameterEnabled("imstat")  )
		{
			otbAppLogINFO("Input image normalization activated.");
			
			// Load input image statistics
			statisticsReader->SetFileName(GetParameterString("imstat"));
			meanMeasurementVector   = statisticsReader->GetStatisticVectorByName("mean");
			stddevMeasurementVector = statisticsReader->GetStatisticVectorByName("stddev");
			
			otbAppLogINFO( "mean used: " << meanMeasurementVector );
			otbAppLogINFO( "standard deviation used: " << stddevMeasurementVector );
			
			// Rescale vector image
			m_Rescaler->SetScale(stddevMeasurementVector);
			m_Rescaler->SetShift(meanMeasurementVector);
			m_Rescaler->SetInput(inImage);

			m_ClassificationFilter->SetInput(m_Rescaler->GetOutput());
		}
		else
		{
			otbAppLogINFO("Input image normalization deactivated.");
			m_ClassificationFilter->SetInput(inImage);
		}

		if(IsParameterEnabled("mask"))
		{
			otbAppLogINFO("Using input mask");
			// Load mask image and cast into LabeledImageType
			MaskImageType::Pointer inMask = GetParameterUInt8Image("mask");

			m_ClassificationFilter->SetInputMask(inMask);
		}

		SetParameterOutputImage<OutputImageType>("out", m_ClassificationFilter->GetOutput());

		// output confidence map
		if (IsParameterEnabled("confmap") && HasValue("confmap"))
		{
			m_ClassificationFilter->SetUseConfidenceMap(true);
			
			if (m_Model->HasConfidenceIndex())
			{
				SetParameterOutputImage<ConfidenceImageType>("confmap",m_ClassificationFilter->GetOutputConfidence());
			}
			else
			{
				otbAppLogWARNING("Confidence map requested but the classifier doesn't support it!");
				this->DisableParameter("confmap");
			}
		}

		//output probabilities map
		if(IsParameterEnabled("proba") && HasValue("proba"))
		{
			m_ClassificationFilter->SetUseProbaMap(true);

			//Si parametre CSV actif initialiser map
			if(IsParameterEnabled("nomen") && HasValue("nomen"))
			{
				std::string in_nomen_file{""};
				in_nomen_file = GetParameterString("nomen");
				std::ifstream csvFile(in_nomen_file);
				
				if(!csvFile)
				{
					otbAppLogWARNING(<<"Could not open file"<<in_nomen_file<<"\n");
					this->DisableParameter("nomen");
				}
				
				//itkGenericExceptionMacro(<<"Could not open file"<<in_nomen_file<<"\n");
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
					maplabels[stoi(chaine[1])]=i;
					//idx++;
				}

				csvFile2.close();
				
				if (maplabels.size() > 1)
				{
					//std::cout << "Actif \n";
					m_ClassificationFilter->SetUseLabelMap(true, maplabels);
				}
				else
				{
					otbAppLogWARNING(<<"Empty file : "<<in_nomen_file<<"\n");
					this->DisableParameter("nomen");
				}
			}
			if(m_Model->HasConfidenceIndex())
			{
				SetParameterOutputImage<ProbaImageType>("proba",m_ClassificationFilter->GetOutputProba());
			}
			else
			{
				otbAppLogWARNING("Probabilities map requested but the classifier doesn't support it!");
				this->DisableParameter("proba");
			}
		}
		// if(IsParameterEnabled("labels"))
		//   {
		// 	std::vector<unsigned int> labels;
		// 	m_Model->GetLabels(&labels);
		// 	std::ofstream outFile;
		// 	outFile.open(this->GetParameterString("labels").c_str());
		// 	outFile << std::fixed;
		// 	outFile << "Correspondance Band -> Label \n";
		// 	for(size_t i=0;i<labels.size();i++)
		// 	  outFile << "Band "<< i << ", Label : " << labels[i] << std::endl;
		// 	outFile.close();
		//   }
	}

	ClassificationFilterType::Pointer m_ClassificationFilter;
	ModelPointerType m_Model;
	RescalerType::Pointer m_Rescaler;
};


}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::RFClassification)
