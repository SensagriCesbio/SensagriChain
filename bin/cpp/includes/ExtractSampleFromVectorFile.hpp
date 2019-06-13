/*
 * ExtractSampleFromVectorFile.hpp
 *
 *  Update: Ludo (add lc, code and crop management) 
 *  Created on: 19/1/2015
 *      Author: silvia
 */

#ifndef ExtractSampleFromVectorFile_HPP_
#define ExtractSampleFromVectorFile_HPP_


#include "itkProcessObject.h"
#include "itkListSample.h"
#include "itkPreOrderTreeIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
//#include "itkMersenneTwisterRandomVariateGenerator.h"


template <class TImage, class TVectorData>
class ITK_EXPORT ExtractSampleFromVectorFile :public itk::ProcessObject
{

public:


	typedef ExtractSampleFromVectorFile Self;
	typedef itk::ProcessObject Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(ExtractSampleFromVectorFile, itk::ProcessObject);

	typedef unsigned int ClassLabelType;

	typedef typename TImage::PixelType SampleType;
	typedef itk::Statistics::ListSample<SampleType> ListSampleType;

	typedef itk::FixedArray<ClassLabelType, 1> CropTypeLabel; 
	typedef itk::Statistics::ListSample<CropTypeLabel> ListCropLabelType;

	typedef itk::FixedArray<ClassLabelType, 1> BinaryCropLabel; 
	typedef itk::Statistics::ListSample<BinaryCropLabel> ListBinaryCropLabelType;

  	typedef itk::FixedArray<unsigned int, 2> CoordinateType; 
  	typedef itk::Statistics::ListSample<CoordinateType> CoordinateList;


	// Build the outputs
	typedef itk::DataObject::Pointer DataObjectPointer;
	virtual DataObjectPointer MakeOutput(unsigned int idx);

	ExtractSampleFromVectorFile();

	void SetInputImage(TImage * InputImage);

	const TImage* GetInput() const;

	void SetInputVectorData(const TVectorData * VectorDataInput);

	const TVectorData* GetInputVectorData() const;

	void GenerateLists();

	int GetNumberOfClasses( );

	std::string SetCropTypeClassKey(std::string inputstring);
	std::string SetIsCropClassKey(std::string inputstring);
	std::string SetLabelClassKey(std::string inputstring);

	std::map<ClassLabelType, std::string>*  GetInfoLabelsOfClasses ( );

	std::map<ClassLabelType, unsigned int>*  GetInfoSamplesOfClasses ( );

	void PrintInformation( );

	/** Returns the Training ListSample as a data object */
	ListSampleType * GetListSample();

	/** Returns the list of crop type labels data object */
	ListCropLabelType * GetCropTypeListLabel();

	/** Returns the list of binary crop labels data object */
	ListBinaryCropLabelType * GetBinaryCropTypeListLabel();

	/** Returns the image coordinates sample list data object */
  	CoordinateList* GetCoordinatesList();


private:

	typedef typename TVectorData::DataNodeType DataNodeType;
	typedef typename DataNodeType::PolygonType PolygonType;
	typedef typename DataNodeType::PolygonPointerType PolygonPointerType;
	typedef typename DataNodeType::PolygonListType PolygonListType;
	typedef typename DataNodeType::PolygonListPointerType PolygonListPointerType;
	typedef itk::PreOrderTreeIterator<typename TVectorData::DataTreeType> TreeIteratorType;

	double GetPolygonAreaPixels(DataNodeType* polygonDataNode, TImage* image);

	int _NumberOfClasses;

	unsigned int _LimitationofMaximumTrainingSamples;

	std::string _CropTypeClassKey;
	std::string _IsCropClassKey;
	std::string _LabelClassKey;

	std::map<ClassLabelType, std::string> _ClassesLabels;

	std::map<ClassLabelType, unsigned int> _ClassesSize;

};


#endif /* ExtractSampleFromVectorFile_HPP_ */
