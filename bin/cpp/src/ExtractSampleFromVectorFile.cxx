/*
 * ExtractSampleFromVectorFile.cxx
 *
 *  Created on: 19/1/2015
 *      Author: silvia
 */

#include "otbVectorImage.h"
#include "otbVectorData.h"
#include "ExtractSampleFromVectorFile.hpp"
#include <iostream>


template <class TImage, class TVectorData>
ExtractSampleFromVectorFile<TImage, TVectorData>::ExtractSampleFromVectorFile ( ) :
_CropTypeClassKey("CODE"), _LabelClassKey("TYPE"),_IsCropClassKey("CROP")
{

	//Para que sirve?
	this->SetNumberOfRequiredInputs(2);
	this->SetNumberOfRequiredOutputs(4);

	// Register the outputs
	this->itk::ProcessObject::SetNthOutput(0, this->MakeOutput(0).GetPointer());
	this->itk::ProcessObject::SetNthOutput(1, this->MakeOutput(1).GetPointer());
	this->itk::ProcessObject::SetNthOutput(2, this->MakeOutput(2).GetPointer());
	this->itk::ProcessObject::SetNthOutput(3, this->MakeOutput(3).GetPointer());
   

};


template <class TImage, class TVectorData>
typename ExtractSampleFromVectorFile<TImage, TVectorData>::DataObjectPointer
ExtractSampleFromVectorFile<TImage, TVectorData>::MakeOutput(unsigned int idx)
{
	DataObjectPointer output;
	switch (idx)
	{
		case 0:
			output = static_cast<itk::DataObject*>(ListSampleType::New().GetPointer());
			break;
		case 1:
			output = static_cast<itk::DataObject*>(ListCropLabelType::New().GetPointer());
			break;
		case 2:
			output = static_cast<itk::DataObject*>(ListBinaryCropLabelType::New().GetPointer());
			break;
		case 3:
			output = static_cast<itk::DataObject*>(CoordinateList::New().GetPointer());
			break;
		default:
			output = static_cast<itk::DataObject*>(CoordinateList::New().GetPointer());
			break;
	}
	return output;
}


template <class TImage, class TVectorData>
void ExtractSampleFromVectorFile<TImage, TVectorData>::SetInputImage(TImage * InputImage)
{

	this->ProcessObject::SetNthInput( 0, const_cast<TImage *>(InputImage) );

};



template <class TImage, class TVectorData>
const TImage* ExtractSampleFromVectorFile<TImage, TVectorData>::GetInput() const
{
	if (this->GetNumberOfInputs() < 1)
	{
		return 0;
	}
  return static_cast<const TImage *>(this->ProcessObject::GetInput(0));
}




template <class TImage, class TVectorData>
void ExtractSampleFromVectorFile<TImage, TVectorData>::SetInputVectorData(const TVectorData * VectorDataInput)
{

	this->ProcessObject::SetNthInput( 1, const_cast<TVectorData *>(VectorDataInput) );
};


template <class TImage, class TVectorData>
const TVectorData* ExtractSampleFromVectorFile<TImage, TVectorData>::GetInputVectorData() const
{

	if (this->GetNumberOfInputs() < 2)
	{
		return 0;
	}
  return static_cast<const TVectorData *>(this->ProcessObject::GetInput(1));
}



template <class TImage, class TVectorData>
double ExtractSampleFromVectorFile<TImage, TVectorData>::GetPolygonAreaPixels(DataNodeType*  datanode, TImage* image)
{
	auto area = 0;
	typename TImage::RegionType imageLargestRegion = image->GetLargestPossibleRegion();

	PolygonPointerType exteriorRing = datanode->GetPolygonExteriorRing();

	typename  TImage::RegionType polygonRegion =
	otb::TransformPhysicalRegionToIndexRegion(exteriorRing->GetBoundingRegion(),image);
	 const bool hasIntersection = polygonRegion.Crop(imageLargestRegion);

	if (hasIntersection)
     	{

	 image->SetRequestedRegion(polygonRegion);
	 image->PropagateRequestedRegion();
	 image->UpdateOutputData();


	 typedef itk::ImageRegionConstIteratorWithIndex<TImage> IteratorType;
	IteratorType it(image, polygonRegion);

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		itk::ContinuousIndex<double, 2> point;
		image->TransformIndexToPhysicalPoint(it.GetIndex(), point);
		bool isInterior = false;
	   if ( exteriorRing->IsInside(point) )
	   {
		PolygonListPointerType interiorRings = datanode->GetPolygonInteriorRings();
		for (typename PolygonListType::Iterator interiorRing = interiorRings->Begin();interiorRing != interiorRings->End();++interiorRing)
		{
			if ( interiorRing.Get()->IsInside(point) )
			{
	  			isInterior = true;	  
               		}
		}

	
		if(!isInterior ){
			++area;						 
		}
	     }
	         
	     }
         
	}

  return area;
}


template <class TImage, class TVectorData>
void ExtractSampleFromVectorFile<TImage, TVectorData>::GenerateLists()
{
	_ClassesSize.clear();

	TImage* image = const_cast<TImage*> (this->GetInput());
	const TVectorData* vectorData = const_cast<TVectorData*>(this->GetInputVectorData());


	typename TImage::RegionType imageLargestRegion = image->GetLargestPossibleRegion();

	typename ListSampleType::Pointer ListSample   = this->GetListSample();
	typename ListCropLabelType::Pointer  ListCropTypeLabel    = this->GetCropTypeListLabel();
	typename ListBinaryCropLabelType::Pointer  ListBinaryCropLabel    = this->GetBinaryCropTypeListLabel();
	typename CoordinateList::Pointer  CoordinatesList  = this->GetCoordinatesList();

		
	ListSample->Clear();
	ListSample->SetMeasurementVectorSize(image->GetNumberOfComponentsPerPixel());

	ListCropTypeLabel->Clear();
	ListCropTypeLabel->SetMeasurementVectorSize(1);

	ListBinaryCropLabel->Clear();
	ListBinaryCropLabel->SetMeasurementVectorSize(1);

	CoordinatesList ->Clear();
	CoordinatesList->SetMeasurementVectorSize(2);
 
	// Compute cumulative area of all polygons of each class
	TreeIteratorType itVector(vectorData->GetDataTree());
	
	for (itVector.GoToBegin(); !itVector.IsAtEnd(); ++itVector)
	{
		DataNodeType* datanode = itVector.Get();
		if (datanode->IsPolygonFeature())
		{
		  PolygonPointerType exteriorRing = datanode->GetPolygonExteriorRing();
	      typename  TImage::RegionType polygonRegion = otb::TransformPhysicalRegionToIndexRegion(exteriorRing->GetBoundingRegion(),image);

		  const bool hasIntersection = polygonRegion.Crop(imageLargestRegion);

		  if (hasIntersection)
		  {
			  double area = GetPolygonAreaPixels(datanode, image);
			 
			  if (area > 0)
			  {

				_ClassesSize[datanode->GetFieldAsInt(_CropTypeClassKey)] += std::floor(area);

					image->SetRequestedRegion(polygonRegion);
					image->PropagateRequestedRegion();
					image->UpdateOutputData();

				

					typedef itk::ImageRegionConstIteratorWithIndex<TImage> IteratorType;
					IteratorType it(image, polygonRegion);
					


					for (it.GoToBegin(); !it.IsAtEnd(); ++it)
					{

						itk::ContinuousIndex<double, 2> point;
						image->TransformIndexToPhysicalPoint(it.GetIndex(), point);
						bool isInterior = false;
						if ( exteriorRing->IsInside(point))
						{
 							PolygonListPointerType interiorRings = datanode->GetPolygonInteriorRings();
						
							for (typename PolygonListType::Iterator interiorRing = interiorRings->Begin();interiorRing != interiorRings->End();
							++interiorRing)
							{
								if ( interiorRing.Get()->IsInside(point)  )
								{
									isInterior = true;
								}
							}

						
							

							if(!isInterior ){

								ListSample->PushBack(it.Get());
								ListCropTypeLabel->PushBack(datanode->GetFieldAsInt(_CropTypeClassKey));
								ListBinaryCropLabel->PushBack(datanode->GetFieldAsInt(_IsCropClassKey));
								CoordinateType coor;
						  		coor[0] = it.GetIndex()[0];
						  		coor[1] = it.GetIndex()[1];

						  		CoordinatesList ->PushBack(coor);
								if ( _ClassesLabels.find(datanode->GetFieldAsInt(_CropTypeClassKey)) == _ClassesLabels.end() )
							    {
									// not found
									ClassLabelType label = datanode->GetFieldAsInt(_CropTypeClassKey);
									std::string labelName = datanode->GetFieldAsString(_LabelClassKey);
									_ClassesLabels.insert ( std::pair<ClassLabelType, std::string>(label,labelName) );
							     }


							}


						}
					}



			}
		  }
		}
	}
}


template <class TImage, class TVectorData>
int ExtractSampleFromVectorFile<TImage, TVectorData>::GetNumberOfClasses ( )
{
	_NumberOfClasses = 	_ClassesLabels.size();
  return (this->_NumberOfClasses);
};


template <class TImage, class TVectorData>
std::map< unsigned, std::string>*  ExtractSampleFromVectorFile<TImage, TVectorData>::GetInfoLabelsOfClasses ( )
{
  return (&(this->_ClassesLabels));
};


template <class TImage, class TVectorData>
std::map< unsigned int, unsigned int>*  ExtractSampleFromVectorFile<TImage, TVectorData>::GetInfoSamplesOfClasses ( )
{
  return (&(this->_ClassesSize));
};


template <class TImage, class TVectorData>
void ExtractSampleFromVectorFile<TImage, TVectorData>::PrintInformation( )
{


	 std::cout<< " The description of both classes is : " <<std::endl;
	 std::map<ClassLabelType, unsigned int>::const_iterator itmap;

	 for (itmap = _ClassesSize.begin(); itmap != _ClassesSize.end(); ++itmap)
	 {
		 std::cout<< " 		Class "<< _ClassesLabels[itmap->first] << " with Label  " << itmap->first << " has " << itmap->second << " Samples"<<std::endl;
	 }

};

// Get the sample list
template <class TImage, class TVectorData>
typename ExtractSampleFromVectorFile<TImage, TVectorData>::ListSampleType*
ExtractSampleFromVectorFile<TImage, TVectorData>
::GetListSample()
{
  return dynamic_cast<ListSampleType*>(this->itk::ProcessObject::GetOutput(0));
}
// Get the crop type labels
template <class TImage, class TVectorData>
typename ExtractSampleFromVectorFile<TImage, TVectorData>::ListCropLabelType*
ExtractSampleFromVectorFile<TImage, TVectorData>
::GetCropTypeListLabel()
{
  return dynamic_cast<ListCropLabelType*>(this->itk::ProcessObject::GetOutput(1));
}

// Get the binary crop labels
template <class TImage, class TVectorData>
typename ExtractSampleFromVectorFile<TImage, TVectorData>::ListBinaryCropLabelType*
ExtractSampleFromVectorFile<TImage, TVectorData>
::GetBinaryCropTypeListLabel()
{
  return dynamic_cast<ListBinaryCropLabelType*>(this->itk::ProcessObject::GetOutput(2));
}


// Get the coordinates
template <class TImage, class TVectorData>
typename ExtractSampleFromVectorFile<TImage, TVectorData>::CoordinateList*
ExtractSampleFromVectorFile<TImage, TVectorData>
::GetCoordinatesList()
{
  return dynamic_cast<CoordinateList*>(this->itk::ProcessObject::GetOutput(3));
}

//The explicit instantiation part
template class ExtractSampleFromVectorFile <otb::VectorImage<float, 2>, otb::VectorData<unsigned int,2> >;



