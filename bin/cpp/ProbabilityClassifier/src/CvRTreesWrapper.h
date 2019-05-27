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
#ifndef __CvRTreesWrapper_h
#define __CvRTreesWrapper_h

#include "otbOpenCVUtils.h"
#include <vector>
#include "otbVectorImage.h"
#include "itkVariableLengthVector.h"
namespace otb
{

/** \class CvRTreesWrapper
 * \brief Wrapper for OpenCV Random Trees
 *
 * \ingroup OTBSupervised
 */

class CV_EXPORTS_W myCvRTreesWrapper : public CvRTrees
{

public:

// typedef std::vector<unsigned int>            VotesVectorType;
typedef std::vector<float>              		VotesVectorType;
typedef std::vector<float>                    	ProbaVectorType;
typedef std::vector<float>                     	OOBeVectorType;
typedef std::vector< std::vector<float> > 		OOBVotesVectorType;
typedef itk::VariableLengthVector<double>      	PixelType;
typedef std::map<unsigned int, unsigned int>   	MapType;

myCvRTreesWrapper();
virtual ~myCvRTreesWrapper();

/** Get the number of classes from the model saved */
unsigned int GetNclasses();

/** Get the Labels from the model saved */
void GetLabels(std::vector<unsigned int>& labels);

unsigned int nbClasses;
MapType MapLabel;
OOBeVectorType oobErrors;
OOBVotesVectorType oobVotesVector;

/** Compute the number of votes for each class. */
void get_votes(
	const cv::Mat& sample, 
	const cv::Mat& missing,
	VotesVectorType& vote_count) const;

/** Predict the confidence of the classifcation by computing the proportion
  of trees which voted for the majority class.
*/
float predict_confidence(
	const cv::Mat& sample, 
	const cv::Mat& missing = cv::Mat()) const;

/** Predict the confidence margin of the classifcation by computing the 
  difference in votes between the first and second most voted classes.
  This measure is preferred to the proportion of votes of the majority
  class, since it provides information about the conflict between the
  most likely classes.
*/
float predict_margin(
	const cv::Mat& sample, 
	const cv::Mat& missing = cv::Mat()) const;

PixelType predict_proba(
	const cv::Mat& sample,
	myCvRTreesWrapper::PixelType& probaOut,
	const cv::Mat& missing = cv::Mat()) const;

PixelType predict_proba(
	const cv::Mat& sample,
	myCvRTreesWrapper::PixelType& probaOut,
	MapType LabelMap,
	const cv::Mat& missing = cv::Mat()) const;


bool oob_train(
	const cv::Mat& trainData, 
	int tflag,
	const cv::Mat& responses, 
	const cv::Mat& varIdx = cv::Mat(),
	const cv::Mat& sampleIdx = cv::Mat(), 
	const cv::Mat& varType = cv::Mat(),
	const cv::Mat& missingDataMask = cv::Mat(),
	CvRTParams params = CvRTParams() );  

bool oob_train( 
	const CvMat* _train_data, 
	int _tflag,
	const CvMat* 
	_responses, 
	const CvMat* _var_idx,
	const CvMat* _sample_idx, 
	const CvMat* _var_type,
	const CvMat* 
	_missing_mask, 
	CvRTParams params );

bool oob_grow_forest( const CvTermCriteria term_crit );

void SetOOBerror(const OOBeVectorType& obb_array)
{
	this->oobErrors = obb_array;
}

OOBeVectorType GetOOBerror() const
{
	return this->oobErrors;
}

void SetOOBvotes(const OOBVotesVectorType& obb_array)
{
	this->oobVotesVector = obb_array;
}

OOBVotesVectorType GetOOBvotes() const
{
	return this->oobVotesVector;
}

};

}

#endif
