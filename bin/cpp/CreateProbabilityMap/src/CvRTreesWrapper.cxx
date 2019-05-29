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
#include "CvRTreesWrapper.h"
#include <algorithm>


namespace otb
{
  
myCvRTreesWrapper::myCvRTreesWrapper(){}

myCvRTreesWrapper::~myCvRTreesWrapper(){}
unsigned int myCvRTreesWrapper::GetNclasses()
{
  
  return nclasses;
}

void myCvRTreesWrapper::GetLabels(std::vector<unsigned int>& labels)
{
  labels.resize(nclasses);
  for (size_t i=0;i<nclasses;i++){
    labels[i] = data->cat_map->data.i[i];
  }
  //std::cout << data->cat_map->data.i[15]  << "\n";
}

void myCvRTreesWrapper::get_votes(const cv::Mat& sample, 
                                const cv::Mat& missing,
                                myCvRTreesWrapper::VotesVectorType& vote_count) const
{
  vote_count.resize(nclasses);
  for( int k = 0; k < ntrees; k++ )
    {
    CvDTreeNode* predicted_node = trees[k]->predict( sample, missing );
    int class_idx = predicted_node->class_idx;
    CV_Assert( 0 <= class_idx && class_idx < nclasses );
    ++vote_count[class_idx];
    }
}

float myCvRTreesWrapper::predict_margin(const cv::Mat& sample, 
                                      const cv::Mat& missing) const
{
  // Sanity check (division by ntrees later on)
  if(ntrees == 0)
    {
    return 0.;
    }
  std::vector<unsigned int> classVotes;
  this->get_votes(sample, missing, classVotes);
// We only sort the 2 greatest elements
  std::nth_element(classVotes.begin(), classVotes.begin(), 
                   classVotes.end(), std::greater<unsigned int>());
  float margin = static_cast<float>(classVotes[0]-classVotes[1])/ntrees;
  return margin;
}

float myCvRTreesWrapper::predict_confidence(const cv::Mat& sample, 
                                  const cv::Mat& missing) const
{
  // Sanity check (division by ntrees later on)
  if(ntrees == 0)
    {
    return 0.;
    }
  std::vector<unsigned int> classVotes;
  this->get_votes(sample, missing, classVotes);
  unsigned int max_votes = *(std::max_element(classVotes.begin(), 
                                              classVotes.end()));
  float confidence = static_cast<float>(max_votes)/ntrees;
  return confidence*1000;
}

myCvRTreesWrapper::PixelType   myCvRTreesWrapper::predict_proba(
                                  const cv::Mat& sample,
				  myCvRTreesWrapper::PixelType& probaOut,
				  const cv::Mat& missing) const
{
 
  myCvRTreesWrapper::ProbaVectorType probaVector(nclasses,0.0) ;
  if(ntrees == 0)
    {
      
      return probaOut;
    }
  std::vector<unsigned int> classVotes;
  this->get_votes(sample, missing, classVotes);
  // Division de chaque element par le nombre d'arbre'
  std::transform(classVotes.begin(), classVotes.end(),probaVector.begin(),std::bind2nd(std::divides<double>(),ntrees));
  
  size_t idx{0};
  for(unsigned i = 0; i < nclasses; i++)
    {
      probaOut[idx] = probaVector[i]*1000;
      idx++;
    }
  return probaOut;
}

myCvRTreesWrapper::PixelType   myCvRTreesWrapper::predict_proba(
                                  const cv::Mat& sample,
				  myCvRTreesWrapper::PixelType& probaOut,
				  myCvRTreesWrapper::MapType LabelMap,
				  const cv::Mat& missing) const
{
  
  myCvRTreesWrapper::ProbaVectorType probaVector(LabelMap.size(),0.0) ;
  if(ntrees == 0)
    {
      return probaOut;
    }
  std::vector<unsigned int> classVotes;
  this->get_votes(sample, missing, classVotes);
  // Division de chaque element par le nombre d'arbre'
  std::transform(classVotes.begin(), classVotes.end(),probaVector.begin(),std::bind2nd(std::divides<double>(),ntrees));
  
  size_t idx{0};
  unsigned int label;
  for(unsigned int i = 0; i<nclasses; i++)
    {
      label = data->cat_map->data.i[i];
      //std::cout << label <<"  " << LabelMap[label] << "\n";
      idx = LabelMap[label];
      //std::cout << idx <<"\n";
      probaOut[idx] = probaVector[i]*1000;
      idx++;
    }
  return probaOut;
}
}

