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
#ifndef __MachineLearningModel_txx
#define __MachineLearningModel_txx

#include "MachineLearningModel.h"

namespace otb
{

  template <class TInputValue, class TOutputValue, class TConfidenceValue>
  myMachineLearningModel<TInputValue,TOutputValue,TConfidenceValue>
::myMachineLearningModel() :
  m_RegressionMode(false),
  m_IsRegressionSupported(false),
  m_ConfidenceIndex(false),
  m_ComputeProba(false),
  m_UseLabelMap(false)
{}


  template <class TInputValue, class TOutputValue, class TConfidenceValue>
myMachineLearningModel<TInputValue,TOutputValue,TConfidenceValue>
::~myMachineLearningModel()
{}

template <class TInputValue, class TOutputValue, class TConfidenceValue>
void
myMachineLearningModel<TInputValue,TOutputValue,TConfidenceValue>
::SetRegressionMode(bool flag)
{
  if (flag && !m_IsRegressionSupported)
    {
    itkGenericExceptionMacro(<< "Regression mode not implemented.");
    }
  if (m_RegressionMode != flag)
    {
    m_RegressionMode = flag;
    this->Modified();
    }
}

template <class TInputValue, class TOutputValue, class TConfidenceValue>
void
myMachineLearningModel<TInputValue,TOutputValue,TConfidenceValue>
::PredictAll()
{
  TargetListSampleType * targets = this->GetTargetListSample();
  targets->Clear();

  for(typename InputListSampleType::ConstIterator sIt = this->GetInputListSample()->Begin();
      sIt!=this->GetInputListSample()->End(); ++sIt)
    {
    targets->PushBack(this->Predict(sIt.GetMeasurementVector()));
    }
}

template <class TInputValue, class TOutputValue, class TConfidenceValue>
void
myMachineLearningModel<TInputValue,TOutputValue,TConfidenceValue>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  // Call superclass implementation
  Superclass::PrintSelf(os,indent);
}
}

#endif
