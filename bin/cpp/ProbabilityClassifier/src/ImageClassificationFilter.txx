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

#ifndef __ImageClassificationFilter_txx
#define __ImageClassificationFilter_txx

#include "ImageClassificationFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace otb
{
    /**
    * Constructor
    */
    template <class TInputImage, class TOutputImage, class TMaskImage>
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::myImageClassificationFilter()
    {
        this->SetNumberOfIndexedInputs(2);
        this->SetNumberOfRequiredInputs(1);
        m_DefaultLabel = itk::NumericTraits<LabelType>::ZeroValue();

        this->SetNumberOfRequiredOutputs(3);
        this->SetNthOutput(0,TOutputImage::New());
        this->SetNthOutput(1,ConfidenceImageType::New());
        this->SetNthOutput(2,ProbaImageType::New());
        m_UseConfidenceMap = false;
        m_UseProbaMap = false;
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    void
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::SetInputMask(const MaskImageType * mask)
    {
        this->itk::ProcessObject::SetNthInput(1, const_cast<MaskImageType *>(mask));
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    const typename myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::MaskImageType *
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::GetInputMask()
    {
        if (this->GetNumberOfInputs() < 2)
        {
            return 0;
        }

        return static_cast<const MaskImageType *>(this->itk::ProcessObject::GetInput(1));
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    typename myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::ConfidenceImageType *
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::GetOutputConfidence()
    {
        if (this->GetNumberOfOutputs() < 3)
        {
            return 0;
        }
        return static_cast<ConfidenceImageType *>(this->itk::ProcessObject::GetOutput(1));
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    typename myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::ProbaImageType *
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::GetOutputProba()
    {
        if (this->GetNumberOfOutputs() < 3)
        {
            return 0;
        }

        return static_cast<ProbaImageType *>(this->itk::ProcessObject::GetOutput(2));
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    void
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::BeforeThreadedGenerateData()
    {
        if (!m_Model)
        {
            itkGenericExceptionMacro(<< "No model for classification");
        }
    }

    template <class TInputImage, class TOutputImage, class TMaskImage>
    void
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
    {

        // Get the input pointers
        InputImageConstPointerType inputPtr     = this->GetInput();
        MaskImageConstPointerType  inputMaskPtr  = this->GetInputMask();
        OutputImagePointerType     outputPtr    = this->GetOutput(0);
        ConfidenceImagePointerType confidencePtr = this->GetOutputConfidence();
        ProbaImagePointerType      probaPtr = this->GetOutputProba();


        // Progress reporting
        itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

        // Define iterators
        typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
        typedef itk::ImageRegionConstIterator<MaskImageType>  MaskIteratorType;
        typedef itk::ImageRegionIterator<OutputImageType>     OutputIteratorType;
        typedef itk::ImageRegionIterator<ConfidenceImageType> ConfidenceMapIteratorType;
        typedef itk::ImageRegionIterator<ProbaImageType>      ProbaMapIteratorType;

        InputIteratorType inIt(inputPtr, outputRegionForThread);
        OutputIteratorType outIt(outputPtr, outputRegionForThread);

        // Eventually iterate on masks
        MaskIteratorType maskIt;
        
        if (inputMaskPtr)
        {
            maskIt = MaskIteratorType(inputMaskPtr, outputRegionForThread);
            maskIt.GoToBegin();
        }

        // setup iterator for confidence map
        bool computeConfidenceMap(m_UseConfidenceMap && m_Model->HasConfidenceIndex() && !m_Model->GetRegressionMode());
        ConfidenceMapIteratorType confidenceIt;
        
        if (computeConfidenceMap)
        {
            confidenceIt = ConfidenceMapIteratorType(confidencePtr,outputRegionForThread);
            confidenceIt.GoToBegin();
        }

        // setup iterator for proba map
        bool computeProbaMap(m_UseProbaMap && m_Model->HasConfidenceIndex() && !m_Model->GetRegressionMode());

        ProbaMapIteratorType probaIt;

        if(computeProbaMap)
        {
            probaIt = ProbaMapIteratorType(probaPtr,outputRegionForThread);
            probaIt.GoToBegin();
        }

        bool validPoint = true;
        double confidenceIndex = 0.0;

        // Taille du vecteur sinon erreur de segmentation en multithread
        itk::VariableLengthVector<double> probaVector{m_Model->GetNumberOfClasses()};
        probaVector.Fill(0.0);
        
        // Walk the part of the image
        for (inIt.GoToBegin(), outIt.GoToBegin(); !inIt.IsAtEnd() && !outIt.IsAtEnd(); ++inIt, ++outIt)
        {

            // Check pixel validity
            if (inputMaskPtr)
            {
                validPoint = maskIt.Get() > 0;
                ++maskIt;
            }
            
            // If point is valid
            if (validPoint)
            {
                // Classifify
                if (computeConfidenceMap && !computeProbaMap)
                {
                    outIt.Set(m_Model->Predict(inIt.Get(),&confidenceIndex)[0]);
                }
                else if (computeProbaMap)
                {
                    outIt.Set(m_Model->Predict(inIt.Get(), &confidenceIndex, &probaVector)[0]);
                }
                else
                {
                    outIt.Set(m_Model->Predict(inIt.Get())[0]);
                }
            }
            else
            {
                // else, set default value
                outIt.Set(m_DefaultLabel);
                confidenceIndex = 0.0;
            }

            if (computeConfidenceMap)
            {
                confidenceIt.Set(confidenceIndex);
                ++confidenceIt;
            }

            if (computeProbaMap)
            {
                probaIt.Set(probaVector);
                ++probaIt;
            }

            progress.CompletedPixel();
        }

    }

    /**
    * PrintSelf Method
    */
    template <class TInputImage, class TOutputImage, class TMaskImage>
    void
    myImageClassificationFilter<TInputImage, TOutputImage, TMaskImage>
    ::PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
} // End namespace otb
#endif
