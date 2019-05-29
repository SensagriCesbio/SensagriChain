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
    // Prototype //
    typedef itk::VariableLengthVector<float> SampleType;
    SampleType  CalculatePrimitives(SampleType Profile, int NbDates, int type);
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
        
        SampleType pixel;
         
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
                //std::cout << m_Model.InputValueType << std::endl; 
                // Classifify
                
		// Get pixel with good primitive /// 
		pixel =  CalculatePrimitives(inIt.Get(), m_NbDates, m_Primitives);                    

                if (computeConfidenceMap && !computeProbaMap)
                {
                    outIt.Set(m_Model->Predict(pixel, &confidenceIndex)[0]);
                }
                else if (computeProbaMap)
                {
                    // Add primitives features - Modif Ludo//
                    //std::cout << pixel << "|" << pixel.Size() << std::endl;
                    //std::cout << m_NbDates << "," << m_Primitives << std::endl;
                    //std::cout << "**********************************************" << std::endl;
                    outIt.Set(m_Model->Predict(pixel, &confidenceIndex, &probaVector)[0]);
                }
                else
                {
                    outIt.Set(m_Model->Predict(pixel)[0]);
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


  SampleType CalculatePrimitives(SampleType Profile, int NbDates, int type)
  {
    SampleType newProfile;
    float A;
    float B;
    float C;
    if(type==0)
    {
      // 10 Optical bands //
      int NbBands = 10;
      newProfile.SetSize(NbDates*NbBands);
      
      for(int i = 0; i < newProfile.Size(); i++)
      {
        newProfile[i] = Profile[i];
      } 
    }
    if(type==1)
    {
      // OSO case //
      int NbBands = 13;
      newProfile.SetSize(NbDates*NbBands);
      int j = 0;  
      for (int i = 0; i < NbDates*10; i++)
      {   
        newProfile[j] =  static_cast<float>(Profile[i]);
        j++;
        if( (i+1)%10 == 0)
        {
          A = static_cast<float>(Profile[i - 11 + 8]); //B8
          B = static_cast<float>(Profile[i - 11 + 4]); //B4
          if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDVI
          else newProfile[j] = (A-B)/0.00000001;
          j++;
         
          A = static_cast<float>(Profile[i - 11 + 3]); //B3
          B = static_cast<float>(Profile[i - 11 + 8]); //B8
          if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDWI
          else newProfile[j] = (A-B)/0.00000001;
          j++; 
        
          A = 0.0;
          for (int k = 0; k<10; k++) A = A + (static_cast<float>(Profile[i - k])*static_cast<float>(Profile[i - k]));
          newProfile[j] = std::sqrt(A);                 //Brightness
          j++; 
        }
      }
    }

    if(type==2)
    {
      // RE case //
      int NbBands = 12;
      newProfile.SetSize(NbDates*NbBands);
      int j = 0;  
      for (int i = 0; i < NbDates*10; i++)
      {   
        newProfile[j] = static_cast<float>(Profile[i]);
        j++;
        if( (i+1)%10 == 0)
        {
          A = static_cast<float>(Profile[i - 11 + 8]); //B8
          B = static_cast<float>(Profile[i - 11 + 4]); //B4
          C = static_cast<float>(Profile[i - 11 + 2]); //B2
          if ( B != 0) newProfile[j] = (B-C)/A; //PSRI_NIR
          else newProfile[j] = (B-C)/0.00000001;
          j++;
         
          A = static_cast<float>(Profile[i - 11 + 5]); //B5 
          B = static_cast<float>(Profile[i - 11 + 8]); //B8
          if ( B != 0) newProfile[j] = A/B; //CHL_RE
          else newProfile[j] = A/0.00000001;
          j++; 
        }
      }
    }
  
    if(type==3)
    {
      // OSO + RE case //
      int NbBands = 15;
      newProfile.SetSize(NbDates*NbBands);
      int j = 0;  
      for (int i = 0; i < NbDates*10; i++)
      {   
        newProfile[j] = static_cast<float>(Profile[i]);
        j++;
        if( (i+1)%10 == 0)
        {
          A = static_cast<float>(Profile[i - 11 + 8]); //B8
          B = static_cast<float>(Profile[i - 11 + 4]); //B4
          if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDVI
          else newProfile[j] = (A-B)/0.00000001;
          j++;
         
          A = static_cast<float>(Profile[i - 11 + 3]); //B3
          B = static_cast<float>(Profile[i - 11 + 8]); //B8
          if ( (A+B) != 0) newProfile[j] = (A-B)/(A+B); //NDWI
          else newProfile[j] = (A-B)/0.00000001;
          j++; 
        
          A = 0.0;
          for (int k = 0; k<10; k++) A = A + (static_cast<float>(Profile[i - k])*static_cast<float>(Profile[i - k]));
          newProfile[j] = std::sqrt(A);                 //Brightness
          j++; 
  
          A = static_cast<float>(Profile[i - 11 + 8]); //B8
          B = static_cast<float>(Profile[i - 11 + 4]); //B4
          C = static_cast<float>(Profile[i - 11 + 2]); //B2
          if ( B != 0) newProfile[j] = (B-C)/A; //PSRI_NIR
          else newProfile[j] = (B-C)/0.00000001;
          j++;
         
          A = static_cast<float>(Profile[i - 11 + 5]); //B5 
          B = static_cast<float>(Profile[i - 11 + 8]); //B8
          if ( B != 0) newProfile[j] = A/B; //CHL_RE
          else newProfile[j] = A/0.00000001;
          j++; 
        }
      } 
    }

    if(type==4)
    {
      // VH and VV Radar Polarizations //
      int NbBands = 2;
      newProfile.SetSize(NbDates*NbBands);
      
      for(int i = 0; i < newProfile.Size(); i++)
      {
        newProfile[i] = Profile[i];
      } 
    }

    if(type==5)
    {
      // 26 primitive S2AGRI //
      int NbBands = 26;
      newProfile.SetSize(NbBands);
      
      for(int i = 0; i < newProfile.Size(); i++)
      {
        newProfile[i] = Profile[i];
      } 
    }
 
    if(type==6)
    {
      // 26 primitive S2AGRI //
      int NbBands = 26;
      newProfile.SetSize(NbBands);
      
      for(int i = 0; i < newProfile.Size(); i++)
      {
        newProfile[i] = Profile[i];
      } 
    }
    return newProfile;
  }

} // End namespace otb
#endif
