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
#ifndef __ImageClassificationFilter_h
#define __ImageClassificationFilter_h

#include "itkImageToImageFilter.h"
#include "MachineLearningModel.h"
#include "otbImage.h"

namespace otb
{
    /** \class ImageClassificationFilter
    *  \brief This filter performs the classification of a VectorImage using a Model.
    *
    *  This filter is streamed and threaded, allowing to classify huge images
    *  while fully using several core.
    *
    * \sa Classifier
    * \ingroup Streamed
    * \ingroup Threaded
    *
    * \ingroup OTBSupervised
    */
    template <class TInputImage, class TOutputImage, class TMaskImage = TOutputImage>
    class ITK_EXPORT myImageClassificationFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage>
    {
    
    public:
        /** Standard typedefs */
        typedef myImageClassificationFilter                       Self;
        typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
        typedef itk::SmartPointer<Self>                            Pointer;
        typedef itk::SmartPointer<const Self>                      ConstPointer;

        /** Type macro */
        itkNewMacro(Self);

        /** Creation through object factory macro */
        itkTypeMacro(ImageClassificationFilter, ImageToImageFilter);

        typedef TInputImage                                     InputImageType;
        typedef typename InputImageType::ConstPointer           InputImageConstPointerType;
        typedef typename InputImageType::InternalPixelType      ValueType;

        typedef TMaskImage                                      MaskImageType;
        typedef typename MaskImageType::ConstPointer            MaskImageConstPointerType;
        typedef typename MaskImageType::Pointer                 MaskImagePointerType;

        typedef TOutputImage                                    OutputImageType;
        typedef typename OutputImageType::Pointer               OutputImagePointerType;
        typedef typename OutputImageType::RegionType            OutputImageRegionType;
        typedef typename OutputImageType::PixelType             LabelType;

        typedef myMachineLearningModel<ValueType, LabelType>    ModelType;
        typedef typename ModelType::Pointer                     ModelPointerType;

        typedef otb::Image<double>                              ConfidenceImageType;
        typedef typename ConfidenceImageType::Pointer           ConfidenceImagePointerType;

        /** Output type for Proba */
        typedef otb::VectorImage<double>                        ProbaImageType;
        typedef typename ProbaImageType::Pointer                ProbaImagePointerType;
        typedef std::map<unsigned int, unsigned int>            MapType;
        
        /** Set/Get the model */
        itkSetObjectMacro(Model, ModelType);
        itkGetObjectMacro(Model, ModelType);

        /** Set/Get the default label */
        itkSetMacro(DefaultLabel, LabelType);
        itkGetMacro(DefaultLabel, LabelType);

        /** Set/Get the confidence map flag */
        itkSetMacro(UseConfidenceMap, bool);
        itkGetMacro(UseConfidenceMap, bool);

        /** Set/Get the confidence map flag */
        //itkSetMacro(UseProbaMap, bool);
        itkGetMacro(UseProbaMap, bool);

        void SetUseProbaMap(bool p)
        {
          m_UseProbaMap = p;
          m_Model->SetComputeProba(p);
        }

        void SetUseLabelMap(bool p, MapType LabelMap)
        {
            m_UseLabelMap = p;
            m_Model->SetUseLabelMap(p);
            m_Model->SetLabelMap(LabelMap);
        }
        
        /** Accessors for the number of bands*/
        itkSetMacro(NumberOfOutputBands, unsigned int);
        itkGetConstMacro(NumberOfOutputBands, unsigned int);
        /**
        * If set, only pixels within the mask will be classified.
        * All pixels with a value greater than 0 in the mask, will be classified.
        * \param mask The input mask.
        */
        void SetInputMask(const MaskImageType * mask);

        /**
        * Get the input mask.
        * \return The mask.
        */
        const MaskImageType * GetInputMask(void);

        /**
        * Get the output confidence map
        */
        ConfidenceImageType * GetOutputConfidence(void);

        /**
        * Get the ouput proba map
        */
        ProbaImageType * GetOutputProba(void);

    protected:
        
        /** Constructor */
        myImageClassificationFilter();
        /** Destructor */
        virtual ~myImageClassificationFilter() {}

        /** Threaded generate data */
        virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId);
        /** Before threaded generate data */
        virtual void BeforeThreadedGenerateData();
        /**PrintSelf method */
        virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

        void GenerateOutputInformation()
        {
            Superclass::GenerateOutputInformation();
            // Define the number of output bands
            this->GetOutputProba()->SetNumberOfComponentsPerPixel(m_Model->GetNumberOfClasses());
        }

    private:

        myImageClassificationFilter(const Self &); //purposely not implemented
        void operator =(const Self&); //purposely not implemented

        /** The model used for classification */
        ModelPointerType m_Model;
        
        /** Default label for invalid pixels (when using a mask) */
        LabelType m_DefaultLabel;
        
        /** Flag to produce the confidence map (if the model supports it) */
        bool m_UseConfidenceMap;
        bool m_UseProbaMap;
        unsigned int m_NumberOfOutputBands;
        bool m_UseLabelMap;
    };
} // End namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "ImageClassificationFilter.txx"
#endif

#endif
