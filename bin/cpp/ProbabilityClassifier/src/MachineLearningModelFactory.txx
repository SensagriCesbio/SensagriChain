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
#include "MachineLearningModelFactory.h"
#include "otbConfigure.h"

#ifdef OTB_USE_OPENCV
//#include "otbKNearestNeighborsMachineLearningModelFactory.h" //saute
#include "RandomForestsMachineLearningModelFactory.h" // change
//#include "otbSVMMachineLearningModelFactory.h" // saute
//#include "otbBoostMachineLearningModelFactory.h" //saute/
//#include "otbNeuralNetworkMachineLearningModelFactory.h" //saute
//#include "otbNormalBayesMachineLearningModelFactory.h" // saute
//#include "otbDecisionTreeMachineLearningModelFactory.h" // saute
//#include "otbGradientBoostedTreeMachineLearningModelFactory.h" // saute
#endif
//#ifdef OTB_USE_LIBSVM
//#include "otbLibSVMMachineLearningModelFactory.h" // saute
//#endif

#include "itkMutexLockHolder.h"


namespace otb
{
template <class TInputValue, class TOutputValue>
typename myMachineLearningModel<TInputValue,TOutputValue>::Pointer
myMachineLearningModelFactory<TInputValue,TOutputValue>
::CreateMachineLearningModel(const std::string& path, FileModeType mode)
{
  RegisterBuiltInFactories();

  std::list<MachineLearningModelTypePointer> possiblemyMachineLearningModel;
  std::list<LightObject::Pointer> allobjects =
    itk::ObjectFactoryBase::CreateAllInstance("myMachineLearningModel");
  for(std::list<LightObject::Pointer>::iterator i = allobjects.begin();
      i != allobjects.end(); ++i)
    {
    myMachineLearningModel<TInputValue,TOutputValue> * io = dynamic_cast<myMachineLearningModel<TInputValue,TOutputValue>*>(i->GetPointer());
    if(io)
      {
      possiblemyMachineLearningModel.push_back(io);
      }
    else
      {
      std::cerr << "Error myMachineLearningModel Factory did not return an myMachineLearningModel: "
                << (*i)->GetNameOfClass()
                << std::endl;
      }
    }
for(typename std::list<MachineLearningModelTypePointer>::iterator k = possiblemyMachineLearningModel.begin();
      k != possiblemyMachineLearningModel.end(); ++k)
    {
      if( mode == ReadMode )
      {
      if((*k)->CanReadFile(path))
        {
        return *k;
        }
      }
    else if( mode == WriteMode )
      {
      if((*k)->CanWriteFile(path))
        {
        return *k;
        }

      }
    }
  return 0;
}

template <class TInputValue, class TOutputValue>
void
myMachineLearningModelFactory<TInputValue,TOutputValue>
::RegisterBuiltInFactories()
{
  itk::MutexLockHolder<itk::SimpleMutexLock> lockHolder(mutex);
  
  //#ifdef OTB_USE_LIBSVM
  //RegisterFactory(LibSVMmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //#endif

#ifdef OTB_USE_OPENCV
  RegisterFactory(myRandomForestsMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(SVMmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(BoostmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(NeuralNetworkmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(NormalBayesmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(DecisionTreemyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(GradientBoostedTreemyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
  //RegisterFactory(KNearestNeighborsmyMachineLearningModelFactory<TInputValue,TOutputValue>::New());
#endif

}

template <class TInputValue, class TOutputValue>
void
myMachineLearningModelFactory<TInputValue,TOutputValue>
::RegisterFactory(itk::ObjectFactoryBase * factory)
{
  // Unregister any previously registered factory of the same class
  // Might be more intensive but static bool is not an option due to
  // ld error.
  itk::ObjectFactoryBase::UnRegisterFactory(factory);
  itk::ObjectFactoryBase::RegisterFactory(factory);
}

} // end namespace otb
