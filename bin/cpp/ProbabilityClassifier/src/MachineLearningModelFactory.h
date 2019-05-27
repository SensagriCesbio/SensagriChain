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
#ifndef __MachineLearningModelFactory_h
#define __MachineLearningModelFactory_h

#include "MachineLearningModel.h"
#include "otbMachineLearningModelFactoryBase.h"

namespace otb
{
/** \class MachineLearningModelFactory
 * \brief Creation of object instance using object factory.
 *
 * \ingroup OTBSupervised
 */
template <class TInputValue, class TOutputValue>
class ITK_EXPORT myMachineLearningModelFactory : public MachineLearningModelFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef myMachineLearningModelFactory                Self;
  typedef itk::Object           Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Class Methods used to interface with the registered factories */

  /** Run-time type information (and related methods). */
  itkTypeMacro(myMachineLearningModelFactory, itk::Object);

  /** Convenient typedefs. */
  typedef myMachineLearningModel<TInputValue,TOutputValue> MachineLearningModelType;
  typedef typename MachineLearningModelType::Pointer MachineLearningModelTypePointer;

  /** Mode in which the files is intended to be used */
  typedef enum { ReadMode, WriteMode } FileModeType;

  /** Create the appropriate MachineLearningModel depending on the particulars of the file. */
  static MachineLearningModelTypePointer CreateMachineLearningModel(const std::string& path, FileModeType mode);

protected:
  myMachineLearningModelFactory();
  ~myMachineLearningModelFactory();

private:
  myMachineLearningModelFactory(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  /** Register Built-in factories */
  static void RegisterBuiltInFactories();

  /** Register a single factory, ensuring it has not been registered
    * twice */
  static void RegisterFactory(itk::ObjectFactoryBase * factory);

};

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "MachineLearningModelFactory.txx"
#endif

#endif
