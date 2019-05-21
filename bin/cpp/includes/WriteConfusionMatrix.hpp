#ifndef CONFUSION_MATRIX_HPP_
#define CONFUSION_MATRIX_HPP_

//#include "itkVariableLengthVector.h"
#include <iostream>
//#include <vector>
#include <fstream>

typedef unsigned int LabelType;
typedef itk::Statistics::ListSample<itk::FixedArray<LabelType, 1>> ListLabelSampleType;
typedef otb::ConfusionMatrixCalculator<ListLabelSampleType,ListLabelSampleType> ConfusionMatrixCalculatorType;
typedef ConfusionMatrixCalculatorType::ConfusionMatrixType ConfusionMatrixType;
typedef ConfusionMatrixCalculatorType::MapOfIndicesType MapOfIndicesType;
typedef ConfusionMatrixCalculatorType::ClassLabelType ClassLabelType;


void writeConfusionMatrix( std::string filename,
                                      ConfusionMatrixCalculatorType::Pointer confMatCalc){

    
    MapOfIndicesType::iterator itMapOfIndicesValid, itMapOfIndicesPred;
    ClassLabelType labelValid = 0;
    
    ConfusionMatrixType confusionMatrix = confMatCalc->GetConfusionMatrix();
    MapOfIndicesType mapOfIndicesValid = confMatCalc->GetMapOfIndices();
    
    unsigned int nbClassesPred = mapOfIndicesValid.size();
    
    /////////////////////////////////////////////
    // Filling the 2 headers for the output file
    const std::string commentValidStr = "#Reference labels (rows):";
    const std::string commentPredStr = "#Produced labels (columns):";
    const char separatorChar = ',';
    std::ostringstream ossHeaderValidLabels, ossHeaderPredLabels;
    
    // Filling ossHeaderValidLabels and ossHeaderPredLabels for the output file
    ossHeaderValidLabels << commentValidStr;
    ossHeaderPredLabels << commentPredStr;
    
    itMapOfIndicesValid = mapOfIndicesValid.begin();
    
    while (itMapOfIndicesValid != mapOfIndicesValid.end())
    {
        // labels labelValid of mapOfIndicesValid are already sorted in otbConfusionMatrixCalculator
        labelValid = itMapOfIndicesValid->second;
    
        
        ossHeaderValidLabels << labelValid;
        ossHeaderPredLabels << labelValid;
        
        ++itMapOfIndicesValid;
        
        if (itMapOfIndicesValid != mapOfIndicesValid.end())
        {
            ossHeaderValidLabels << separatorChar;
            ossHeaderPredLabels << separatorChar;
        }
        else
        {
            ossHeaderValidLabels << std::endl;
            ossHeaderPredLabels << std::endl;
        }
    }
    
    std::ofstream  outFile;
    outFile.open(filename.c_str());
    outFile << std::fixed;
    outFile.precision(10);
    
    /////////////////////////////////////
    // Writing the 2 headers
    outFile << ossHeaderValidLabels.str();
    outFile << ossHeaderPredLabels.str();
    /////////////////////////////////////
    
    unsigned int indexLabelValid = 0, indexLabelPred = 0;
    
    for (itMapOfIndicesValid = mapOfIndicesValid.begin(); itMapOfIndicesValid != mapOfIndicesValid.end(); ++itMapOfIndicesValid)
    {
        indexLabelPred = 0;
        
        for (itMapOfIndicesPred = mapOfIndicesValid.begin(); itMapOfIndicesPred != mapOfIndicesValid.end(); ++itMapOfIndicesPred)
        {
            // Writing the confusion matrix (sorted in otbConfusionMatrixCalculator) in the output file
            outFile << confusionMatrix(indexLabelValid, indexLabelPred);
            if (indexLabelPred < (nbClassesPred - 1))
            {
                outFile << separatorChar;
            }
            else
            {
                outFile << std::endl;
            }
            ++indexLabelPred;
        }
        
        ++indexLabelValid;
    }
    
    outFile.close();
    
}



#endif /* IO_CONFUSION_MATRIX_HPP_ */
