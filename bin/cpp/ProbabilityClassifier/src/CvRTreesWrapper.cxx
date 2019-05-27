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
#include <fstream>


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
	
	for (size_t i=0;i<nclasses;i++)
	{
		labels[i] = data->cat_map->data.i[i];
	}
	//std::cout << data->cat_map->data.i[15]  << "\n";
}

void myCvRTreesWrapper::get_votes(
	const cv::Mat& sample, 
	const cv::Mat& missing,
	myCvRTreesWrapper::VotesVectorType& vote_count) const
{
	vote_count.resize(nclasses);

    myCvRTreesWrapper::OOBeVectorType oob_array;
    // myCvRTreesWrapper::OOBVotesVectorType oob_votes(ntrees, std::vector<float>(nclasses, 0));

    oob_array = this->GetOOBerror();
    // oob_votes = this->GetOOBvotes();
	
	for( int k = 0; k < ntrees; k++ )
	{
		CvDTreeNode* predicted_node = trees[k]->predict( sample, missing );
		
		int class_idx = predicted_node->class_idx;
		CV_Assert( 0 <= class_idx && class_idx < nclasses );
		
		// ++vote_count[class_idx];
		vote_count[class_idx] = vote_count[class_idx] + (1 - oob_array[k]);
        // std::cout << "OOB Vote class: " << oob_votes[k][class_idx] <<"\n";
        // std::cout << "Vote count: " << vote_count[class_idx] <<"\n";
        // vote_count[class_idx] = vote_count[class_idx] + oob_votes[k][class_idx];
	}
}

float myCvRTreesWrapper::predict_margin(
	const cv::Mat& sample, 
	const cv::Mat& missing) const
{
	// Sanity check (division by ntrees later on)
	if(ntrees == 0)
	{
		return 0.;
	}
	
	myCvRTreesWrapper::VotesVectorType classVotes;
	this->get_votes(sample, missing, classVotes);
	
	// We only sort the 2 greatest elements
	std::nth_element(classVotes.begin(), classVotes.begin(), 
		classVotes.end(), std::greater<unsigned int>());

	float margin = static_cast<float>(classVotes[0]-classVotes[1])/ntrees;
	
	return margin;
}

float myCvRTreesWrapper::predict_confidence(
	const cv::Mat& sample, 
	const cv::Mat& missing) const
{
	// Sanity check (division by ntrees later on)
	if(ntrees == 0)
	{
		return 0.;
	}

	myCvRTreesWrapper::VotesVectorType classVotes;

    // std::cout << "Before get_votes (Conf)." << "\n";

	this->get_votes(sample, missing, classVotes);

    // std::cout << "After get_votes (Conf)." << "\n";
	
	float max_votes = *(std::max_element(classVotes.begin(), classVotes.end()));
	float confidence = static_cast<float>(max_votes)/static_cast<float>(ntrees);
	
	return confidence * 1000;
}

myCvRTreesWrapper::PixelType myCvRTreesWrapper::predict_proba(
	const cv::Mat& sample,
	myCvRTreesWrapper::PixelType& probaOut,
	const cv::Mat& missing) const
{
 
	myCvRTreesWrapper::ProbaVectorType probaVector(nclasses, 0.0) ;
	
	if(ntrees == 0)
	{
		return probaOut;
	}

	myCvRTreesWrapper::VotesVectorType classVotes;
	this->get_votes(sample, missing, classVotes);
	
	// Division de chaque element par le nombre d'arbre'
	std::transform(
        classVotes.begin(), 
        classVotes.end(),
		probaVector.begin(),
        std::bind2nd(std::divides<float>(), ntrees)
        );

	size_t idx{0};
	
	for(unsigned i = 0; i < nclasses; i++)
	{
		probaOut[idx] = probaVector[i] * 1000;
		idx++;
	}

	return probaOut;
}

myCvRTreesWrapper::PixelType myCvRTreesWrapper::predict_proba(
	const cv::Mat& sample,
	myCvRTreesWrapper::PixelType& probaOut,
	myCvRTreesWrapper::MapType LabelMap,
	const cv::Mat& missing) const
{
  
	myCvRTreesWrapper::ProbaVectorType probaVector(LabelMap.size(), 0.0) ;
	
	if(ntrees == 0)
	{
		return probaOut;
	}
	
    // std::cout << "Before get_votes (Proba)." << "\n";

	myCvRTreesWrapper::VotesVectorType classVotes;
	this->get_votes(sample, missing, classVotes);

    // std::cout << "After get_votes (Proba)." << "\n";
	
	// Division de chaque element par le nombre d'arbre'
	std::transform(
        classVotes.begin(), 
        classVotes.end(),
        probaVector.begin(),
        std::bind2nd(std::divides<double>(),ntrees));

	size_t idx{0};
	unsigned int label;
	
	for(unsigned int i = 0; i<nclasses; i++)
	{
		label = data->cat_map->data.i[i];
		//std::cout << label <<"  " << LabelMap[label] << "\n";
		idx = LabelMap[label];
		//std::cout << idx <<"\n";
		probaOut[idx] = probaVector[i] * 1000;
		idx++;
	}
	
	return probaOut;
}

bool myCvRTreesWrapper::oob_train( 
	const CvMat* _train_data, 
	int _tflag,
    const CvMat* 
    _responses, 
    const CvMat* _var_idx,
    const CvMat* _sample_idx, 
    const CvMat* _var_type,
    const CvMat* 
    _missing_mask, 
    CvRTParams params )
{
    // std::cout << "Method Train 2" << std::endl;

	clear();

	CvDTreeParams tree_params(
		params.max_depth, 
		params.min_sample_count,
		params.regression_accuracy, 
		params.use_surrogates, 
		params.max_categories,
		params.cv_folds, 
		params.use_1se_rule, 
		false, 
		params.priors);

	data = new CvDTreeTrainData();

    data->set_data( 
    	_train_data, 
    	_tflag, 
    	_responses, 
    	_var_idx,
        _sample_idx,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        _var_type, 
        _missing_mask, 
        tree_params, 
        true);

	int var_count = data->var_count;

	if( params.nactive_vars > var_count )
		params.nactive_vars = var_count;
	else if( params.nactive_vars == 0 )
		params.nactive_vars = (int)sqrt((double)var_count);
	else if( params.nactive_vars < 0 )
		CV_Error( CV_StsBadArg, "<nactive_vars> must be non-negative" );

	// Create mask of active variables at the tree nodes
	active_var_mask = cvCreateMat( 1, var_count, CV_8UC1 );
	
	if( params.calc_var_importance )
	{
		var_importance  = cvCreateMat( 1, var_count, CV_32FC1 );
		cvZero(var_importance);
	}
	{ // initialize active variables mask
		CvMat submask1, submask2;
		CV_Assert( (active_var_mask->cols >= 1) && (params.nactive_vars > 0) && (params.nactive_vars <= active_var_mask->cols) );
		cvGetCols( active_var_mask, &submask1, 0, params.nactive_vars );
		cvSet( &submask1, cvScalar(1) );
		if( params.nactive_vars < active_var_mask->cols )
		{
			cvGetCols( active_var_mask, &submask2, params.nactive_vars, var_count );
			cvZero( &submask2 );
		}
	}

    std::cout << "Grow Forest." << std::endl;

	return this->oob_grow_forest( params.term_crit );
}

bool myCvRTreesWrapper::oob_grow_forest( const CvTermCriteria term_crit )
{
    CvMat* sample_idx_mask_for_tree = 0;
    CvMat* sample_idx_for_tree      = 0;

    const int max_ntrees = term_crit.max_iter;
    const double max_oob_err = term_crit.epsilon;

    std::cout << "max_oob_err: "        << max_oob_err << std::endl;
    std::cout << "term_crit.type: "     << term_crit.type << std::endl;
    std::cout << "CV_TERMCRIT_ITER: "   << CV_TERMCRIT_ITER << std::endl;
    std::cout << "var_importance: "     << var_importance << std::endl;
    std::cout << "max_ntrees: "         << max_ntrees << std::endl;

    const int dims = data->var_count;
    float maximal_response = 0;
    float oob_oa = 0;

    CvMat* oob_sample_votes    = 0;
    CvMat* oob_responses       = 0;

    float* oob_samples_perm_ptr= 0;

    float* samples_ptr     = 0;
    uchar* missing_ptr     = 0;
    float* true_resp_ptr   = 0;

    bool is_oob_or_vimportance = (max_oob_err > 0 && term_crit.type != CV_TERMCRIT_ITER) || var_importance;

    // oob_predictions_sum[i] = sum of predicted values for the i-th sample
    // oob_num_of_predictions[i] = number of summands
    //                            (number of predictions for the i-th sample)
    // initialize these variable to avoid warning C4701
    CvMat oob_predictions_sum = cvMat( 1, 1, CV_32FC1 );
    CvMat oob_num_of_predictions = cvMat( 1, 1, CV_32FC1 );

    nsamples = data->sample_count;
    nclasses = data->get_num_classes();

    std::cout << "#Classes: " << nclasses << std::endl;
    std::cout << "#Samples: " << nsamples << std::endl;

    if ( is_oob_or_vimportance )
    {
        if( data->is_classifier )
        {
            oob_sample_votes = cvCreateMat( nsamples, nclasses, CV_32SC1 );
            cvZero(oob_sample_votes);
        }
        else
        {
            // oob_responses[0,i] = oob_predictions_sum[i]
            //    = sum of predicted values for the i-th sample
            // oob_responses[1,i] = oob_num_of_predictions[i]
            //    = number of summands (number of predictions for the i-th sample)
               
            oob_responses = cvCreateMat( 2, nsamples, CV_32FC1 );
            cvZero(oob_responses);
            cvGetRow( oob_responses, &oob_predictions_sum, 0 );
            cvGetRow( oob_responses, &oob_num_of_predictions, 1 );
        }

        // std::cout << "Before heavy allocation." << std::endl;

        oob_samples_perm_ptr     = (float*)cvAlloc( sizeof(float)*nsamples*dims );
        samples_ptr              = (float*)cvAlloc( sizeof(float)*nsamples*dims );
        missing_ptr              = (uchar*)cvAlloc( sizeof(uchar)*nsamples*dims );
        true_resp_ptr            = (float*)cvAlloc( sizeof(float)*nsamples );

        data->get_vectors( 0, samples_ptr, missing_ptr, true_resp_ptr );

        double minval, maxval;
        CvMat responses = cvMat(1, nsamples, CV_32FC1, true_resp_ptr);
        cvMinMaxLoc( &responses, &minval, &maxval );
        maximal_response = (float)MAX( MAX( fabs(minval), fabs(maxval) ), 0 );
    }

    trees = (CvForestTree**)cvAlloc( sizeof(trees[0])*max_ntrees );
    memset( trees, 0, sizeof(trees[0])*max_ntrees );

    sample_idx_mask_for_tree = cvCreateMat( 1, nsamples, CV_8UC1 );
    sample_idx_for_tree      = cvCreateMat( 1, nsamples, CV_32SC1 );

    /* OBB errors array, one per each tree */

    std::cout << "Allocating OOB Error array ..." << "\n";
    myCvRTreesWrapper::OOBeVectorType oob_array;
    oob_array.resize(max_ntrees);

    std::cout << "Allocating OOB Votes array ..." << "\n";
    myCvRTreesWrapper::OOBVotesVectorType oob_votes(max_ntrees, std::vector<float>(nclasses, 0));

    itk::VariableLengthVector<float> oobVotesOK{(unsigned int)nclasses};
    itk::VariableLengthVector<float> oobVotesAC{(unsigned int)nclasses};
    itk::VariableLengthVector<float> oobNumSamplesPerClass{(unsigned int)nclasses};

    std::ofstream oob_error_file;
    std::ofstream oob_acc_class_file;
    std::ofstream oob_oa_file;

    oob_error_file.open("oob_error.txt");
    oob_acc_class_file.open("oob_accuracy_per_class.txt");
    oob_oa_file.open("oob_oa.txt");

    ntrees = 0;
    while( ntrees < max_ntrees )
    {
        oobVotesOK.Fill(0.0);
        oobVotesAC.Fill(0.0);
        oobNumSamplesPerClass.Fill(0.0);
        // std::cout << "Tree: " << ntrees << "\n";

        int i, oob_samples_count = 0;
        double ncorrect_responses = 0; // used for estimation of variable importance
        CvForestTree* tree = 0;

        cvZero( sample_idx_mask_for_tree );
        
        for(i = 0; i < nsamples; i++ ) //form sample for creation one tree
        {
            int idx = (*rng)(nsamples);

            sample_idx_for_tree->data.i[i] = idx;
            sample_idx_mask_for_tree->data.ptr[idx] = 0xFF;
        }

        // std::cout << "Start training Tree_" << ntrees << "\n";

        trees[ntrees] = new CvForestTree();
        tree = trees[ntrees];
        tree->train( data, sample_idx_for_tree, this );

        // std::cout << "End training Tree_" << ntrees << "\n";

        if ( is_oob_or_vimportance )
        {
            CvMat sample, missing;

            // form array of OOB samples indices and get these samples
            sample   = cvMat( 1, dims, CV_32FC1, samples_ptr );
            missing  = cvMat( 1, dims, CV_8UC1,  missing_ptr );

            oob_error = 0;
            oob_oa = 0;

            int rep = 0;
            for( i = 0; i < nsamples; i++, sample.data.fl += dims, missing.data.ptr += dims )
            {
                CvDTreeNode* predicted_node = 0;
                // check if the sample is OOB
                if( sample_idx_mask_for_tree->data.ptr[i] )
                    continue;

                rep++;

                // predict oob samples
                if( !predicted_node )
                {
                    // std::cout << "Start prediction oob samples ..." << "\n";
                    predicted_node = tree->predict(&sample, &missing, true);
                    // std::cout << "End prediction oob samples ..." << "\n";
                }

                if( !data->is_classifier ) //regression
                {

                    double avg_resp, resp = predicted_node->value;
                    oob_predictions_sum.data.fl[i] += (float)resp;
                    oob_num_of_predictions.data.fl[i] += 1;

                    // compute oob error
                    avg_resp = oob_predictions_sum.data.fl[i]/oob_num_of_predictions.data.fl[i];
                    avg_resp -= true_resp_ptr[i];
                    oob_error += avg_resp*avg_resp;
                    resp = (resp - true_resp_ptr[i])/maximal_response;
                    ncorrect_responses += exp( -resp*resp );
                }
                else //classification
                {
                    double prdct_resp;
                    CvPoint max_loc;
                    CvMat votes;

                    cvGetRow(oob_sample_votes, &votes, i);
                    votes.data.i[predicted_node->class_idx]++;

                    if (predicted_node->value == true_resp_ptr[i])
                    {
                        oobVotesOK[predicted_node->class_idx]++;
                    }

                    // oobNumSamplesPerClass[predicted_node->class_idx]++;

                    for (auto j = 0; j < nclasses; j++ )
                    {
                        if (true_resp_ptr[i] == data->cat_map->data.i[j])
                        {
                            oobNumSamplesPerClass[j]++;
                        }
                    } 

                    // compute oob error
                    cvMinMaxLoc( &votes, 0, 0, 0, &max_loc );

                    prdct_resp = data->cat_map->data.i[max_loc.x];
                    oob_error += (fabs(prdct_resp - true_resp_ptr[i]) < FLT_EPSILON) ? 0 : 1;

                    ncorrect_responses += cvRound(predicted_node->value - true_resp_ptr[i]) == 0;
                }

                oob_samples_count++;     
            }

            // std::cout << "Tree[ " << ntrees << " ]: Counter = " << rep << "\n";

            
            if( oob_samples_count > 0 )
            {
                oob_error /= (double)oob_samples_count;
                oob_oa = (double)(ncorrect_responses/oob_samples_count);
            }
                

            for (auto k = 0; k < nclasses; k++ )
            {
                // oobVotesAC[k] = oobVotesOK[k] / oobNumSamplesPerClass[k];
                oobVotesAC[k] = oobNumSamplesPerClass[k];
                oob_votes[ntrees][k] = oobVotesAC[k];
            } 

            oob_array[ntrees] = oob_error;

            // estimate variable importance
            if( var_importance && oob_samples_count > 0 )
            {
                int m;

                memcpy( oob_samples_perm_ptr, samples_ptr, dims*nsamples*sizeof(float));
                for( m = 0; m < dims; m++ )
                {
                    double ncorrect_responses_permuted = 0;
                    // randomly permute values of the m-th variable in the oob samples
                    float* mth_var_ptr = oob_samples_perm_ptr + m;

                    for( i = 0; i < nsamples; i++ )
                    {
                        int i1, i2;
                        float temp;

                        if( sample_idx_mask_for_tree->data.ptr[i] ) //the sample is not OOB
                            continue;
                        i1 = (*rng)(nsamples);
                        i2 = (*rng)(nsamples);
                        CV_SWAP( mth_var_ptr[i1*dims], mth_var_ptr[i2*dims], temp );

                        // turn values of (m-1)-th variable, that were permuted
                        // at the previous iteration, untouched
                        if( m > 1 )
                            oob_samples_perm_ptr[i*dims+m-1] = samples_ptr[i*dims+m-1];
                    }

                    // predict "permuted" cases and calculate the number of votes for the
                    // correct class in the variable-m-permuted oob data
                    sample  = cvMat( 1, dims, CV_32FC1, oob_samples_perm_ptr );
                    missing = cvMat( 1, dims, CV_8UC1, missing_ptr );
                    for( i = 0; i < nsamples; i++,
                        sample.data.fl += dims, missing.data.ptr += dims )
                    {
                        double predct_resp, true_resp;

                        if( sample_idx_mask_for_tree->data.ptr[i] ) //the sample is not OOB
                            continue;

                        predct_resp = tree->predict(&sample, &missing, true)->value;
                        true_resp   = true_resp_ptr[i];
                        if( data->is_classifier )
                            ncorrect_responses_permuted += cvRound(true_resp - predct_resp) == 0;
                        else
                        {
                            true_resp = (true_resp - predct_resp)/maximal_response;
                            ncorrect_responses_permuted += exp( -true_resp*true_resp );
                        }
                    }
                    var_importance->data.fl[m] += (float)(ncorrect_responses
                        - ncorrect_responses_permuted);
                }
            }
        }
        ntrees++;

        if( term_crit.type != CV_TERMCRIT_ITER && oob_error < max_oob_err )
            break;


        oob_error_file << oob_error << std::endl;
        oob_oa_file << oob_oa << std::endl;
   
        for ( auto i = 0; i < nclasses; ++i )
        {
            // if (i == nclasses - 1)
            // {
            //     oob_acc_class_file << (unsigned int)(oobVotesAC[i] * 1000) << std::endl;
            // }
            // else
            // {
            //     oob_acc_class_file << (unsigned int)(oobVotesAC[i] * 1000) << ",";
            // }
            if (i == nclasses - 1)
            {
                oob_acc_class_file << (unsigned int)(oobVotesAC[i]) << std::endl;
            }
            else
            {
                oob_acc_class_file << (unsigned int)(oobVotesAC[i]) << ",";
            }
        }

    }

    oob_error_file.close();
    oob_oa_file.close();
    oob_acc_class_file.close();

    std::cout << "Setting OOB Error array ..." << "\n";
    this->SetOOBerror(oob_array);
    this->SetOOBvotes(oob_votes);

    if( var_importance )
    {
        for ( int vi = 0; vi < var_importance->cols; vi++ )
                var_importance->data.fl[vi] = ( var_importance->data.fl[vi] > 0 ) ?
                    var_importance->data.fl[vi] : 0;
        cvNormalize( var_importance, var_importance, 1., 0, CV_L1 );
    }

    cvFree( &oob_samples_perm_ptr );
    cvFree( &samples_ptr );
    cvFree( &missing_ptr );
    cvFree( &true_resp_ptr );

    cvReleaseMat( &sample_idx_mask_for_tree );
    cvReleaseMat( &sample_idx_for_tree );

    cvReleaseMat( &oob_sample_votes );
    cvReleaseMat( &oob_responses );

    return true;
}

bool myCvRTreesWrapper::oob_train( 
	const cv::Mat& _train_data, 
	int _tflag,
    const cv::Mat& _responses, 
    const cv::Mat& _var_idx,
    const cv::Mat& _sample_idx, 
    const cv::Mat& _var_type,
    const cv::Mat& _missing_mask, 
    CvRTParams _params )
{
    CvMat tdata = _train_data, responses = _responses, vidx = _var_idx,
    sidx = _sample_idx, vtype = _var_type, mmask = _missing_mask;

    // std::cout << "Method Train 1" << std::endl;
    
    return this->oob_train(
    	&tdata, 
    	_tflag, 
    	&responses, 
    	vidx.data.ptr ? &vidx : 0,
    	sidx.data.ptr ? &sidx : 0, 
    	vtype.data.ptr ? &vtype : 0,
        mmask.data.ptr ? &mmask : 0, 
        _params);
}

}

