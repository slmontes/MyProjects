/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "StochasticTargetProportionBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "Exception.hpp"

StochasticTargetProportionBasedCellCycleModel::StochasticTargetProportionBasedCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mTargetProportion(DOUBLE_UNSET),
	  mCC_Scale(DOUBLE_UNSET)
{
}

StochasticTargetProportionBasedCellCycleModel::~StochasticTargetProportionBasedCellCycleModel()
{
}

void StochasticTargetProportionBasedCellCycleModel::SetTargetProportion(double targetProportion)
{
	mTargetProportion = targetProportion;
}

double StochasticTargetProportionBasedCellCycleModel::GetTargetProportion()
{
	return mTargetProportion;
}

AbstractCellCycleModel* StochasticTargetProportionBasedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
	StochasticTargetProportionBasedCellCycleModel* p_model = new StochasticTargetProportionBasedCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetTargetProportion(mTargetProportion);
    p_model->SetCellCycleLengthScale(mCC_Scale);

    return p_model;
}

void StochasticTargetProportionBasedCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double cc_scaling = GetCellCycleLengthScale();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 2*p_gen->ranf(); // U[0,2]
    }
    //If we have a transit cell, i.e. stem or Paneth cell
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (mpCell->GetMutationState()->IsType<PanethCellMutationState>()))
    {
        mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf())*cc_scaling; // U[0,2]
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && (!mpCell->GetMutationState()->IsType<PanethCellMutationState>()))
	{
		mG1Duration = (GetTransitCellG1Duration() + 2*p_gen->ranf()); // U[0,2]
	}
    //if we have a non-epithelial cell
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void StochasticTargetProportionBasedCellCycleModel::InitialiseDaughterCell()
{

    /*
     * We set the daughter cell to have a probability mTargetProportion of being a stem cell or
     * a paneth cell
     */

	//Get the target proportion probability
	double target_proportion = GetTargetProportion();

	//Get random number
	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
	double uniform_random_number = p_gen->ranf();

	//Set daughter cell to be a transit cell
	boost::shared_ptr<AbstractCellProperty> p_transit_type =
			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
	mpCell->SetCellProliferativeType(p_transit_type);

	if (uniform_random_number < target_proportion) // set probability of becoming a Paneth cell
	{
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
		mpCell->SetMutationState(p_wildtype_state);
	}
	else
	{
		boost::shared_ptr<AbstractCellProperty> p_paneth_state =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
		mpCell->SetMutationState(p_paneth_state);
	}

	AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

void StochasticTargetProportionBasedCellCycleModel::SetCellCycleLengthScale(double cc_scale){

	mCC_Scale = cc_scale;

}

double StochasticTargetProportionBasedCellCycleModel::GetCellCycleLengthScale(){

	return mCC_Scale;

}


void StochasticTargetProportionBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<TargetProportion>" << mTargetProportion << "</TargetProportion>\n";

    // Nothing to output, so just call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticTargetProportionBasedCellCycleModel)
