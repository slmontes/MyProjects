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

#ifndef STOCHASTICTARGETPROPORTIONBASEDCELLCYCLEMODEL_HPP_
#define STOCHASTICTARGETPROPORTIONBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * This class contains all the things common to simple generation-based cell cycle
 * models, i.e. models in which the length of cell cycle phases are determined
 * when the cell-cycle model is created, rather than evaluated 'on the fly'
 * by ODEs and suchlike, and in which each cell has a 'generation'.
 *
 * N.B. Whether or not the cell should actually divide may depend on
 * Wnt / Oxygen etc. in subclasses.
 */
class StochasticTargetProportionBasedCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mTargetProportion;
        archive & mCC_Scale;
    }


    /*
     * @param mTargetProportion
     */
    double mTargetProportion;
    double mCC_Scale;

public:

    /**
     * Default constructor - creates an AbstractSimplePhaseBasedCellCycleModel.
     */
    StochasticTargetProportionBasedCellCycleModel();

    virtual double GetAverageTransitCellCycleTime()
    {
        return 1.0;
    }

    virtual double GetAverageStemCellCycleTime()
    {
        return 1.0;
    }

    virtual void SetCellCycleDuration()
    {

    }


    /**
     * Destructor.
     */
    virtual ~StochasticTargetProportionBasedCellCycleModel();

    /*
     * @return target proportion
     */
    double GetTargetProportion();

    /*
     * Method to set target proportion
     */
    void SetTargetProportion(double targetProportion);

    /**
     * Overridden SetG1Duration Method to add stochastic cell cycle times
     */
    void SetG1Duration();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();


    /**
     * Set the new cell's type after division based on probability.
     */
    void InitialiseDaughterCell();

    void SetCellCycleLengthScale(double cc_scale);

    double GetCellCycleLengthScale();
    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StochasticTargetProportionBasedCellCycleModel)

#endif /*STOCHASTICTARGETPROPORTIONBASEDCYCLEMODEL_HPP_*/
