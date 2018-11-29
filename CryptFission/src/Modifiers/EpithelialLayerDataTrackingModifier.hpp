#ifndef EPITHELIALLAYERDATATRACKINGMODIFIER_HPP_
#define EPITHELIALLAYERDATATRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"


template<unsigned DIM>
class EpithelialLayerDataTrackingModifier : public AbstractCellBasedSimulationModifier<DIM>
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
    }

    //Output file for data
    out_stream mpEpithelialLayerDataFile;

protected:


    /** An output stream for writing data. */
    out_stream mpOutStream;

public:

    /**
     * Constructor.
     * @param a vector containing the domain's specifications, length/width/height etc
     */
    EpithelialLayerDataTrackingModifier();

    /**
     * Destructor.
     */
    virtual ~EpithelialLayerDataTrackingModifier();


    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the volume of each cell in the population and store these in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /* Overridden UpdateAtEndOfSolve() method */
    virtual void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /* Get the total volume of the box of cells */
    double CalculateTotalVolume(AbstractCellPopulation<DIM>& rCellPopulation);

    /* Get the total volume of the ring of cells */
    double CalculateRingVolume(AbstractCellPopulation<DIM>& rCellPopulation);

    /* Get the width of the cell population, accounting for ghost nodes */
    c_vector<double,DIM> CalculateCellPopulationWidth(AbstractCellPopulation<DIM>& rCellPopulation);

    /*Get the perimeter of the ring of cells */
    double CalculateRingPerimeter(AbstractCellPopulation<DIM>& rCellPopulation);

    /* Count the cell proliferative types*/
    std::vector<unsigned> CountCellProliferativeTypes(AbstractCellPopulation<DIM>& rCellPopulation);

    /* Count the mutation states*/
    std::vector<unsigned> CountCellMutationStates(AbstractCellPopulation<DIM>& rCellPopulation);

    std::vector<unsigned> GetCellsInRingInOrder(AbstractCellPopulation<DIM>& rCellPopulation);

    /*Method for compiling all the data */
    void CalculateModifierData(AbstractCellPopulation<DIM>& rCellPopulation);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EpithelialLayerDataTrackingModifier)

#endif /*EPITHELIALLAYERDATATRACKINGMODIFIER_HPP_*/
