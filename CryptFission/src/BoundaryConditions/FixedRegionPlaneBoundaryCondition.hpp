#ifndef FIXEDREGIONPLANEBOUNDARYCONDITION_HPP_
#define FIXEDREGIONPLANEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A plane cell boundary coundition with a fixed region, derived from
 * PlaneBoundaryCondition.
 */
template<unsigned DIM>
class FixedRegionPlaneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /**
     * A point on the boundary plane.
     */
    c_vector<double, DIM> mPointOnPlane;

    /**
     * The outward-facing unit normal vector to the boundary plane.
     */
    c_vector<double, DIM> mNormalToPlane;

    /**
     * Whether to jiggle the cells on the bottom surface, initialised to false
     * in the constructor.
     */
    bool mUseJiggledNodesOnPlane;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mUseJiggledNodesOnPlane;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     */
    FixedRegionPlaneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                           c_vector<double, DIM> point,
                           c_vector<double, DIM> normal);

    /**
     * @return #mPointOnPlane.
     */
    const c_vector<double, DIM>& rGetPointOnPlane() const;

    /**
     * @return #mNormalToPlane.
     */
    const c_vector<double, DIM>& rGetNormalToPlane() const;

    /**
     * Set method for mUseJiggledNodesOnPlane
     *
     * @param useJiggledNodesOnPlane whether to jiggle the nodes on the surface of the plane, can help stop overcrowding on plane.
     */
    void SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane);

    /** @return #mUseJiggledNodesOnPlane. */
    bool GetUseJiggledNodesOnPlane();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedRegionPlaneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const FixedRegionPlaneBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> point = t->rGetPointOnPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << point[i];
    }
    c_vector<double, DIM> normal = t->rGetNormalToPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << normal[i];
    }
}

/**
 * De-serialize constructor parameters and initialize a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, FixedRegionPlaneBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> point;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> point[i];
    }
    c_vector<double, DIM> normal;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> normal[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedRegionPlaneBoundaryCondition<DIM>(p_cell_population, point, normal);
}
}
} // namespace ...

#endif /*FIXEDREGIONPLANEBOUNDARYCONDITION_HPP_*/
