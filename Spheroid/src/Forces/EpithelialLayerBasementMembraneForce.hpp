#ifndef EPITHELIALLAYERBASEMENTMEMBRANEFORCE_HPP_
#define EPITHELIALLAYERBASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "NodeBasedCellPopulation.hpp"

#include <cmath>
#include <list>
#include <fstream>
#include <algorithm> //in order to use .erase
#include <math.h> //in order to use acos

/*
 * Author: Sandra M.
 * Created on: 15/03/2019
 * Last modified: 26/03/2019
 */

/**
 * A force class that defines the force due to the basement membrane.
 */

class EpithelialLayerBasementMembraneForce : public AbstractForce<3>
{
private :

    /** Parameter that defines the stiffness of the basement membrane */
    double mBasementMembraneParameter;

    /** Target curvature for the layer of cells */
    double mTargetCurvature;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<3> >(*this);
        archive & mBasementMembraneParameter;
        archive & mTargetCurvature;
    }

public :

    /**
     * Constructor.
     */
	EpithelialLayerBasementMembraneForce();

    /**
     * Destructor.
     */
    ~EpithelialLayerBasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of Target Curvature in epithelial layer */
    void SetTargetCurvature(double targetCurvature = 0.0);

    /* Get method for Target Curvature
     *
     */
    double GetTargetCurvature();

    /* Finding the epithelial nodes in the tissue
     */
    std::vector<unsigned> GetEpithelialNodes(NodeBasedCellPopulation<3>& rCellPopulation);

    /* Finding the neighbours of each epithelial node (only epithelial neighbours)
     */
    std::vector<c_vector<unsigned, 3> > GetNodeNeighbours(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index);

    /* Finding the vectors from the p_node to its epithelial neighbours
     */
    std::vector<c_vector<unsigned, 3> > GetNodetoNeighbourVectors(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index);

    /* Finding the vectors from each epithelial neighbour to its immediate next
     * epithelial neighbour (from the pool of neighbours found before)
     */
    std::vector<c_vector<unsigned, 3> > GetNeighbourtoNeighbourVectors(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index);

    /* Finding the angles between each p_node to its epithelial neighbours vectors
     */
    std::vector<c_vector<unsigned, 3> > GetNodetoNeighbourAngles(NodeBasedCellPopulation<3>& rCellPopulation);

    /* Finding the area of the triangles formed by the network of p_node and
     * its neighbours
     */
    std::vector<c_vector<unsigned, 3> > GetNodetoNeighbourAreas(NodeBasedCellPopulation<3>& rCellPopulation);

    /* Takes the set of angles and areas to calculate the Gaussian curvature
     * at p_node
     */
    double GetGaussianCurvature(NodeBasedCellPopulation<3>& rCellPopulation);

    /* Takes all the ephitelial nodes and adds the BM force that acts in each node
     * Overridden AddForceContribution method
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(NodeBasedCellPopulation<3>& rCellPopulation);

    /*
     * Outputs force Parameters to file
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(EpithelialLayerBasementMembraneForce)

#endif /*EPITHELIALLAYERBASEMENTMEMBRANEFORCE_HPP_*/
