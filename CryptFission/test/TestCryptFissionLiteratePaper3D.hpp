/*
 * TestCryptFissionLiteratePaper3D.hpp
 *
 * Created on: 31/01/2019
 * Last modified: 31/01/2019
 * 		Author: Sandra Montes
 */


#ifndef TESTCRYPTFISSIONLITERATEPAPER3D_HPP_
#define TESTCRYPTFISSIONLITERATEPAPER3D_HPP_
/*
 * = Simulation of a Node-Based Organoid =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this test we will generate a node-based cell mesh in the form of an spheroid.
 *
 * == Including header files ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files. The first ones are common to all cell_based Chaste simulations
 which are used in the TestRunningNodeBasedSimulationsTutorial
 */

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)

#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "PetscSetupAndFinalize.hpp" //Not sure it needed but just in case. In the case of Langlands paper they used:
//#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
//#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "UniformCellCycleModel.hpp" //In langlands we used the following line instead
//#include "StochasticTargetProportionBasedCellCycleModel.hpp" //Asymmetric-division-based cell cycle model
#include "GeneralisedLinearSpringForce.hpp" // We can also use the following line
//#include "EpithelialLayerLinearSpringForce3D.hpp" //Spring force law to account for different cell type pairs
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "NodesOnlyMesh.hpp" //Defines the class storing the spatial info of cells
#include "NodeBasedCellPopulation.hpp" //Defines a node-based cell population class
//#include "SphereGeometryBoundaryCondition.hpp" //In case we would like to set a spheroid boundary condition to the system

#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations


/* This set of classes are not part of the core code and were made for the Langlands et al and Almet et al model. */
/*
#include "StochasticTargetProportionBasedCellCycleModel.hpp" //Asymmetric-division-based cell cycle model
#include "PanethCellMutationState.hpp" //Mutation class that defines Paneth cells

#include "EpithelialLayerBasementMembraneForce3D.hpp" //Basement membrane force, as defined in Dunn et al. (2012)
#include "EpithelialLayerAnoikisCellKiller3D.hpp" //Cell killer to remove proliferative cells that have fallen into the lumen

#include "EpithelialLayerDataTrackingModifier3D.hpp" //Modifier for all the necessary data

#include "FixedRegionPlaneBoundaryCondition.hpp" //Boundary condition that fixes cells past a given plane
*/

#include <boost/lexical_cast.hpp>
#include "CellAgesWriter.hpp"

//Additions for Langlands et al update Chaste_2018.1
//#include "WildTypeCellMutationState.hpp"

/*
 * Define the Chaste simulation as a test class. This is how all simulations
 * in Chaste are defined.
 */
class TestCryptFissionLiteratePaper3D : public AbstractCellBasedTestSuite
{

public:

    void TestSpheroidNodeBasedOrganoid() throw(Exception)
    {
	    EXIT_IF_PARALLEL;
//		/* We first set all the simulation parameters. */
//
//		//Simulation time parameters
		double dt = 0.005; //Set dt
		double end_time = 100; //250.0; //Set end time
		double sampling_timestep = 0.5/dt; //Set sampling timestep
//        unsigned start_sim = 1;
//        unsigned num_sims = 1; //It was set to 100 but took to long to run*
//
//		for(unsigned index=start_sim; index < start_sim + num_sims; index++)
//		{
//
//			// Seed the random number generator for each simulation (affects initial placement of cells)
//			RandomNumberGenerator::Instance()->Reseed(index);
//
//			/* Generate the initial mesh of cells. */
			std::vector<Node<3>*> nodes;
			nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
            nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);

            boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
            boost::shared_ptr<AbstractCellProperty> p_paneth_state = CellPropertyRegistry::Instance()->Get<PanethCellMutationState>();

            /*Create a vector of cells and define cell types
             * *How do I include more cell types?
             */

            std::vector<CellPtr> cells;
//            MAKE_PTR(TransitCellProliferativeType, p_stem_type);
            CellsGenerator<UniformCellCycleModel, 3> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);

            //Create cell population
            NodeBasedCellPopulation<3> cell_population(mesh, cells);

			/*
			 * Randomly assign cells in the layer to be Paneth cells.
			 */

			std::vector<unsigned> cells_in_layer; //Initialise vector

			//Obtain the proliferative cells
			for (NodeBasedCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End();
					++cell_iter)
			{
//			    cell_iter->SetCellProliferativeType(p_stem_type);
			    unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
				cells_in_layer.push_back(node_index);
			 }

		    //For each cell in the ring, we draw a random number and assign cells to be stem cells with
			// a probability equal to the target proportion, as defined above.

			for (unsigned i = 0; i < cells_in_layer.size(); i++)
			{
				unsigned node_index = cells_in_layer[i];

				CellPtr cell = cell_population.GetCellUsingLocationIndex(node_index);

				//Randomly generate number
				double random_number = RandomNumberGenerator::Instance()->ranf();

				if(random_number >= 0.4) //Assign cells to be Paneth with 1 - target_proportion
				{
					cell->SetMutationState(p_paneth_state);
				}
			}

//
			/* Define the simulation class. */
			OffLatticeSimulation<3> simulator(cell_population);
//
//			//Set output directory
			simulator.SetOutputDirectory("Node_based_Organoid");
//			simulator.SetSamplingTimestepMultiple(12);
//          simulator.SetEndTime(50.0);
			simulator.SetDt(dt); //Set the timestep dt for force volution
			simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set the sampling timestep multiple for animations
			simulator.SetEndTime(end_time); //Set the number of hours to run the simulation to
//
//			/* Add linear spring force
//			 */
			MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
			simulator.AddForce(p_force);
//
//			/* Add the basement membrane force.
//			MAKE_PTR(EpithelialLayerBasementMembraneForce3D, p_bm_force);
//			p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
//			p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
//			simulator.AddForce(p_bm_force);
//			*/
//
//			/* Add an anoikis-based cell killer.
//			MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller3D, p_anoikis_killer, (&cell_population));
//			simulator.AddCellKiller(p_anoikis_killer);
//			*/
//
//			// Add cell age writer
//			//cell_population.AddCellWriter<CellAgesWriter>();
//
//			/* Run the simulation. */
			simulator.Solve();
//
//
//			SimulationTime::Destroy();
//			SimulationTime::Instance()->SetStartTime(0.0);
//			//RandomNumberGenerator::Destroy();
//
    }
    void TestSpheroid()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedSpheroid");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
	}

/*
 * To visualize the results, open a new terminal and {{{cd}}} to the Chaste directory, then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CryptFissionLiteratePaper/results_from_time_0}}}.
 * You may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the java executable. You should also select the axes equal option.
 */

};

#endif /* TESTCRYPTFISSIONLITERATEPAPER3D_HPP_ */
