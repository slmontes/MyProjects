/*
 * TestSingleOrganoidSimulation.hpp
 *
 *  Created on: 13 Feb 2018
 *      Author: chastebox
 */

#ifndef PROJECTS_CRYPTFISSIONPLOS2016_TEST_TESTSINGLEORGANOIDSIMULATION_HPP_
#define PROJECTS_CRYPTFISSIONPLOS2016_TEST_TESTSINGLEORGANOIDSIMULATION_HPP_

/*
 * = Simulation of fission in epithelial layer =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this test we show how Chaste can be used to simulate a buckling layer of epithelial cells.
 * Details of the computational model can be found in
 * Langlands et al (2016) "Paneth cell-rich regions separated by a cluster of Lgr5+ cells initiate
 * fission in the intestinal stem cell niche".
 *
 * == Including header files ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files. The first ones are common to all cell_based Chaste simulations
 */

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing

/* The next set of classes are needed specifically for the simulation, which can be found in the core code. */

#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "TransitCellProliferativeType.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel

/* This set of classes are not part of the core code and were made for this model. */

#include "StochasticTargetProportionBasedCellCycleModel.hpp" //Asymmetric-division-based cell cycle model
#include "PanethCellMutationState.hpp" //Mutation class that defines Paneth cells

#include "EpithelialLayerLinearSpringForce.hpp" //Spring force law to account for different cell type pairs
#include "EpithelialLayerBasementMembraneForce.hpp" //Basement membrane force, as defined in Dunn et al. (2012)
#include "EpithelialLayerAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have fallen into the lumen

#include "EpithelialLayerDataTrackingModifier.hpp" //Modifier for all the necessary data

#include "FixedRegionPlaneBoundaryCondition.hpp" //Boundary condition that fixes cells past a given plane
#include <boost/lexical_cast.hpp>
#include "CellAgesWriter.hpp"

#include "WildTypeCellMutationState.hpp"

/*
 * Define the Chaste simulation as a test class. This is how all simulations
 * in Chaste are defined.
 */
class TestSingleOrganoidSimulation : public AbstractCellBasedTestSuite
{
public:
	void TestEpithelialLayerUndergoingFission() throw(Exception)
	{
		/* We first set all the simulation parameters. */

		//Simulation time parameters
		double dt = 0.005; //Set dt
		double end_time = 250.0; //Set end time
		double sampling_timestep = 0.5/dt; //Set sampling timestep


        // THIS IS A HACK WHICH HIJACKS PETSC TO GET ARGUMENTS INTO A TEST!
		CommandLineArguments* p_args = CommandLineArguments::Instance();
		unsigned argc = *(p_args->p_argc); // has the number of arguments, and
		char **argv = *(p_args->p_argv); // is a char** of them.
		std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

			   if (argc != 3)
			   {
				   std::cerr << "TestGenerateCryptSteadyState::Please input two arguments\n"
								"* Cell-type ratio & cell-cycle length scaling\n" << std::flush;
				   return;
			   }


	    /* Set the target proportion for stem cells. Note that the target proportion and
		 * the actual proportion of stem cells will be different, due to cell death and
		 * mechanical heterogeneities.
		 */
		//double target_proportion = 0.3; //This corresponds to the 20% stem cell (80% Paneth cell) case
		// double target_proportion = 0.8; //This corresponds to the 60% stem cell (40% Paneth cell) case.
		double target_proportion = atof(argv[1]);
		double cc_scale = atof(argv[2]);


		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 15.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections

		//Set the stiffness ratio for Paneth cells to stem cells. This is the
		double stiffness_ratio = 4.5;



		//Set the BM force parameters
		double bm_force = 10.0; //Set the basement membrane stiffness
		double target_curvature = 0.2; //Set the target curvature, i.e. how circular the layer wants to be

		//Matrigel domain setup (usually square 24 accross 27 up)
		unsigned cells_across = 36; //Desired width + a few more layers, to avoid the box collapsing
		unsigned cells_up = 40.5; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) = desired height
		unsigned ghosts = 4; //Define a sufficient layer of ghost nodes to avoid infinite tessellations and hence excessively large forces

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a sufficient layer of fixed cells around the boundary.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.5;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = -1.5;

		/* Define the initially circular lumen by centre and radius */
		c_vector<double,2> circle_centre;
		circle_centre(0) = 15.75;
		circle_centre(1) = 15.0;

		double circle_radius = 2.5; //Size of hole
		assert(circle_radius > 0); //Just in case someone does something crazy.

		double ring_radius = circle_radius + 2.0; //Radius of the ring of cells. This isn't the actual radius, just has to be large enough for later
		assert((ring_radius <= cells_across)&&(ring_radius <= cells_up)); //Again, just in case.

		// Set output directory
		std::string arg1_string(argv[1]);
		std::string arg2_string(argv[2]);
		std::string simParams = arg1_string + "_" + arg2_string;
		std::string output_directory = "CryptFissionSingleSimulation_" + simParams;

			// Seed the random number generator for each simulation (affects initial placement of cells)
			RandomNumberGenerator::Instance()->Reseed(time(NULL));

			/* Generate the initial mesh of cells. */
			HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
			MutableMesh<2,2>* p_mesh = generator.GetMesh();

			//Translate mesh appropriately
			p_mesh->Translate(translate_left);
			p_mesh->Translate(translate_down);

			/* Define the lumen as an inner region of ghost nodes. */

			std::vector<unsigned> initial_real_indices = generator.GetCellLocationIndices(); //Obtain the locations of real nodes

			std::vector<unsigned> real_indices; //Vector used to define the locations of non-ghost nodes

			//Sweep over the initial real indices
			for (unsigned i = 0; i < initial_real_indices.size(); i++)
			{
				unsigned cell_index = initial_real_indices[i];
				double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
				double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

				// If the location of the node falls inside the defined lumen region, then it becomes a ghost node.
				if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) > pow(circle_radius,2))
				{
					real_indices.push_back(cell_index);
				}
			}

			/* Define cell types: non-epithelial cells are differentiated, all proliferative epithelial cells
			 * will be transit cells (to define cell cycle duration). In addition, stem cells are assigned a wildtype mutation state,
			 * while Paneth cells are assigned a Paneth cell mutation state.
			 */

			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_paneth_state = CellPropertyRegistry::Instance()->Get<PanethCellMutationState>();
			boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

			//Create vector of cells
			std::vector<CellPtr> cells;

			/*
			 * Create asymmetric-division-based cell cycle for each cell. However, we initially set all cells
			 * to be non-epithelial cells before defining our layer of epithelial cells.
			 */

			for (unsigned i = 0; i<real_indices.size(); i++)
			{
				//Set cell cycle
				StochasticTargetProportionBasedCellCycleModel* p_cycle_model = new StochasticTargetProportionBasedCellCycleModel();
				p_cycle_model->SetTargetProportion(target_proportion); //Set the division parameter
				p_cycle_model->SetCellCycleLengthScale(cc_scale);
				p_cycle_model->SetDimension(2);

				//To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
				// ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf();
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_diff_type); //Set the cell to be differentiated and hence non-epithelial
				p_cell->InitialiseCellCycleModel();

				cells.push_back(p_cell);
			}

			//Create cell population
			MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

			//Create layer of proliferative cells
			for (unsigned i = 0; i < real_indices.size(); i++)
			{
				unsigned cell_index = real_indices[i];
				CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
				double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
				double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];


				/* We 'un-differentiate' any cells adjacent to the lumen into epithelial cells.
				 * We only consider cells within the pre-defined ring radius. Not only does this narrow down
				 * our search, but it also means we won't accidentally consider any cells on the outside and
				 * turn them into epithelial cells.
				 */

				if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(ring_radius,2))
				{
					Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);


					//Iterate over all possible neighbours of the node
					for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
							iter != p_node->ContainingElementsEnd();
							++iter)
					{
						bool element_contains_ghost_nodes = false;

						// Get a pointer to the element
						Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

						// Check whether it's triangulation contains a ghost node
						for (unsigned local_index=0; local_index<3; local_index++)
						{
							unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

							if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
							{
								element_contains_ghost_nodes = true;
								break; 				// This should break out of the inner for loop
							}
						}

						//If a cell has a ghost node as a neighbour, we make it an epithelial cells.
						if(element_contains_ghost_nodes)
						{
							cell_iter->SetCellProliferativeType(p_stem_type);
						}
					}
				}
			}

			/* Iterate again and check that proliferative cells are also attached to non-epithelial
			 * cells. If they are not, remove them from the simulation.
			 */

			for (unsigned i = 0; i < real_indices.size(); i++)
			{
				unsigned cell_index = real_indices[i];
				CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
				double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
				double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

				//Only consider this inside the pre-defined ring to narrow down our search

				if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(ring_radius,2))
				{
					Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);


					//Only iterate over the initial layer of transit cells
					if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false)
					{
						bool element_contains_gel_nodes = false;

						//Iterate over elements (triangles) containing the node
						for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
								iter != p_node->ContainingElementsEnd();
								++iter)
						{
							// Get a pointer to the element (triangle)
							Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

							// Check if its triangulation contains a gel node
							for (unsigned local_index=0; local_index<3; local_index++)
							{
								unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);
								bool is_ghost_node = cell_population.IsGhostNode(nodeGlobalIndex);

								if (is_ghost_node == false) //Make sure we're not dealing with ghost nodes (otherwise this stuff will fail)
								{
									CellPtr p_local_cell = cell_population.GetCellUsingLocationIndex(nodeGlobalIndex);
									if (p_local_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()==true)
									{
										element_contains_gel_nodes = true;
										break; 				// This should break out of the inner for loop
									}
								}
							}
						}

						if(element_contains_gel_nodes == false)
						{
							cell_iter->Kill();
						}
					}
				}
			}

			/*
			 * Randomly assign cells in the layer to be Paneth cells.
			 */

			std::vector<unsigned> cells_in_layer; //Initialise vector

			//Obtain the proliferative cells
			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End();
					++cell_iter)
			{
				unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

				//If the cell is an epithelial cell
				if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
				{
					cells_in_layer.push_back(node_index); //Add the angle and node index
				}
			}

			/* For each cell in the ring, we draw a random number and assign cells to be stem cells with
			 * a probability equal to the target proportion, as defined above.
			 */

			for (unsigned i = 0; i < cells_in_layer.size(); i++)
			{
				unsigned node_index = cells_in_layer[i];

				CellPtr cell = cell_population.GetCellUsingLocationIndex(node_index);

				//Randomly generate number
				double random_number = RandomNumberGenerator::Instance()->ranf();

				if(random_number >= target_proportion) //Assign cells to be Paneth with 1 - target_proportion
				{
					cell->SetMutationState(p_paneth_state);
				}
			}

			//Allow output in Paraview, a program that can be used to visualise Chaste simulations
			cell_population.AddPopulationWriter<VoronoiDataWriter>();

			/* Define the simulation class. */
			OffLatticeSimulation<2> simulator(cell_population);

			//Set output directory
			simulator.SetOutputDirectory(output_directory);

			simulator.SetDt(dt); //Set the timestep dt for force volution
			simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set the sampling timestep multiple for animations
			simulator.SetEndTime(end_time); //Set the number of hours to run the simulation to

			/* We add a modifier class to track relevant cell population numbers and shape measurements.*/
			MAKE_PTR(EpithelialLayerDataTrackingModifier<2>, p_data_tracking_modifier);
			simulator.AddSimulationModifier(p_data_tracking_modifier);

			/* Add linear spring force which has different spring stiffness constants, depending
			 * on the pair of cells it is connecting.
			 */
			MAKE_PTR(EpithelialLayerLinearSpringForce<2>, p_spring_force);
			p_spring_force->SetCutOffLength(1.5);
			//Set the spring stiffnesses
			p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
			p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
			p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
			p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
			simulator.AddForce(p_spring_force);

			/* Add the basement membrane force. */
			MAKE_PTR(EpithelialLayerBasementMembraneForce, p_bm_force);
			p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
			p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
			simulator.AddForce(p_bm_force);

			/* Add an anoikis-based cell killer. */
			MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller, p_anoikis_killer, (&cell_population));
			simulator.AddCellKiller(p_anoikis_killer);

			// Add cell age writer
			cell_population.AddCellWriter<CellAgesWriter>();

			/* We fix all cells outside of the 20 x 20 box. */

			c_vector<double,2> point = zero_vector<double>(2);
			c_vector<double,2> normal = zero_vector<double>(2);

			//Fix cells in the region x < 0
			normal(0) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc1);

			//Fix cells in the region x > 20
			point(0) = 32.0;
			normal(0) = 1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc2);

			//Fix cells in the region y < 0
			point(0) = 0.0;
			point(1) = 0.0;
			normal(0) = 0.0;
			normal(1) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc3);

			//Fix cells in the region y > 20
			point(1) = 32.0;
			normal(1) = 1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc4);

			/* Run the simulation. */
			simulator.Solve();

			SimulationTime::Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);
			RandomNumberGenerator::Destroy();



	}
};

#endif /* PROJECTS_CRYPTFISSIONPLOS2016_TEST_TESTSINGLEORGANOIDSIMULATION_HPP_ */
