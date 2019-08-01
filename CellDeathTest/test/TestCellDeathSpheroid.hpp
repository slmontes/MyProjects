/*
 * Author: Sandra Montes
 * Date: 08/Feb/3019
 *
 */

#ifndef TESTCELLDEATHSPHEROID_HPP_
#define TESTCELLDEATHSPHEROID_HPP_

/*
 * = Examples showing how to create, run and visualize 2D node-based simulations =
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * Necessary header files.
 */

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is usually included in all cell-based test suites. It enables us to write tests where the {{{SimulationTime}}} is handled automatically and simplifies the tests.*/
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. We encountered some of these header files in
 * UserTutorials/RunningMeshBasedSimulations. */
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "EpithelialLayerLinearSpringForce3D.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeMap.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimulationTime.hpp"
#include "TimeStepper.hpp"
#include "SimTimer.hpp"
#include "SmartPointers.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "TransitCellProliferativeType.hpp"
#include "AbstractCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"

/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
#include "NodeAttributes.hpp"
/* The next header file defines a boundary condition which can be found in the source folder*/
#include "SphereBoundaryCondition.hpp"
// Cell population writers
#include "VoronoiDataWriter.hpp"

//Include oxigen dependant cell cycle
#include "SimpleOxygenBasedCellCycleModel.hpp"
//Include cell types
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
//Include a way to kill the apoptotic cells (due to lack of oxygen)
#include "ApoptoticCellKiller.hpp"
#include "ApoptoticCellProperty.hpp"
#include <cmath>
using namespace std;
/* Next, we define the test class.
 */
class TestCellDeathSpheroid : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test - a node-based simulation on a restricted spheroid geometry ==
     *
     * EMPTYLINE
     *
     */
    void TestCellDeathOnSphere()
    {
        EXIT_IF_PARALLEL;

        // Prepare initial state
        double n_cells = 100;
        double dist_to_nb = 0.75;
        double lumen_radius = pow(n_cells / 0.64, 1. / 3) * dist_to_nb / 2;
        double mean_le = 5;
        double sd_le = 1;
        double simulation_duration = 10;
        double sample_each_time = 2.5;
        double cell_radius = 0.5;
        // Set spheriod of cells with centre cell for lumen then outer monolayer
        std::vector<Node<3>*> nodes;
        for (int i = 0; i < n_cells; i++)
        {
            auto r = lumen_radius;
            auto phi = rand() / (RAND_MAX + 1.) * 2 * M_PI;
            auto theta = acos(2. * rand() / (RAND_MAX + 1.) - 1);
            nodes.push_back(new Node<3>(i, false, r * sin(theta) * cos(phi),
                                        r * sin(theta) * sin(phi), r * cos(theta)));
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5); // Distance cut-off
        std::vector<unsigned> real_indices = mesh.GetAllNodeIndices();
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Create population with monolayer cells default radius at 0.5 and centred
        // around (0,0,0)
        NodeBasedCellPopulation<3> cell_population(mesh, cells, real_indices);


        // Create "spring" forces between cells for interactions
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeDivisionRestingSpringLength(cell_radius);

        // The life expectancy of the epithelial cells will be gaussian
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{mean_le,sd_le};
        vector<double> life_expectancies(cell_population.GetNumNodes());
        for (unsigned i=0; i < cell_population.GetNumNodes(); i++)
        {
            life_expectancies[i] = d(gen);
        }

        //Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Run through each timestep in the simulation

        // ** Here is where I am struggling, when I try to call my SimTimer class
        // the test fails. We believe it may be because of the simulation time
        // function being used ** \\

        SimTimer counter;
        int sim_ongoing = 1;
        // Runs this loop while the timesteps have not taken the simulation beyond
        // its duration
        while (sim_ongoing=1)
        {
              int cell_index = 0;
              double apoptosis_count = 0;
              vector<int> dead_cells_ind = {};
              // Kills the cells thats ages are greater than their life expectancies
              // and counts how many are killed each time.
              for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
                   cell_iter != r_cells.end();
                   ++cell_iter)
              {
                    double age = (*cell_iter)->GetAge();
                    double life_expectancy = life_expectancies[cell_index];
                    if (!(*cell_iter)->IsDead() &&
                    life_expectancy<=age)
                    {
                        (*cell_iter)->Kill();
                        apoptosis_count++;
                        dead_cells_ind.insert(dead_cells_ind.begin(),cell_index);
                    }
                    cell_index++;

              }

              cell_population.RemoveDeadCells();
              for (int i : dead_cells_ind)
              {
                life_expectancies.erase(life_expectancies.begin()+i);
              }
              // NodeMap map(mesh.GetNumNodes());
              // mesh.ReMesh(map);

              // The increase of the radius of the organoid is proportional to the number of
              // cells dying in each timestep. The cells are then pushed back.
              lumen_radius = pow(pow(lumen_radius,3)+apoptosis_count*pow(cell_radius,3), 1./3);
              for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
                   cell_iter != r_cells.end();
                   ++cell_iter)
              {
                    c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
                    double previous_radius = norm_2(cell_location);
                    if (previous_radius<lumen_radius)
                    {
                            c_vector<double, 3> location_on_sphere = lumen_radius*(cell_location)/previous_radius;
                            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                            Node<3>* node = cell_population.GetNode(node_index);
                            node->rGetModifiableLocation() = location_on_sphere;
                    }
              }

              // ** Here is where it seems to break down ** \\ 
              counter.MoveOneTimeStep();
              counter.IsSimFinished(sim_ongoing);

        }

        // Setup simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestCellDeathSpheroid");
        simulator.SetSamplingTimestepMultiple(sample_each_time);
        simulator.SetEndTime(simulation_duration);
        simulator.SetNoBirth(true);
        simulator.AddForce(p_spring_force);

        // Run simulation
        simulator.Solve();

    }

};

#endif /* TESTCELLDEATHSPHEROID_HPP_ */
