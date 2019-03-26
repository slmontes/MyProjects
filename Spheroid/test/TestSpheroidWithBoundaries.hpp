/*

Copyright (c) 2005-2018, University of Oxford.
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
/*
 * Author: Sandra Montes
 * Date: 08/Feb/2019
 *
 */

#ifndef TESTSPHEROIDWITHBOUNDARIES_HPP_
#define TESTSPHEROIDWITHBOUNDARIES_HPP_

/*
 * = Examples showing how to create, run and visualize node-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize node-based simulations.
 * Full details of the mechanical model can be found in Pathamathan et al "A computational study of
 * discrete mechanical tissue models", Physical Biology. Vol. 6. No. 3. 2009.. DOI (10.1088/1478-3975/6/3/036001).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials (UserTutorials/RunningMeshBasedSimulations), we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is usually included in all cell-based test suites. It enables us to write tests where the {{{SimulationTime}}} is handled automatically and simplifies the tests.*/
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. We encountered some of these header files in
 * UserTutorials/RunningMeshBasedSimulations. */
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp" //Mutation class that defines Paneth cells
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "EpithelialLayerLinearSpringForce.hpp" //Spring force law to account for different cell type pairs

#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
/* The next header file defines a boundary condition which can be found in the source folder*/
#include "SphereBoundaryCondition.hpp"
/* The next headers are included for the 2nd test, trying 3D mesh*/
#include "WildTypeCellMutationState.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
// Cell population writers
//#include "CellProliferativeTypesCountWriter.hpp"
//#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
//TestCellPopulation3dTRY
#include "ApoptoticCellProperty.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"


/* Next, we define the test class.
 */
class TestSpheroidWithBoundaries : public AbstractCellBasedTestSuite
{

public:
    /* EMPTYLINE
     *
     * == Test - a node-based simulation on a restricted spheroid geometry ==
     *
     * EMPTYLINE
     *
     */
     void TestCube()
     {
         /** The next line is needed because we cannot currently run node based simulations in parallel. */
         EXIT_IF_PARALLEL;

         //Set the stiffness ratio for Paneth cells to stem cells. This is the
         double stiffness_ratio = 1.0;

         //Set the relative adhesiveness of hard cells. The values used in the paper were 1.0, 1.3 and 2.0
		     double hard_cell_drag_multiplier = 1.0;

         //Set all the spring stiffness variables
    		 double epithelial_epithelial_stiffness = 15.0; //Epithelial-epithelial spring connections
    		 double epithelial_nonepithelial_stiffness = 30.0; //Epithelial-non-epithelial spring connections
    		 double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections

         std::vector<Node<3>*> nodes;
         //Nodes of a Cube of 3x3x3

         nodes.push_back(new Node<3>(0u,  false,  -0.25, -0.25, -0.25));
         nodes.push_back(new Node<3>(1u,  false,  0.0, -0.25, -0.25));
         nodes.push_back(new Node<3>(2u,  false,  0.25, -0.25, -0.25));
         nodes.push_back(new Node<3>(3u,  false,  -0.25, 0, -0.25));
         nodes.push_back(new Node<3>(4u,  false,  0, 0, -0.25));
         nodes.push_back(new Node<3>(5u,  false,  0.25, 0, -0.25));
         nodes.push_back(new Node<3>(6u,  false,  -0.25, 0.25, -0.25));
         nodes.push_back(new Node<3>(7u,  false,  0, 0.25, -0.25));
         nodes.push_back(new Node<3>(8u,  false,  0.25, 0.25, -0.25));

         nodes.push_back(new Node<3>(9u,  false,  -0.25, -0.25, 0));
         nodes.push_back(new Node<3>(10u,  false,  0, -0.25, 0));
         nodes.push_back(new Node<3>(11u,  false,  0.25, -0.25, 0));
         nodes.push_back(new Node<3>(12u,  false,  -0.25, 0, 0));
         nodes.push_back(new Node<3>(13u,  false,  0, 0, 0));
         nodes.push_back(new Node<3>(14u,  false,  0.25, 0, 0));
         nodes.push_back(new Node<3>(15u,  false,  -0.25, 0.25, 0));
         nodes.push_back(new Node<3>(16u,  false,  0, 0.25, 0));
         nodes.push_back(new Node<3>(17u,  false,  0.25, 0.25, 0));

         nodes.push_back(new Node<3>(18u,  false,  -0.25, -0.25, 0.25));
         nodes.push_back(new Node<3>(19u,  false,  0, -0.25, 0.25));
         nodes.push_back(new Node<3>(20u,  false,  0.25, -0.25, 0.25));
         nodes.push_back(new Node<3>(21u,  false,  -0.25, 0, 0.25));
         nodes.push_back(new Node<3>(22u,  false,  0, 0, 0.25));
         nodes.push_back(new Node<3>(23u,  false,  0.25, 0, 0.25));
         nodes.push_back(new Node<3>(24u,  false,  -0.25, 0.25, 0.25));
         nodes.push_back(new Node<3>(25u,  false,  0, 0.25, 0.25));
         nodes.push_back(new Node<3>(26u,  false,  0.25, 0.25, 0.25));

         NodesOnlyMesh<3> p_mesh;
         p_mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        //  //Translate mesh 1.5 units left and 1.5 units down, so that we have a sufficient layer of fixed cells around the boundary.
     		// c_vector<double, 3> translate_x = zero_vector<double>(3);
     		// translate_x(0) = -1.5;
        //
     		// c_vector<double, 3> translate_y = zero_vector<double>(3);
     		// translate_y(1) = -1.5;
        //
        //  c_vector<double, 3> translate_z = zero_vector<double>(3);
        //  translate_z(2) = -1.5;
        //
        //  //Translate mesh appropriately
        //  p_mesh.Translate(translate_x);
        //  p_mesh.Translate(translate_y);
        //  p_mesh.Translate(translate_z);


         /* Define the initially spheroid lumen by centre and radius */
         c_vector<double,3> spheroid_centre;
         spheroid_centre(0) = 0;
         spheroid_centre(1) = 0;
         spheroid_centre(2) = 0;

         std::vector<unsigned> real_indices = p_mesh.GetAllNodeIndices(); //Obtain the locations of real nodes

         //std::vector<unsigned> real_indices; //Vector used to define the locations of non-ghost nodes

         // //Sweep over the initial real indices
         //   for (unsigned i = 0; i < real_indices.size(); i++)
         //   {
         //     unsigned cell_index = initial_real_indices[i];
         //     double x = p_mesh.GetNode(cell_index)->rGetLocation()[0];
         //     double y = p_mesh.GetNode(cell_index)->rGetLocation()[1];
         //     double z = p_mesh.GetNode(cell_index)->rGetLocation()[2];
         //
         //     // If the location of the node falls inside the defined lumen region, then it becomes a ghost node.
         //     //if (pow(x-spheroid_centre[0],3) + pow(y-spheroid_centre[1],3) + pow(z-spheroid_centre[2],3) > pow(spheroid_radius,3))
         //     //if ((x!=spheroid_centre[0]) && (y!=spheroid_centre[1]) && (z!=spheroid_centre[2]))
         //     if ((x,y,z)!=(spheroid_centre(0),spheroid_centre(1),spheroid_centre(2)))
         //     {
         //       real_indices.push_back(cell_index);
         //     }
         //   }

         boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
         boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
         boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
         boost::shared_ptr<AbstractCellProperty> p_paneth_state = CellPropertyRegistry::Instance()->Get<PanethCellMutationState>();

         std::vector<CellPtr> cells;
         CellsGenerator<UniformCellCycleModel, 3> cells_generator;
         cells_generator.GenerateBasicRandom(cells, p_mesh.GetNumNodes(), p_diff_type);

         /*Create cell population*/
         NodeBasedCellPopulation<3> cell_population(p_mesh, cells, real_indices);
         //MeshBasedCellPopulationWithGhostNodes<3> cell_population(p_mesh, cells, real_indices);

         //Create layer of proliferative cells
   			for (unsigned i = 0; i < real_indices.size(); i++)
   			{
   				unsigned cell_index = real_indices[i];
   				CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
   				double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
   				double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];
          double z = cell_population.GetLocationOfCellCentre(cell_iter)[2];

          //If the cell is located in the spheroid_centre, then it will be a differentiated cell (in the future we will declare this a Dead cell type)
          if ((x==spheroid_centre(0))&&(y==spheroid_centre(1))&&(z==spheroid_centre(2)))
          {
            cell_iter->SetCellProliferativeType(p_diff_type);
          }

          else
   				{
						cell_iter->SetCellProliferativeType(p_stem_type);
   				}
   			}

        //Obtain the proliferative cells
        std::vector<unsigned> cells_in_layer; //Initialise vector

    		for (NodeBasedCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
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



         OffLatticeSimulation<3> simulator(cell_population);
         simulator.SetOutputDirectory("NodeBasedCube_2");
         simulator.SetSamplingTimestepMultiple(12);
         simulator.SetEndTime(10.0);

         // MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
         // simulator.AddForce(p_force);

         /* Add linear spring force which has different spring stiffness constants, depending
          * on the pair of cells it is connecting.
          */
         MAKE_PTR(EpithelialLayerLinearSpringForce<3>, p_spring_force);
         p_spring_force->SetCutOffLength(1.5);
         //Set the spring stiffnesses
         p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
         p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
         p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
         p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
         simulator.AddForce(p_spring_force);

         simulator.Solve();

         // for (unsigned i=0; i<nodes.size(); i++)
         // {
         //     delete nodes[i];
         // }
       }


    // void TestOnSurfaceOfSphere()
    // {
    //     /** The next line is needed because we cannot currently run node based simulations in parallel. */
    //     EXIT_IF_PARALLEL;
    //
    //     /*
    //      * We begin with exactly the same code as the previous test: we create a cell population
    //      * from a mesh and vector of cells, and use this in turn to create
    //      * a simulation object.
    //      */
    //
    //     std::vector<Node<3>*> nodes;
    //     nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
    //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
    //     NodesOnlyMesh<3> mesh;
    //     /* To run node-based simulations you need to define a cut off length (second argument in
    //      * {{{ConstructNodesWithoutMesh}}}), which defines the connectivity of the nodes by defining
    //      * a radius of interaction. */
    //     mesh.ConstructNodesWithoutMesh(nodes, 1.5);
    //
    //     std::vector<CellPtr> cells;
    //     MAKE_PTR(TransitCellProliferativeType, p_transit_type);
    //     CellsGenerator<UniformCellCycleModel, 3> cells_generator;
    //     cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
    //
    //     NodeBasedCellPopulation<3> cell_population(mesh, cells);
    //
    //     OffLatticeSimulation<3> simulator(cell_population);
    //     simulator.SetOutputDirectory("NodeBasedSpheroidWithBoundaries");
    //     simulator.SetSamplingTimestepMultiple(12);
    //     simulator.SetEndTime(100.0);
    //
    //     /* As before, we create a linear spring force and pass it to the simulation object. */
    //     MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
    //     simulator.AddForce(p_force);
    //
    //     /*
    //      * This time we create a {{{CellPopulationBoundaryCondition}}} and pass this to
    //      * the {{{OffLatticeSimulation}}}. Here we use a {{{SphereBoundaryCondition}}}
    //      * which restricts cells to lie on a sphere (in 3D).
    //      *
    //      * For a list of possible boundary conditions see subclasses of {{{AbstractCellPopulationBoundaryCondition}}}.
    //      * These can be found in the inheritance diagram, here, [class:AbstractCellPopulationBoundaryCondition AbstractCellPopulationBoundaryCondition].
    //      * Note that some of these boundary conditions are not compatible with node-based simulations see the specific class documentation for details,
    //      * if you try to use an incompatible class then you will receive a warning.
    //      *
    //      * First we set the centre (0,0,1) and radius of the sphere (1).
    //      */
    //     c_vector<double,3> centre = zero_vector<double>(3);
    //     centre(2) = 1.0;
    //     double radius = 5.0;
    //     /* We then make a pointer to the boundary condition using the MAKE_PTR_ARGS macro, and pass
    //      * it to the {{{OffLatticeSimulation}}}. */
    //     MAKE_PTR_ARGS(SphereBoundaryCondition<3>, p_boundary_condition, (&cell_population, centre, radius));
    //     simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);
    //
    //     /* To run the simulation, we call {{{Solve()}}}. */
    //     simulator.Solve();
    //
    //     /* The next two lines are for test purposes only and are not part of this tutorial.
    //      */
    //     //TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
    //     //TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    //
    //     /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.
    //     for (unsigned i=0; i<nodes.size(); i++)
    //     {
    //         delete nodes[i];
    //     }
    //     */
    // }

//     void TestGhostNodesSpheroidSimulation3D()
//     {
//         unsigned width = 3;
//         unsigned height = 3;
//         unsigned depth = 3;
//
//         MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
//         p_mesh->ConstructCuboid(width, height, depth);
//
//         c_vector<double, 3> spheroid_centre;
//         spheroid_centre[0] = 0.5*((double) width);
//         spheroid_centre[1] = 0.5*((double) height);
//         spheroid_centre[2] = 0.5*((double) depth);
//
//         // Set up cells by iterating through the mesh nodes
//         unsigned num_nodes = p_mesh->GetNumAllNodes();
//
//         std::vector<CellPtr> cells;
//         std::vector<CellPtr> cells2;
//         std::vector<unsigned> location_indices;
//
//         MAKE_PTR(WildTypeCellMutationState, p_state);
//         MAKE_PTR(StemCellProliferativeType, p_stem_type);
//
//         for (unsigned i=0; i<num_nodes; i++)
//         {
//             c_vector<double, 3> node_location = p_mesh->GetNode(i)->rGetLocation();
//
//             unsigned min_spatial_dimension;
//             if (width <= height && width <= depth)
//             {
//                 min_spatial_dimension = width;
//             }
//             else
//             {
//                 if (height <= depth)
//                 {
//                     min_spatial_dimension = height;
//                 }
//                 else
//                 {
//                     min_spatial_dimension = depth;
//                 }
//             }
//
//             FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
//             p_model->SetGeneration(0);
//             CellPtr p_cell(new Cell(p_state, p_model));
//             p_cell->SetCellProliferativeType(p_stem_type);
//             p_cell->SetBirthTime(-RandomNumberGenerator::Instance()->ranf()*
//                                  (p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration()));
//
//             cells2.push_back(p_cell);
//
//             if (norm_2(node_location - spheroid_centre) <= 0.5*sqrt(3.0)*1.01*((double) min_spatial_dimension)/3.0)
//             {
//                 location_indices.push_back(i);
//                 cells.push_back(p_cell);
//             }
//         }
//
// //        // Test Save() with a MeshBasedCellPopulationWithGhostNodes
//         MeshBasedCellPopulationWithGhostNodes<3> cell_population(*p_mesh, cells, location_indices);
// //        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
//         cell_population.AddPopulationWriter<VoronoiDataWriter>();
// //        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
// //        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
// //        cell_population.AddCellWriter<CellVolumesWriter>();
//
//         OffLatticeSimulation<3> simulator(cell_population);
//         simulator.SetOutputDirectory("TestSpheroid3DMesh");
//         simulator.SetEndTime(25);
//
//         // Create a force law and pass it to the simulation
//         MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
//         p_linear_force->SetCutOffLength(1.0);
//         simulator.AddForce(p_linear_force);
//
//         simulator.Solve();
//     }

    // void TestMesh3D()
    // {
    //     /** The next line is needed because we cannot currently run node based simulations in parallel. */
    //     EXIT_IF_PARALLEL;
    //
    //     // Create a simple 3D mesh
    //     std::vector<Node<3>*> nodes;
    //     nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
    //     nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
    //     nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
    //     nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
    //     MutableMesh<3,3> mesh(nodes);
    //     /*
    //     std::vector<Node<3>*> nodes;
    //     nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
    //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
    //     NodesOnlyMesh<3> mesh;
    //     mesh.ConstructNodesWithoutMesh(nodes, 1.5);
    //     */
    //
    //     // Set up cells
    //     // std::vector<CellPtr> cells;
    //     // CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
    //     // cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
    //
    //     std::vector<CellPtr> cells;
    //     MAKE_PTR(TransitCellProliferativeType, p_transit_type);
    //     CellsGenerator<UniformCellCycleModel, 3> cells_generator;
    //     cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
    //
    //     // Create cell population
    //     MeshBasedCellPopulation<3> cell_population(mesh, cells);
    //     cell_population.CreateVoronoiTessellation();
    //     VertexMesh<3,3>* p_tessellation = cell_population.GetVoronoiTessellation();
    //     //NodeBasedCellPopulation<3> cell_population(mesh, cells);
    //
    //     // Create Voronoi tessellation
    //     cell_population.CreateVoronoiTessellation();
    //
    //     OffLatticeSimulation<3> simulator(cell_population);
    //     simulator.SetOutputDirectory("Mesh3D");
    //     simulator.SetSamplingTimestepMultiple(12);
    //     simulator.SetEndTime(100.0);
    //
    //     MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
    //     simulator.AddForce(p_force);
    //
    //     simulator.Solve();
    //     /*
    //     TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
    //     TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    //
    //     for (unsigned i=0; i<nodes.size(); i++)
    //     {
    //         delete nodes[i];
    //     }*/
    //  }


};

#endif /* TESTSPHEROIDWITHBOUNDARIES_HPP_ */
