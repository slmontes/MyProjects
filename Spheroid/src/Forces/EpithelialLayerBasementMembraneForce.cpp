#include "EpithelialLayerBasementMembraneForce.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include <algorithm> //in order to use .erase
#include <math.h> //in order to use acos

/*
 * Author: Sandra M.
 * Created on: 15/03/2019
 * Last modified: 30/04/2019
 */

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
EpithelialLayerBasementMembraneForce::EpithelialLayerBasementMembraneForce()
   :  AbstractForce<3>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mTargetCurvature(DOUBLE_UNSET)
{
}

EpithelialLayerBasementMembraneForce::~EpithelialLayerBasementMembraneForce()
{

}


void EpithelialLayerBasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double EpithelialLayerBasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}


void EpithelialLayerBasementMembraneForce::SetTargetCurvature(double targetCurvature)
{
	mTargetCurvature = targetCurvature;
}


double EpithelialLayerBasementMembraneForce::GetTargetCurvature()
{
	return mTargetCurvature;
}

/*
 * A method to find all the epithelial nodes
 */
 std::vector<unsigned> EpithelialLayerBasementMembraneForce::GetEpithelialNodes(NodeBasedCellPopulation<3>& rCellPopulation)
 {
 	 //Get cell population
 	 NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation);

 	 // Create a vector to record the nodes corresponding only to epithelial nodes
   std::vector<unsigned> epithelial_nodes;

   // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
   for (NodeBasedCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
   {
     boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

   	 // Only epithelial cell nodes
   	 if ( (p_type->IsType<DifferentiatedCellProliferativeType>()==false) && (!cell_iter->IsDead()) )
      {
 	      Node<3>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
   	    unsigned node_index = p_node->GetIndex();
        //unsigned p_node_location = p_node->rGetLocation();

         epithelial_nodes.push_back(node_index);
       }
     }

   return epithelial_nodes;

 }


/*
 * A method to find all the neighbouring cells of each epithelial node
 */
std::vector<c_vector<unsigned, 3> > EpithelialLayerBasementMembraneForce::GetNodeNeighbours(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index)
{
	//Get cell population
	NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation);
  //Get the current node using the node_index
  //unsigned p_node= p_tissue->GetNode(node_index);
  //Create a vector to record the current node neighbours
  std::set<unsigned> current_node_neighbours;
	// // Create a vector to record all the neighbouring nodes corresponding to epithelial nodes
  // c_vector<unsigned, 3> node_neighbours;
  //
  // // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
  // for (NodeBasedCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
  //      cell_iter != rCellPopulation.End();
  //      ++cell_iter)
  // {
    boost::shared_ptr<AbstractCellProperty> p_type = p_node->GetCellProliferativeType();

  	// Need these to not be gel or dead cells
  	if ( (p_type->IsType<DifferentiatedCellProliferativeType>()==false) && (!cell_iter->IsDead()) )	// an epithelial cell
    {
    	// unsigned node_index = p_node->GetIndex();
    	// unsigned p_node_location = p_node->rGetLocation();

  		/* Find the neighbours of this node
       * It is important to note that, different from the previous definition
       * of BM force in Langlands et al (2016) and Almet et al (2017), in the
       * NodeBasedCellPopulation there are no elements or mesh
       */
       // Get the set of node indices corresponding to this cell's neighbours
       current_node_neigbours = p_tissue->GetNeighbouringNodeIndices(node_index);

       //Remove neighbours that are not epithelial Nodes
       current_node_neigbours.erase( std::remove_if(current_node_neigbours.begin(), current_node_neigbours.end(),
                                    (p_type->IsType<DifferentiatedCellProliferativeType>()), current_node_neigbours.end() );

       // node_neighbours.push_back(p_node->rGetNeighbours());

	     //Remove neighbours that are not epithelial Nodes
	     // node_neighbours.erase( std::remove_if(v.begin(), v.end(), (p_type->IsType<DifferentiatedCellProliferativeType>()), v.end() );
    }
  // }

    return current_node_neigbours;

}

/*
 * A method to find all the vectors that connect an epithelial cell with
 * its neighbouring cells
 */

std::vector<c_vector<unsigned, 3> > EpithelialLayerBasementMembraneForce::GetNodetoNeighbourVectors(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index)
{
  //Get cell population
	NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation);
  //Get p_node
  unsigned p_node= p_tissue->GetNode(node_index);
  //Get the p_node location with node_index
  unsigned p_node_location = p_node->rGetLocation();
  // Start by identifying the node neighbours
  std::vector<c_vector<unsigned, 3> > current_node_neigbours = GetNodeNeighbours(rCellPopulation, node_index);
  //Create a vector where the node to neighbour vectors will be stored
  c_vector<unsigned, 3> node_to_neighbour_vectors;  //in our diagram these will be either b or c

  // Iterate among the neighbours
  // for (std::vector<unsigned>::iterator iter = current_node_neigbours.begin();
  //      iter != current_node_neigbours.end();
  //      ++iter)
  for (unsigned i=0; i<current_node_neigbours.size(); i++)
  {
    Node<3>* p_node_j = p_tissue->GetNode(current_node_neigbours[i]);

    // Get the location of this node
    unsigned p_node_j_location = p_node_j->rGetLocation();

    // Get the vector the two nodes (using GetVectorFromAtoB)
    c_vector<double, 3> node_to_node_j_vector = p_tissue->rGetMesh().GetVectorFromAtoB(p_node_location, p_node_j_location);

    //Vector with all the vectors from p_node to its neighbours
    node_to_neighbour_vectors.push_back(node_to_node_j_vector);
  }

  return node_to_neighbour_vectors;

}

/*
 * A method to find all the vectors that connect one neighbouring cell with the
 * immediate next neighbouring cell (side a of our diagram triangle)
 */
std::vector<c_vector<unsigned, 3> > EpithelialLayerBasementMembraneForce::GetNeighbourtoNeighbourVectors(NodeBasedCellPopulation<3>& rCellPopulation, unsigned node_index)
{
  //Get cell population
	NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation);
  //Get p_node
  unsigned p_node= p_tissue->GetNode(node_index);
  // Start by identifying the node neighbours
  std::vector<c_vector<unsigned, 3> > current_node_neigbours = GetNodeNeighbours(rCellPopulation, node_index);
  //Create a vector where the neighbour to neighbour vectors will be stored
  c_vector<unsigned, 3> neighbour_to_neighbour_vectors; //in our diagram it will be side a of our triangle

  // Iterate among the neighbours
  for (std::vector<unsigned>::iterator iter = current_node_neigbours.begin();
       iter != current_node_neigbours.end();
       ++iter)
  {
    //Get a pair of neighbours to calculate the vector between them
    Node<3>* p_node_n1 = rCellPopulation.GetNode(*iter);

    if (((*iter)+1) != current_node_neigbours.end())
    {
      Node<3>* p_node_n2 = rCellPopulation.GetNode(*iter);
    }
    //But if we are at the end of our vector, the last value should be compared with the first one
    else if (((*iter)+1) == current_node_neigbours.end())
    {
      Node<3>* p_node_n2 = rCellPopulation.GetNode(current_node_neigbours.begin());
    }
    // Get the location of these nodes
    c_vector<double, 3> p_node_n1_location = p_node_n1->rGetLocation();
    c_vector<double, 3> p_node_n2_location = p_node_n2->rGetLocation();

    // Get the vector the between the two neighbour nodes (using GetVectorFromAtoB)
    c_vector<double, 3> node_n1_to_node_n2_vector = p_tissue->rGetMesh().GetVectorFromAtoB(p_node_n1_location, p_node_n2_location);

    neighbour_to_neighbour_vectors.push_back(node_n1_to_node_n2_vector);

  }

  return neighbour_to_neighbour_vectors;

}

/*
 * A method to find all the angles that connect p_node with its neighbouring cells
 * We can calculate this using the Law of Cosines
 */
std::vector<c_vector<unsigned, 3> > EpithelialLayerBasementMembraneForce::GetNodetoNeighbourAngles(NodeBasedCellPopulation<3>& rCellPopulation)
{
  // Start by identifying the node_to_neighbour_vectors
  std::vector<c_vector<unsigned, 3> > node_to_neighbour_vectors = GetNodetoNeighbourVectors(rCellPopulation, node_index);
  // Start by identifying the neighbour_to_neighbour_vectors
  std::vector<c_vector<unsigned, 3> > neighbour_to_neighbour_vectors = GetNeighbourtoNeighbourVectors(rCellPopulation, node_index);
  //Create a vectore to store all the calculated angles
  c_vector<unsigned, 3> angles_vector; //in our diagram it will be side a of our triangle

  //Iterate among the node_to_neighbour_vectors
  // for (std::vector<unsigned>::iterator iter = node_to_neighbour_vectors.begin();
  //     iter != node_to_neighbour_vectors.end();
  //     ++iter)
  for (unsigned i=0; i<node_to_neighbour_vectors.size(); i++)
  {
    unsigned Side_B = norm_2(node_to_neighbour_vectors[i]);
    unsigned Side_C = norm_2(node_to_neighbour_vectors[i+1]);
    unsigned Side_A = norm_2(neighbour_to_neighbour_vectors[i]);

    // unsigned Side_B = node_to_neighbour_vectors[i];
    // unsigned Side_C = node_to_neighbour_vectors[i+1];
    // unsigned Side_A = neighbour_to_neighbour_vectors[i];

    unsigned angle = (acos)*(((pow(Side_B,2))+(pow(Side_C,2))+(pow(Side_A,2)))/(2*Side_B*Side_C));

    angles_vector.push_back(angle);
  }

 return angles_vector;

}

/*
* A method to find all the Areas of the triangles that connect p_node with
* its neighbouring cells. As we can calculate the sides of the triangle, we
* use Heron's formula to calculate the area.
*/
std::vector<c_vector<unsigned, 3> > EpithelialLayerBasementMembraneForce::GetNodetoNeighbourAreas(NodeBasedCellPopulation<3>& rCellPopulation)
{
  // Start by identifying the node_to_neighbour_vectors
  std::vector<c_vector<unsigned, 3> > node_to_neighbour_vectors = GetNodetoNeighbourVectors(rCellPopulation, node_index);
  // Start by identifying the neighbour_to_neighbour_vectors
  std::vector<c_vector<unsigned, 3> > neighbour_to_neighbour_vectors = GetNeighbourtoNeighbourVectors(rCellPopulation, node_index);

  //Create a vector to store the triangles' areas
  c_vector<unsigned, 3> Areas_vector; //in our diagram it will be side a of our triangle

  // for (std::vector<unsigned>::iterator iter = node_to_neighbour_vectors.begin();
  //      iter != node_to_neighbour_vectors.end();
  //      ++iter) ***************************************************************************Fix the end of the loop!!!!!
  for (unsigned i=0; i<node_to_neighbour_vectors.size(); i++)
  {
    unsigned Side_B = norm_2(node_to_neighbour_vectors[i]);
    unsigned Side_C = norm_2(node_to_neighbour_vectors[i+1]);
    unsigned Side_A = norm_2(neighbour_to_neighbour_vectors[i]);

    // unsigned Side_B = node_to_neighbour_vectors[i];
    // unsigned Side_C = node_to_neighbour_vectors[i+1];
    // unsigned Side_A = neighbour_to_neighbour_vectors[i];

    unsigned semiperimeter = (Side_A+Side_B+Side_C)/2;

    unsigned Triangle_area = sqrt(semiperimeter*(semiperimeter-Side_A)*(semiperimeter-Side_B)*(semiperimeter-Side_C));

    Areas_vector.push_back(Triangle_area);
  }

  return Areas_vector;

}

/*
 * A method to calculate the Gaussian Curvature.
 * Which references? Please Add
 */
 double EpithelialLayerBasementMembraneForce::GetGaussianCurvature(NodeBasedCellPopulation<3>& rCellPopulation)
 {
   // Start by identifying the angles_vector
   c_vector<unsigned, 3> angles_vector = GetNodetoNeighbourAngles(rCellPopulation);
   // Start by identifying the Areas_vector
   c_vector<unsigned, 3> Areas_vector = GetNodetoNeighbourAreas(rCellPopulation);

   //Define Pi value
   const double pi = 3.1415926535897;
   //Calculate the sum of all angles obtained from p_node neighbours
   double sum_of_angles_vector = std::accumulate(angles_vector.begin(), angles_vector.end(), 0.0);
   //Calculate the sum of all areas of triangles obtained from p_node neighbours
   double sum_of_Areas_vector = std::accumulate(Areas_vector.begin(), Areas_vector.end(), 0.0);

   double GaussianCurvature = ((2*pi)-sum_of_angles_vector)/((1/3)*sum_of_Areas_vector);

   //Return the GaussianCurvature calculated at p_node
   return GaussianCurvature;

 }

/*
* A function to calculate the force applied to p_node due to the BM force.
* It overrides the virtual method for AbstractForce and adds the BM force
* to the other forces acting on the epithelial node.
*/
void EpithelialLayerBasementMembraneForce::AddForceContribution(NodeBasedCellPopulation<3>& rCellPopulation)
{
  //Get cell population
	NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation);
  // Start by identifying all the epithelial nodes index
  std::vector< c_vector<unsigned, 3> > epithelial_nodes = GetEpithelialNodes(rCellPopulation);

  //Define the centre of the spheriod to help define the direction of the force
  //(Will this also work when the fingerlike structures are generated?***)
  //Help define the initial conditions centre
  c_vector<double,3> spheroid_centre;
  spheroid_centre(0) = 0;
  spheroid_centre(1) = 0;
  spheroid_centre(2) = 0;

  // We iterate over all epithelial cells in the tissue, and find the force
  // acting on them and its direction.
  for (unsigned i=0; i<epithelial_nodes.size(); i++)
  {
    //Get the current epithelial node index
    unsigned epithelial_node_index = epithelial_nodes[i];
    //Get the location of each epithelial cell
    c_vector<double, 3> epithelial_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();

    // Get the vector from the node -> spheroid_centre to define the direction
    // of the BM force.
    c_vector<double, 3> curvature_force_direction = p_tissue->rGetMesh().GetVectorFromAtoB(epithelial_location, spheroid_centre);

    //Calculate the unit vector in the direction of the given node -> spheroid_centre vector
    double distance_between_node_spheroidcentre = norm_2(curvature_force_direction);
    assert(distance_between_node_spheroidcentre > 0);
		assert(!std::isnan(distance_between_node_spheroidcentre));

    curvature_force_direction /= distance_between_node_spheroidcentre;

    double GaussianCurvature = GetGaussianCurvature(rCellPopulation);

    //Get the BM parameter fromthe main script
    double basement_membrane_parameter = GetBasementMembraneParameter();

    //Calculate the force that the BM is exerting to the epithelial node
    c_vector<double, 3> force_due_to_basement_membrane = basement_membrane_parameter*GaussianCurvature*curvature_force_direction;

    //Add this force to the other forces acting on that nodes
    rCellPopulation.GetNode(epithelial_node_index)->AddAppliedForceContribution(force_due_to_basement_membrane);
  }

}

void EpithelialLayerBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n";
	*rParamsFile <<  "\t\t\t<TargetCurvature>"<< mTargetCurvature << "</TargetCurvature> \n";

	// Call direct parent class
	AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EpithelialLayerBasementMembraneForce)
