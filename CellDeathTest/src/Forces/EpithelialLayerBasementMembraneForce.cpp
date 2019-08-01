#include "EpithelialLayerBasementMembraneForce.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

/*
 * Created on: 21/12/2014
 * Last modified: 02/10/2015
 */

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
EpithelialLayerBasementMembraneForce::EpithelialLayerBasementMembraneForce()
   :  AbstractForce<2>(),
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

void EpithelialLayerBasementMembraneForce::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * A method to find all the pairs of connections between healthy epithelial cells and labelled gel cells.
 * Returns a vector of node pairings, without repeats. The first of each pair is the epithelial node index,
 * and the second is the gel node index. Updating so that it also returns mutant-labelled cell pairs.
 */

std::vector<c_vector<unsigned, 2> > EpithelialLayerBasementMembraneForce::GetEpithelialGelPairs(AbstractCellPopulation<2>& rCellPopulation)
{

    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Create a vector to record the pairs of nodes corresponding to *joined* epithelial and gel nodes
    std::vector<c_vector<unsigned, 2> > node_pairs;
    c_vector<double, 2> pair;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

    	// Need these to not be stromal cells (and not dead)
    	if ( (p_type->IsType<DifferentiatedCellProliferativeType>()==false) && (!cell_iter->IsDead()) )	// an epithelial cell
    	{
    		Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
    		unsigned node_index = p_node->GetIndex();

    		assert(!(p_tissue->IsGhostNode(node_index)));  // bit unnecessary at this stage but paranoia demands it

    		// ITERATE OVER CONTAINING ELEMENTS and only work with those that DO NOT contain ghost nodes

    		std::vector<unsigned> gel_nodes;

    		for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
		         iter != p_node->ContainingElementsEnd();
		         ++iter)
    		{
    			bool element_contains_ghost_nodes = false;

    			// Get a pointer to the element
    			Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(*iter);

    			// ITERATE OVER NODES owned by this element
    			for (unsigned local_index=0; local_index<3; local_index++)
    			{
    				unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

    				if (p_tissue->IsGhostNode(nodeBGlobalIndex) == true)
    				{
    					element_contains_ghost_nodes = true;
    					break; 				// This should break out of the inner for loop
    				}
    			}

				if (element_contains_ghost_nodes==false)
				{
                    // ITERATE OVER NODES owned by this element
                    for (unsigned local_index=0; local_index<3; local_index++)
                    {
                        unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

						if (p_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType	>()==true)
						{
							// Store the index of each gel node that is attached to the epithelial
							// node. There will be repetitions due to iterating over neighbouring elements
							gel_nodes.push_back(nodeBGlobalIndex);
						}
                    }
				}
			}

    		// Remove any nodes that have been found twice
    		RemoveDuplicates1D(gel_nodes);

    		// Now construct the vector of node pairs
    		for (unsigned i=0; i<gel_nodes.size(); i++)
    		{
    			pair[0] = node_index;
    			pair[1] = gel_nodes[i];
    			node_pairs.push_back(pair);

    			// Check that these node share a common element
				bool has_common_element = false;

				// The elements that contain this epithelial node:
				std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(node_index)->rGetContainingElementIndices();
				assert(epithelial_elements.size() != 0);

				// The elements that contain the gel node:
				std::set<unsigned> gel_elements = rCellPopulation.GetNode(gel_nodes[i])->rGetContainingElementIndices();
				assert(gel_elements.size() != 0);

				// Loop over all elements that contain the gel node
				for (Node<2>::ContainingElementIterator elt_it = rCellPopulation.GetNode(gel_nodes[i])->ContainingElementsBegin();
				         elt_it != rCellPopulation.GetNode(gel_nodes[i])->ContainingElementsEnd();
				         ++elt_it)
				{
					unsigned elt_index = *elt_it;

					bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);

					// Keep only those elements that also contain the epithelial node, but do not have ghost nodes
					if ( (elt_contains_ghost_nodes == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
					{
						// Common element
						has_common_element = true;
						break;
					}
				}

				if (!has_common_element)
				{
					TRACE("No common element between:");
					PRINT_2_VARIABLES(node_index,gel_nodes[i]);
				}
				assert(has_common_element);
    		}
    	}
    }

	return node_pairs;
}

/*
 * Method to determine whether an element contains ghost nodes
 */

bool EpithelialLayerBasementMembraneForce::DoesElementContainGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<3; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;

		}
	}

	return element_contains_ghost_nodes;
}

/*
 * Using the vector of node pairs found in GetEpithelialGelPairs to determine the curvature
 * of the curve that passes through the midpoints of the neighbouring springs, given a epithelial-gel
 * node pairing. (Note - it ignores the end pairs because one of the common elements will contain ghost nodes, but this
 * should only crop up if you don't have periodic boundary conditions)
 * Updating this so that it will still find the curvature if one of the epithelial cells is a mutant cell, eg. apc2 hit
 */

double EpithelialLayerBasementMembraneForce::GetCurvatureFromNodePair(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialNodeIndex,
																unsigned gelNodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Seeking the common elements that contain both the epithelial node and the gel node
	// Note: Don't want to count any elements that have ghost nodes

	std::vector<unsigned> common_elements;	// Initialising

	// The elements that contain this epithelial node:
    std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(epithelialNodeIndex)->rGetContainingElementIndices();
    assert(epithelial_elements.size() != 0);

    assert(rCellPopulation.GetNode(gelNodeIndex)->GetNumContainingElements() != 0);

    // Loop over all elements that contain the gel node
    for (Node<2>::ContainingElementIterator elt_it = rCellPopulation.GetNode(gelNodeIndex)->ContainingElementsBegin();
         elt_it != rCellPopulation.GetNode(gelNodeIndex)->ContainingElementsEnd();
         ++elt_it)
    {
    	unsigned elt_index = *elt_it;

    	bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);

    	// Keep only those elements that also contain the epithelial node, but do not have ghost nodes

        if ( (elt_contains_ghost_nodes == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
        {
        	// Common element
           	common_elements.push_back(elt_index);
        }
    }

	assert(common_elements.size() != 0);		// This is bad - the nodes should be connected in the first place...

	// We iterate over these common elements to find the midpoints of the springs that
	// connect epithelial and  gel nodes

	c_vector<double, 2> spring_midpoint_a, spring_midpoint_b, spring_midpoint_c;

	spring_midpoint_b = p_tissue->GetNode(epithelialNodeIndex)->rGetLocation() + 0.5*(p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(epithelialNodeIndex)->rGetLocation(), p_tissue->GetNode(gelNodeIndex)->rGetLocation()));

    // If there is only one such common element, then this epithelial node will be at either end and so we don't
    // consider the force along the very first / very last spring (only happens if you don't use a cylindrical
	// mesh!)
    if (common_elements.size() == 1)
    {
    	double curvature = 0.0;
    	return curvature;
    }

    else
    {
    	assert(common_elements.size() == 2);		// Should only be two common elements

    	for (std::vector<unsigned>::iterator elem_iter = common_elements.begin();
    									   	 elem_iter != common_elements.end();
    									   	 ++elem_iter)
    	{
    		// Need to determine the cell type of each local node in the element
    		// Want to find the midpoint between epithelial and gel pairs

    		// Looking at the three nodes which form the vertices of the element
    		unsigned global_index_A = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(0);
	   		unsigned global_index_B = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(1);
	   		unsigned global_index_C = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(2);

	   		unsigned E=UINT_MAX;  // E - the epithelial node we're interested in,
	   		unsigned G=UINT_MAX;  // G - the tissue node it's connected to, and
	   		unsigned P=UINT_MAX;  // P - the other node in the element.

	   	    assert(!(p_tissue->IsGhostNode(global_index_A)));
	   	    assert(!(p_tissue->IsGhostNode(global_index_B)));
	   	    assert(!(p_tissue->IsGhostNode(global_index_C)));

	   		if (global_index_A == epithelialNodeIndex)
	   		{
	   			E = global_index_A;

	   			if (global_index_B == gelNodeIndex)
	   			{
	   				G = global_index_B;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_B;
	   				G = global_index_C;
	   			}
	   		}
	   		else if (global_index_B == epithelialNodeIndex)
	   		{
	   			E = global_index_B;

	   			if (global_index_A == gelNodeIndex)
	   			{
	   				G = global_index_A;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				G = global_index_C;
	   			}
	   		}
	   		else if (global_index_C == epithelialNodeIndex)
	   		{
	   			E = global_index_C;

	   			if (global_index_A == gelNodeIndex)
	   			{
	   				G = global_index_A;
	   				P = global_index_B;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				G = global_index_B;
	   			}
	   		}

	   		assert(E<UINT_MAX);
	   		assert(G<UINT_MAX);
	   		assert(P<UINT_MAX);

	   		// Check that we have assigned the epithelial node correctly
	   		boost::shared_ptr<AbstractCellProperty> p_E_type = p_tissue->GetCellUsingLocationIndex(E)->GetCellProliferativeType();
	   		assert(!p_E_type->IsType<DifferentiatedCellProliferativeType>());
			assert(E == epithelialNodeIndex);

	   		// Check that we have assigned the gel node correctly
	   		CellPtr p_G_cell = rCellPopulation.GetCellUsingLocationIndex(G);
	   		boost::shared_ptr<AbstractCellProperty> p_G_type = p_tissue->GetCellUsingLocationIndex(G)->GetCellProliferativeType();
	   		assert(p_G_type->IsType<DifferentiatedCellProliferativeType>());
	   		assert(G == gelNodeIndex);

	   		/*
	   		 * Now we work with E (the epithelial node), G (the tissue node it's connected to) and P, the other node,
	   		 * which we will now assign as either P1 (det < 0) or P2 (det > 0).
	   		 *
	   		 * We also need to determine whether P1 and P2 are epithelial or tissue, as this will affect which vector we
	   		 * choose to take the spring midpoint from.
	   		 */
	   		c_vector<double, 2> vector_G_to_E = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(G)->rGetLocation(),p_tissue->GetNode(E)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_G = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(E)->rGetLocation(),p_tissue->GetNode(G)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_P = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(E)->rGetLocation(),p_tissue->GetNode(P)->rGetLocation());
	   		c_vector<double, 2> vector_G_to_P = vector_G_to_E + vector_E_to_P;
	   		c_vector<double, 2> vector_P_to_G = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(P)->rGetLocation(),p_tissue->GetNode(G)->rGetLocation());

	   		// Now calculate the determinant (GE x EP)

	   		double det = vector_G_to_E[0]*vector_E_to_P[1] - vector_G_to_E[1]*vector_E_to_P[0];
	  		assert(!std::isnan(det));
	   		assert(det != 0.0);

	   		/*
	   		 * If det < 0 then P = P1 and we can assign spring_midpoint_c
	   		 * If det > 0 then P = P2 and we can assign spring_midpoint_a
	   		 * Also need to take into account whether P is a epithelial or tissue node
	   		 * to choose the right spring
	   		 */
	   		boost::shared_ptr<AbstractCellProperty> p_type = p_tissue->GetCellUsingLocationIndex(P)->GetCellProliferativeType();
	   		CellPtr p_cell_P = p_tissue->GetCellUsingLocationIndex(P);

	   		if ( (det < 0) && (p_type->IsType<DifferentiatedCellProliferativeType>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_c = p_tissue->GetNode(E)->rGetLocation() + vector_E_to_P + 0.5*vector_P_to_G;
	   		}
	   		else if ((det < 0) && (p_type->IsType<DifferentiatedCellProliferativeType>()==true))	// P = Gel
	   		{
	   			spring_midpoint_c = p_tissue->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}

	   		else if ( (det > 0) && (p_type->IsType<DifferentiatedCellProliferativeType>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_a = p_tissue->GetNode(E)->rGetLocation() + vector_E_to_G + 0.5*vector_G_to_P;
	   		}

	   		else if ((det > 0) && (p_type->IsType<DifferentiatedCellProliferativeType>()==true))	// P = Gel
	   		{
	   			spring_midpoint_a = p_tissue->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}
    	}

    	double curvature = FindParametricCurvature(rCellPopulation, spring_midpoint_a, spring_midpoint_b, spring_midpoint_c);

    	// Need to subtract the target curvature if the epithelial node lies in the crypt base
    	c_vector<double, 2> epithelial_location = p_tissue->GetNode(epithelialNodeIndex)->rGetLocation();

    	//Subtract the target curvature
    	curvature -= mTargetCurvature;

    	assert(!std::isnan(curvature));
    	return curvature;
    }
}

/*
 * Function to return the curvature between three midpoints parametrically - in this case, we find the normal
 * to the vector joining the left and right midpoints, and then find the perpendicular distance of the centre midpoint
 * from the left->right vector
 */

double EpithelialLayerBasementMembraneForce::GetCurvatureFromMidpoints(AbstractCellPopulation<2>& rCellPopulation,
																c_vector<double, 2> leftMidpoint,
																c_vector<double, 2> centreMidpoint,
																c_vector<double, 2> rightMidpoint)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Firstly find the normal to the vector joining the left and right midpoints
	c_vector<double, 2>	vector_left_to_right = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, rightMidpoint);
	c_vector<double, 2> normal_vector;
	normal_vector[0] = vector_left_to_right[1];
	normal_vector[1] = -vector_left_to_right[0];

	c_vector<double, 2> vector_left_to_centre = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, centreMidpoint);

	double curvature = normal_vector[0]*vector_left_to_centre[0] + normal_vector[1]*vector_left_to_centre[1];

	return curvature;
}

/*
* Function to return the curvature between three points parametrically - the midpoints of the springs connecting the
* transit cells to the differentiated cells. NB. The input arguments need to be in order from either left to right
* or right to left. If they are wrongly arranged (eg. middle, left, right) then you get a different curvature,
* but left->right = -(right-> left).
*/

double EpithelialLayerBasementMembraneForce::FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
															c_vector<double, 2> leftMidpoint,
															c_vector<double, 2> centreMidpoint,
															c_vector<double, 2> rightMidpoint)
{
	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftMidpoint, centreMidpoint);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centreMidpoint, rightMidpoint);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftMidpoint, rightMidpoint);

	// Firstly find the parametric intervals
	double left_s = sqrt(pow(left_to_centre[0],2) + pow(left_to_centre[1],2));
	double right_s = sqrt(pow(centre_to_right[0],2) + pow(centre_to_right[1],2));

	double sum_intervals = left_s + right_s;

	//Calculate finite difference of first derivatives
	double x_prime = (left_to_right[0])/sum_intervals;
	double y_prime = (left_to_right[1])/sum_intervals;

	//Calculate finite difference of second derivatives
	double x_double_prime = 2*(left_s*centre_to_right[0] - right_s*left_to_centre[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*centre_to_right[1] - right_s*left_to_centre[1])/(left_s*right_s*sum_intervals);

	//Calculate curvature using formula
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)),3/2);

	return curvature;
}


/*
 * A method to return the number of elements that contain a particular node,
 * excluding those elements that have ghost nodes
 */

unsigned EpithelialLayerBasementMembraneForce::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Get pointer to the node
    Node<2>* p_node = p_tissue->GetNode(nodeIndex);
    assert(!(p_tissue->IsGhostNode(nodeIndex)));

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd(); ++iter)
    {
        bool element_contains_ghost_nodes = false;
        Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(*iter);

        // Iterate over nodes owned by this element
        for (unsigned local_index=0; local_index<3; local_index++)
        {
            if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
            {
                element_contains_ghost_nodes = true;
                break; // I think this should break out of the inner for loop
            }
        }

        if (element_contains_ghost_nodes==false)
        {
            // This element contains no ghost nodes
            num_elements_with_no_ghost_nodes++;
        }
    }

    return num_elements_with_no_ghost_nodes;
}

/*
 * Method to return the nodes connected to a particular node via the Delaunay
 * triangulation, excluding ghost nodes.
 */

std::set<unsigned> EpithelialLayerBasementMembraneForce::GetNeighbouringNodeIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	assert(!(p_tissue->IsGhostNode(nodeIndex)));

	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
	std::set<unsigned> containing_elem_indices = p_tissue->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Get all the nodes contained in this element
        // Don't want to include the current node
        unsigned neighbour_global_index;

        for (unsigned local_index=0; local_index<3; local_index++)
        {
            neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

            if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
        }
    }

    return neighbouring_node_indices;
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */

bool EpithelialLayerBasementMembraneForce::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool has_cell_detached = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, nodeIndex);

   	unsigned num_gel_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of gel cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
   		boost::shared_ptr<AbstractCellProperty> p_type = p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType();
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_type->IsType<DifferentiatedCellProliferativeType>()==true))
   		{
			num_gel_neighbours += 1;
		}
   	}

   	if(num_gel_neighbours < 1)
   	{
   		has_cell_detached = true;
   	}

	return has_cell_detached;
}
//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void EpithelialLayerBasementMembraneForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// First determine the force acting on each epithelial cell due to the basement membrane
	// Start by identifying the epithelial-gel node pairs (now also returns any apc2hit-gel pairs)
	std::vector<c_vector<unsigned, 2> > node_pairs = GetEpithelialGelPairs(rCellPopulation);

	// We loop over the epithelial-gel node pairs to find the force acting on that
	// epithelial node, and the direction in which it acts
	for (unsigned i=0; i<node_pairs.size(); i++)
	{
		unsigned epithelial_node_index = node_pairs[i][0];
		unsigned gel_node_index = node_pairs[i][1];

   		CellPtr p_cell_epithelial = p_tissue->GetCellUsingLocationIndex(epithelial_node_index);
   		assert(p_cell_epithelial->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false);

		CellPtr p_cell_gel = p_tissue->GetCellUsingLocationIndex(gel_node_index);
		assert(p_cell_gel->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == true);

		c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();
		c_vector<double, 2> gel_location = rCellPopulation.GetNode(gel_node_index)->rGetLocation();

		// The force due to the basal lamina acts along the spring connecting the epithelial and gel nodes, G->E direction
		c_vector<double, 2> curvature_force_direction = p_tissue->rGetMesh().GetVectorFromAtoB(gel_location, epithelial_location);

		double distance_between_nodes = norm_2(curvature_force_direction);
		assert(distance_between_nodes > 0);
		assert(!std::isnan(distance_between_nodes));

		curvature_force_direction /= distance_between_nodes;

		double curvature = GetCurvatureFromNodePair(rCellPopulation, epithelial_node_index, gel_node_index);

		double basement_membrane_parameter = GetBasementMembraneParameter();

		c_vector<double, 2> force_due_to_basement_membrane = basement_membrane_parameter*curvature*curvature_force_direction;

		// Add the force due to the basal lamina to the forces acting on that epithelial node
		rCellPopulation.GetNode(epithelial_node_index)->AddAppliedForceContribution(force_due_to_basement_membrane);
	}

}

void EpithelialLayerBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n";
	*rParamsFile <<  "\t\t\t<TargetCurvature>"<< mTargetCurvature << "</TargetCurvature> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EpithelialLayerBasementMembraneForce)
