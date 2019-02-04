#include "PanethCellMutationState.hpp"
#include "EpithelialLayerDataTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


template<unsigned DIM>
EpithelialLayerDataTrackingModifier<DIM>::EpithelialLayerDataTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()

{
}

template<unsigned DIM>
EpithelialLayerDataTrackingModifier<DIM>::~EpithelialLayerDataTrackingModifier()
{
}


template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
//    UpdateCellData(rCellPopulation);
//
//    CalculateModifierData(rCellPopulation);
//
//    *mpEpithelialLayerDataFile << "\n";

}

template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);

    // Create output file
	OutputFileHandler output_file_handler( outputDirectory + "/", false);
    mpEpithelialLayerDataFile = output_file_handler.OpenOutputFile("EpithelialLayerdata.dat");

    //Initialise method
    CalculateModifierData(rCellPopulation);

    *mpEpithelialLayerDataFile << "\n";
}

template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

}

//Overridden method to close file at the end of solve
template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);

	CalculateModifierData(rCellPopulation);

	// Close output file.
    mpEpithelialLayerDataFile->close();
}

//Method to calculate the volume of all the cells and Matrigel
template<unsigned DIM>
double EpithelialLayerDataTrackingModifier<DIM>::CalculateTotalVolume(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	//Find the first node that's a Matrigel node, and record its volume. (We're assuming that the box maintains a constant area)
	double volume_of_tissue = 1.0;
//	double unit_gel_volume = 0.0;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Get the domain width
		c_vector<double,DIM> domain_width = CalculateCellPopulationWidth(rCellPopulation);

		//Multiply the gel node volume by the dimensions of the area
		for (unsigned i = 0; i<DIM; i++)
		{
			volume_of_tissue *= domain_width[i];
		}

	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Get the domain width
		c_vector<double,DIM> domain_width = CalculateCellPopulationWidth(rCellPopulation);

		for (unsigned i = 0; i<DIM; i++)
		{
			volume_of_tissue *= domain_width[i];
		}

	}

	return volume_of_tissue;
}

//Method to calculate the volume of the organoid ring, based on the formula
// A = 0.5*|\sum^{n-1}_{i=0}x_iy_i+1 - x_i+1y_i|, where (x_i, y_i) are the vertices of
// the polygon spanning the ring
template<unsigned DIM>
double EpithelialLayerDataTrackingModifier<DIM>::CalculateRingVolume(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	//Initialise area
	double area_of_ring = 0.0;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Get the cells in the ring, in order

		std::vector<unsigned> cells_in_ring = GetCellsInRingInOrder(rCellPopulation);

		for (unsigned i = 0; i < (cells_in_ring.size() - 1); i++)
		{
			//Get the nodes
			unsigned first_node_index = cells_in_ring[i];
			unsigned second_node_index = cells_in_ring[i+1];

			//Get their locations
			c_vector<double,DIM> first_cell_location = rCellPopulation.GetNode(first_node_index)->rGetLocation();
			c_vector<double,DIM> second_cell_location = rCellPopulation.GetNode(second_node_index)->rGetLocation();

			//Get the x and y coordinates
			double x_i = first_cell_location[0];
			double y_i = first_cell_location[1];
			double x_ip1 = second_cell_location[0];
			double y_ip1 = second_cell_location[1];

			area_of_ring += (x_i*y_ip1 - x_ip1*y_i);
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Get the cells in the ring, in order

		std::vector<unsigned> cells_in_ring = GetCellsInRingInOrder(rCellPopulation);

		for (unsigned i = 0; i < (cells_in_ring.size() - 1); i++)
		{
			//Get the nodes
			unsigned first_node_index = cells_in_ring[i];
			unsigned second_node_index = cells_in_ring[i+1];

			//Get their locations
			c_vector<double,DIM> first_cell_location = rCellPopulation.GetNode(first_node_index)->rGetLocation();
			c_vector<double,DIM> second_cell_location = rCellPopulation.GetNode(second_node_index)->rGetLocation();

			//Get the x and y coordinates
			double x_i = first_cell_location[0];
			double y_i = first_cell_location[1];
			double x_ip1 = second_cell_location[0];
			double y_ip1 = second_cell_location[1];

			area_of_ring += (x_i*y_ip1 - x_ip1*y_i);
		}
	}

	//Halve and ensure area is positive (according to formula)
	area_of_ring = 0.5*fabs(area_of_ring);

	return area_of_ring;
}

// Method to calculate the width of the cell population in the 0th, 1st and 2nd dimension.
// This accounts for the fact that we have ghost nodes, which GetWidth() apparently does
// not do.
template<unsigned DIM>
c_vector<double,DIM> EpithelialLayerDataTrackingModifier<DIM>::CalculateCellPopulationWidth(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	//Initialise vectors
	c_vector<double,DIM> domain_width;
	c_vector<double,DIM> minimum_points;
	c_vector<double,DIM> maximum_points;

	minimum_points = scalar_vector<double>(DIM, DBL_MAX);
	maximum_points = scalar_vector<double>(DIM, -DBL_MAX);

	//Get the minimum and maximum value of points in each dimension
	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //Get the location index corresponding to cell
			Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

			if(!p_cell_population->IsGhostNode(node_index)) //Make sure it's a cell or gel
			{
				c_vector<double,DIM> node_location = p_node->rGetLocation();

				for (unsigned i = 0; i < DIM; i++)
				{
					if (node_location[i] < minimum_points[i])
					{
						minimum_points[i] = node_location[i];
					}
					else if (node_location[i] > maximum_points[i])
					{
						maximum_points[i] = node_location[i];
					}
				}
			}
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //Get the location index corresponding to cell
			Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

			if(!p_cell_population->IsGhostNode(node_index)) //Make sure it's a cell or gel
			{
				c_vector<double,DIM> node_location = p_node->rGetLocation();

				for (unsigned i = 0; i < DIM; i++)
				{
					if (node_location[i] < minimum_points[i])
					{
						minimum_points[i] = node_location[i];
					}
					else if (node_location[i] > maximum_points[i])
					{
						maximum_points[i] = node_location[i];
					}
				}
			}
		}
	}

	//Get the width in all dimensions
	for (unsigned i = 0; i < DIM; i++)
	{
		domain_width[i] = fabs(maximum_points[i] - minimum_points[i]);
	}

	return domain_width;
}

//Method to calculate the perimeter of the organoid ring. Essentially we scan for all the epithelial pairs and
//Calculate the length of the spring attaching them.
template<unsigned DIM>
double EpithelialLayerDataTrackingModifier<DIM>::CalculateRingPerimeter(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	//Initialise perimeter
	double ring_perimeter = 0.0;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Iterate over all the springs
		for (typename MeshBasedCellPopulationWithGhostNodes<DIM>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
				spring_iter != p_cell_population->SpringsEnd();
				++spring_iter)
		{
			unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
			unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

			//If they're both real nodes
			if ( (!p_cell_population->IsGhostNode(nodeA_global_index))&&(!p_cell_population->IsGhostNode(nodeB_global_index)) )
			{
				//Get the cells corresponding to them
				CellPtr p_cell_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
				CellPtr p_cell_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index);
				//If they're both epithelial cells

				boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_A->GetCellProliferativeType();
				boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_B->GetCellProliferativeType();
				if ( (!p_type_A->IsType<DifferentiatedCellProliferativeType>())&&(!p_type_B->IsType<DifferentiatedCellProliferativeType>()) )
				{
					double spring_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index); //Get spring length
					ring_perimeter += spring_length;
				}
			}
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

		//Iterate over all the springs
		for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
				spring_iter != p_cell_population->SpringsEnd();
				++spring_iter)
		{
			unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
			unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

			//If they're both real nodes (don't really need this, but paranoia demands it)
			if ( (!p_cell_population->IsGhostNode(nodeA_global_index))&&(!p_cell_population->IsGhostNode(nodeB_global_index)) )
			{
				//Get the cells corresponding to them
				CellPtr p_cell_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
				CellPtr p_cell_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index);
				//If they're both epithelial cells

				boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_A->GetCellProliferativeType();
				boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_B->GetCellProliferativeType();
				if ( (!p_type_A->IsType<DifferentiatedCellProliferativeType>())&&(!p_type_B->IsType<DifferentiatedCellProliferativeType>()) )
				{
					double spring_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index); //Get spring length
					ring_perimeter += spring_length;
				}
			}
		}
	}

	return ring_perimeter;
}

//Method to count the proliferative cell types
template <unsigned DIM>
std::vector<unsigned> EpithelialLayerDataTrackingModifier<DIM>::CountCellProliferativeTypes(AbstractCellPopulation<DIM>& rCellPopulation)
{
	std::vector<unsigned> proliferative_cell_count = rCellPopulation.GetCellProliferativeTypeCount();

	return proliferative_cell_count;
}


//Method to count the mutant states (don't need at the moment, but implemented for future purposes)
template<unsigned DIM>
std::vector<unsigned> EpithelialLayerDataTrackingModifier<DIM>::CountCellMutationStates(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	std::vector<unsigned> mutant_cell_count = rCellPopulation.GetCellMutationStateCount();

	//Need to count Paneth cells as well now
	unsigned num_paneth_cells = 0;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			if (!p_cell_population->IsGhostNode(node_index))
			{
				boost::shared_ptr<AbstractCellProperty> p_state = cell_iter->GetMutationState();

				if(p_state->IsType<PanethCellMutationState>()) //If we have a paneth cell
				{
					num_paneth_cells += 1;
				}
			}
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);


		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			if (!p_cell_population->IsGhostNode(node_index))
			{
				boost::shared_ptr<AbstractCellProperty> p_state = cell_iter->GetMutationState();

				if(p_state->IsType<PanethCellMutationState>()) //If we have a paneth cell
				{
					num_paneth_cells += 1;
				}
			}
		}
	}

	mutant_cell_count.push_back(num_paneth_cells);

	return mutant_cell_count;
}

//Method to get the cells in the epithelial ring in anti-clockwise order
template<unsigned DIM>
std::vector<unsigned> EpithelialLayerDataTrackingModifier<DIM>::GetCellsInRingInOrder(AbstractCellPopulation<DIM>& rCellPopulation)
{
	//First we obtain the 'centre of mass' of the ring.
	c_vector<double, 2> centre_of_mass = zero_vector<double>(2);
	unsigned num_epithelial_cells = 0;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			if (!p_cell_population->IsGhostNode(node_index)) //Paranoia
			{
				boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

				if(!p_type->IsType<DifferentiatedCellProliferativeType>()) //If it's an epithelial node, get its location
				{
					c_vector<double,DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

					//Add to the pre-existing CoM location vector
					centre_of_mass += epithelial_location;
					num_epithelial_cells += 1;
				}
			}
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			if(!p_type->IsType<DifferentiatedCellProliferativeType>()) //If it's an epithelial node, get its location
			{
				c_vector<double,DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

				//Add to the pre-existing CoM location vector
				centre_of_mass += epithelial_location;
				num_epithelial_cells += 1;
			}
		}
	}

	//Average CoM vector
	centre_of_mass /= num_epithelial_cells;

	std::vector<std::pair<double, unsigned> > angles_and_cells; //Initialise vector

	//Obtain the proliferative cells

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);
		// Iterate over cell population
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			//Get location of cell
			double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
			double y = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			if (!p_cell_population->IsGhostNode(node_index))
			{
				if (!cell_iter->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>() )
				{
					std::pair<double, unsigned> angle_cell;

					//Get point relative to the centre of mass
					double rel_x = x - centre_of_mass[0];
					double rel_y = y - centre_of_mass[1];

					double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

					if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
					{
						circle_angle += M_PI;
					}
					else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
					{
						circle_angle += 2*M_PI;
					}

					angle_cell = std::make_pair(circle_angle, node_index);

					angles_and_cells.push_back(angle_cell); //Add the angle and node index
				}
			}
		}
	}
	if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
		// Iterate over cell population
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			//Get location of cell
			double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
			double y = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			if (!p_cell_population->IsGhostNode(node_index))
			{
				if (!cell_iter->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>() )
				{
					std::pair<double, unsigned> angle_cell;

					//Get point relative to the centre of mass
					double rel_x = x - centre_of_mass[0];
					double rel_y = y - centre_of_mass[1];

					double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

					if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
					{
						circle_angle += M_PI;
					}
					else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
					{
						circle_angle += 2*M_PI;
					}

					angle_cell = std::make_pair(circle_angle, node_index);

					angles_and_cells.push_back(angle_cell); //Add the angle and node index
				}
			}
		}
	}

	//Sort the vector by the angle
	std::sort(angles_and_cells.begin(), angles_and_cells.end());

	//Create vector
	std::vector<unsigned> cells_in_ring;

	for (unsigned i = 0; i < angles_and_cells.size(); i++)
	{
		//Get angle and cell pair
		std::pair<double, unsigned> angle_and_cell = angles_and_cells[i];

		//Get node index
		unsigned node_index = angle_and_cell.second;

		cells_in_ring.push_back(node_index);
	}

	//We need to 'close' the ring, so add the first node at the end
	cells_in_ring.push_back(cells_in_ring[0]);

	return cells_in_ring;
}

//Method to compile all the data together as one vector
template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::CalculateModifierData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	std::vector<double> time_and_data;

	//Add the time
	double time = SimulationTime::Instance()->GetTime();
	time_and_data.push_back(time);

	//Add the ring volume
	double ring_volume = CalculateRingVolume(rCellPopulation);
	time_and_data.push_back(ring_volume);

	//Add the total volume
	double total_volume = CalculateTotalVolume(rCellPopulation);
	time_and_data.push_back(total_volume);

	//Add the ring perimeter
	double ring_perimeter = CalculateRingPerimeter(rCellPopulation);
	time_and_data.push_back(ring_perimeter);

	//Add the cell proliferative type counts
	std::vector<unsigned> proliferative_cells = CountCellProliferativeTypes(rCellPopulation);

	for (unsigned i = 0; i < proliferative_cells.size(); i++)
	{
		time_and_data.push_back(proliferative_cells[i]);
	}

	//Add the cell mutation state counts
	std::vector<unsigned> mutation_states = CountCellMutationStates(rCellPopulation);

	for (unsigned i = 0; i < mutation_states.size(); i++)
	{
		time_and_data.push_back(mutation_states[i]);
	}


	for (unsigned i = 0; i < time_and_data.size(); i++)
	{
		*mpEpithelialLayerDataFile << time_and_data[i] << "\t";
	}
}


template<unsigned DIM>
void EpithelialLayerDataTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class EpithelialLayerDataTrackingModifier<1>;
template class EpithelialLayerDataTrackingModifier<2>;
//template class EpithelialLayerDataTrackingModifier<3>; //3D won't quite work regarding perimeter etc.

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EpithelialLayerDataTrackingModifier)
