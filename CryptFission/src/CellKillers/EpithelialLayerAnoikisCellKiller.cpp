/*
 * LAST MODIFIED: 02/10/2015
 * Anoikis cell killer created for epithelial layer model. Removes any epithelial cells that have detached from
 * the non-epithelial region and entered the lumen.
 *
 * Created on: Dec 21 2014
 * Last modified:
 * 			Author: Axel Almet
 */

#include "EpithelialLayerAnoikisCellKiller.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PanethCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

EpithelialLayerAnoikisCellKiller::EpithelialLayerAnoikisCellKiller(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellKiller<2>(pCellPopulation),
    mCellsRemovedByAnoikis(0),
    mCutOffRadius(1.5)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

EpithelialLayerAnoikisCellKiller::~EpithelialLayerAnoikisCellKiller()
{
//    mAnoikisOutputFile->close();
}

//Method to get mCutOffRadius
double EpithelialLayerAnoikisCellKiller::GetCutOffRadius()
{
	return mCutOffRadius;
}

//Method to set mCutOffRadius
void EpithelialLayerAnoikisCellKiller::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

std::set<unsigned> EpithelialLayerAnoikisCellKiller::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		// Need access to the mesh but can't get to it because the cell killer only owns a
		// pointer to an AbstractCellPopulation
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

		// Find the indices of the elements owned by this node
		std::set<unsigned> containing_elem_indices = p_tissue->rGetMesh().GetNode(nodeIndex)->rGetContainingElementIndices();

		// Iterate over these elements
		for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
				elem_iter != containing_elem_indices.end();
				++elem_iter)
		{
			// Get all the nodes contained in this element
			unsigned neighbour_global_index;

			for (unsigned local_index=0; local_index<3; local_index++)
			{
				neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
				// Don't want to include the original node or ghost nodes
				if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
				{
					neighbouring_node_indices.insert(neighbour_global_index);
				}
			}
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		// Need access to the mesh but can't get to it because the cell killer only owns a
		// pointer to an AbstractCellPopulation
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		//Update cell population
		p_tissue->Update();

		double radius = GetCutOffRadius();

		neighbouring_node_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);
	}

    return neighbouring_node_indices;
}

/** Method to determine if an epithelial cell has lost all contacts with the gel cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool EpithelialLayerAnoikisCellKiller::HasCellPoppedUp(unsigned nodeIndex)
{
	bool has_cell_popped_up = false;	// Initialising

	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

		std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

		unsigned num_gel_neighbours = 0;

		// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

		for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
				neighbour_iter != neighbours.end();
				++neighbour_iter)
		{
			if ( (!p_tissue->IsGhostNode(*neighbour_iter))&&(p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) )
			{
				num_gel_neighbours += 1;
			}
		}

		if(num_gel_neighbours < 1)
		{
			has_cell_popped_up = true;
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

		unsigned num_gel_neighbours = 0;

		// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

		for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
				neighbour_iter != neighbours.end();
				++neighbour_iter)
		{
			if (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
			{
				num_gel_neighbours += 1;
			}
		}

		if(num_gel_neighbours < 1)
		{
			has_cell_popped_up = true;
		}
	}

	return has_cell_popped_up;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by compression-driven apoptosis
 */
std::vector<c_vector<unsigned,2> > EpithelialLayerAnoikisCellKiller::RemoveByAnoikis()
{

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
    {
    	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
    	//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    	c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
    		assert((!p_tissue->IsGhostNode(node_index)));

    		// Initialise
    		individual_node_information[0] = node_index;
    		individual_node_information[1] = 0;

    		// Examine each epithelial node to see if it should be removed by anoikis and then if it
    		// should be removed by compression-driven apoptosis
    		if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    		{
    			// Determining whether to remove this cell by anoikis

    			if(this->HasCellPoppedUp(node_index))
    			{
    				individual_node_information[1] = 1;
    			}
    		}

    		cells_to_remove.push_back(individual_node_information);
    	}
    }
    else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
    {
    	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

    	c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();

    		// Initialise
    		individual_node_information[0] = node_index;
    		individual_node_information[1] = 0;

    		// Examine each epithelial node to see if it should be removed by anoikis and then if it
    		// should be removed by compression-driven apoptosis
    		if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    		{
    			// Determining whether to remove this cell by anoikis

    			if(this->HasCellPoppedUp(node_index))
    			{
    				individual_node_information[1] = 1;
    			}
    		}

    		cells_to_remove.push_back(individual_node_information);
    	}
    }

	return cells_to_remove;
}


/*
 * Cell Killer that kills healthy cells that pop outwards and become detached from
 * the labelled tissue cells, i.e. removal by anoikis
 *
 * Also will remove differentiated cells at the orifice if mSloughOrifice is true
 */
void EpithelialLayerAnoikisCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
		//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

		// Get the information at this timestep for each node index that says whether to remove by anoikis or random apoptosis
		std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

		// Keep a record of how many cells have been removed at this timestep
		this->SetNumberCellsRemoved(cells_to_remove);
		this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

		// Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
		// Loop over these vectors individually and kill any cells that they tell you to

		for (unsigned i=0; i<cells_to_remove.size(); i++)
		{
			if (cells_to_remove[i][1] == 1)
			{
				// Get cell associated to this node
				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
				p_cell->Kill();
			}
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		// Get the information at this timestep for each node index that says whether to remove by anoikis or random apoptosis
		std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

		// Keep a record of how many cells have been removed at this timestep
		this->SetNumberCellsRemoved(cells_to_remove);
		this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

		// Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
		// Loop over these vectors individually and kill any cells that they tell you to

		for (unsigned i=0; i<cells_to_remove.size(); i++)
		{
			if (cells_to_remove[i][1] == 1)
			{
				// Get cell associated to this node
				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
				p_cell->Kill();
			}
		}
	}
}

void EpithelialLayerAnoikisCellKiller::SetNumberCellsRemoved(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
}

unsigned EpithelialLayerAnoikisCellKiller::GetNumberCellsRemoved()
{
	return mCellsRemovedByAnoikis;
}

void EpithelialLayerAnoikisCellKiller::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
		double x_location, y_location;
		c_vector<double, 3> time_and_location;

		// Need to use the node indices to store the locations of where cells are removed
		for (unsigned i=0; i<cellsRemoved.size(); i++)
		{
			if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
			{
				time_and_location[0] = SimulationTime::Instance()->GetTime();

				unsigned node_index = cellsRemoved[i][0];

				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
				x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
				y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

				time_and_location[1] = x_location;
				time_and_location[2] = y_location;

				mLocationsOfAnoikisCells.push_back(time_and_location);
			}
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);
		double x_location, y_location;
		c_vector<double, 3> time_and_location;

		// Need to use the node indices to store the locations of where cells are removed
		for (unsigned i=0; i<cellsRemoved.size(); i++)
		{
			if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
			{
				time_and_location[0] = SimulationTime::Instance()->GetTime();

				unsigned node_index = cellsRemoved[i][0];

				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
				x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
				y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

				time_and_location[1] = x_location;
				time_and_location[2] = y_location;

				mLocationsOfAnoikisCells.push_back(time_and_location);
			}
		}
	}
}

std::vector<c_vector<double,3> > EpithelialLayerAnoikisCellKiller::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}

void EpithelialLayerAnoikisCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
    *rParamsFile << "\t\t\t<CutOffRadius>" << mCutOffRadius << "</CutOffRadius> \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}




#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EpithelialLayerAnoikisCellKiller)
