#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

/* Class of boundary conditions that fixes cells that sit outside of the plane defined by
 * the point and normal vectors, rather than pushing them back in. Also labels the fixed cells,
 * for visualisation purposes.
 */
template<unsigned DIM>
FixedRegionPlaneBoundaryCondition<DIM>::FixedRegionPlaneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                    c_vector<double, DIM> point,
                                                    c_vector<double, DIM> normal)
        : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
          mPointOnPlane(point),
          mUseJiggledNodesOnPlane(false)
{
    assert(norm_2(normal) > 0.0);
    mNormalToPlane = normal/norm_2(normal);
}

template<unsigned DIM>
const c_vector<double, DIM>& FixedRegionPlaneBoundaryCondition<DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template<unsigned DIM>
const c_vector<double, DIM>& FixedRegionPlaneBoundaryCondition<DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}


template<unsigned DIM>
void FixedRegionPlaneBoundaryCondition<DIM>::SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane)
{
    mUseJiggledNodesOnPlane = useJiggledNodesOnPlane;
}

template<unsigned DIM>
bool FixedRegionPlaneBoundaryCondition<DIM>::GetUseJiggledNodesOnPlane()
{
    return mUseJiggledNodesOnPlane;
}

template<unsigned DIM>
void FixedRegionPlaneBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(this->mpCellPopulation)==NULL)
    {
        EXCEPTION("PlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
            || (dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation))
            || (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation)) );

    // THis is a magic number.
    double max_jiggle = 1e-4;

    if (DIM != 1)
    {
        if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

                c_vector<double, DIM> current_location = p_node->rGetLocation();

                //Get the previous location as well
    			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
    			c_vector<double, DIM> previous_location = it->second;

                c_vector<double, DIM> nearest_point;

                //Get the vector from the point on the plane to the current and previous location
                c_vector<double, DIM> plane_to_current_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, current_location);
                c_vector<double, DIM> plane_to_previous_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, previous_location);

                double current_signed_distance = inner_prod(plane_to_current_point, mNormalToPlane);
                double previous_signed_distance = inner_prod(plane_to_previous_point, mNormalToPlane);

                if (current_signed_distance > 0.0)
                {
                	//If the cell outside of the boundary has not been given a CellLabel
                	if (!cell_iter->template HasCellProperty<CellLabel>())
                	{
                			if (mUseJiggledNodesOnPlane)
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = current_location;

                				//Differentiate cell, otherwise something weird will happen with proliferation
                				boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
                				cell_iter->SetCellProliferativeType(p_diff_type);

                				//Label the fixed cell for visualisation
                				boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
                				cell_iter->AddCellProperty(p_label);
                			}
                			else
                			{
                				//Assign nearest point, which is just the current location
                				nearest_point = current_location;

                				//Differentiate cell, otherwise something weird will happen with proliferation
                				boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
                				cell_iter->SetCellProliferativeType(p_diff_type);

                				//Label the fixed cell for visualisation
                				boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
                				cell_iter->AddCellProperty(p_label);
                			}
                			p_node->rGetModifiableLocation() = nearest_point;
                	}
                	else //Otherwise it will already have one, so don't re-label it
                	{
                		if (previous_signed_distance > 0.0)
                		{
                			if (mUseJiggledNodesOnPlane)
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = previous_location;
                			}
                			else
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = previous_location;
                			}
                			p_node->rGetModifiableLocation() = nearest_point;
                		}
                		else
                		{
                			if (mUseJiggledNodesOnPlane)
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = current_location;
                			}
                			else
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = current_location;
                			}
                			p_node->rGetModifiableLocation() = nearest_point;
                		}
                	}
                }
                else if ( (current_signed_distance < 0.0 )&&(cell_iter->template HasCellProperty<CellLabel>()) )
                {
            		if (mUseJiggledNodesOnPlane)
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		else
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
        else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation))
        {
        	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
        			cell_iter != this->mpCellPopulation->End();
        			++cell_iter)
        	{
        		unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        		Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

                c_vector<double, DIM> current_location = p_node->rGetLocation();

                //Get the previous location as well
    			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
    			c_vector<double, DIM> previous_location = it->second;

                c_vector<double, DIM> nearest_point;

                //Get the vector from the point on the plane to the current and previous location
                c_vector<double, DIM> plane_to_current_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, current_location);
                c_vector<double, DIM> plane_to_previous_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, previous_location);

                double current_signed_distance = inner_prod(plane_to_current_point, mNormalToPlane);
                double previous_signed_distance = inner_prod(plane_to_previous_point, mNormalToPlane);

        		if (current_signed_distance > 0.0)
        		{
        			//If the cell outside of the boundary has not been given a CellLabel
        			if (!cell_iter->template HasCellProperty<CellLabel>())
        			{
        				if (mUseJiggledNodesOnPlane)
        				{
        					//Assign nearest point, which is just the current location
        					nearest_point = current_location;

        					//Differentiate cell, otherwise something weird will happen with proliferation
        					boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        					cell_iter->SetCellProliferativeType(p_diff_type);

        					//Label the fixed cell for visualisation
        					boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
        					cell_iter->AddCellProperty(p_label);
        				}
        				else
        				{
        					//Assign nearest point, which is just the urrnt
        					nearest_point = current_location;

        					//Differentiate cell, otherwise something weird will happen with proliferation
        					boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        					cell_iter->SetCellProliferativeType(p_diff_type);

        					//Label the fixed cell for visualisation
        					boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
        					cell_iter->AddCellProperty(p_label);
        				}
        				p_node->rGetModifiableLocation() = nearest_point;
        			}
        			else //Otherwise it will already have one, so don't re-label it
        			{
        				if (previous_signed_distance > 0.0)
        				{
        					if (mUseJiggledNodesOnPlane)
        					{
        						//Assign nearest point, which is just the old location
        						nearest_point = previous_location;
        					}
        					else
        					{
        						//Assign nearest point, which is just the old location
        						nearest_point = previous_location;
        					}
        					p_node->rGetModifiableLocation() = nearest_point;
        				}
        				else
        				{
           					if (mUseJiggledNodesOnPlane)
            					{
            						//Assign nearest point, which is just the old location
            						nearest_point = current_location;
            					}
            					else
            					{
            						//Assign nearest point, which is just the old location
            						nearest_point = current_location;
            					}
            					p_node->rGetModifiableLocation() = nearest_point;
        				}
        			}
        		}
                else if ( (current_signed_distance < 0.0 )&&(cell_iter->template HasCellProperty<CellLabel>()) )
                {
            		if (mUseJiggledNodesOnPlane)
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		else
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		p_node->rGetModifiableLocation() = nearest_point;
                }
        	}
        }
        else
        {
            assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation));

            VertexBasedCellPopulation<DIM>* pStaticCastCellPopulation = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);

            // Iterate over all nodes and update their positions according to the boundary conditions
            unsigned num_nodes = pStaticCastCellPopulation->GetNumNodes();
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<DIM>* p_node = pStaticCastCellPopulation->GetNode(node_index);
        		c_vector<double, DIM> nearest_point;

                CellPtr cell_iter = pStaticCastCellPopulation->GetCellUsingLocationIndex(node_index);

                c_vector<double, DIM> current_location = p_node->rGetLocation();

                //Get the previous location as well
    			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
    			c_vector<double, DIM> previous_location = it->second;

                //Get the vector from the point on the plane to the current and previous location
                c_vector<double, DIM> plane_to_current_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, current_location);
                c_vector<double, DIM> plane_to_previous_point = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, previous_location);

                double current_signed_distance = inner_prod(plane_to_current_point, mNormalToPlane);
                double previous_signed_distance = inner_prod(plane_to_previous_point, mNormalToPlane);

                if (current_signed_distance > 0.0)
                {
                	//If the cell outside of the boundary has not been given a CellLabel
                	if (!cell_iter->template HasCellProperty<CellLabel>())
                	{
                		if (mUseJiggledNodesOnPlane)
                		{
                			//Assign nearest point, which is just the current
                			nearest_point = current_location;

        					//Differentiate cell, otherwise something weird will happen with proliferation
        					boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        					cell_iter->SetCellProliferativeType(p_diff_type);

                			//Label the fixed cell for visualisation
                			boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
                			cell_iter->AddCellProperty(p_label);
                		}
                		else
                		{
                			//Assign nearest point, which is just the old location
                			nearest_point = current_location;

        					//Differentiate cell, otherwise something weird will happen with proliferation
        					boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        					cell_iter->SetCellProliferativeType(p_diff_type);

                			//Label the fixed cell for visualisation
                			boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
                			cell_iter->AddCellProperty(p_label);
                		}
                		p_node->rGetModifiableLocation() = nearest_point;
                	}
                	else //Otherwise it already has a CellLabel and thus we do not relabel it.
                	{
                		if (previous_signed_distance > 0.0)
                		{
                			if (mUseJiggledNodesOnPlane)
                			{
                				//Assign nearest point, which is just the old location
                				typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
                				nearest_point = it->second;
                			}
                			else
                			{
                				//Assign nearest point, which is just the old location
                				typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
                				nearest_point = it->second;
                			}
                			p_node->rGetModifiableLocation() = nearest_point;
                		}
                		else
                		{
                			if (mUseJiggledNodesOnPlane)
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = current_location;
                			}
                			else
                			{
                				//Assign nearest point, which is just the old location
                				nearest_point = current_location;
                			}
                			p_node->rGetModifiableLocation() = nearest_point;
                		}
                	}
                }
                else if ( (current_signed_distance < 0.0 )&&(cell_iter->template HasCellProperty<CellLabel>()) )
                {
            		if (mUseJiggledNodesOnPlane)
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		else
            		{
            			//Assign nearest point, which is just the old location
            			typename std::map<Node<DIM>*, c_vector<double, DIM> >::const_iterator it = rOldLocations.find(p_node);
            			nearest_point = it->second;
            		}
            		p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        //PlaneBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template<unsigned DIM>
bool FixedRegionPlaneBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (DIM == 1)
    {
        EXCEPTION("PlaneBoundaryCondition is not implemented in 1D");
    }
    else
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

            c_vector<double, DIM> plane_to_cell_location = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(mPointOnPlane, cell_location);

            if ( (inner_prod(plane_to_cell_location, mNormalToPlane) > 0.0)&&(!cell_iter->template HasCellProperty<CellLabel>()) )
            {
                condition_satisfied = false;
                break;
            }
        }
    }

    return condition_satisfied;
}

template<unsigned DIM>
void FixedRegionPlaneBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[DIM-1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mNormalToPlane[index] << ",";
    }
    *rParamsFile << mNormalToPlane[DIM-1] << "</NormalToPlane>\n";
    *rParamsFile << "\t\t\t<UseJiggledNodesOnPlane>" << mUseJiggledNodesOnPlane << "</UseJiggledNodesOnPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class FixedRegionPlaneBoundaryCondition<1>;
template class FixedRegionPlaneBoundaryCondition<2>;
template class FixedRegionPlaneBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedRegionPlaneBoundaryCondition)
