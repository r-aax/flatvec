#include "stdafx.h"

#include "control_graph.h"

namespace fv
{
	/// <summary>
	/// One control graph for all cases.
	/// </summary>
	ControlGraph CG;

	/// <summary>
	/// Default constructor.
	/// </summary>
	ControlGraph::ControlGraph()
		: free_zmm_id(0)
	{
	}

	/// <summary>
	/// Get free ZMM id and shift.
	/// </summary>
	void ControlGraph::clear()
	{
		free_zmm_id = 0;
	}

	/// <summary>
	/// Get free zmm id and shift id.
	/// </summary>
	int ControlGraph::zmm_id()
	{
		return free_zmm_id++;
	}

	/// <summary>
	/// Add link.
	/// </summary>
	/// <param name="from">Index of from register.</param>
	/// <param name="to">Index of to register.</param>
	void ControlGraph::add_link(int from, int to)
	{
		std::cout << from << " -> " << to << std::endl;
	}

	/// <summary>
	/// Print info.
	/// </summary>
	void ControlGraph::print()
	{
		std::cout << "Control graph:" << std::endl;

		std::cout << "\tlinks:" << std::endl;
	}
}
