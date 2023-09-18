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
		: is_active(false),
		  free_zmm_id(0)
	{
	}

	/// <summary>
	/// Switch on.
	/// </summary>
	void ControlGraph::switch_on()
	{
		is_active = true;
		clear();
	}

	/// <summary>
	/// Switch off.
	/// </summary>
	void ControlGraph::switch_off()
	{
		is_active = false;
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
		if (is_active)
		{
			return free_zmm_id++;
		}
		else
		{
			return -1;
		}
	}

	/// <summary>
	/// Add link.
	/// </summary>
	/// <param name="from">Index of from register.</param>
	/// <param name="to">Index of to register.</param>
	void ControlGraph::add_link(int from, int to)
	{
		if (is_active)
		{
			links.push_back(std::vector<int> { from, to });
		}
	}

	/// <summary>
	/// Print info.
	/// </summary>
	void ControlGraph::print()
	{
		if (!is_active)
		{
			// Print information only in active mode.
			return;
		}

		std::cout << "Control graph:" << std::endl;

		std::cout << "\tlinks:" << std::endl;
		for (int i = 0; i < links.size(); ++i)
		{
			std::cout << "\t\t[" << links[i][0] << " -> " << links[i][1] << "]" << std::endl;
		}
	}
}
