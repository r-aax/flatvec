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
		srcs.clear();
		links.clear();
		nodes.clear();
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
	/// Register new ZMM.
	/// </summary>
	/// <param name="i">ZMM number.</param>
	/// <param name="src">Source.</param>
	void ControlGraph::register_zmm(int i, std::string src)
	{
		if (i >= 0)
		{
			while (srcs.size() <= i)
			{
				srcs.push_back("");
			}

			if (srcs[i] == "")
			{
				srcs[i] = src;
			}
			else
			{
				srcs[i] = srcs[i] + " | " + src;
			}
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

		std::cout << "\tsrcs:" << std::endl;
		for (int i = 0; i < srcs.size(); ++i)
		{
			std::cout << "\t\t" << std::setw(5) << i << " : " << srcs[i] << std::endl;
		}

		std::cout << "\tlinks:" << std::endl;
		for (int i = 0; i < links.size(); ++i)
		{
			int from = links[i][0];
			std::string from_str = (from == -1) ? "GLOB" : std::to_string(from);
			int to = links[i][1];

			std::cout << "\t\t[" << std::setw(5) << from_str << " -> " << std::setw(5) << to
				      << "] // " << srcs[to] << std::endl;
		}

		std::cout << "\tgraph:" << std::endl;
		for (int i = 0; i < nodes.size(); ++i)
		{
			std::cout << "\t\t" << std::setw(5) << nodes[i].get_id() << ") " << nodes[i].get_name() << std::endl;
			std::cout << "\t\t\tpred : ";

			for (int j = 0; j < nodes[i].pred.size(); ++j)
			{
				std::cout << nodes[i].pred[j] << " ";
			}

			std::cout << std::endl;
			std::cout << "\t\t\tsucc : ";

			for (int j = 0; j < nodes[i].succ.size(); ++j)
			{
				std::cout << nodes[i].succ[j] << " ";
			}
				
			std::cout << std::endl;
		}
	}

	/// <summary>
	/// Construct graph.
	/// </summary>
	void ControlGraph::construct_graph()
	{
		for (int i = 0; i < srcs.size(); ++i)
		{
			nodes.push_back(NetNode(i, srcs[i]));
		}

		for (int i = 0; i < links.size(); ++i)
		{
			int from = links[i][0];
			int to = links[i][1];

			if (from != -1)
			{
				// Node can refer to global variable,
				// it must not be included in the graph.
				nodes[from].succ.push_back(to);
			}

			nodes[to].pred.push_back(from);
		}
	}

	/// <summary>
	/// Analyze control graph.
	/// </summary>
	void ControlGraph::analyze()
	{
		if (!is_active)
		{
			// Analyze only in active mode.
			return;
		}

		construct_graph();
		print();
	}
}
