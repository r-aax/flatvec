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
		links.clear();
		nodes.clear();
		inputs.clear();
		outputs.clear();
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
	/// <param name="act">Act.</param>
	void ControlGraph::register_zmm(int i, std::string act)
	{
		if (i >= 0)
		{
			while (nodes.size() <= i)
			{
				nodes.push_back(NetNode(i));
			}

			nodes[i].acts.push_back(act);
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

		std::cout << "Control graph :" << std::endl;

		std::cout << "\tlinks :" << std::endl;
		for (int i = 0; i < links.size(); ++i)
		{
			int from = links[i][0];
			std::string from_str = (from == -1) ? "GLOB" : std::to_string(from);
			int to = links[i][1];

			std::cout << "\t\t[" << std::setw(5) << from_str << " -> " << std::setw(5) << to << "]" << std::endl;
		}

		std::cout << "\tgraph :" << std::endl;
		std::cout << "\tinputs : ";

		for (int i = 0; i < inputs.size(); ++i)
		{
			std::cout << inputs[i].get_id() << " ";
		}

		std::cout << std::endl;
		std::cout << "\toutputs : ";

		for (int i = 0; i < outputs.size(); ++i)
		{
			std::cout << outputs[i].get_id() << " ";
		}

		std::cout << std::endl;

		for (int i = 0; i < nodes.size(); ++i)
		{
			std::cout << "\t\t" << std::setw(5) << nodes[i].get_id() << ") ";

			for (int j = 0; j < nodes[i].acts.size(); ++j)
			{
				std::cout << "[" << nodes[i].acts[j] << "] ";
			}

			std::cout << std::endl;
			std::cout << "\t\t\tpreds : ";

			for (int j = 0; j < nodes[i].preds.size(); ++j)
			{
				std::cout << nodes[i].preds[j] << " ";
			}

			std::cout << std::endl;
			std::cout << "\t\t\tsuccs : ";

			for (int j = 0; j < nodes[i].succs.size(); ++j)
			{
				std::cout << nodes[i].succs[j] << " ";
			}
				
			std::cout << std::endl;
		}
	}

	/// <summary>
	/// Construct graph.
	/// </summary>
	void ControlGraph::construct_graph()
	{
		for (int i = 0; i < links.size(); ++i)
		{
			int from = links[i][0];
			int to = links[i][1];

			if (from != -1)
			{
				// Node can refer to global variable,
				// it must not be included in the graph.
				nodes[from].succs.push_back(to);
			}

			nodes[to].preds.push_back(from);
		}

		for (int i = 0; i < nodes.size(); ++i)
		{
			NetNode& n = nodes[i];

			if (n.preds.size() == 0)
			{
				if (n.succs.size() == 0)
				{
					// It may be register that was created and rewritten in the end.
					// Example:
					// _mm512 a;            // created
					// void fun(_mm512 &r); // declaration with reference
					// fun(a);              // call function by reference
					// r = <result>;        // write information into register (rewrite it with move semantic).
					if ((n.acts.size() == 2)
						&& (n.acts[0] == "new")
						&& (n.acts[1].find("rewrite", 0) == 0))
					{
						; // ok
					}
					else
					{
						throw std::runtime_error("unexpected node");
					}
				}
				else
				{
					// Input.
					inputs.push_back(n);
				}
			}
			else if (n.succs.size() == 0)
			{
				if (n.acts.back() == "store")
				{
					// Output.
					outputs.push_back(n);
				}
				else
				{
					throw std::runtime_error("output node must has store as its last act");
				}
			}
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
