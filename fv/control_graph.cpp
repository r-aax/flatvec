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
		  free_id(0)
	{
	}

	/// <summary>
	/// Get new identificator.
	/// </summary>
	/// <returns>New register identifier.</returns>
	int ControlGraph::new_id()
	{
		if (is_active)
		{
			return free_id++;
		}
		else
		{
			return -1;
		}
	}

	/// <summary>
	/// Register new register.
	/// </summary>
	/// <param name="act">Act.</param>
	/// <param name="i">Identifier.</param>
	void ControlGraph::reg(std::string act,
						   int i)
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
	void ControlGraph::link(int from1,
							int to)
	{
		if (is_active)
		{
			links.push_back(std::vector<int> { from1, to });
		}
	}

	/// <summary>
	/// Add links for operation from 2 sources.
	/// </summary>
	/// <param name="from1">First source index.</param>
	/// <param name="from2">Second source index.</param>
	/// <param name="to">Destination index.</param>
	void ControlGraph::link(int from1,
							int from2,
							int to)
	{
		link(from1, to);
		link(from2, to);
	}

	/// <summary>
	/// Add links for operation from 3 sources.
	/// </summary>
	/// <param name="from1">First source index.</param>
	/// <param name="from2">Second source index.</param>
	/// <param name="from3">Third source index.</param>
	/// <param name="to">Destination index.</param>
	void ControlGraph::link(int from1,
							int from2,
							int from3,
							int to)
	{
		link(from1, to);
		link(from2, to);
		link(from3, to);
	}

	/// <summary>
	/// Add links for operation from 4 sources.
	/// </summary>
	/// <param name="from1">First source index.</param>
	/// <param name="from2">Second source index.</param>
	/// <param name="from3">Third source index.</param>
	/// <param name="from4">4th source index.</param>
	/// <param name="to">Destination index.</param>
	void ControlGraph::link(int from1,
							int from2,
							int from3,
							int from4,
							int to)
	{
		link(from1, to);
		link(from2, to);
		link(from3, to);
		link(from4, to);
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
		free_id = 0;
		links.clear();
		nodes.clear();
		inputs.clear();
		outputs.clear();
	}

	/// <summary>
	/// Get free zmm id and shift id.
	/// </summary>
	/// <returns>New ZMM identifier.</returns>
	int ControlGraph::zmm_id()
	{
		return new_id();
	}

	/// <summary>
	/// Get new mask identifier.
	/// </summary>
	/// <returns>New mask identifier.</returns>
	int ControlGraph::mask_id()
	{
		return new_id();
	}

	/// <summary>
	/// Register without link.
	/// </summary>
	/// <param name="act">Action.</param>
	/// <param name="to">Destination id.</param>
	void ControlGraph::reglink(std::string act,
							   int to)
	{
		reg(act, to);
	}

	/// <summary>
	/// Register with 1 link.
	/// </summary>
	/// <param name="act">Action.</param>
	/// <param name="from1">First source id.</param>
	/// <param name="to">Destination id.</param>
	void ControlGraph::reglink(std::string act,
							   int from1,
							   int to)
	{
		reg(act, to);
		link(from1, to);
	}

	/// <summary>
	/// Register with 2 links.
	/// </summary>
	/// <param name="act">Action.</param>
	/// <param name="from1">First source id.</param>
	/// <param name="from2">Second source id.</param>
	/// <param name="to">Destination id.</param>
	void ControlGraph::reglink(std::string act,
							   int from1,
							   int from2,
							   int to)
	{
		reg(act, to);
		link(from1, from2, to);
	}

	/// <summary>
	/// Register with 3 links.
	/// </summary>
	/// <param name="act">Action.</param>
	/// <param name="from1">First source id.</param>
	/// <param name="from2">Second source id.</param>
	/// <param name="from3">Third source id.</param>
	/// <param name="to">Destination id.</param>
	void ControlGraph::reglink(std::string act,
							   int from1,
							   int from2,
							   int from3,
							   int to)
	{
		reg(act, to);
		link(from1, from2, from3, to);
	}

	/// <summary>
	/// Register with 4 links.
	/// </summary>
	/// <param name="act">Action.</param>
	/// <param name="from1">First source id.</param>
	/// <param name="from2">Second source id.</param>
	/// <param name="from3">Third source id.</param>
	/// <param name="from4">4th source id.</param>
	/// <param name="to">Destination id.</param>
	void ControlGraph::reglink(std::string act,
							   int from1,
							   int from2,
							   int from3,
							   int from4,
							   int to)
	{
		reg(act, to);
		link(from1, from2, from3, from4, to);
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
						// Do not throw exception yet.

						// throw std::runtime_error("unexpected node");
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
			}
		}
	}

	/// <summary>
	/// Unmark all nodes.
	/// </summary>
	void ControlGraph::unmark_nodes()
	{
		for (NetNode& n : nodes)
		{
			n.unmark();
		}
	}

	/// <summary>
	/// Mark node with recursive mark of all successors.
	/// </summary>
	/// <param name="n">Node.</param>
	/// <returns>
	/// true - process has completed successfully, no conflicts detected,
	/// false - process was interrupted with conflict detection.
	/// </returns>
	bool ControlGraph::mark_node(NetNode& n)
	{
		if (n.is_marked())
		{
			return false;
		}

		for (int pi : n.preds)
		{
			if (pi != -1)
			{
				if (!mark_node(nodes[pi]))
				{
					return false;
				}
			}
		}

		n.mark();

		return true;
	}

	/// <summary>
	/// Count vector opers.
	/// </summary>
	/// <param name="n">Node.</param>
	/// <returns>Count of vector opers.</returns>
	int ControlGraph::count_vector_opers(NetNode& n)
	{
		int c = 0;

		if (!n.is_marked())
		{
			c += n.vector_opers();

			n.mark();

			for (int pi : n.preds)
			{
				if (pi != -1)
				{
					c += count_vector_opers(nodes[pi]);
				}
			}
		}

		return c;
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

		bool is_finished = false;

		// Warning for hanging nodes.
		for (std::vector<NetNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
			NetNode& n = *it;
			std::string act = n.acts.back();
			bool is_correct_last_use = (act == "store") || (act == "use mask") || (act.find("rewrite") == 0);

			if (n.succs.empty() && !is_correct_last_use)
			{
				// Hanging node detected.
				std::cout << "! Warning ! : hanging node " << n.get_id() << " detected. "
						  << "Last action : " << n.acts.back() << std::endl;

				is_finished = true;
			}
		}

		if (is_finished)
		{
			return;
		}

		// Warning for multiple rewriting.
		for (std::vector<NetNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
			NetNode& n = *it;
			int cnt = 0;

			for (std::vector<std::string>::iterator si = n.acts.begin(); si != n.acts.end(); ++si)
			{
				std::string s = *si;

				if (s.find("rewrite") == 0)
				{
					++cnt;
				}
			}

			if (cnt > 1)
			{
				// Multiple rewrite.
				std::cout << "! Warning ! : mupliple rewrite for " << n.get_id() << " detected. "
						  << "Last action : " << n.acts.back() << std::endl;

				is_finished = true;
			}
			else if (cnt == 1)
			{
				if ((n.acts.size() == 2) && (n.acts[0].find("new") == 0))
				{
					; // ok
				}
				else
				{
					// More complex template can take place.
					// std::cout << "! Warning ! : wrong rewrite template for " << n.get_id() << " detected. "
					//           << "Last action : " << n.acts.back() << std::endl;

					//is_finished = true;
				}
			}
		}

		if (is_finished)
		{
			return;
		}

		//
		// Analyze all roots.
		//
		/*
		for (NetNode& n : nodes)
		{
			if (n.is_root())
			{
				unmark_nodes();

				if (mark_node(n))
				{
					; // ok
				}
				else
				{
					// Loop is detected.

					std::cout << "! Warning ! : loop for final node " << n.get_id() << " is detected." << std::endl;

					is_finished = true;
				}
			}

			if (is_finished)
			{
				return;
			}
		}
		*/

		//
		// Count opers with trees.
		//

		int vector_opers_count = 0;

		unmark_nodes();

		for (NetNode& n : nodes)
		{

			if (n.is_root())
			{
				vector_opers_count += count_vector_opers(n);
			}
		}

		std::cout << "OPERS COUNT : " << vector_opers_count << std::endl;
	}
}
