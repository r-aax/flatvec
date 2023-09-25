#include "stdafx.h"

#ifndef CONTROL_GRAPH_H
#define CONTROL_GRAPH_H

#include "zmm.h"
#include "net_node.h"

namespace fv
{
	/// <summary>
	/// Control graph of flat loop body.
	/// </summary>
	class ControlGraph
	{

	private:

		/// <summary>
		/// Active flag.
		/// </summary>
		bool is_active;

		/// <summary>
		/// Free register id (ZMM or Mask).
		/// </summary>
		int free_id;

		/// <summary>
		/// Links between vector registers.
		/// </summary>
		std::vector<std::vector<int>> links;

		/// <summary>
		/// Nodes of graph.
		/// </summary>
		std::vector<NetNode> nodes;

		/// <summary>
		/// Inputs vector.
		/// </summary>
		std::vector<NetNode> inputs;

		/// <summary>
		/// Outputs vector.
		/// </summary>
		std::vector<NetNode> outputs;

		// Get new id.
		int new_id();

	public:

		// Constructor.
		ControlGraph();

		// Switch on for create graph.
		void switch_on();

		// Switch off.
		void switch_off();

		// Clear.
		void clear();

		// Get ZMM id.
		int zmm_id();

		// Get mask id.
		int mask_id();

		// Register new ZMM.
		void reg(int i, std::string act);

		// Add link.
		void link(int from, int to);

		// Add links for 2 sources.
		void link2(int from1, int from2, int to);

		// Add links for 3 sources.
		void link3(int from1, int from2, int from3, int to);

		// Add links for 4 sources.
		void link4(int from1, int from2, int from3, int from4, int to);

		// Print.
		void print();

		// Construct graph.
		void construct_graph();

		// Analyze graph.
		void analyze();
	};

	/// <summary>
	/// Extern declaration.
	/// </summary>
	extern ControlGraph CG;
}

#endif
