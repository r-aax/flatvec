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

		// Register new register.
		void reg(std::string act,
				 int i);

		// Add link.
		void link(int from1,
				  int to);

		// Add links for 2 sources.
		void link(int from1,
				  int from2,
				  int to);

		// Add links for 3 sources.
		void link(int from1,
				  int from2,
				  int from3,
				  int to);

		// Add links for 4 sources.
		void link(int from1,
				  int from2,
				  int from3,
				  int from4,
				  int to);

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

		// Register without link.
		void reglink(std::string act,
					 int to);

		// Register with 1 link.
		void reglink(std::string act,
					 int from1,
					 int to);

		// Register with 2 links.
		void reglink(std::string act,
					 int from1,
					 int from2,
					 int to);

		// Register with 3 links.
		void reglink(std::string act,
					 int from1,
					 int from2,
					 int from3,
					 int to);

		// Register with 4 links.
		void reglink(std::string act,
					 int from1,
					 int from2,
					 int from3,
					 int from4,
					 int to);

		// Print.
		void print();

		// Construct graph.
		void construct_graph();

		// Unmark nodes.
		void unmark_nodes();

		// Mark node.
		bool mark_node(NetNode& n);

		// Count vector opers.
		double count_vector_opers(NetNode& n);

		// Count vector opers with blend reduce.
		double count_vector_opers_with_blend_reduce(NetNode& n);

		// Analyze graph.
		void analyze();
	};

	/// <summary>
	/// Extern declaration.
	/// </summary>
	extern ControlGraph CG;
}

#endif
