#include "stdafx.h"

#ifndef CONTROL_GRAPH_H
#define CONTROL_GRAPH_H

#include "zmm.h"

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
		/// Free ZMM register id.
		/// </summary>
		int free_zmm_id;

		/// <summary>
		/// Sources of zmm registers.
		/// </summary>
		std::vector<std::string> srcs;

		/// <summary>
		/// Links between vector registers.
		/// </summary>
		std::vector<std::vector<int>> links;

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

		// Register new ZMM.
		void register_zmm(int i, std::string src);

		// Add link.
		void add_link(int from, int to);

		// Print.
		void print();
	};

	/// <summary>
	/// Extern declaration.
	/// </summary>
	extern ControlGraph CG;
}

#endif
