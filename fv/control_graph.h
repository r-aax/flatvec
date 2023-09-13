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
		/// Free ZMM register id.
		/// </summary>
		int free_zmm_id;

	public:

		// Constructor.
		ControlGraph();

		// Clear.
		void clear();

		// Get ZMM id.
		int zmm_id();

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
