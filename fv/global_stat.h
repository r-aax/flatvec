#include "stdafx.h"

#ifndef GLOBAL_STAT_H
#define GLOBAL_STAT_H

namespace fv
{
	/// <summary>
	/// Global statistics.
	/// </summary>
	class GlobalStat
	{

	private:

		// Vectors operations.
		int vector_opers_count;

		// Mask operations.
		int mask_opers_count;

	public:

		// Constructor.
		GlobalStat();

		// Clean.
		void clean();

		// Append vector operation information.
		void append_vector_oper();

		// Append mask operation information.
		void append_mask_oper();

		// Print statistics.
		void print();
	};

	// Extern declaration.
	extern GlobalStat GS;
}

#endif
