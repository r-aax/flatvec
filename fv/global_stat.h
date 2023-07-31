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

		// Vector operations.
		int vector_opers_count;

		// Scalar operations.
		int scalar_opers_count;

		// Mask operations.
		int mask_opers_count;

		// Vector operations mask total density.
		double vector_opers_masks_total_density;

	public:

		// Constructor.
		GlobalStat();

		// Clean.
		void clean();

		// Append vector operation information.
		void append_vector_oper(int width, int scalar_opers);

		// Append mask operation information.
		void append_mask_oper();

		// Mean masks density.
		double mean_masks_density() const;

		// Print statistics.
		void print();
	};

	// Extern declaration.
	extern GlobalStat GS;
}

#endif
