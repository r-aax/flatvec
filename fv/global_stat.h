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

		// Vectors created.
		int vectors_created;

		// Vectors copied.
		int vectors_copied;

		// Vectors assigned.
		int vectors_assigned;

		// Vector operations.
		int vector_opers_count;

		// Scalar operations.
		int scalar_opers_count;

		// Mask operations.
		int mask_opers_count;

		// Vector operations mask total density.
		double vector_opers_masks_total_density;

		// Time point before calc.
		std::chrono::steady_clock::time_point time_before;

		// Time point in the middle.
		std::chrono::steady_clock::time_point time_middle;

		// Time point after.
		std::chrono::steady_clock::time_point time_after;

	public:

		// Constructor.
		GlobalStat();

		// Clean.
		void clean();

		// Create vector.
		void create_vector();

		// Copy vecror.
		void copy_vector();

		// Assign vector.
		void assign_vector();

		// Append vector operation information.
		void append_vector_oper(int width, int scalar_opers);

		// Append mask operation information.
		void append_mask_oper();

		// Mean masks density.
		double mean_masks_density() const;

		// Fix time before.
		void fix_time_before();

		// Fix time in the middle.
		void fix_time_middle();

		// Fix time in the end.
		void fix_time_after();

		// Time acceleration.
		double real_time_acceleration() const;

		// Print statistics.
		void print();
	};

	// Extern declaration.
	extern GlobalStat GS;
}

#endif
