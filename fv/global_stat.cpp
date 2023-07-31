#include "stdafx.h"

#include "global_stat.h"

namespace fv
{
	// Global variable.
	GlobalStat GS;

	/// <summary>
	/// Default constructor.
	/// </summary>
	GlobalStat::GlobalStat()
		: vector_opers_count(0),
		  scalar_opers_count(0),
		  mask_opers_count(0),
		  vector_opers_masks_total_density(0.0)
	{
	}

	/// <summary>
	/// Clean stat.
	/// </summary>
	void GlobalStat::clean()
	{
		vector_opers_count = 0;
		scalar_opers_count = 0;
		mask_opers_count = 0;
		vector_opers_masks_total_density = 0.0;
	}

	/// <summary>
	/// Append vector oper information.
	/// </summary>
	/// <param name="width"></param>
	/// <param name="scalar_opers"></param>
	void GlobalStat::append_vector_oper(int width,
										int scalar_opers)
	{
		vector_opers_count++;
		scalar_opers_count += scalar_opers;
		vector_opers_masks_total_density += (static_cast<double>(scalar_opers) / static_cast<double>(width));
	}

	/// <summary>
	/// Append mask oper information.
	/// </summary>
	void GlobalStat::append_mask_oper()
	{
		mask_opers_count++;
	}

	/// <summary>
	/// Mean masks density.
	/// </summary>
	/// <returns>Mean masks density</returns>
	double GlobalStat::mean_masks_density() const
	{
		return vector_opers_masks_total_density / static_cast<double>(vector_opers_count);
	}

	/// <summary>
	/// Print statistics.
	/// </summary>
	void GlobalStat::print()
	{
		std::cout << "Global stat:" << std::endl;
		std::cout << "  vector opers           : " << vector_opers_count << std::endl;
		std::cout << "  scalar opers           : " << scalar_opers_count << std::endl;
		std::cout << "  mask opers             : " << mask_opers_count << std::endl;
		std::cout << "  effective vector width : " << (static_cast<double>(scalar_opers_count) / static_cast<double>(vector_opers_count)) << std::endl;
		std::cout << "  mean masks density     : " << mean_masks_density() << std::endl;
	}
}
