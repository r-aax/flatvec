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
		  mask_opers_count(0)
	{
	}

	/// <summary>
	/// Clean stat.
	/// </summary>
	void GlobalStat::clean()
	{
		vector_opers_count = 0;
		mask_opers_count = 0;
	}

	/// <summary>
	/// Append vector oper information.
	/// </summary>
	void GlobalStat::append_vector_oper()
	{
		vector_opers_count++;
	}

	/// <summary>
	/// Append mask oper information.
	/// </summary>
	void GlobalStat::append_mask_oper()
	{
		mask_opers_count++;
	}

	/// <summary>
	/// Print statistics.
	/// </summary>
	void GlobalStat::print()
	{
		std::cout << "Global stat:" << std::endl;
		std::cout << "  vector opers : " << vector_opers_count << std::endl;
		std::cout << "  mask_opers   : " << mask_opers_count << std::endl;
	}
}
