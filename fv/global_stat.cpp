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
		: vectors_created(0),
		  vectors_copied(0),
		  vectors_moved(0),
		  vector_opers_count(0),
		  scalar_opers_count(0),
		  mask_opers_count(0),
		  vector_opers_masks_total_density(0.0),
		  blend_reduce_coefficient(0.0)
	{
	}

	/// <summary>
	/// Clean stat.
	/// </summary>
	void GlobalStat::clear()
	{
		vectors_created = 0;
		vectors_copied = 0;
		vectors_moved = 0;
		vector_opers_count = 0;
		scalar_opers_count = 0;
		mask_opers_count = 0;
		vector_opers_masks_total_density = 0.0;
		blend_reduce_coefficient = 0.0;
	}

	/// <summary>
	/// Create vector.
	/// </summary>
	void GlobalStat::create_vector()
	{
		++vectors_created;
	}

	/// <summary>
	/// Copy vector.
	/// </summary>
	void GlobalStat::copy_vector()
	{
		++vectors_copied;
	}

	/// <summary>
	/// Move vector.
	/// </summary>
	void GlobalStat::move_vector()
	{
		++vectors_moved;
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
		vector_opers_masks_total_density += (static_cast<double>(scalar_opers)
											 / static_cast<double>(width));
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
		return vector_opers_masks_total_density
			   / static_cast<double>(vector_opers_count);
	}

	/// <summary>
	/// Fix time before.
	/// </summary>
	void GlobalStat::fix_time_before()
	{
		time_before = std::chrono::steady_clock::now();
	}

	/// <summary>
	/// Fix time in the middle.
	/// </summary>
	void GlobalStat::fix_time_middle()
	{
		time_middle = std::chrono::steady_clock::now();
	}

	/// <summary>
	/// Fix time in the end.
	/// </summary>
	void GlobalStat::fix_time_after()
	{
		time_after = std::chrono::steady_clock::now();
	}

	/// <summary>
	/// Time acceleration.
	/// </summary>
	/// <returns>Time acceleration (ms). </returns>
	double GlobalStat::real_time_acceleration() const
	{
		std::chrono::duration<double> s = time_middle - time_before;
		std::chrono::duration<double> v = time_after - time_middle;

		return s.count() / v.count();
	}

	/// <summary>
	/// Print statistics.
	/// </summary>
	void GlobalStat::print()
	{
		double th_a = (static_cast<double>(scalar_opers_count)
					  / static_cast<double>(vector_opers_count));

		std::cout << "Global stat:" << std::endl;
		std::cout << "\tvectors created          : " << vectors_created << std::endl;
		std::cout << "\tvectors copied           : " << vectors_copied << std::endl;
		std::cout << "\tvectors moved            : " << vectors_moved << std::endl;
		std::cout << "\tvector opers             : "
			      << vector_opers_count << std::endl;
		std::cout << "\tscalar opers             : "
			      << scalar_opers_count << std::endl;
		std::cout << "\tmask opers               : "
			      << mask_opers_count << std::endl;
		std::cout << "\ttheoretical acceleration : " << th_a << std::endl;
		
		if (blend_reduce_coefficient > 0.0)
		{
			std::cout << "\t      (with blend corr.) : "
				      << (th_a * blend_reduce_coefficient) << std::endl;
		}
		else
		{
			std::cout << "\t      (with blend corr.) : "
				      << "not available" << std::endl;
		}

		std::cout << "\tmean masks density       : "
			      << mean_masks_density() << std::endl;

#ifdef LINUX_ICC_BUILD
		std::cout << "\treal time acceleration   : "
			      << std::fixed << real_time_acceleration() << std::endl;
#endif

		std::cout << "--------------------------" << std::endl;
	}
}
