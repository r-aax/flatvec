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
		  vector_opers_count(0),
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
		vectors_created = 0;
		vectors_copied = 0;
		vector_opers_count = 0;
		scalar_opers_count = 0;
		mask_opers_count = 0;
		vector_opers_masks_total_density = 0.0;
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
		std::cout << "Global stat:" << std::endl;
		std::cout << "  vectors created          : " << vectors_created << std::endl;
		std::cout << "  vectors copied           : " << vectors_copied << std::endl;
		std::cout << "  vector opers             : " << vector_opers_count << std::endl;
		std::cout << "  scalar opers             : " << scalar_opers_count << std::endl;
		std::cout << "  mask opers               : " << mask_opers_count << std::endl;
		std::cout << "  theoretical acceleration : " << (static_cast<double>(scalar_opers_count)
													  / static_cast<double>(vector_opers_count)) << std::endl;
		std::cout << "  mean masks density       : " << mean_masks_density() << std::endl;

#ifdef LINUX_ICC_BUILD
		std::cout << "  real time acceleration   : " << std::fixed << real_time_acceleration() << std::endl;
#endif

		std::cout << "--------------------------" << std::endl;
	}
}
