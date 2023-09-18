#include "stdafx.h"

#ifndef NET_NODE_H
#define NET_NODE_H

namespace fv
{
	/// <summary>
	/// Node of the net.
	/// </summary>
	class NetNode
	{

	private:

		/// <summary>
		/// Identifier.
		/// </summary>
		int id;

		/// <summary>
		/// Name.
		/// </summary>
		std::string name;

	public:

		/// <summary>
		/// Predeccessors vector.
		/// </summary>
		std::vector<int> pred;

		/// <summary>
		/// Successors vector.
		/// </summary>
		std::vector<int> succ;

	public:

		NetNode(int id_par,
				std::string name_par);

		/// <summary>
		/// Get identifier.
		/// </summary>
		/// <returns>Identifier.</returns>
		inline int get_id() const
		{
			return id;
		}

		/// <summary>
		/// Get name.
		/// </summary>
		/// <returns>Name.</returns>
		inline std::string get_name() const
		{
			return name;
		}
	};
}

#endif
