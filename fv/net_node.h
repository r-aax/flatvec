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

	public:

		/// <summary>
		/// Vector of actions.
		/// </summary>
		std::vector<std::string> acts;

		/// <summary>
		/// Predeccessors vector.
		/// </summary>
		std::vector<int> preds;

		/// <summary>
		/// Successors vector.
		/// </summary>
		std::vector<int> succs;

	public:

		NetNode(int id_par);

		/// <summary>
		/// Get identifier.
		/// </summary>
		/// <returns>Identifier.</returns>
		inline int get_id() const
		{
			return id;
		}
	};
}

#endif
