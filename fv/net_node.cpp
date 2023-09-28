#include "stdafx.h"

#include "net_node.h"

namespace fv
{
	/// <summary>
	/// Contstructor.
	/// </summary>
	/// <param name="id_par">Identifier.</param>
	NetNode::NetNode(int id_par)
		: id(id_par)
	{
	}

	/// <summary>
	/// Check if node is root.
	/// </summary>
	/// <returns>
	/// true - if the node is a root,
	/// false - if the node is not a root.
	/// </returns>
	bool NetNode::is_root() const
	{
		return succs.empty();
	}

	/// <summary>
	/// Vector opers count.
	/// </summary>
	/// <returns>Vector opers count.</returns>
	int NetNode::vector_opers() const
	{
		int c = 0;

		for (std::string s : acts)
		{
			if ((s == "load") || (s == "store")
				|| (s == "arith1") || (s == "arith2") || (s == "arith3")
				|| (s == "blend")
				|| (s == "cmp"))
			{
				++c;
			}
		}

		return c;
	}

}
