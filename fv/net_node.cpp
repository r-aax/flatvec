#include "stdafx.h"

#include "net_node.h"

namespace fv
{
	/// <summary>
	/// Contstructor.
	/// </summary>
	/// <param name="id_par">Identifier.</param>
	/// <param name="name_par">Name.</param>
	NetNode::NetNode(int id_par,
					 std::string name_par)
		: id(id_par),
		  name(name_par)
	{
	}
}
