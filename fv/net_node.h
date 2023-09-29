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
		/// Mask flag.
		/// </summary>
		bool marked;

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

		bool is_root() const;

		/// <summary>
		/// Set mark flag.
		/// </summary>
		inline void mark()
		{
			marked = true;
		}

		/// <summary>
		/// Clear mark flag.
		/// </summary>
		inline void unmark()
		{
			marked = false;
		}

		/// <summary>
		/// Check mark flag.
		/// </summary>
		/// <returns>
		/// true - if node is marked,
		/// false - if node is not marked.
		/// </returns>
		inline bool is_marked() const
		{
			return marked;
		}

		// Vector opers count.
		int vector_opers() const;

		// Check if is blend.
		bool is_blend() const;
	};
}

#endif
