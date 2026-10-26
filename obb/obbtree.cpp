#include <algorithm>
#include <numeric>
#include "obbtree.h"

extern "C" {
	int tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3],
						  const float U0[3], const float U1[3], const float U2[3]);
}

namespace obb
{
	static bool g_Diagnostics = false;

	void set_diagnostics(bool v)
	{
		g_Diagnostics = v;
	}

	bool OBBNode::verify() const
	{
		Vector3f error;
		for (const auto& i : pointIndices)
		{
			const auto& p = points->at(i);
			if (!contains(box, p, error))
				return false;
		}
		if (leaf)
			return true;
		return left->verify() && right->verify();
	}

	OBBTree::OBBTree()
		: m_Points(std::make_shared<PointList>())
	{
		m_Root = std::make_shared<OBBNode>(m_Points, 1);
	}

	void OBBTree::add_triangle(const Triangle& tri)
	{
		m_Points->push_back(tri.a);
		m_Points->push_back(tri.b);
		m_Points->push_back(tri.c);
	}

	static void add_triangle_points(std::vector<PointIndex>& indices, PointIndex index)
	{
		// each triangle has 3 points.  index can have any 3 modulo value
		index -= (index % 3);
		indices.push_back(index);
		indices.push_back(index + 1);
		indices.push_back(index + 2);
	}

	size_t OBBTree::split_node(OBBNodePtr node,
							   size_t max_points_per_leaf)
	{
		if (node->pointIndices.size() <= max_points_per_leaf)
		{
			node->leaf = true;
			if (g_Diagnostics)
			{
				if (!node->verify())
					throw std::runtime_error("OBBTree build verification failed on leaf node");
			}
			return node->depth; // leaf node
		}
		Eigen::Index split_axis;
		node->box.half_size.maxCoeff(&split_axis);
		node->left = std::make_shared<OBBNode>(node->points, node->depth + 1);
		node->right = std::make_shared<OBBNode>(node->points, node->depth + 1);
		for(size_t j=0;j<node->pointIndices.size();j+=3)
		{
			bool to_left = false;
			bool to_right = false;
			for(size_t k=0;k<3;++k)
			{
				size_t i = node->pointIndices[j+k];
				const Vector3f& p = node->points->at(i);
				Vector3f local = to_local(node->box, p, false);
				if (local[split_axis] <= 0.0f) to_left = true;
				else to_right = true;
			}
			if (to_left) add_triangle_points(node->left->pointIndices, node->pointIndices[j]);
			else add_triangle_points(node->right->pointIndices, node->pointIndices[j]);
		}
		if (node->left->pointIndices.empty() || node->right->pointIndices.empty() || 
			node->left->pointIndices.size()==node->pointIndices.size() ||
			node->right->pointIndices.size() == node->pointIndices.size())
		{
			// failed to split, make leaf
			node->leaf = true;
			node->left = nullptr;
			node->right = nullptr;
			return node->depth;
		}
		node->left->box = compute_obb(*(node->points), node->left->pointIndices);
		node->right->box = compute_obb(*(node->points), node->right->pointIndices);
		size_t rc = std::max(split_node(node->left, max_points_per_leaf),
						split_node(node->right, max_points_per_leaf));
		if (g_Diagnostics)
		{
			if (!node->left->verify())
				throw std::runtime_error("OBBTree build verification failed on left node");
			if (!node->right->verify())
				throw std::runtime_error("OBBTree build verification failed on right node");
		}
		return rc;
	}

	size_t OBBTree::build(size_t max_triangles_per_leaf)
	{
		size_t max_points_per_leaf = max_triangles_per_leaf * 3;
		m_Root->box = compute_obb(*m_Points);
		m_Root->pointIndices.resize(m_Points->size());
		std::iota(m_Root->pointIndices.begin(), m_Root->pointIndices.end(), 0);
		return split_node(m_Root, max_points_per_leaf);
	}

	IntersectionResult OBBTree::bf_check_intersection(const OBBTree& other) const
	{
		IntersectionResult result;
		result.hit = false;
		for (size_t i = 0; i < m_Points->size(); i += 3)
		{
			for(size_t j = 0; j < other.m_Points->size(); j += 3)
			{
				const Vector3f& A0 = m_DynamicTransform.apply(m_Points->at(i + 0));
				const Vector3f& A1 = m_DynamicTransform.apply(m_Points->at(i + 1));
				const Vector3f& A2 = m_DynamicTransform.apply(m_Points->at(i + 2));
				const Vector3f& B0 = other.m_DynamicTransform.apply(other.m_Points->at(j + 0));
				const Vector3f& B1 = other.m_DynamicTransform.apply(other.m_Points->at(j + 1));
				const Vector3f& B2 = other.m_DynamicTransform.apply(other.m_Points->at(j + 2));
				if (tri_tri_intersect(
					A0.data(), A1.data(), A2.data(),
					B0.data(), B1.data(), B2.data()))
				{
					result.hit = true;
					result.triangle_index1 = i / 3;
					result.triangle_index2 = j / 3;
					return result;
				}
			}
		}
		return result;
	}

	bool OBBTree::check_intersection_leaves(const OBBNodePtr& node,
											const Transform& node_transform,
											const OBBNodePtr& other_node,
											const Transform& other_transform,
											IntersectionResult& result)
	{
		for(size_t i=0;i<node->pointIndices.size();i+=3)
		{
			for(size_t j=0;j<other_node->pointIndices.size();j+=3)
			{
				const Vector3f& A0 = node_transform.apply(node->points->at( node->pointIndices[i + 0] ));
				const Vector3f& A1 = node_transform.apply(node->points->at( node->pointIndices[i + 1] ));
				const Vector3f& A2 = node_transform.apply(node->points->at( node->pointIndices[i + 2] ));
				const Vector3f& B0 = other_transform.apply(other_node->points->at( other_node->pointIndices[j + 0] ));
				const Vector3f& B1 = other_transform.apply(other_node->points->at( other_node->pointIndices[j + 1] ));
				const Vector3f& B2 = other_transform.apply(other_node->points->at( other_node->pointIndices[j + 2] ));
				if (tri_tri_intersect(
					A0.data(), A1.data(), A2.data(),
					B0.data(), B1.data(), B2.data()))
				{
					result.hit = true;
					result.triangle_index1 = node->pointIndices[i] / 3;
					result.triangle_index2 = other_node->pointIndices[j] / 3;
					if (!result.find_all)
						return true;
					result.pairs.emplace_back(result.triangle_index1, result.triangle_index2);
				}
			}
		}
		return false;
	}

	bool OBBTree::check_intersection_recursive(
		const OBBNodePtr& node,
		const Transform& node_transform,
		const OBBNodePtr& other_node,
		const Transform& other_transform,
		IntersectionResult& result)
	{
		++result.checks;
		if (!intersects(node->box, node_transform, other_node->box, other_transform))
			return false;
		if (node->leaf)
		{
			if (!other_node->leaf)
			{
				return
					check_intersection_recursive(node, node_transform, other_node->left, other_transform, result) ||
					check_intersection_recursive(node, node_transform, other_node->right, other_transform, result);
			}
			return check_intersection_leaves(node, node_transform, other_node, other_transform, result);
		}
		if (other_node->leaf)
		{
			return
				check_intersection_recursive(node->left, node_transform, other_node, other_transform, result) ||
				check_intersection_recursive(node->right, node_transform, other_node, other_transform, result);
		}
		if (node->box.half_size.maxCoeff() > other_node->box.half_size.maxCoeff())
		{
			return
				check_intersection_recursive(node->left, node_transform, other_node, other_transform, result) ||
				check_intersection_recursive(node->right, node_transform, other_node, other_transform, result);
		}
		else
		{
			return
				check_intersection_recursive(node, node_transform, other_node->left, other_transform, result) ||
				check_intersection_recursive(node, node_transform, other_node->right, other_transform, result);
		}
	}

	IntersectionResult OBBTree::check_intersection(const OBBTree& other) const
	{
		IntersectionResult result;
		check_intersection_recursive(m_Root, m_DynamicTransform, other.m_Root, other.m_DynamicTransform, result);
		return result;
	}

	IntersectionResult OBBTree::find_all_intersections(const OBBTree& other) const
	{
		IntersectionResult result;
		result.find_all = true;
		check_intersection_recursive(m_Root, m_DynamicTransform, other.m_Root, other.m_DynamicTransform, result);
		return result;
	}

	bool OBBTree::verify() const
	{
		return m_Root->verify();
	}
}
