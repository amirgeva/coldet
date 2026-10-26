#pragma once

#include <memory>
#include "obb.h"

namespace obb
{

	using OBBNodePtr = std::shared_ptr<struct OBBNode>;
	using PointList = std::vector<Vector3f>;
	using PointListPtr = std::shared_ptr<PointList>;
	using PointIndex = uint32_t;
	using Point = Vector3f;

	struct Triangle
	{
		Point a, b, c;
	};

	struct OBBNode
	{
		OBBNode(PointListPtr points, size_t depth) 
			: points(points)
			, depth(depth)
		{}
		OBB	box;
		std::vector<PointIndex>	pointIndices;
		PointListPtr points;
		OBBNodePtr	left, right;
		size_t depth = 0;
		bool leaf = false;

		bool verify() const;
	};

	using TrianglePair = std::pair<size_t, size_t>;

	struct IntersectionResult
	{
		bool	hit = false;
		size_t  triangle_index1, triangle_index2;
		size_t	checks = 0;
		bool	find_all = false;			// keep searching after the first hit
		std::vector<TrianglePair> pairs;	// filled when find_all is set
	};

	class OBBTree
	{
		OBBNodePtr m_Root;
		PointListPtr m_Points;
		Transform m_DynamicTransform;

		static bool check_intersection_recursive(const OBBNodePtr& node,
												 const Transform& node_transform,
												 const OBBNodePtr& other_node,
												 const Transform& other_transform,
												 IntersectionResult& result);
		static bool check_intersection_leaves(const OBBNodePtr& node,
											  const Transform& node_transform,
											  const OBBNodePtr& other_node,
											  const Transform& other_transform,
											  IntersectionResult& result);
		static size_t split_node(OBBNodePtr node, 
								 size_t max_points_per_leaf);
	public:
		OBBTree();
		bool verify() const;
		void add_triangle(const Triangle& tri);
		size_t build(size_t max_triangles_per_leaf);
		void set_dynamic_transform(const Transform& t) { m_DynamicTransform = t; }

		IntersectionResult check_intersection(const OBBTree& other) const;
		IntersectionResult find_all_intersections(const OBBTree& other) const;
		IntersectionResult bf_check_intersection(const OBBTree& other) const;
	};

	void set_diagnostics(bool v);
}
