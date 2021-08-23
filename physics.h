#pragma once
#include "engine.h"
#include "math.h"

class Collisions {
public:
	static bool AABB(const Mesh& m1, const Mesh& m2, Engine3D& engine) {

		if (m1.mass + m2.mass == 0.0f) return false;

		const vf3d& pos1 = engine.GetWorldPoint(m1.pos);
		const vf3d& pos2 = engine.GetWorldPoint(m2.pos);

		return (
			pos1.x < pos2.x + m2.size.x && pos2.x < pos1.x + m1.size.x &&
			pos1.y < pos2.y + m2.size.y && pos2.y < pos1.y + m1.size.y &&
			pos1.z < pos2.z + m2.size.z && pos2.z < pos1.z + m1.size.z
		);
	}

	static bool CircleCircle(const Mesh& m1, const Mesh& m2, Engine3D& engine) {

		if (m1.mass + m2.mass == 0.0f) return false;

		const vf3d& pos1 = engine.GetWorldPoint(m1.pos);
		const vf3d& pos2 = engine.GetWorldPoint(m2.pos);

		return ((pos2 - pos1).mag2() <= (m1.radius + m2.radius) * (m1.radius + m2.radius));
	}

	static bool CircleRect(const Mesh& circle, const Mesh& rect, Engine3D& engine, float& depth) {

		if (circle.mass + rect.mass == 0.0f) return false;

		const vf3d& circle_pos = engine.GetWorldPoint(circle.pos);
		const vf3d& rect_pos = engine.GetWorldPoint(rect.pos);
		float r = circle.radius;

		vf3d direction = circle_pos - rect_pos;
		const vf3d& rect_extents = rect.size / 2.0f;
		vf3d clamped_direction;
		clamped_direction.x = Clamp(direction.x, -rect_extents.x, rect_extents.x);
		clamped_direction.y = Clamp(direction.y, -rect_extents.y, rect_extents.y);
		clamped_direction.z = Clamp(direction.z, -rect_extents.z, rect_extents.z);

		const vf3d& abs_direction = direction.abs();

		vf3d normal;
		if (abs_direction.min() == abs_direction.x) { normal.x = Sign(direction.x); }
		if (abs_direction.min() == abs_direction.y) { normal.y = Sign(direction.y); }
		if (abs_direction.min() == abs_direction.z) { normal.z = Sign(direction.z); }

		vf3d res = direction - clamped_direction;
		float d = res.mag2();
		if (d > r * r || d == 0.0f) return false;

		d = std::sqrt(d);
		depth = r - d;

		return true;
	}
};

class Resolutions {
public:
	static void AABB(Mesh& m1, Mesh& m2, Engine3D& engine) {
		const vf3d& pos1 = engine.GetWorldPoint(m1.pos);
		const vf3d& pos2 = engine.GetWorldPoint(m2.pos);

		const vf3d& dir = (pos2 + m2.size / 2.0f) - (pos1 + m1.size / 2.0f);
		const vf3d& total_extents = (m1.size + m2.size) / 2.0f;

		const vf3d& depth = total_extents - dir.abs();

		vf3d normal, min_depth;
		if (depth.x > 0.0f && depth.min() == depth.x) { normal.x = Sign(dir.x); min_depth.x = depth.x; }
		if (depth.y > 0.0f && depth.min() == depth.y) { normal.y = Sign(dir.y); min_depth.y = depth.y; }
		if (depth.z > 0.0f && depth.min() == depth.z) { normal.z = Sign(dir.z); min_depth.z = depth.z; }

		/*if (normal == vf3d::ZERO) {
			normal.x = Sign(Random(-1.0f, 1.0f));
			min_depth = total_extents;
		}*/

		float total_mass = m1.mass + m2.mass;

		m1.pos -= normal * min_depth * m1.mass / total_mass;
		m2.pos += normal * min_depth * m2.mass / total_mass;
	}

	static void CircleCircle(Mesh& m1, Mesh& m2, Engine3D& engine) {
		const vf3d& pos1 = engine.GetWorldPoint(m1.pos);
		const vf3d& pos2 = engine.GetWorldPoint(m2.pos);
		float r1 = m1.radius, r2 = m2.radius;

		vf3d normal = pos2 - pos1;
		float len = normal.mag();
		float depth = (r1 + r2) - len;
		if (len == 0.0f) return;
		normal /= len;

		float total_mass = m1.mass + m2.mass;

		m1.pos -= normal * depth * m1.mass / total_mass;
		m2.pos += normal * depth * m2.mass / total_mass;
	}

	static void CircleRect(Mesh& circle, Mesh& rect, Engine3D& engine, float depth) {
		vf3d normal = (circle.pos - rect.pos).norm();
		
		float total_mass = circle.mass + rect.mass;
		
		rect.pos -= depth * normal * rect.mass / total_mass;
		circle.pos += depth * normal * circle.mass / total_mass;
	}
};