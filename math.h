#pragma once

struct vf3d {
	float x = 0.0f, y = 0.0f, z = 0.0f, w = 1.0f;
	static const vf3d ZERO, ONE;

	float dot(const vf3d& other) const { return x * other.x + y * other.y + z * other.z; }
	
	float mag2() const { return this->dot(*this); }
	float mag() const { return std::sqrtf(mag2()); }
	float min() const { return x < y && x < z ? x : (y < x && y < z ? y : z); }
	float max() const { return x > y && x > z ? x : (y > x && y > z ? y : z); }
	vf3d abs() const { return { std::fabsf(x), std::fabsf(y), std::fabsf(z) }; }
	
	vf3d cross(const vf3d& other) const { return { y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x }; }
	vf3d norm() const { return (*this) / mag(); }
	static vf3d Lerp(const vf3d& v1, const vf3d& v2, float t) { return v1 + (v2 - v1) * t; }

	vf3d operator+(const vf3d& other) const { return { x + other.x, y + other.y, z + other.z, w }; }
	void operator+=(const vf3d& other) { (*this) = (*this) + other; }

	vf3d operator-(const vf3d& other) const { return { x - other.x, y - other.y, z - other.z, w }; }
	void operator-=(const vf3d& other) { (*this) = (*this) - other; }

	vf3d operator-() const { return -1.0f * (*this); }

	vf3d operator*(const vf3d& other) const { return { x * other.x, y * other.y, z * other.z, w }; }
	vf3d operator*(float other) const { return { x * other, y * other, z * other, w }; }
	friend vf3d operator*(float other, const vf3d& v) { return v * other; }
	void operator*=(const vf3d& other) { (*this) = (*this) * other; }
	void operator*=(float other) { (*this) = (*this) * other; }

	vf3d operator/(const vf3d& other) const { return { x / other.x, y / other.y, z / other.z, w }; }
	vf3d operator/(float other) const { return { x / other, y / other, z / other, w }; }
	void operator/=(const vf3d& other) { (*this) = (*this) / other; }
	void operator/=(float other) { (*this) = (*this) / other; }

	bool operator==(const vf3d& other) { return (x == other.x && y == other.y && z == other.z); }
	bool operator!=(const vf3d& other) { return (x != other.x && y != other.y && z != other.z); }
};

vf3d const vf3d::ZERO = { 0.0f, 0.0f, 0.0f, 1.0f };
vf3d const vf3d::ONE  = { 1.0f, 1.0f, 1.0f, 1.0f };

struct Tex2D {
	float u = 0.0f, v = 0.0f;
	float w = 1.0f; // Depth

	static Tex2D Lerp(const Tex2D& a, const Tex2D& b, float t) {
		return {
			(1.0f - t) * a.u + t * b.u,
			(1.0f - t) * a.v + t * b.v,
			(1.0f - t) * a.w + t * b.w,
		};
	}

	friend Tex2D operator/(const Tex2D& t, float other) {
		return {
			t.u / other,
			t.v / other,
			t.w
		};
	}
};

struct Triangle { 
	vf3d v[3];
	Tex2D tex[3];
	olc::Pixel shade;
	float t = 0.0f;

	Triangle operator+(const Triangle& other) { return { v[0] + other.v[0], v[1] + other.v[1], v[2] + other.v[2] }; }
	Triangle operator+(const vf3d& other) { return { v[0] + other, v[1] + other, v[2] + other }; }
	void operator+=(const Triangle& other) { (*this) = (*this) + other; }
	void operator+=(const vf3d& other) { (*this) = (*this) + other; }

	Triangle operator-(const Triangle& other) { return { v[0] - other.v[0], v[1] - other.v[1], v[2] - other.v[2] }; }
	Triangle operator-(const vf3d& other) { return { v[0] - other, v[1] - other, v[2] - other }; }
	void operator-=(const Triangle& other) { (*this) = (*this) - other; }
	void operator-=(const vf3d& other) { (*this) = (*this) - other; }

	Triangle operator*(const Triangle& other) { return { v[0] * other.v[0], v[1] * other.v[1], v[2] * other.v[2] }; }
	Triangle operator*(const vf3d& other) { return { v[0] * other, v[1] * other, v[2] * other }; }
	void operator*=(const Triangle& other) { (*this) = (*this) * other; }
	void operator*=(const vf3d& other) { (*this) = (*this) * other; }

	Triangle operator/(const Triangle& other) { return { v[0] / other.v[0], v[1] / other.v[1], v[2] / other.v[2] }; }
	Triangle operator/(const vf3d& other) { return { v[0] / other, v[1] / other, v[2] / other }; }
	void operator/=(const Triangle& other) { (*this) = (*this) / other; }
	void operator/=(const vf3d& other) { (*this) = (*this) / other; }
};

struct Plane {
	vf3d point, normal; // Point in the plane and its normalized normal

	vf3d LinePlaneIntersection(const vf3d& a, const vf3d& b, float& t) {
		t = (point - a).dot(normal) / (b - a).dot(normal);
		return vf3d::Lerp(a, b, t);
	}

	int TrianglePlaneClip(Triangle& input, Triangle* output) {
		auto Distance = [&](const vf3d& p) -> float {
			/* 
				Equation of a plane = Ax + By + Cz - D = 0
				=> (p-point).dot(normal) as D = point.dot(normal)
			*/
			return (p - point).dot(normal);
		};
	
		vf3d* inside_points[3]; int n_inside = 0;
		vf3d* outside_points[3]; int n_outside = 0;
		Tex2D* inside_tex[3]; int n_inside_tex = 0;
		Tex2D* outside_tex[3]; int n_outside_tex = 0;

		float d0 = Distance(input.v[0]);
		float d1 = Distance(input.v[1]);
		float d2 = Distance(input.v[2]);

		if (d0 >= 0.0f) { inside_points[n_inside++] = &input.v[0]; inside_tex[n_inside_tex++] = &input.tex[0]; }
		else { outside_points[n_outside++] = &input.v[0]; outside_tex[n_outside_tex++] = &input.tex[0]; }
		if (d1 >= 0.0f) { inside_points[n_inside++] = &input.v[1]; inside_tex[n_inside_tex++] = &input.tex[1]; }
		else { outside_points[n_outside++] = &input.v[1]; outside_tex[n_outside_tex++] = &input.tex[1]; }
		if (d2 >= 0.0f) { inside_points[n_inside++] = &input.v[2]; inside_tex[n_inside_tex++] = &input.tex[2]; }
		else { outside_points[n_outside++] = &input.v[2]; outside_tex[n_outside_tex++] = &input.tex[2]; }
	
		if (n_outside == 3) {
			// All points are in the negative plane
			return 0; 
		}
		if (n_inside == 3) {
			// All points are in the positive plane
			output[0] = input;
			return 1;
		}
		if (n_inside == 1 && n_outside == 2) {
			// One point in the positive plane and two points in the negative plane
			// A smaller sub-triangle forms

			//output[0].shade = olc::RED;
			output[0].shade = input.shade;
			output[0].t = input.t;

			float t = 0.0f;

			output[0].v[0] = *inside_points[0];
			output[0].tex[0] = *inside_tex[0];

			output[0].v[1] = LinePlaneIntersection(*inside_points[0], *outside_points[0], t);
			output[0].tex[1] = Tex2D::Lerp(*inside_tex[0], *outside_tex[0], t);

			output[0].v[2] = LinePlaneIntersection(*inside_points[0], *outside_points[1], t);
			output[0].tex[2] = Tex2D::Lerp(*inside_tex[0], *outside_tex[1], t);
		
			return 1;
		}
		if (n_inside == 2 && n_outside == 1) {
			// Two points in the positive plane and one point in the negative plane
			// A quad forms -> Two sub-triangles

			//output[0].shade = olc::GREEN;
			//output[1].shade = olc::BLUE;

			output[0].shade = input.shade;
			output[0].t = input.t;
			output[1].shade = input.shade;
			output[1].t = input.t;

			output[0].v[0] = *inside_points[0];
			output[0].tex[0] = *inside_tex[0];

			output[0].v[1] = *inside_points[1];
			output[0].tex[1] = *inside_tex[1];

			float t = 0.0f;

			output[0].v[2] = LinePlaneIntersection(*inside_points[0], *outside_points[0], t);
			output[0].tex[2] = Tex2D::Lerp(*inside_tex[0], *outside_tex[0], t);

			output[1].v[0] = *inside_points[1];
			output[1].tex[0] = *inside_tex[1];

			output[1].v[1] = output[0].v[2];
			output[1].tex[1] = output[0].tex[2];

			output[1].v[2] = LinePlaneIntersection(*inside_points[1], *outside_points[0], t);
			output[1].tex[2] = Tex2D::Lerp(*inside_tex[1], *outside_tex[0], t);

			return 2;
		}

		return 0;
	}
};

struct Matrix4x4 {
	float m[4][4]{ 0.0f };

	vf3d Multiply(const vf3d& v) const {
		float x = v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0] + v.w * m[3][0];
		float y = v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1] + v.w * m[3][1];
		float z = v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2] + v.w * m[3][2];
		float w = v.x * m[0][3] + v.y * m[1][3] + v.z * m[2][3] + v.w * m[3][3];

		vf3d m = { x, y, z, w };

		//if (w) m /= w;

		return m;
	}

	Triangle MultiplyTriangle(const Triangle& t) const {
		const vf3d& v0 = Multiply(t.v[0]);
		const vf3d& v1 = Multiply(t.v[1]);
		const vf3d& v2 = Multiply(t.v[2]);

		return { v0, v1, v2 };
	}

	void Identity() {
		memset(m, 0, 4 * 4 * sizeof(float));
		m[0][0] = 1.0f;
		m[1][1] = 1.0f;
		m[2][2] = 1.0f;
		m[3][3] = 1.0f;
	}

	void SetProjectionMatrix(int32_t w, int32_t h, float FOV, float z_near, float z_far) {
		memset(m, 0, 4 * 4 * sizeof(float));

		float aspect_ratio = (float)h / (float)w;

		m[0][0] = aspect_ratio;
		m[1][1] = 1.0f / tanf(FOV / 2.0f);
		m[2][2] = z_far / (z_far - z_near);
		m[3][2] = -z_near * (z_far / (z_far - z_near));
		m[2][3] = 1.0f;
	}

	static Matrix4x4 MakeTranslation(const vf3d& pos) {
		Matrix4x4 translation;
		translation.m[0][0] = 1.0f;
		translation.m[1][1] = 1.0f;
		translation.m[2][2] = 1.0f;
		translation.m[3][3] = 1.0f;
		translation.m[3][0] = pos.x;
		translation.m[3][1] = pos.y;
		translation.m[3][2] = pos.z;
		return translation;
	}

	static Matrix4x4 GetRotationMatrixX(float angle) {
		Matrix4x4 rotation_x;
		rotation_x.m[0][0] = 1.0f;
		rotation_x.m[1][1] = cosf(angle);
		rotation_x.m[1][2] = -sinf(angle);
		rotation_x.m[2][1] = sinf(angle);
		rotation_x.m[2][2] = cosf(angle);
		rotation_x.m[3][3] = 1.0f;

		return rotation_x;
	}

	static Matrix4x4 GetRotationMatrixY(float angle) {
		Matrix4x4 rotation_y;
		rotation_y.m[0][0] = cosf(angle);
		rotation_y.m[0][2] = sinf(angle);
		rotation_y.m[1][1] = 1.0f;
		rotation_y.m[2][0] = -sinf(angle);
		rotation_y.m[2][2] = cosf(angle);
		rotation_y.m[3][3] = 1.0f;

		return rotation_y;
	}

	static Matrix4x4 GetRotationMatrixZ(float angle) {
		Matrix4x4 rotation_z;
		rotation_z.m[0][0] = cosf(angle);
		rotation_z.m[0][1] = -sinf(angle);
		rotation_z.m[1][0] = sinf(angle);
		rotation_z.m[1][1] = cosf(angle);
		rotation_z.m[2][2] = 1.0f;
		rotation_z.m[3][3] = 1.0f;

		return rotation_z;
	}

	static Matrix4x4 PointAtMatrix(const vf3d& pos, const vf3d& target, const vf3d& up) {
		const vf3d& forward_dir = (target - pos).norm();
		const vf3d& up_dir = (up - forward_dir * forward_dir.dot(up)).norm();
		const vf3d& right_dir = up_dir.cross(forward_dir);

		Matrix4x4 matrix;
		matrix.m[0][0] = right_dir.x;	matrix.m[0][1] = right_dir.y;	matrix.m[0][2] = right_dir.z;	matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = up_dir.x;		matrix.m[1][1] = up_dir.y;		matrix.m[1][2] = up_dir.z;		matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = forward_dir.x; matrix.m[2][1] = forward_dir.y; matrix.m[2][2] = forward_dir.z; matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
		return matrix;
	}

	static Matrix4x4 LookAtMatrix(const Matrix4x4& point_at) {
		// Inverse of PointAtMatrix
		Matrix4x4 matrix;
		matrix.m[0][0] = point_at.m[0][0]; matrix.m[0][1] = point_at.m[1][0]; matrix.m[0][2] = point_at.m[2][0]; matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = point_at.m[0][1]; matrix.m[1][1] = point_at.m[1][1]; matrix.m[1][2] = point_at.m[2][1]; matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = point_at.m[0][2]; matrix.m[2][1] = point_at.m[1][2]; matrix.m[2][2] = point_at.m[2][2]; matrix.m[2][3] = 0.0f;

		matrix.m[3][0] = -(point_at.m[3][0] * matrix.m[0][0] + point_at.m[3][1] * matrix.m[1][0] + point_at.m[3][2] * matrix.m[2][0]); 
		matrix.m[3][1] = -(point_at.m[3][0] * matrix.m[0][1] + point_at.m[3][1] * matrix.m[1][1] + point_at.m[3][2] * matrix.m[2][1]);
		matrix.m[3][2] = -(point_at.m[3][0] * matrix.m[0][2] + point_at.m[3][1] * matrix.m[1][2] + point_at.m[3][2] * matrix.m[2][2]);
		matrix.m[3][3] = 1.0f;

		return matrix;
	}

	static Matrix4x4 MatrixMultiply(const Matrix4x4& m1, const Matrix4x4& m2) {
		Matrix4x4 p;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				p.m[j][i] = m1.m[j][0] * m2.m[0][i] + m1.m[j][1] * m2.m[1][i] + m1.m[j][2] * m2.m[2][i] + m1.m[j][3] * m2.m[3][i];
			}
		}
		return p;
	}
};

// Global functions
vf3d RotateX(const vf3d& v, float angle) {
	return Matrix4x4::GetRotationMatrixX(angle).Multiply(v);
}

Triangle RotateTriangleX(const Triangle& t, float angle) {
	return {
		RotateX(t.v[0], angle),
		RotateX(t.v[1], angle),
		RotateX(t.v[2], angle)
	};
}

vf3d RotateY(const vf3d& v, float angle) {
	return Matrix4x4::GetRotationMatrixY(angle).Multiply(v);
}

Triangle RotateTriangleY(const Triangle& t, float angle) {
	return {
		RotateY(t.v[0], angle),
		RotateY(t.v[1], angle),
		RotateY(t.v[2], angle)
	};
}

vf3d RotateZ(const vf3d& v, float angle) {
	return Matrix4x4::GetRotationMatrixZ(angle).Multiply(v);
}

Triangle RotateTriangleZ(const Triangle& t, float angle) {
	return {
		RotateZ(t.v[0], angle),
		RotateZ(t.v[1], angle),
		RotateZ(t.v[2], angle)
	};
}