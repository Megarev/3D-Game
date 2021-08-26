#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include "math.h"

constexpr float PI = 3.1415926f;
constexpr float Sign(float value) { return value == 0 ? 0.0f : (value > 0.0f ? 1.0f : -1.0f); }
static int Random(int a, int b) {

	if (a > b) std::swap(a, b);

	std::random_device rd;
	static std::mt19937 m(rd());
	std::uniform_int_distribution<> dist(a, b);

	return dist(m);
}

static float Random(float a, float b) {

	if (a > b) std::swap(a, b);

	std::random_device rd;
	static std::mt19937 m(rd());
	std::uniform_real_distribution<> dist(a, b);

	return dist(m);
}

constexpr float Clamp(float value, float a, float b) {
	if (a > b) std::swap(a, b);
	return value <= a ? a : (value >= b ? b : value);
}


enum class ShapeType {
	CUBE,
	SPHERE
};

struct Mesh {
	std::vector<Triangle> t;
	vf3d pos, velocity;
	vf3d size = vf3d::ONE; // Cubes and other shapes
	float radius = 0.5f; // Spheres

	olc::Pixel init_color, color;
	bool is_render = true;
	float mass = 1.0f, inv_mass = 1.0f;
	ShapeType type = ShapeType::CUBE;

	Mesh() {}
	Mesh(const std::string& filename) { LoadObj(filename); }

	bool LoadObj(const std::string& filename, bool is_default_texels = true) {
		t.clear();

		std::ifstream reader(filename);
		if (!reader.is_open()) return false;

		std::vector<vf3d> vertices;

		std::string line;
		while (std::getline(reader, line)) {
			std::istringstream iss(line);
			if (line[0] == 'v') {
				vf3d v;
				char type;
				iss >> type >> v.x >> v.y >> v.z;
				vertices.push_back(v);
			}


			if (line[0] == 'f') {
				int f[3]{ 0 };
				char type;
				iss >> type >> f[0] >> f[1] >> f[2];
				t.push_back(Triangle{ vertices[f[0] - 1], vertices[f[1] - 1], vertices[f[2] - 1] });
			}
		}

		if (is_default_texels) {
			const Tex2D tex1[3] = { 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f };
			const Tex2D tex2[3] = { 0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f };
			//const Tex2D tex2[3] = { 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f };
			//const Tex2D tex2[3] = { 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f };

			int n = 0;
			for (auto& tri : t) {
				const Tex2D* texels = (n % 2) ? tex1 : tex2;
				
				tri.tex[0] = texels[0];
				tri.tex[1] = texels[1];
				tri.tex[2] = texels[2];

				n++;
			}
		}

		return true;
	}

	void Update(float dt) {
		pos += velocity * dt;
	}

	void Scale(float s) {
		size *= s;
		radius *= s;

		for (auto& tri : t) {
			for (auto& v : tri.v) v *= s;
		}
	}

	void Scale(float x, float y, float z) {
		size *= vf3d{ x, y, z };
		radius = 0.0f;
		
		for (auto& tri : t) {
			for (auto& v : tri.v) v *= size;
		}
	}
};

class Engine3D {
public:
	Matrix4x4 projection, world, cam, view;
	float render_distance = 100.0f, z_translation = 0.0f, yaw = 0.0f, xaw = 0.0f;
	float angle = 0.0f;
	vf3d look_dir, cam_pos, light_dir;
	
	int32_t width = 0, height = 0; // Screen dimensions

	// Cache
	olc::vf2d prev_m_pos; // Previous mouse position

	float* depth = nullptr; // Depth buffer
private:
	void RasterizeMesh(const Mesh& mesh, std::vector<Triangle>& raster) {
	
		if (!mesh.is_render) return;

		auto GetShade = [](float p, int n, const olc::Pixel& ref = olc::WHITE) -> olc::Pixel {
			int luminosity = p * n;
			uint8_t lum = 25;

			if (luminosity <= 0) return olc::Pixel{ lum, lum, lum };
			return olc::Pixel{ (uint8_t)(ref.r / n * luminosity + lum), (uint8_t)(ref.g / n * luminosity + lum), (uint8_t)(ref.b / n * luminosity + lum) };
		};
		
		for (const auto& v : mesh.t) {
			Triangle translated = world.MultiplyTriangle(v) + mesh.pos;
			translated.tex[0] = v.tex[0];
			translated.tex[1] = v.tex[1];
			translated.tex[2] = v.tex[2];

			const vf3d& line_a = translated.v[1] - translated.v[0];
			const vf3d& line_b = translated.v[2] - translated.v[0];
			const vf3d& normal = line_a.cross(line_b).norm();

			if ((translated.v[0] - cam_pos).mag2() > render_distance * render_distance) continue;

			if (normal.dot(translated.v[0] - cam_pos) < 0.0f) {

				float p = (normal.dot(light_dir.norm()));
				const olc::Pixel& shade = GetShade(p, 10, mesh.color);

				// Transform world space to cam space
				Triangle view_triangle = view.MultiplyTriangle(translated);
				view_triangle.shade = shade;
				view_triangle.t = p;
				view_triangle.tex[0] = translated.tex[0];
				view_triangle.tex[1] = translated.tex[1];
				view_triangle.tex[2] = translated.tex[2];

				// Clipping
				Plane forward_plane{ { 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f } };
				Triangle clipped[2];
				int n_clip_triangles = forward_plane.TrianglePlaneClip(view_triangle, clipped);

				for (int n = 0; n < n_clip_triangles; n++) {
					// Project triangle onto screen (3D -> 2D)
					Triangle proj = projection.MultiplyTriangle(clipped[n]);
					proj.shade = clipped[n].shade;
					proj.t = clipped[n].t;
					proj.tex[0] = clipped[n].tex[0];
					proj.tex[1] = clipped[n].tex[1];
					proj.tex[2] = clipped[n].tex[2];

					proj.v[0] /= proj.v[0].w;
					proj.v[1] /= proj.v[1].w;
					proj.v[2] /= proj.v[2].w;

					proj.tex[0] = proj.tex[0] / proj.v[0].w;
					proj.tex[1] = proj.tex[1] / proj.v[1].w;
					proj.tex[2] = proj.tex[2] / proj.v[2].w;
					proj.tex[0].w = 1.0f / proj.v[0].w;
					proj.tex[1].w = 1.0f / proj.v[1].w;
					proj.tex[2].w = 1.0f / proj.v[2].w;

					// Invert x/y axis
					proj.v[0].x *= -1.0f; proj.v[1].x *= -1.0f; proj.v[2].x *= -1.0f;
					proj.v[0].y *= -1.0f; proj.v[1].y *= -1.0f; proj.v[2].y *= -1.0f;

					// Translate to center
					proj.v[0].x += 1.0f; proj.v[1].x += 1.0f; proj.v[2].x += 1.0f;
					proj.v[0].y += 1.0f; proj.v[1].y += 1.0f; proj.v[2].y += 1.0f;

					// Scale to pixel space
					proj.v[0].x *= 0.5f * width; proj.v[1].x *= 0.5f * width; proj.v[2].x *= 0.5f * width;
					proj.v[0].y *= 0.5f * height; proj.v[1].y *= 0.5f * height; proj.v[2].y *= 0.5f * height;

					raster.push_back(proj);
				}
			}
		}
	}
	
	void UpdateWorldMatrix(float z = 16.0f) {

		z_translation = z;

		// XZ rotation
		world = Matrix4x4::MatrixMultiply(Matrix4x4::GetRotationMatrixZ(0.5f * angle), Matrix4x4::GetRotationMatrixX(angle));
		world = Matrix4x4::MatrixMultiply(world, Matrix4x4::MakeTranslation({ 0.0f, 0.0f, z }));
	}

	void UpdateViewMatrix(const vf3d& forward = { 0.0f, 0.0f, 1.0f }, const vf3d & up = { 0.0f, 1.0f, 0.0f }) {

		Matrix4x4 rotation = Matrix4x4::MatrixMultiply(Matrix4x4::GetRotationMatrixX(xaw), Matrix4x4::GetRotationMatrixY(yaw));

		look_dir = rotation.Multiply(forward);

		cam = Matrix4x4::PointAtMatrix(cam_pos, cam_pos + look_dir, up);
		view = Matrix4x4::LookAtMatrix(cam);
	}
	
	void UpdateMatrix() {
		UpdateWorldMatrix(4.0f);
		UpdateViewMatrix();
	}

	void SetProjectionMatrix(int32_t w, int32_t h, float FOV, float z_near, float z_far) {
		width = w;
		height = h;
		projection.SetProjectionMatrix(w, h, FOV, z_near, z_far);
	}
private:
	void DrawTexture(olc::PixelGameEngine* pge, std::vector<Triangle>& raster, olc::Sprite* sprite, bool is_wire_frame = false) {
		for (auto& t : raster) {
			Triangle clipped[2];
			std::list<Triangle> clipped_triangles = { t };
			int n_new_triangles = 1;

			for (int n = 0; n < 4; n++) {
				while (n_new_triangles > 0) {
					Triangle test = clipped_triangles.front();
					clipped_triangles.pop_front();
					n_new_triangles--;

					Plane plane;
					switch (n) {
					case 0: // Top plane
						plane = Plane{ { 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f } };
						break;
					case 1: // Bottom plane
						plane = Plane{ { 0.0f, height - 1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f } };
						break;
					case 2: // Left plane
						plane = Plane{ { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f } };
						break;
					case 3:
						plane = Plane{ { width - 1.0f, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f } };
						break;
					}

					int n_clip_triangles = plane.TrianglePlaneClip(test, clipped);

					for (int i = 0; i < n_clip_triangles; i++) {
						clipped_triangles.push_back(clipped[i]);
					}
				}
				n_new_triangles = clipped_triangles.size();
			}

			for (auto& t_clipped : clipped_triangles) {

				SpriteData s1 = { (int)t_clipped.v[0].x, (int)t_clipped.v[0].y, t_clipped.tex[0] };
				SpriteData s2 = { (int)t_clipped.v[1].x, (int)t_clipped.v[1].y, t_clipped.tex[1] };
				SpriteData s3 = { (int)t_clipped.v[2].x, (int)t_clipped.v[2].y, t_clipped.tex[2] };

				DrawTexturedTriangle(pge, s1, s2, s3, sprite);

				if (is_wire_frame) {
					pge->DrawTriangle((int)t_clipped.v[0].x, (int)t_clipped.v[0].y,
						(int)t_clipped.v[1].x, (int)t_clipped.v[1].y,
						(int)t_clipped.v[2].x, (int)t_clipped.v[2].y,
						olc::WHITE);
				}
			}
		}
	}
private:
	struct SpriteData {
		float x, y;
		Tex2D tex;

		SpriteData operator+(const SpriteData& other) {
			return {
				x + other.x,
				y + other.y,
				tex.u + other.tex.u,
				tex.v + other.tex.v,
				tex.w + other.tex.w
			};
		}

		SpriteData operator-(const SpriteData& other) {
			return {
				x - other.x,
				y - other.y,
				tex.u - other.tex.u,
				tex.v - other.tex.v,
				tex.w - other.tex.w
			};
		}

		SpriteData operator*(float other) {
			return {
				x * other,
				y * other,
				tex.u * other,
				tex.v * other,
				tex.w * other
			};
		}

		friend SpriteData operator/(const SpriteData& s, float other) {
			return {
				s.x / other,
				s.y / other,
				s.tex.u / other,
				s.tex.v / other,
				s.tex.w / other
			};
		}

		friend SpriteData operator*(const SpriteData& s, float other) {
			return s * other;
		}
	};
public:
	Engine3D() {
		light_dir = { 0.0f, 0.0f, -1.0f };
	}
	~Engine3D() {
		delete[] depth;
	}

	void Initialize(int32_t w, int32_t h, float FOV, float z_near, float z_far) {
		SetProjectionMatrix(w, h, FOV, z_near, z_far);
		
		depth = new float[w * h]{ 0.0f };
	}

	vf3d GetWorldPoint(const vf3d& point) const { return world.Multiply(point); }
	vf3d GetCamPoint(const vf3d& point) const { return cam.Multiply(point); }
	vf3d GetViewPoint(const vf3d& point) const { return view.Multiply(point); }
	vf3d GetScreenPoint(const vf3d& point) const { return point - vf3d{ 0.0f, 0.0f, z_translation }; }
	vf3d GetNearPoint(const vf3d& point) const { return point + vf3d{ 0.0f, 0.0f, z_translation }; }
	vf3d ToWorld(const olc::vf2d& point, float z_offset = 10.0f) const {
		const vf3d& pos = vf3d{ -2.0f * (point.x / width - 0.5f), -2.0f * (point.y / height - 0.5f), 1.0f, 0.0f };
		const vf3d& center = pos / vf3d{ projection.m[0][0], 1.0f, 1.0f, 1.0f };
		const vf3d& cam_point = GetScreenPoint(cam_pos) + z_offset * GetCamPoint(center);

		return cam_point;
	}

	void Input(olc::PixelGameEngine* pge, float dt, float move_speed = 8.0f, float rotation_sensitivity = 2.0f) {
		// Input
		const olc::vf2d& m_pos = olc::vf2d(pge->GetMousePos());

		float speed = move_speed;

		const vf3d& up = { 0.0f, 1.0f, 0.0f };
		if (pge->GetKey(olc::SPACE).bHeld) cam_pos += up * speed * dt;
		else if (pge->GetKey(olc::SHIFT).bHeld) cam_pos -= up * speed * dt;

		const vf3d& right = (up - look_dir * look_dir.dot(up)).norm().cross(look_dir);
		if (pge->GetKey(olc::A).bHeld) cam_pos += right * speed * dt;
		else if (pge->GetKey(olc::D).bHeld) cam_pos -= right * speed * dt;

		const vf3d& forward = look_dir * speed;
		if (pge->GetKey(olc::W).bHeld) cam_pos += forward * dt;
		else if (pge->GetKey(olc::S).bHeld) cam_pos -= forward * dt;
		
		if (pge->GetMouse(0).bHeld) {
			float rotation_x = Sign(m_pos.x - prev_m_pos.x);
			yaw += rotation_sensitivity * rotation_x * dt;

			float rotation_y = Sign(m_pos.y - prev_m_pos.y);
			// Constrain xaw to prevent jumpy camera rotation
			xaw = Clamp(xaw - rotation_sensitivity * rotation_y * dt, -PI / 2.0f + 0.2f, PI / 2.0f - 0.2f);
		}
		prev_m_pos = m_pos;
	}
	
	void Update(float dt = 0.0f) {
		UpdateMatrix();
		memset(depth, 0.0f, width * height * sizeof(float));
	}

	void RenderTriangles(olc::PixelGameEngine* pge, std::vector<Triangle>& raster, bool is_wire_frame = false, const olc::Pixel& wire_frame_color = olc::WHITE) {

		for (auto& t : raster) {
			Triangle clipped[2];
			std::list<Triangle> clipped_triangles = { t };
			int n_new_triangles = 1;

			for (int n = 0; n < 4; n++) {
				while (n_new_triangles > 0) {
					Triangle test = clipped_triangles.front();
					clipped_triangles.pop_front();
					n_new_triangles--;

					Plane plane;
					switch (n) {
					case 0: // Top plane
						plane = Plane{ { 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f } };
						break;
					case 1: // Bottom plane
						plane = Plane{ { 0.0f, height - 1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f } };
						break;
					case 2: // Left plane
						plane = Plane{ { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f } };
						break;
					case 3: // Right plane
						plane = Plane{ { width - 1.0f, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f } };
						break;
					}

					int n_clip_triangles = plane.TrianglePlaneClip(test, clipped);

					for (int i = 0; i < n_clip_triangles; i++) {
						clipped_triangles.push_back(clipped[i]);
					}
				}
				n_new_triangles = clipped_triangles.size();
			}

			for (auto& t_clipped : clipped_triangles) {

				const olc::vi2d& p1 = { (int)t_clipped.v[0].x, (int)t_clipped.v[0].y };
				float w1 = t_clipped.tex[0].w;
				const olc::vi2d& p2 = { (int)t_clipped.v[1].x, (int)t_clipped.v[1].y };
				float w2 = t_clipped.tex[1].w;
				const olc::vi2d& p3 = { (int)t_clipped.v[2].x, (int)t_clipped.v[2].y };
				float w3 = t_clipped.tex[2].w;


				DrawColoredTriangle(pge, p1, w1, p2, w2, p3, w3, t_clipped.shade);
				if (is_wire_frame) {
					pge->DrawTriangle((int)t_clipped.v[0].x, (int)t_clipped.v[0].y,
						(int)t_clipped.v[1].x, (int)t_clipped.v[1].y,
						(int)t_clipped.v[2].x, (int)t_clipped.v[2].y,
						wire_frame_color);
				}
			}
		}
	}

	void DrawTexturedTriangle(olc::PixelGameEngine* pge, SpriteData s1, SpriteData s2, SpriteData s3, olc::Sprite* sprite = nullptr) {
		// y-sort the points
		if (s2.y < s1.y) std::swap(s1, s2);
		if (s3.y < s1.y) std::swap(s1, s3);
		if (s3.y < s2.y) std::swap(s2, s3);

		// Offsets between position and texture UV coordinates
		SpriteData offset1 = s2 - s1;
		SpriteData offset2 = s3 - s1;

		// Steps along the lines
		SpriteData move_step1, move_step2;
		
		// Calculate gradients for upper triangle
		if (offset1.y) move_step1 = offset1 / std::fabsf(offset1.y);
		if (offset2.y) move_step2 = offset2 / std::fabsf(offset2.y);

		// Draw upper triangle
		if (offset1.y) {
			for (int i = s1.y; i <= s2.y; i++) {
				SpriteData a = s1 + move_step1 * (i - s1.y);
				SpriteData b = s1 + move_step2 * (i - s1.y);
			
				if (a.x > b.x) std::swap(a, b);

				float step = 1.0f / (float)((int)b.x - (int)a.x);
				float t = 0.0f;

				for (int j = a.x; j < b.x; j++) {
					Tex2D texel = Tex2D::Lerp(a.tex, b.tex, t);
					if (texel.w > depth[i * width + j]) {
						pge->Draw(j, i, sprite->Sample(texel.u / texel.w, texel.v / texel.w));
						depth[i * width + j] = texel.w;
					}
					t += step;
				}
			}
		}

		// Calculate gradients for lower triangle
		offset1 = s3 - s2;
		move_step1 = SpriteData{}; // Reset move_step1
		if (offset1.y) move_step1 = offset1 / std::fabsf(offset1.y);

		// Draw lower triangle
		if (offset1.y) {
			for (int i = s2.y; i <= s3.y; i++) {

				SpriteData a = s2 + move_step1 * (i - s2.y);
				SpriteData b = s1 + move_step2 * (i - s1.y);

				if (a.x > b.x) std::swap(a, b);

				float step = 1.0f / (float)((int)b.x - (int)a.x);
				float t = 0.0f;

				for (int j = a.x; j < b.x; j++) {
					Tex2D texel = Tex2D::Lerp(a.tex, b.tex, t);
					if (texel.w > depth[i * width + j]) {
						pge->Draw(j, i, sprite->Sample(texel.u / texel.w, texel.v / texel.w));
						depth[i * width + j] = texel.w;
					}
					t += step;
				}
			}
		}
	}

	void DrawColoredTriangle(olc::PixelGameEngine* pge, olc::vi2d p1, float w1, olc::vi2d p2, float w2, olc::vi2d p3, float w3, const olc::Pixel& color = olc::WHITE) {
		// y-sort the points
		if (p2.y < p1.y) { std::swap(p1, p2); std::swap(w1, w2); }
		if (p3.y < p1.y) { std::swap(p1, p3); std::swap(w1, w3); }
		if (p3.y < p2.y) { std::swap(p2, p3); std::swap(w2, w3); }

		// Offsets between position
		olc::vf2d offset1 = olc::vf2d(p2 - p1);
		olc::vf2d offset2 = olc::vf2d(p3 - p1);

		// Steps along the lines
		olc::vf2d move_step1, move_step2;
		float move_step_w1 = 0.0f, move_step_w2 = 0.0f;

		// Calculate gradients for upper triangle
		if (offset1.y) {
			move_step1 = offset1 / std::fabsf(offset1.y);
			move_step_w1 = (w2 - w1) / std::fabsf(offset1.y);
		}
		if (offset2.y) {
			move_step2 = offset2 / std::fabsf(offset2.y);
			move_step_w1 = (w3 - w1) / std::fabsf(offset2.y);
		}

		// Draw upper triangle
		if (offset1.y) {
			for (int i = p1.y; i <= p2.y; i++) {
				olc::vi2d a = p1 + move_step1 * (i - p1.y);
				float aw = w1 + (i - p1.y) * move_step_w1;
				olc::vi2d b = p1 + move_step2 * (i - p1.y);
				float bw = w1 + (i - p1.y) * move_step_w2;

				if (a.x > b.x) { std::swap(a, b); std::swap(aw, bw); }

				float step = 1.0f / float((int)b.x - (int)a.x);
				float t = 0.0f;

				for (int j = a.x; j < b.x; j++) {
					float w = aw + t * (bw - aw);
					if (w > depth[i * width + j]) {
						pge->Draw(j, i, color);
						depth[i * width + j] = w;
					}
					t += step;
				}
			}
		}

		// Calculate gradients for lower triangle
		offset1 = olc::vf2d(p3 - p2);

		// Reset move step
		move_step1 = olc::vf2d{};
		move_step_w1 = 0.0f;

		if (offset1.y) {
			move_step1 = offset1 / std::fabsf(offset1.y);
			move_step_w1 = (w3 - w2) / std::fabsf(offset1.y);
		}

		// Draw lower triangle
		if (offset1.y) {
			for (int i = p2.y; i <= p3.y; i++) {
				olc::vi2d a = p2 + move_step1 * (i - p2.y);
				float aw = w1 + (i - p2.y) * move_step_w1;
				olc::vi2d b = p1 + move_step2 * (i - p1.y);
				float bw = w1 + (i - p1.y) * move_step_w2;

				if (a.x > b.x) { std::swap(a, b); std::swap(aw, bw); }

				float step = 1.0f / float((int)b.x - (int)a.x);
				float t = 0.0f;

				for (int j = a.x; j < b.x; j++) {
					float w = aw + t * (bw - aw);
					if (w > depth[i * width + j]) {
						pge->Draw(j, i, color);
						depth[i * width + j] = w;
					}
					t += step;
				}
			}
		}
	}

	void DrawSprite(olc::PixelGameEngine* pge, const Mesh& mesh, olc::Sprite* sprite = nullptr, bool is_wire_frame = false) {
		std::vector<Triangle> tris;
		RasterizeMesh(mesh, tris);
		DrawTexture(pge, tris, sprite, is_wire_frame);
	}

	void DrawCube(olc::PixelGameEngine* pge, const vf3d& position, const olc::Pixel& color, bool is_wire_frame = false, const olc::Pixel& wire_frame_color = olc::WHITE) {
		Mesh cube;

		cube.t = {
			{ 0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f },

			{ 1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f },

			{ 1.0f, 0.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f },

			{ 0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f },

			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f },

			{ 0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 1.0f },
			{ 0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f }
		};

		cube.pos = position;
		cube.init_color = color;
		cube.color = color;
	
		std::vector<Triangle> tri_cube;
		RasterizeMesh(cube, tri_cube);
		RenderTriangles(pge, tri_cube, is_wire_frame, wire_frame_color);
	}
};