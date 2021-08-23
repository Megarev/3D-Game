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

	bool LoadObj(const std::string& filename) {
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
	std::unordered_map<std::string, Mesh> mesh_map;
	std::vector<Triangle> triangles_raster;

	// Cache
	olc::vf2d prev_m_pos; // Previous mouse position
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
			Triangle translated = world.MultiplyTriangle(v);
			translated.v[0] += mesh.pos;
			translated.v[1] += mesh.pos;
			translated.v[2] += mesh.pos;

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

				// Clipping
				Plane forward_plane{ { 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f } };
				Triangle clipped[2];
				int n_clip_triangles = forward_plane.TrianglePlaneClip(view_triangle, clipped);

				for (int n = 0; n < n_clip_triangles; n++) {
					// Project triangle onto screen (3D -> 2D)
					Triangle proj = projection.MultiplyTriangle(clipped[n]);
					proj.shade = clipped[n].shade;
					proj.t = clipped[n].t;
					proj.v[0] /= proj.v[0].w;
					proj.v[1] /= proj.v[1].w;
					proj.v[2] /= proj.v[2].w;

					proj *= vf3d{ -1.0f, -1.0f, 1.0f }; // Invert x/y axis
					proj += vf3d{ 1.0f, 1.0f, 0.0f };
					proj *= vf3d{ 0.5f * width, 0.5f * height, 1.0f };

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
public:
	Engine3D() {
		light_dir = { 0.0f, 0.0f, -1.0f };
	}

	vf3d GetWorldPoint(const vf3d& point) const { return world.Multiply(point); }
	vf3d GetCamPoint(const vf3d& point) const { return cam.Multiply(point); }
	vf3d GetViewPoint(const vf3d& point) const { return view.Multiply(point); }
	vf3d Convert2DTo3D(const olc::vf2d& point, float z_offset = 10.0f) const {
		const vf3d& pos = vf3d{ -2.0f * (point.x / width - 0.5f), -2.0f * (point.y / height - 0.5f), 1.0f, 0.0f };
		const vf3d& center_to_mouse = pos / vf3d{ projection.m[0][0], 1.0f, 1.0f, 1.0f };
		const vf3d& cam_mouse = GetCamPoint(center_to_mouse);
		const vf3d& cam_point = GetScreenPoint(cam_pos) + z_offset * cam_mouse;

		return cam_point;
	}

	vf3d GetScreenPoint(const vf3d& point) const { return point - vf3d{ 0.0f, 0.0f, z_translation }; }
	vf3d GetNearPoint(const vf3d& point) const { return point + vf3d{ 0.0f, 0.0f, z_translation }; }

	void SetProjectionMatrix(int32_t w, int32_t h, float FOV, float z_near, float z_far) {
		width = w;
		height = h;
		projection.SetProjectionMatrix(w, h, FOV, z_near, z_far);
	}

	void Input(olc::PixelGameEngine* pge, float dt, float move_speed = 8.0f, float rotation_sensitivity = 1.2f) {
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

		//angle += 1.0f * dt;

		triangles_raster.clear();
		UpdateMatrix();
	
		for (auto& mesh : mesh_map) {
			mesh.second.Update(dt);
			RasterizeMesh(mesh.second, triangles_raster);
		}
	}

	void AddMesh(const std::string& identifier, const Mesh& mesh) {
		mesh_map.insert(std::make_pair(identifier, mesh));
	}

	bool AddMesh(const std::string& identifier, const std::string& filename, const olc::Pixel& color = olc::WHITE) {
		Mesh mesh;
		if (mesh.LoadObj(filename)) {
			mesh.color = color;
			mesh_map.insert(std::make_pair(identifier, mesh));
			return true;
		}

		return false;
	}

	bool RemoveMesh(const std::string& identifier) {
		if (!GetMesh(identifier)) return false;
		mesh_map.erase(identifier);

		return true;
	}

	Mesh* GetMesh(const std::string& identifier) {
		auto m = mesh_map.find(identifier);
		if (m != mesh_map.end()) return &m->second;

		return nullptr;
	}

	Mesh& operator[](const std::string& identifier) {
		return *GetMesh(identifier);
	}

	std::unordered_map<std::string, Mesh>& GetMeshMap() { return mesh_map; }

	size_t GetMeshCount() const { return mesh_map.size(); }

	void RenderTriangles(olc::PixelGameEngine* pge, std::vector<Triangle>& raster) {
		
		std::sort(raster.begin(), raster.end(), [](const Triangle& t1, const Triangle& t2) {
			float z1 = (t1.v[0].z + t1.v[1].z + t1.v[2].z) / 3.0f;
			float z2 = (t2.v[0].z + t2.v[1].z + t2.v[2].z) / 3.0f;

			return z1 > z2;
		});

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
				pge->FillTriangle((int)t_clipped.v[0].x, (int)t_clipped.v[0].y,
					(int)t_clipped.v[1].x, (int)t_clipped.v[1].y,
					(int)t_clipped.v[2].x, (int)t_clipped.v[2].y,
					t_clipped.shade);

				/*DrawTriangle((int)t_clipped.v[0].x, (int)t_clipped.v[0].y,
							 (int)t_clipped.v[1].x, (int)t_clipped.v[1].y,
							 (int)t_clipped.v[2].x, (int)t_clipped.v[2].y,
							 olc::BLACK);*/
			}
		}
	}

	void Render(olc::PixelGameEngine* pge) {
		RenderTriangles(pge, triangles_raster);
	}

	void DrawCube(olc::PixelGameEngine* pge, const vf3d& position, const olc::Pixel& color, bool is_fill = false) {
		Mesh cube;

		cube.t = {
			{ 0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f },
			{ 0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f },

			{ 1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f },
			{ 1.0f, 0.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 1.0f, 1.0f },

			{ 1.0f, 0.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f },
			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f },

			{ 0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f },
			{ 0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f },

			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f },
			{ 1.0f, 0.0f, 1.0f, 1.0f,  0.0f, 0.0f, 0.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f },

			{ 0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f },
			{ 0.0f, 1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f }
		};

		cube.pos = position;
		cube.init_color = color;
		cube.color = color;
	
		std::vector<Triangle> tri_cube;
		RasterizeMesh(cube, tri_cube);
		RenderTriangles(pge, tri_cube);
	}
};