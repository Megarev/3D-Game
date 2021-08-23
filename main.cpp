#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include "engine.h"

class PhysicsGame : public olc::PixelGameEngine {
private:
	Engine3D engine;
	Mesh cube;
public:
	PhysicsGame() {
		sAppName = "Title";
	}

	bool OnUserCreate() override {

		cube.LoadObj("cube.txt");
		engine.SetProjectionMatrix(ScreenWidth(), ScreenHeight(), PI / 2.0f, 0.5f, 1000.0f);

		Mesh box = cube;
		box.color = olc::RED;
		box.type = ShapeType::CUBE;
		engine.AddMesh("box", box);

		return true;
	}

	bool OnUserUpdate(float dt) override {

		const olc::vf2d& m_pos = GetMousePos() * 1.0f;

		// Input
		engine.Input(this, dt, 8.0f, 2.0f);

		// Logic
		engine.light_dir = -engine.look_dir;
		engine.Update(dt);

		// Render
		Clear(olc::BLACK);
		engine.Render(this);

		return true;
	}
};

int main() {

	PhysicsGame game;
	if (game.Construct(400, 400, 1, 1, false, true)) {
		game.Start();
	}

	return 0;
}