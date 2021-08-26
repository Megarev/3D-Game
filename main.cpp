#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include "engine.h"

class Game : public olc::PixelGameEngine {
private:
	Mesh cube;
	Engine3D engine;

	olc::Sprite* box = nullptr;
public:
	Game() {
		sAppName = "Title";
	}

	bool OnUserCreate() override {

		engine.Initialize(ScreenWidth(), ScreenHeight(), PI / 2.0f, 0.5f, 1000.0f);

		box = new olc::Sprite("images/box.png");

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

		//cube.LoadObj("cube.txt");
		cube.color = olc::WHITE;
		cube.type = ShapeType::CUBE;

		engine.AddMesh("cube", cube);

		return true;
	}


	bool OnUserUpdate(float dt) override {

		// Input
		engine.Input(this, dt);

		// Logic
		engine.light_dir = -engine.look_dir;
		engine.Update(dt);

		// Render
		Clear(olc::BLACK);
		//engine.Render(this);
		engine.DrawTexture(this, box);

		return true;
	}
};

int main() {

	Game game;
	if (game.Construct(400, 400, 1, 1, false, true)) {
		game.Start();
	}

	return 0;
}