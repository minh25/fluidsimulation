#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include "fluid.h"

#define TWO_PI 6.2831855

int _SIZE = 256;
int _SCALE = 1;
float _DELTA_TIME = 0.2;
float _DIFFUSION = 0.0000001;
float _VISCOSITY = 0.0000001;

// Override base class with your custom functionality
class FluidSimulation : public olc::PixelGameEngine
{
public:
	FluidSimulation()
	{
		// Name your application
		sAppName = "Fluid simulation";
	}

private:
	Fluid *fluid;
	int pMouseX, pMouseY; // previous position of mouse
	int cx, cy;

public:
	bool OnUserCreate() override
	{
		fluid = new Fluid(_SIZE, _DELTA_TIME, _DIFFUSION, _VISCOSITY);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		cx = GetMouseX();
		cy = GetMouseY();

		if (GetMouse(0).bHeld)
		{
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					fluid->addDensity(cx + i, cy + j, 200);
				}
			}

			fluid->addVelocity(cx, cy, cx - pMouseX, cy - pMouseY);
		}

		if (GetMouse(1).bHeld)
		{
			fluid->addVelocity(cx, cy, cx - pMouseX, cy - pMouseY);
		}

		pMouseX = cx;
		pMouseY = cy;

		fluid->step();
		// Called once per frame, draws random coloured pixels
		for (int x = 0; x < ScreenWidth(); x++)
			for (int y = 0; y < ScreenHeight(); y++)
			{
				auto d = (int)(fluid->getDensity(x, y));
				if (d > 255)
					d = 255;
				Draw(x, y, olc::Pixel(d, d, d));
			}

		return true;
	}

	bool OnUserDestroy() override
	{
		delete fluid;
		return true;
	}
};

int main()
{
	FluidSimulation my_app;
	if (my_app.Construct(_SIZE, _SIZE, _SCALE, _SCALE))
		my_app.Start();
	return 0;
}
