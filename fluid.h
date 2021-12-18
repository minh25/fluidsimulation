#include <math.h>

class Fluid
{
public:
    Fluid(int size, float dt, float diffusion, float viscosity)
    {
        int N = size;

        this->size = N;
        this->dt = dt;
        this->diff = diffusion;
        this->visc = viscosity;

        this->s = new float[N * N];
        this->density = new float[N * N];

        this->Vx = new float[N * N];
        this->Vy = new float[N * N];

        this->Vx0 = new float[N * N];
        this->Vy0 = new float[N * N];
    }

    ~Fluid()
    {
        delete[] s;
        delete[] density;
        delete[] Vx;
        delete[] Vy;
        delete[] Vx0;
        delete[] Vy0;
    }

    void step()
    {
        auto N = this->size;
        auto visc = this->visc;
        auto diff = this->diff;
        auto dt = this->dt;
        auto Vx = this->Vx;
        auto Vy = this->Vy;
        auto Vx0 = this->Vx0;
        auto Vy0 = this->Vy0;
        auto s = this->s;
        auto density = this->density;

        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);

        project(Vx0, Vy0, Vx, Vy);

        advect(1, Vx, Vx0, Vx0, Vy0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, dt);

        project(Vx, Vy, Vx0, Vy0);
        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, dt);
    }

    float getDensity(int x, int y)
    {
        return this->density[IX(x, y)];
    }

    void addDensity(int x, int y, float amount)
    {
        this->density[IX(x, y)] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        auto index = IX(x, y);
        this->Vx[index] += amountX;
        this->Vy[index] += amountY;
    }

    // function to render density
    // void renderD()
    // {
    // 	for (int i = 0; i < N; i++)
    // 	{
    // 		for (int j = 0; j < N; j++)
    // 		{
    // 			auto x = i * SCALE;
    // 			auto y = j * SCALE;
    // 			auto d = this->density[IX(i, j)];
    // 			fill(d);
    // 			noStroke();
    // 			rect(x, y, SCALE, SCALE);
    // 		}
    // 	}
    // }

    // function to render velocity
    // void renderV()
    // {
    // 	for (auto i = 0; i < N; i++)
    // 	{
    // 		for (auto j = 0; j < N; j++)
    // 		{
    // 			auto x = i * SCALE;
    // 			auto y = j * SCALE;
    // 			auto vx = this->Vx[IX(i, j)];
    // 			auto vy = this->Vy[IX(i, j)];
    // 			this->canvas.stroke(0);
    //
    // 			if (!(abs(vx) < 0.1 && abs(vy) <= 0.1))
    // 			{
    // 				line(x, y, x + vx * SCALE, y + vy * SCALE);
    // 			}
    // 		}
    // 	}
    // }

    int IX(int x, int y)
    {
        return x + y * size;
    }

    void linear_solve(int b, float *x, float *x0, float a, float c)
    {
        int N = this->size;

        auto cRecip = 1.0 / c;
        for (auto t = 0; t < iter; t++)
        {
            for (auto j = 1; j < N - 1; j++)
            {
                for (auto i = 1; i < N - 1; i++)
                {
                    x[IX(i, j)] =
                        (x0[IX(i, j)] +
                         a *
                             (x[IX(i + 1, j)] +
                              x[IX(i - 1, j)] +
                              x[IX(i, j + 1)] +
                              x[IX(i, j - 1)])) *
                        cRecip;
                }
            }
            set_bnd(b, x);
        }
    }

    void diffuse(int b, float *x, float *x0, float diff, float dt)
    {
        int N = this->size;

        auto a = dt * diff * (N - 2) * (N - 2);
        linear_solve(b, x, x0, a, 1 + 6 * a);
    }

    /*
    Remember when I said that we're only simulating incompressible fluids? This means that the amount of fluid in each box has to stay constant. That means that the amount of fluid going in has to be exactly equal to the amount of fluid going out. The other operations tend to screw things up so that you get some boxes with a net outflow, and some with a net inflow. This operation runs through all the cells and fixes them up so everything is in equilibrium.
    */
    void project(float *velocX, float *velocY, float *p, float *div)
    {
        int N = this->size;

        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                div[IX(i, j)] =
                    (-0.5 *
                     (velocX[IX(i + 1, j)] -
                      velocX[IX(i - 1, j)] +
                      velocY[IX(i, j + 1)] -
                      velocY[IX(i, j - 1)])) /
                    N;
                p[IX(i, j)] = 0;
            }
        }

        set_bnd(0, div);
        set_bnd(0, p);
        linear_solve(0, p, div, 1, 6);

        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
                velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
            }
        }

        set_bnd(1, velocX);
        set_bnd(2, velocY);
    }

    /*
    Every cell has a set of velocities, and these velocities make things move. This is called advection. As with diffusion, advection applies both to the dye and to the velocities themselves.
    */
    void advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt)
    {
        int N = this->size;

        float i0, i1, j0, j1;

        auto dtx = dt * (N - 2);
        auto dty = dt * (N - 2);

        float s0, s1, t0, t1;
        float tmp1, tmp2, tmp3, x, y;

        auto Nfloat = N - 2;
        float ifloat, jfloat;
        int i, j, k;

        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++)
        {
            for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++)
            {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;

                if (x < 0.5)
                    x = 0.5;
                if (x > Nfloat + 0.5)
                    x = Nfloat + 0.5;
                i0 = std::floor(x);
                i1 = i0 + 1.0;
                if (y < 0.5)
                    y = 0.5;
                if (y > Nfloat + 0.5)
                    y = Nfloat + 0.5;
                j0 = std::floor(y);
                j1 = j0 + 1.0;

                s1 = x - i0;
                s0 = 1.0 - s1;
                t1 = y - j0;
                t0 = 1.0 - t1;

                int i0i = int(i0);
                int i1i = int(i1);
                int j0i = int(j0);
                int j1i = int(j1);

                d[IX(i, j)] =
                    s0 * (t0 * d0[IX(i0i, j0i)] +
                          t1 * d0[IX(i0i, j1i)]) +
                    s1 * (t0 * d0[IX(i1i, j0i)] +
                          t1 * d0[IX(i1i, j1i)]);
            }
        }

        set_bnd(b, d);
    }

    /*
    This is short for "set bounds", and it's a way to keep fluid from leaking out of your box. Not that it could really leak, since it's just a simulation in memory, but not having walls really screws up the simulation code. Walls are added by treating the outer layer of cells as the wall. Basically, every velocity in the layer next to this outer layer is mirrored. So when you have some velocity towards the wall in the next-to-outer layer, the wall gets a velocity that perfectly counters it.    
    b = 1 -> Ox
    b = 2 -> Oy
    */
    void set_bnd(int b, float *x)
    {
        int N = this->size;

        for (auto i = 1; i < N - 1; i++)
        {
            x[IX(i,     0)] = b == 2 ? -x[IX(i,     1)] : x[IX(i,     1)];
            x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
        }
        for (auto j = 1; j < N - 1; j++)
        {
            x[IX(    0, j)] = b == 1 ? -x[IX(    1, j)] : x[IX(    1, j)];
            x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
        }

        x[IX(    0,     0)] = 0.5 * (x[IX(    1,     0)] + x[IX(    0,     1)]);
        x[IX(    0, N - 1)] = 0.5 * (x[IX(    1, N - 1)] + x[IX(    0, N - 2)]);
        x[IX(N - 1,     0)] = 0.5 * (x[IX(N - 2,     0)] + x[IX(N - 1,     1)]);
        x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
    }

private:
    int size;
    int iter = 4;
    float dt, diff, visc;
    float *s, *density;
    float *Vx, *Vy;
    float *Vx0, *Vy0;
};