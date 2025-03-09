#include "core.hpp"

void IcoNS::setBoundaryConditions()
{
    std::shared_ptr<BoundaryFunction> u_func;
    std::shared_ptr<BoundaryFunction> v_func, v_func1;
    std::shared_ptr<BoundaryFunction> w_func, w_func1;

    if (testCase == 1)
    {
        u_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0; });
        v_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0.0; });
        v_func1 = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                              { return 1.0; });
        w_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0; });
        for (int i = 0; i < 6 /*nfaces*/; i++)
        {
            boundary.addFunction(U, u_func);

            boundary.addFunction(W, v_func);
            if (i == RIGHT)
            {
                boundary.addFunction(V, v_func1);
            }
            else
            {
                boundary.addFunction(V, v_func);
            }
        }
    }
    else if (testCase == 2)
    {
        u_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0; });
        v_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0.0; });
        v_func1 = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                              { return 1.0; });
        w_func = std::make_shared<Dirichlet>([&](Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/)
                                             { return 0; });
        for (int i = 0; i < 6 /*nfaces*/; i++)
        {
            boundary.addFunction(U, u_func);

            boundary.addFunction(W, v_func);
            if (i == LEFT)
            {
                boundary.addFunction(V, v_func1);
            }
            else
            {
                boundary.addFunction(V, v_func);
            }
        }
    }
    else
    {
        u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                             { return std::sin(SX + (x + 0.5) * DX) * std::cos(SY + y * DY) * std::sin(SZ + z * DZ) * std::sin(t); });
        v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                             { return std::cos(SX + x * DX) * std::sin(SY + (y + 0.5) * DY) * std::sin(SZ + z * DZ) * std::sin(t); });
        w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                             { return 2 * (std::cos(SZ + x * DX) * std::cos(SY + y * DY) * std::cos(SZ + (z + 0.5) * DZ) * std::sin(t)); });
        for (int i = 0; i < 6 /*nfaces*/; i++)
        {
            boundary.addFunction(U, u_func);
            boundary.addFunction(V, v_func);
            boundary.addFunction(W, w_func);
        }
    }
}
