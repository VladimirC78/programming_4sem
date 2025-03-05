#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;


class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0];
    double dy = x[1];
    values[0] = 500 * exp(-(dx*dx + dy*dy) * 1000);
    values[1] = 2000 * exp(-(dx*dx + dy*dy) * 1000);
  }
};


class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 20 * exp(-1*(x[0]*x[0] + x[1]*x[1]))*sin(x[0]*x[0] + x[1]*x[1])*sin(x[0]*x[0] + x[1]*x[1]);
  }
};


class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    auto x1 = x[0], y1 = x[1];
    return pow((x1 * x1 + y1 * y1 - 0.8), 3) - x1*x1*y1*y1*y1 < 0;
  }
};

int main()
{

  auto mesh = std::make_shared<Mesh>("meshes/domain.xml");
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(2.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);


  File file("snapshots/poisson.pvd");
  file << u;

  return 0;
}
