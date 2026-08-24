#include <math.h>
#include <vector>

static int ciMode = 0;

static dfloat rhoRef = 998.0;  //kg / m^3
static dfloat muRef = 1e-3;    //Pa-s
static dfloat lenScale = 0.01; //m
                               
static dfloat velScale = 1.0;
static dfloat timeScale = 1.0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("CUBATURE POLYNOMIAL DEGREE", std::string("10"));
  options.setArgs("CONSTANT FLOW RATE", "TRUE");
  options.setArgs("CONSTANT FLOW DIRECTION", "X");
  options.setArgs("RESTART FILE NAME", std::string("r.fld+time=0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "-1");
  options.setArgs("FLUID VISCOSITY", std::to_string(1.0 / 43500.0));
  options.setArgs("FLUID DENSITY", std::string("1.0"));
  options.setArgs("NUMBER TIMESTEPS", std::string("-1"));
  options.setArgs("END TIME", std::string("10"));
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("DT", std::string("4e-2"));

  options.setArgs("SUBCYCLING STEPS", std::string("1"));

  options.setArgs("FLUID PRESSURE SOLVER", "GMRES+FLEXIBLE");
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-4"));
  options.setArgs("FLUID PRESSURE PRECONDITIONER", "MULTIGRID");
  options.setArgs("FLUID PRESSURE MULTIGRID COARSE SOLVE", "TRUE");
  options.setArgs("FLUID PRESSURE MULTIGRID SMOOTHER", "FOURTHOPTCHEBYSHEV+ASM");
  options.setArgs("FLUID PRESSURE MULTIGRID CHEBYSHEV DEGREE", std::string("3"));
  options.setArgs("FLUID PRESSURE MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR", std::string("1.1"));
  options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");

  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-6"));
  options.setArgs("FLUID VELOCITY SOLVER", "CG+BLOCK");
  options.setArgs("FLUID VELOCITY INITIAL GUESS", "EXTRAPOLATION");

  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-8"));
  options.setArgs("SCALAR00 INITIAL GUESS", "EXTRAPOLATION");

  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-6"));
  options.setArgs("SCALAR01 INITIAL GUESS", "EXTRAPOLATION");

  if (ciMode == 2) { //dimensional run

    options.setArgs("FLUID DENSITY", to_string_f(rhoRef));
    options.setArgs("FLUID VISCOSITY", to_string_f(muRef));

    const dfloat Re = 43500;
    velScale = Re * muRef / rhoRef / lenScale;
    timeScale = lenScale / velScale;

    options.setArgs("FLOW RATE", to_string_f(velScale));

    dfloat dt;
    options.getArgs("DT", dt);
    options.setArgs("DT", to_string_f(dt * timeScale));

    dfloat endTime;
    options.getArgs("END TIME", endTime);
    options.setArgs("END TIME", to_string_f(endTime * timeScale));

    //adjust tolerances
    options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-3"));
    options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-5"));
  }
}

void ciDimensionalizeRestart()
{
  auto mesh = nrs->meshV;
  platform->linAlg->scale(mesh->Nlocal, lenScale, mesh->o_x);
  platform->linAlg->scale(mesh->Nlocal, lenScale, mesh->o_y);
  platform->linAlg->scale(mesh->Nlocal, lenScale, mesh->o_z);

  platform->linAlg->scale(mesh->dim * nrs->fluid->fieldOffset, velScale, nrs->fluid->o_U);
  platform->linAlg->scale(mesh->Nlocal, velScale * velScale, nrs->scalar->o_solution("k"));
  platform->linAlg->scale(mesh->Nlocal, timeScale, nrs->scalar->o_solution("tau"));
}

