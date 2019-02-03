#ifndef GET_STUFF_H
#define GET_STUFF_H

#include <iostream>
#include <ctime>
#include <sstream>
#include <utility>
#include <Eigen/Core>
//#include <Tudat/Basics/basicTypedefs.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace bislip {
namespace variables {

double getGamma ( const double A, const double B, const double phi );

double getEps_T ( const double &AoA, const Eigen::VectorXd &p, const Eigen::Vector6d &newCoefficients );

Eigen::Vector2d getGravs (
        const double mu,
        const double J2,
        const double J3,
        const double J4,
        const double R_E,
        const double r,
        const double delta_rad );

double getRoot_NewtRaph (
        const double &x_old,
        const double &f,
        const double &fp);

double getRoot_Halley (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p);

double getRoot_Schroder (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p);

bool getSorting(const std::tuple<double, double> &a, const std::tuple<double, double> &b);

}
}
#endif // GET_STUFF_H
