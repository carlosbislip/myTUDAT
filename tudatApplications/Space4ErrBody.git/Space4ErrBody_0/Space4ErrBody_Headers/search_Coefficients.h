#ifndef SEARH_COEFFICIENTS_H
#define SEARH_COEFFICIENTS_H

double getCm (
        const double &AoA,
        const double &M);

double getCL (
        const double &AoA,
        const double &M);

double getCD (
        const double &AoA,
        const double &M,
        const double &h);

#endif // SEARH_COEFFICIENTS_H
