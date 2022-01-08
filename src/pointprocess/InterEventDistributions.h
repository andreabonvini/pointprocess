//
// Created by Andrea Bonvini on 16/04/21.
//

#ifndef POINTPROCESS_INTEREVENTDISTRIBUTIONS_H
#define POINTPROCESS_INTEREVENTDISTRIBUTIONS_H


enum PointProcessDistributions : unsigned  char
{
    // TODO: Pay attention! Do not change the order of the following distributions (there must be a 1 to 1 mapping
    //  between these values and the one defined in the Python enum in pointprocesslib.py
    //  FIXME: Find a way to automatically test the mapping.
    Gaussian = 0,
    InverseGaussian = 1,
    LogNormal = 2,

};

#endif //POINTPROCESS_INTEREVENTDISTRIBUTIONS_H
