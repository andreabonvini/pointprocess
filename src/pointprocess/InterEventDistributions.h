//
// Created by Andrea Bonvini on 16/04/21.
//

#ifndef POINTPROCESS_INTEREVENTDISTRIBUTIONS_H
#define POINTPROCESS_INTEREVENTDISTRIBUTIONS_H


namespace pointprocess{
    enum Distributions : unsigned  char
    {
        Gaussian = 0,
        InverseGaussian = 1,
        LogNormal = 2,

    };
}
#endif //POINTPROCESS_INTEREVENTDISTRIBUTIONS_H
