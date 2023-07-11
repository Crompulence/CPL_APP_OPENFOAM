CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"


/*
volScalarField divU(fvc::div(phi));
CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU,
    pimple
);
*/
