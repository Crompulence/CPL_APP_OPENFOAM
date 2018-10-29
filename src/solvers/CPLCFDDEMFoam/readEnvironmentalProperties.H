Info<< "\nReading environmentalProperties" << endl;

IOdictionary environmentalProperties
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedVector g(environmentalProperties.lookup("g"));

// default max alpha if not defined in dictionary
//bool coupled;
//coupled = environmentalProperties.lookupOrDefault<bool>("coupled", true);
if (!environmentalProperties.found("coupled"))
{
    Info<< "*** coupled set to true by default" << endl;
}
else
{
    Info<< "*** coupled set to " << coupled << endl;
}
