// Create a single cell mesh to hold OpenFOAM thermo
// Note that the actual geometry does not affect anything
Foam::argList args(argc, argv);
Foam::Time runTime(Foam::Time::controlDictName, args);
Foam::Time runTime_helper(Foam::Time::controlDictName, args);
Foam::labelList owner(6, Foam::label(0));
Foam::labelList owner_helper(6, Foam::label(0));
Foam::labelList neighbour(0);
Foam::labelList neighbour_helper(0);

Foam::pointField points(8);
Foam::pointField points_helper(8);
points[0] = Foam::vector(0, 0, 0);
points[1] = Foam::vector(1, 0, 0);
points[2] = Foam::vector(1, 1, 0);
points[3] = Foam::vector(0, 1, 0);
points[4] = Foam::vector(0, 0, 1);
points[5] = Foam::vector(1, 0, 1);
points[6] = Foam::vector(1, 1, 1);
points[7] = Foam::vector(0, 1, 1);
points_helper[0] = Foam::vector(0, 0, 0);
points_helper[1] = Foam::vector(1, 0, 0);
points_helper[2] = Foam::vector(1, 1, 0);
points_helper[3] = Foam::vector(0, 1, 0);
points_helper[4] = Foam::vector(0, 0, 1);
points_helper[5] = Foam::vector(1, 0, 1);
points_helper[6] = Foam::vector(1, 1, 1);
points_helper[7] = Foam::vector(0, 1, 1);

const Foam::cellModel& hexa = *(Foam::cellModeller::lookup("hex"));
const Foam::cellModel& hexa_helper = *(Foam::cellModeller::lookup("hex"));
Foam::faceList faces = hexa.modelFaces();
Foam::faceList faces_helper = hexa_helper.modelFaces();

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::READ_IF_PRESENT
    ),
    Foam::xferMove<Foam::Field<Foam::vector>>(points),
    faces.xfer(),
    owner.xfer(),
    neighbour.xfer()
);

Foam::fvMesh mesh_helper
(
    Foam::IOobject
    (
        fvMesh::defaultRegion,
        runTime_helper.timeName(),
        runTime_helper,
        IOobject::READ_IF_PRESENT
    ),
    Foam::xferMove<Foam::Field<Foam::vector>>(points_helper),
    faces_helper.xfer(),
    owner_helper.xfer(),
    neighbour_helper.xfer()
);

Foam::List<Foam::polyPatch*> patches(1);
Foam::List<Foam::polyPatch*> patches_helper(1);

patches[0] = new Foam::emptyPolyPatch
(
    "boundary",
    6,
    0,
    0,
    mesh.boundaryMesh(),
    Foam::emptyPolyPatch::typeName
);
patches_helper[0] = new Foam::emptyPolyPatch
(
    "boundary",
    6,
    0,
    0,
    mesh_helper.boundaryMesh(),
    Foam::emptyPolyPatch::typeName
);

mesh.addFvPatches(patches);
mesh_helper.addFvPatches(patches);

{
    Foam::volScalarField Ydefault
    (
        Foam::IOobject
        (
            "Ydefault",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("Ydefault", dimless, 1)
    );

    Ydefault.write();

    Foam::volScalarField p
    (
        Foam::IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("p", dimPressure, 101325)
    );

    p.write();

    Foam::volScalarField T
    (
        Foam::IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("T", dimTemperature, 300)
    );

    T.write();
}
{
    Foam::volScalarField Ydefault_helper
    (
        Foam::IOobject
        (
            "Ydefault_h",
            runTime_helper.timeName(),
            mesh_helper,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh_helper,
        Foam::dimensionedScalar("Ydefault", dimless, 1)
    );

    Ydefault_helper.write();

    Foam::volScalarField p_helper
    (
        Foam::IOobject
        (
            "p_h",
            runTime_helper.timeName(),
            mesh_helper,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh_helper,
        Foam::dimensionedScalar("p", dimPressure, 101325)
    );

    p_helper.write();

    Foam::volScalarField T_helper
    (
        Foam::IOobject
        (
            "T_h",
            runTime_helper.timeName(),
            mesh_helper,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh_helper,
        Foam::dimensionedScalar("T", dimTemperature, 300)
    );

    T_helper.write();
}
