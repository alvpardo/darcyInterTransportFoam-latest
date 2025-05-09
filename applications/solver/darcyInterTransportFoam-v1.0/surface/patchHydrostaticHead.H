    const fvPatchList& numPatches(mesh.boundary());

    Info << "Number of patches: " << numPatches.size() << endl;

    // Calculation of the hydrostatic head at the inlet, outlet and riverbed
    // patches of the surface domain (m)
    for (label patchID = 0;patchID < 6;patchID++)
    {
        string patchName = mesh.boundaryMesh()[patchID].name();

        if (patchName == "inletWater" || patchName == "outlet" || patchName == interPatchName)
        {
            fvPatchScalarField piezo = piezoSur.boundaryField()[patchID];

            if (patchName == "inletWater" || patchName == "outlet")
            {
                piezo *= pos0(alpha1.boundaryField()[patchID] - 0.5);

                if (patchName == "inletWater")
                {
                    patchName = "inlet";

                    label subPatchID = meshSub.boundaryMesh().findPatchID(patchName);

                    // Check the patch type to transfer the data to the BC
                    if (isA<surfaceHydrostaticHeadFvPatchScalarField>(hSub.boundaryField()[subPatchID]))
                    {
                        surfaceHydrostaticHeadFvPatchScalarField& temphInSubp =
                            dynamic_cast<surfaceHydrostaticHeadFvPatchScalarField&>(hSub.boundaryFieldRef()[subPatchID]);

                        temphInSubp.hInFromSur() = gMax(piezo);
                    }
                }

                if (patchName == "outlet")
                {
                    label subPatchID = meshSub.boundaryMesh().findPatchID(patchName);

                    // Check the patch type to transfer the data to the BC
                    if (isA<surfaceHydrostaticHeadFvPatchScalarField>(hSub.boundaryField()[subPatchID]))
                    {
                        surfaceHydrostaticHeadFvPatchScalarField& temphOutSubp =
                            dynamic_cast<surfaceHydrostaticHeadFvPatchScalarField&>(hSub.boundaryFieldRef()[subPatchID]);

                        temphOutSubp.hOutFromSur() = gMax(piezo);
                    }
                }

                inOutStatic_hSur.set(patchName,gMax(piezo));
                inOutStatic_hSur.regIOobject::write();
            }
            else
            {
                // Check the patch type to transfer the data to the BC
                if (isA<surfaceHydrostaticHeadFvPatchScalarField>(hSub.boundaryField()[subSurfaceBottomPatchID]))
                {
                    scalarField htmp_sub(hSub.boundaryField()[subSurfaceBottomPatchID].size(),Zero);

                    forAll(faceCtr_subBottom,i)
                    {
                        label minPointi = 0;
                        scalar minDistSqr = magSqr(faceCtr_surBottom[minPointi] - faceCtr_subBottom[i]);

                        for (label j=1;j<faceCtr_surBottom.size();j++)
                        {
                            scalar distSqr = magSqr(faceCtr_surBottom[j] - faceCtr_subBottom[i]);

                            if (distSqr < minDistSqr)
                            {
                                minDistSqr = distSqr;
                                minPointi = j;
                            }
                        }

                        htmp_sub[i] = piezo[minPointi];
                    }

                    surfaceHydrostaticHeadFvPatchScalarField& temphRiverSubp =
                        dynamic_cast<surfaceHydrostaticHeadFvPatchScalarField&>(hSub.boundaryFieldRef()[subSurfaceBottomPatchID]);

                    //temphRiverSubp.hRiverFromSur() = surBottomTosubBottomInter.faceInterpolate<scalar>(piezo);
                    temphRiverSubp.hRiverFromSur() = htmp_sub;
                }
            }
        }
    }
