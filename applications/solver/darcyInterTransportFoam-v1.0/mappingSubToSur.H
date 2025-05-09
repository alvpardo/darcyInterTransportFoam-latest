bool mapSubToSur = false;

// Mapping the subsurface top velocity to the surface bottom patch
// Using Darcy velocity for consistency in fluid flux across 
// the interface. May not work if the upwelling/downwelling is significant!
// Only use the normal component may help because the only thing matters
// is the normal velocity for flux.
if (flow or (flowSur and flowSub))
{
    mapSubToSur = true;

    Info << "Mapping subsurface velocity to surface domain..." << endl;

    U.boundaryFieldRef()[surfaceBottomPatchID] *= 0.0;
}

// Mapping the subsurface top temperature to the surface bottom patch
if (heatTransfer or (heatTransferSur and heatTransferSub))
{
    mapSubToSur = true;

    Info << nl << "Mapping subsurface temperature to surface domain..." << endl;
}

// Mapping the subsurface top solute to the surface bottom patch
if (soluteTransport or (soluteTransportSur and soluteTransportSub))
{
    mapSubToSur = true;

    Info << nl << "Mapping subsurface solute to surface domain..." << endl;
}

// Mapping of values from the subsurface top to the surface bottom patch
if (mapSubToSur)
{
    scalarField Ctmp_sur(CSur.boundaryField()[surfaceBottomPatchID].size(),Zero);
    scalarField Ttmp_sur(TSur.boundaryField()[surfaceBottomPatchID].size(),Zero);

    // Selection of the matching nodes of both meshes at the interface patch
    // for the suitable mapping of values between them
    forAll(faceCtr_surBottom,i)
    {
        label minPointi = 0;
        scalar minDistSqr = magSqr(faceCtr_subTop[minPointi] - faceCtr_surBottom[i]);

        for (label j=1;j<faceCtr_subTop.size();j++)
        {
            scalar distSqr = magSqr(faceCtr_subTop[j] - faceCtr_surBottom[i]);

            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minPointi = j;
            }
        }

        // Assignment of values from the subsurface to the surface domain
        if (flow or (flowSur and flowSub))
        {
            U.boundaryFieldRef()[surfaceBottomPatchID][i] = q.boundaryField()[subSurfaceTopPatchID][minPointi];
        }

        if (heatTransfer or (heatTransferSur and heatTransferSub))
        {
            Ttmp_sur[i] = TSub.boundaryField()[subSurfaceTopPatchID][minPointi];
        }

        if (soluteTransport or (soluteTransportSur and soluteTransportSub))
        {
            Ctmp_sur[i] = CSub.boundaryField()[subSurfaceTopPatchID][minPointi];
        }
    }

    if (heatTransfer or (heatTransferSur and heatTransferSub))
    {
        surInterfaceHeatFvPatchScalarField& tempTSurp =
            dynamic_cast<surInterfaceHeatFvPatchScalarField&>(TSur.boundaryFieldRef()[surfaceBottomPatchID]);

        tempTSurp.TFromSub() = Ttmp_sur;
    }

    if (soluteTransport or (soluteTransportSur and soluteTransportSub))
    {
        surInterfaceSoluteFvPatchScalarField& tempCSurp =
            dynamic_cast<surInterfaceSoluteFvPatchScalarField&>(CSur.boundaryFieldRef()[surfaceBottomPatchID]);

        tempCSurp.CFromSub() = Ctmp_sur;
    }
}

Info << nl << endl;
