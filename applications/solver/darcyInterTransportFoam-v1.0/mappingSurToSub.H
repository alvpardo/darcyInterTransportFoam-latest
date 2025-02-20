bool mapSurToSub = false;

// Mapping the surface bottom pressure to the subsurface top patch
if (flow or (flowSur and flowSub))
{
    mapSurToSub = true;

    Info << "Mapping surface head to subsurface domain..." << endl;

    hSub.boundaryFieldRef()[subSurfaceTopPatchID] *= 0.0;
}

// Mapping the surface bottom temperature to the subsurface top patch
if (heatTransfer or (heatTransferSur and heatTransferSub))
{
    mapSurToSub = true;

    Info << nl << "Mapping surface temperature to subsurface domain..." << endl;
}

// Mapping the surface bottom solute to the subsurface top patch
if (soluteTransport or (soluteTransportSur and soluteTransportSub))
{
    mapSurToSub = true;

    Info << nl << "Mapping surface solute to subsurface domain..." << endl;
}

// Mapping of values from the surface bottom to the subsurface top patch
if (mapSurToSub)
{
    scalarField Ctmp_sub(CSub.boundaryField()[subSurfaceTopPatchID].size(),Zero);
    scalarField Ttmp_sub(TSub.boundaryField()[subSurfaceTopPatchID].size(),Zero);

    // Selection of the matching nodes of both meshes at the interface patch
    // for the suitable mapping of values between them
    forAll(faceCtr_subTop,i)
    {
        label minPointi = 0;
        scalar minDistSqr = magSqr(faceCtr_surBottom[minPointi] - faceCtr_subTop[i]);

        for (label j=1;j<faceCtr_surBottom.size();j++)
        {
            scalar distSqr = magSqr(faceCtr_surBottom[j] - faceCtr_subTop[i]);

            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minPointi = j;
            }
        }

        // Assignment of values from the surface to the subsurface domain
        if (flow or (flowSur and flowSub))
        {
            hSub.boundaryFieldRef()[subSurfaceTopPatchID][i] = hSur.boundaryField()[surfaceBottomPatchID][minPointi];
        }

        if (heatTransfer or (heatTransferSur and heatTransferSub))
        {
            Ttmp_sub[i] = TSur.boundaryField()[surfaceBottomPatchID][minPointi];
        }

        if (soluteTransport or (soluteTransportSur and soluteTransportSub))
        {
            Ctmp_sub[i] = CSur.boundaryField()[surfaceBottomPatchID][minPointi];
        }
    }

    if (heatTransfer or (heatTransferSur and heatTransferSub))
    {
        subInterfaceHeatFvPatchScalarField& tempTSubp =
            dynamic_cast<subInterfaceHeatFvPatchScalarField&>(TSub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempTSubp.TFromSur() = Ttmp_sub;
    }

    if (soluteTransport or (soluteTransportSur and soluteTransportSub))
    {
        subInterfaceSoluteFvPatchScalarField& tempCSubp =
            dynamic_cast<subInterfaceSoluteFvPatchScalarField&>(CSub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempCSubp.CFromSur() = Ctmp_sub;
    }
}

Info << nl << endl;
