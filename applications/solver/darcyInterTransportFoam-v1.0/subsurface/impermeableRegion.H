if (imperRegion and runTime.timeName() <= endTimeUpdateRegion)
{
    forAll(cellCtr_meshSub,i)
    {
        if (imperProperty == "K" or imperProperty == "all")
        {
            K[i] = K_orig[i];
        }

        if (imperProperty == "theta" or imperProperty == "all")
        {
            theta[i] = theta_orig[i];
        }

        if (soluteTransport or soluteTransportSub)
        {
            CSub[i] = CSub_orig[i];

            if (imperProperty == "alpha_L" or imperProperty == "all")
            {
                alpha_L[i] = alpha_L_orig[i];
            }

            if (imperProperty == "alpha_T" or imperProperty == "all")
            {
                alpha_T[i] = alpha_T_orig[i];
            }
        }

        if (heatTransfer or heatTransferSub)
        {
            TSub[i] = TSub_orig[i];
        }

        label minPointi = 0;
        scalar minDistSqr = magSqr(faceCtr_surBottom[minPointi] - cellCtr_meshSub[i]);

        for (label j=1;j<faceCtr_surBottom.size();j++)
        {
            scalar distSqr = magSqr(faceCtr_surBottom[j] - cellCtr_meshSub[i]);

            if (distSqr < minDistSqr)
            {
                minDistSqr = distSqr;
                minPointi = j;
            }
        }

        if (dryBottomSur.boundaryField()[surfaceBottomPatchID][minPointi] == 1.0 and 
            meshSub.cellCentres()[i].component(vector::Z) > mesh.C().boundaryField()[surfaceBottomPatchID][minPointi].component(vector::Z)+d_imper.value())
        {
            if (imperProperty == "K" or imperProperty == "all")
            {
                K[i] = K_imper.value()*tensor(I);
            }

            if (imperProperty == "theta" or imperProperty == "all")
            {
                theta[i] = theta_imper.value();
            }

            if (soluteTransport or soluteTransportSub)
            {
                CSub[i] = CSub_imper.value();

                if (imperProperty == "alpha_L" or imperProperty == "all")
                {
                    alpha_L[i] = alpha_L_imper.value();
                }

                if (imperProperty == "alpha_T" or imperProperty == "all")
                {
                    alpha_T[i] = alpha_T_imper.value();
                }
            }

            if (heatTransfer or heatTransferSub)
            {
                TSub[i] = TSub_imper.value();
            }
        }
    }

    K.correctBoundaryConditions();
    theta.correctBoundaryConditions();
    alpha_L.correctBoundaryConditions();
    alpha_T.correctBoundaryConditions();
    CSub.correctBoundaryConditions();
    TSub.correctBoundaryConditions();
}
