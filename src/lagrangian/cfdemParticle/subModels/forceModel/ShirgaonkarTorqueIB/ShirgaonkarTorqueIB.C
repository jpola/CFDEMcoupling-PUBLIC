/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "ShirgaonkarTorqueIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShirgaonkarTorqueIB, 0);

addToRunTimeSelectionTable
(
        forceModel,
        ShirgaonkarTorqueIB,
        dictionary
        );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
ShirgaonkarTorqueIB::ShirgaonkarTorqueIB
(
        const dictionary &dict,
        cfdemCloud &sm
        )
    :
      forceModel(dict, sm),
      propsDict_(dict.subDict(typeName + "Props")),
      twoDimensional_(false),
      depth_(1),
      velFieldName_(propsDict_.lookup("velFieldName")),
      U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
      pressureFieldName_(propsDict_.lookup("pressureFieldName")),
      p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_))

{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "ShirgaonkarTorqueIB.logDat");
    particleCloud_.probeM().vectorFields_.append("hdtorque"); //first entry must the be the force
    particleCloud_.probeM().writeHeader();


    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        depth_ = readScalar(propsDict_.lookup("depth"));
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
        Info << "depth of domain is assumed to be :" << depth_ << endl;
    }

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0, true); // activate treatExplicit switch

    forceSubM(0).setSwitchesList(3, true); // activate search for verbose switch

    forceSubM(0).setSwitches(0,true);  // enable treatExplicit, otherwise this force would be implicit in slip vel! - IMPORTANT!

    // read those switches defined above, if provided in dict
    for (int iFSub = 0; iFSub < nrForceSubModels(); iFSub++)
    {
        forceSubM(iFSub).readSwitches();
    }

    particleCloud_.checkCG(true);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
ShirgaonkarTorqueIB::~ShirgaonkarTorqueIB()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ShirgaonkarTorqueIB::setForce() const
{
    if(forceSubM(0).verbose())
    {
        Info << "ShirgaonkarTorqueIB setForce" << endl;
    }

    label cellI;
    vector torque;
    vector force;

    volVectorField h = forceSubM(0).IBDragPerV(U_, p_);

    #include "setupProbeModel.H"

    label noParticles = particleCloud_.numberOfParticles();

    for ( int index = 0; index < noParticles; index++)
    {
        torque = vector::zero;

        // it is 2D array but with 1 in second dimm, this is due the dataExchangeModel function?
        for (int subCell = 0;
             subCell < particleCloud_.cellsPerParticle()[index][0];
             subCell++)
        {
            //Get the index of the cell which is inside of a particle with number = index
            cellI = particleCloud_.cellIDs()[index][subCell];

            // Get the reference to cell centres locations in our mesh;
            const volVectorField &cellCentres = h.mesh().C();

            if (cellI > -1) // particle is found
            {
                // Calculate the vector from the center of the particle to the current cell
                vector pos = particleCloud_.position(index);
                vector dist(0, 0, 0);

                //go trough x, y, z components;
                for (int i = 0; i < 3; i++)
                {
                    dist[i] = cellCentres[index][i] - pos[i];
                }
                //drag force equation see ShigaonkarIB.C
                force = h[cellI] * h.mesh().V()[cellI];

                vector localTorque = dist ^ force; // r x F

                torque += localTorque;
            }
        }

        if(forceSubM(0).verbose())
        {
            Info << "Particle " << index <<" torque: " << torque << endl;
        }

        for (int i = 0; i < 3; i++)
        {
            DEMTorques()[index][i] += torque[i];
        }


    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// ************************************************************************* //

