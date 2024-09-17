/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "truncated3D.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

//#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(truncated3D, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    truncated3D,
    externalWaveForcing
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


truncated3D::truncated3D
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh),

    waveProps_(io.db().lookupObject<IOdictionary>("waveProperties")),
    coeffDict_(waveProps_.subDict("externalForcingCoeffs")),
    
    seaLevel_(readScalar(waveProps_.lookup("seaLevel"))),
    
    NSurface_(readScalar(waveProps_.lookup("NSurface"))),
    timeUandEta_(readScalar(waveProps_.lookup("timeUandEta"))),
    Nsample(readScalar(waveProps_.lookup("Nsample"))),
    Nsets(readScalar(waveProps_.lookup("Nsets"))),
    
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),
    g_( uniformDimensionedVectorField
        (
            mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
        ).value() ),
    direction_( g_/mag(g_) )
{
    {
        IOdictionary transProp
        (
            IOobject
            (
                "transportProperties",
                "constant",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        dictionary sD(transProp.subDict(Foam::waves2Foam::waterPhase()));
        rhoWater_ = (dimensionedScalar(sD.lookup("rho"))).value();
        IFstream dataStream1("kinematics/surfaceElevation.dat");
	for (label i=0; i<NSurface_; i++) //for (label i=0; i<2002; i++)
        {
            dataStream1 >> ttEta[i] >> yyEta[i] >> yyEta3 >> yyEta3 >> yyEta3;
        }
        IFstream dataStream2("kinematics/set.dat");
	for (label i=0; i<Nsets; i++)
        {
            dataStream2 >> readTime[i];
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void truncated3D::step()
{
    // Nothing to be done
}


scalar truncated3D::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar OFtime = time + timeUandEta_; //scalar OFtime = time + 10;  // 32 can be user input in waveproperties say.
    scalar yy2 = 0;

    for (label i=0; i<NSurface_-1; i++) //for (label i=0; i<2002-1; i++)
    {
        if ( ttEta[i] >= OFtime )
        {
            yy2 = yyEta[i-1] + ((yyEta[i] - yyEta[i-1])*(OFtime - ttEta[i-1])/(ttEta[i] - ttEta[i-1]));
            break;
        }
    }

    return yy2;
}


scalar truncated3D::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    // SAME IMPLEMENTATION IS OF P BUT THIS IS MEANT FOR PRESSURE GRADIENT
    scalar Z(returnZ(x));

    label index=0;
    scalar OFtime = time;

    for (label i=0; i<Nsets-1; i++)
    {
        if ( readTime[i] <= (OFtime + timeUandEta_) && readTime[i+1] >= (OFtime + timeUandEta_) )
            {
                index = i;
            }
    }

    std::ostringstream os;
    os << readTime[index];

    std::ostringstream os1;
    os1 << readTime[index+1];

    scalar xx[999], yy[999], zz[999], P[999];

	IFstream dataStream1("kinematics/sets/"+ os.str() + "/lineX1_dPdX.xy");
    	for (label i=0; i<Nsample; i++)
    	{
        	dataStream1 >> xx[i] >> yy[i] >> zz[i] >> P[i];
    	}

    scalar xx1[999], yy1[999], zz1[999], P1[999];

	IFstream dataStream2("kinematics/sets/"+ os1.str() + "/lineX1_dPdX.xy");
	for (label i=0; i<Nsample; i++)
    	{
        	dataStream2 >> xx1[i] >> yy1[i] >> zz1[i] >> P1[i];
    	}

    scalar Pp = 0; scalar Pp1 = 0; scalar Pp2 = 0;

    for (label i=0; i<Nsample-1; i++)
    {
        if ( yy[i] <= Z && yy[i+1] >= Z )
        {
            // do linear interpolation for Uxx, Uyy and Uzz, later on do interpolation in time if necessary
            Pp = P[i] + ((P[i+1] - P[i])*(Z - yy[i])/(yy[i+1] - yy[i]));

            Pp1 = P1[i] + ((P1[i+1] - P1[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
        }
    }

    Pp2 = Pp + ((Pp1 - Pp)*(OFtime+timeUandEta_ - readTime[index])/(readTime[index+1] - readTime[index]));

    return Pp2;
}


scalar truncated3D::pExcess
(
    const point& x,
    const scalar& time
) const
{
    // THIS PART IS NOT BEING READ BY EXTERNALWAVEFORCING
	
    scalar Z(returnZ(x));

    label index=0;
    scalar OFtime = time;

    for (label i=0; i<Nsets-1; i++)
    {
        if ( readTime[i] <= (OFtime + timeUandEta_) && readTime[i+1] >= (OFtime + timeUandEta_) )
            {
                index = i;
            }
    }

    std::ostringstream os;
    os << readTime[index];

    std::ostringstream os1;
    os1 << readTime[index+1];

    scalar xx[999], yy[999], zz[999], P[999];

	IFstream dataStream1("kinematics/sets/"+ os.str() + "/lineX1_P.xy");
    	for (label i=0; i<Nsample; i++)
    	{
        	dataStream1 >> xx[i] >> yy[i] >> zz[i] >> P[i];
    	}

    scalar xx1[999], yy1[999], zz1[999], P1[999];

	IFstream dataStream2("kinematics/sets/"+ os1.str() + "/lineX1_P.xy");
	for (label i=0; i<Nsample; i++)
    	{
        	dataStream2 >> xx1[i] >> yy1[i] >> zz1[i] >> P1[i];
    	}

    scalar Pp = 0; scalar Pp1 = 0; scalar Pp2 = 0;

    for (label i=0; i<Nsample-1; i++)
    {
        if ( yy[i] <= Z && yy[i+1] >= Z )
        {
            // do linear interpolation for Uxx, Uyy and Uzz, later on do interpolation in time if necessary
            Pp = P[i] + ((P[i+1] - P[i])*(Z - yy[i])/(yy[i+1] - yy[i]));

            Pp1 = P1[i] + ((P1[i+1] - P1[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
        }
    }

    Pp2 = Pp + ((Pp1 - Pp)*(OFtime+timeUandEta_ - readTime[index])/(readTime[index+1] - readTime[index]));

    return Pp2;
}


vector truncated3D::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    label index=0;
    scalar OFtime = time;

    for (label i=0; i<Nsets-1; i++)
    {
        if ( readTime[i] <= (OFtime + timeUandEta_) && readTime[i+1] >= (OFtime + timeUandEta_) )
            {
                index = i;
            }
    }

    std::ostringstream os;
    os << readTime[index];

    std::ostringstream os1;
    os1 << readTime[index+1];

    scalar xx[999], yy[999], zz[999], Uxx[999], Uyy[999], Uzz[999];

		IFstream dataStream1("kinematics/sets/"+ os.str() + "/lineX1_U.xy");
    		for (label i=0; i<Nsample; i++)
    			{
        			dataStream1 >> xx[i] >> yy[i] >> zz[i] >> Uxx[i] >> Uyy[i] >> Uzz[i];
    			}

    scalar xx1[999], yy1[999], zz1[999], Uxx1[999], Uyy1[999], Uzz1[999];

		IFstream dataStream2("kinematics/sets/"+ os1.str() + "/lineX1_U.xy");
		for (label i=0; i<Nsample; i++)
    			{
        			dataStream2 >> xx1[i] >> yy1[i] >> zz1[i] >> Uxx1[i] >> Uyy1[i] >> Uzz1[i];
    			}

    scalar UUxx = 0; scalar UUyy = 0; scalar UUzz = 0;
    scalar UUxx1 = 0; scalar UUyy1 = 0; scalar UUzz1 = 0;
    scalar UUxx2 = 0; scalar UUyy2 = 0; scalar UUzz2 = 0;
    for (label i=0; i<Nsample-1; i++)
    {
        if ( yy[i] <= Z && yy[i+1] >= Z )
        {
            UUxx = Uxx[i] + ((Uxx[i+1] - Uxx[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
            UUyy = Uyy[i] + ((Uyy[i+1] - Uyy[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
            UUzz = Uzz[i] + ((Uzz[i+1] - Uzz[i])*(Z - yy[i])/(yy[i+1] - yy[i]));

            UUxx1 = Uxx1[i] + ((Uxx1[i+1] - Uxx1[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
            UUyy1 = Uyy1[i] + ((Uyy1[i+1] - Uyy1[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
            UUzz1 = Uzz1[i] + ((Uzz1[i+1] - Uzz1[i])*(Z - yy[i])/(yy[i+1] - yy[i]));
        }
    }

    UUxx2 = UUxx + ((UUxx1 - UUxx)*(OFtime+timeUandEta_ - readTime[index])/(readTime[index+1] - readTime[index]));
    UUyy2 = UUyy + ((UUyy1 - UUyy)*(OFtime+timeUandEta_ - readTime[index])/(readTime[index+1] - readTime[index]));
    UUzz2 = UUzz + ((UUzz1 - UUzz)*(OFtime+timeUandEta_ - readTime[index])/(readTime[index+1] - readTime[index]));

    return UUxx2*k_/K_ - UUyy2*direction_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
