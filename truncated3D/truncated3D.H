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

Class
    Foam::waveTheories::truncated3D

Description


SourceFiles
    truncated3D.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.
    Bo Terp Paulsen, Deltares, All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef truncated3D_H
#define truncated3D_H

#include "externalWaveForcing.H"
#include "waveTheory.H"

//#include "relaxationShape.H"
//#include "crossVersionCompatibility.H"
//#include "setWaveField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class truncated3D Declaration
\*---------------------------------------------------------------------------*/

class truncated3D
:
    public externalWaveForcing
{

protected:

    // Protected data
        const IOdictionary& waveProps_;
        const dictionary& coeffDict_;

        scalar seaLevel_;

	scalar NSurface_;
	scalar timeUandEta_;
	scalar Nsample;
 	scalar Nsets;

        scalar rhoWater_;
        //scalar Tsoft_;
        vector g_;
        vector direction_;
        vector wind_;
        vector k_;
        scalar K_;

	scalar ttEta[999999], yyEta[999999], ttEta1[999999], yyEta3, readTime[999999];

    // Protected member functions
        inline scalar returnZ(const point& x) const
        {
            return (-(direction_ & x) - seaLevel_);
        }


public:

    //- Runtime type information
    TypeName("truncated3D");

    // Constructors

        //- from components
        truncated3D
        (
            IOobject,
            Time&,
            const fvMesh&
        );


    // Destructor

        ~truncated3D()
        {}


    // Member Functions
        virtual void step();

        virtual inline void close()
        {
            // Nothing to be done
        };


        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;


        virtual scalar ddxPd
        (
            const point&,
            const scalar&,
            const vector&
        ) const;


        virtual vector U
        (
            const point&,
            const scalar&
        ) const;


        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;
        
        
        bool implementPressure() const
        {
            return false;
        };


        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
