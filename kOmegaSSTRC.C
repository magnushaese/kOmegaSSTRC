/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTRC.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSSTRC, 0);
addToRunTimeSelectionTable(RASModel, kOmegaSSTRC, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOmegaSSTRC::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> kOmegaSSTRC::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


tmp<volScalarField> kOmegaSSTRC::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*nu()/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


tmp<volScalarField> kOmegaSSTRC::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTRC::kOmegaSSTRC
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    alphaK1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "b1",
            coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            coeffDict_,
            false
        )
    ),
    Cr1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cr1",
            coeffDict_,
            1.0
        )
    ),
    Cr2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cr2",
            coeffDict_,
            2.0
        )
    ),
    Cr3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cr3",
            coeffDict_,
            1.0
        )
    ),

    spalartShurCorrection_
    (
        Switch::lookupOrAddToDict
	(
		"spalartShurCorrection", 
		coeffDict_,
		false
	)
    ),
    Arollabifurcation_
    (
        Switch::lookupOrAddToDict
        (
                "Arollabifurcation",
                coeffDict_,
                false
        )
    ),
    DhakalAndWalters_
    (
        Switch::lookupOrAddToDict
        (
                "DhakalAndWalters",
                coeffDict_,
                false
        )
    ),
    y_(mesh_),

    fr1_
    (
        IOobject
        (
            "fr1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
	dimensionedScalar("fr1",dimless,1)
    ),
    Cmu_
    (
        IOobject
        (
            "Cmu",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Cmu",dimless,1)
    ),
    eta1_
    (
        IOobject
        (
            "eta1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eta1",dimless,0)
    ),
    eta2_
    (
        IOobject
        (
            "eta2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eta2",dimless,0)
    ),
    eta_
    (
        IOobject
        (
            "eta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eta",dimless,1)
    ), 
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_, U_.db())
    ),
    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("v2", mesh_, U_.db())
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_, U_.db())
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_, U_.db())
    )
{
    bound(k_, k0_);
    bound(omega_, omega0_);

    nut_ =
    (
        a1_*k_/
        max
        (
            a1_*omega_,
            b1_*F23()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );

    if(Arollabifurcation_)
    {

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volSymmTensorField S(symm(tgradU()));
    volTensorField W(skew(tgradU()));

	volSymmTensorField divS =
	    (
		fvc::ddt(S)
	       +fvc::div
		(
		    phi_, S
		)
	    );
	volSymmTensorField DSDtTMP(divS);
        volTensorField OmegaSS(((S & DSDtTMP)- (DSDtTMP & S))/(2*magSqr(S) + dimensionedScalar("small", dimensionSet(0, 0, -2, 0, 0), 1E-14)));

        volTensorField OmegaMod = W; // - OmegaSS;
        volScalarField T1(1/(betaStar_*omega_));
        volScalarField T2(scalar(6)*sqrt(nu()/(betaStar_*k_*omega_)));
        volScalarField T3(pow(pow(T1,1.625)*T2,1/2.625));
        volScalarField T = max(T1,T3); 

        eta1_ = magSqr(S)*pow(T,2.0);
        eta2_.storePrevIter();
	eta2_ = magSqr(OmegaMod)*pow(T,2.0);
	eta2_.relax();	

        volScalarField eta3(eta1_ - eta2_);      

	Cmu_.storePrevIter();
        Cmu_ = min
          (
            pow((0.04645*(mag(eta3) - (eta3)) + sqrt(1-min(0.25*(eta3), 0.99))),-1),
            2.5
          );
	Cmu_.relax();
        nut_ = Cmu_*k_/omega_;	
    }
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaSSTRC::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kOmegaSSTRC::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kOmegaSSTRC::divDevReff() const
{
    return
    (
      - fvm::laplacian(nuEff(), U_)
      - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
    );
}


bool kOmegaSSTRC::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        b1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        F3_.readIfPresent("F3", coeffDict());
	Cr1_.readIfPresent(coeffDict());
        Cr2_.readIfPresent(coeffDict());
        Cr3_.readIfPresent(coeffDict());
        spalartShurCorrection_.readIfPresent
        (
            "spalartShurCorrection", coeffDict()
        );
        Arollabifurcation_.readIfPresent
        (
            "Arollabifurcation", coeffDict()
        );
        DhakalAndWalters_.readIfPresent
        (
            "DhakalAndWalters", coeffDict()
        );
        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaSSTRC::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(omega_, omega0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField G("RASModel::G", nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));


////////////////Spalart Shur Menter 2009/////////////
  if (spalartShurCorrection_)

	{
    tmp<volTensorField> tgradU = fvc::grad(U_);
    tmp<volTensorField> tSkew = skew(tgradU()); // tmp<volTensorField> tSkew = skew(fvc::grad(U)); 
tmp<volSymmTensorField> tSymm = symm(tgradU()); // tmp<volSymmTensorField> tSymm = symm(fvc::grad(U));
// Compute rStar
    volScalarField symInnerProduct(2.0*tSymm() && tSymm()); 
    volScalarField asymInnerProduct
    (
        max(2.0*tSkew() && tSkew(),
        dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0.0))
    );
    volScalarField w
    (
        atan(dimensionedScalar("4",dimensionSet(0,0,2,0,0),1.0e-02)*asymInnerProduct)*2.0/scalar(3.14159)*(asymInnerProduct-symInnerProduct)+symInnerProduct
    );
    volScalarField rStar(sqrt(symInnerProduct/max(w, dimensionedScalar("minw", w.dimensions(), SMALL)))); //avoiding dividing by zero    volScalarField rStar(sqrt(symInnerProduct/w));
     
    // Compute rTilda 
    volScalarField D(sqrt(max(symInnerProduct, 0.09*omega_*omega_)));
    tmp<volSymmTensorField> divS =
    (
        fvc::ddt(tSymm())
       +fvc::div
        (
            phi_, tSymm()
        )
    );
    volScalarField rT((tSkew().T() & tSymm) && divS);

    divS.clear();
    tSkew.clear();
    tSymm.clear();
    
    volScalarField w2
    (
        atan(dimensionedScalar("1",dimensionSet(0,0,2,0,0),1.0e-2)*asymInnerProduct)*2.0/scalar(3.14159)*(sqrt(asymInnerProduct)-D)+D //T
    );
    volScalarField rTilda(2.0*rT/w2/D/D/D);

    (
fr1_ =	max
    	(
		min
		(
        		(scalar(1.0)+Cr1_)*2.0*rStar/(scalar(1)+rStar)*(scalar(1.0)
			-Cr3_*atan(Cr2_*rTilda))-Cr1_,
        		scalar(1.25)
		),
			scalar(0.0)
	)
);
/*
        gamma1_ = 0.35705; //0.5555; 
        gamma2_ = 0.4403546667;
        alphaOmega1_ = 0.65; //0.5;

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volSymmTensorField S(symm(tgradU()));
    volTensorField W(skew(tgradU()));

	volSymmTensorField divS =
	    (
		fvc::ddt(S)
	       +fvc::div
		(
		    phi_, S
		)
	    );

	volSymmTensorField DSDtTMP(divS);
        volTensorField OmegaMod = W;

        volScalarField sqrS(2.0*magSqr(S));
        volScalarField sqrOmegaMod(2.0*magSqr(OmegaMod));  
        volScalarField sqrD(max(sqrS, 0.09*sqr(omega_)));           
    

        volScalarField rStar(sqrt(sqrS)/sqrt(sqrOmegaMod));

    volScalarField rT(
			(OmegaMod.T() & S) && DSDtTMP			
		     );

    volScalarField rTilda(2.0*rT/((sqrD)*sqrt(sqrD)*sqrt(sqrOmegaMod)));
             
        fr1_ = max
        (
            min
            (
                (1 + Cr1_)*2.0*rStar/(1.0 + rStar)*(1.0 - Cr3_*atan(Cr2_*rTilda)) - Cr1_,
                1.25
            ),
            0.0
        );  
*/  
    }

/////////////////////////////////////////////////////
/////////Dhakal and Walters 2011 /////////////////////////
 if (DhakalAndWalters_)
	{
	 //alphaOmega1_ = 0.5; //0.65;
	 gamma1_ = 0.35705; //beta1_/betaStar_-alphaOmega1_*kappa_*kappa_/sqrt(betaStar_)
         gamma2_ = 0.4403546667; //beta2_/betaStar_-alphaOmega2_*kappa_*kappa_/sqrt(betaStar_)
	 a1_ = 0.31;

            alphaK1_=0.5;
            alphaK2_=1.0;
            alphaOmega1_=0.5;
            alphaOmega2_=0.856;
            gamma1_=0.55555; //5.0/9.0;
            gamma2_=0.44; //0.44; //0.7;
            beta1_=0.0750;
            beta2_=0.0828;

            betaStar_=0.09;


    tmp<volTensorField> tgradU = fvc::grad(U_);
    volSymmTensorField S(symm(tgradU()));
    volTensorField W(skew(tgradU()));


	volSymmTensorField divS =
	    (
		fvc::ddt(S)
	       +fvc::div
		(
		    phi_, S
		)
	    );

	volSymmTensorField DSDtTMP(divS);
        volTensorField OmegaSS(((S & DSDtTMP)- (DSDtTMP & S))/(2*magSqr(S) + dimensionedScalar("small", dimensionSet(0, 0, -2, 0, 0), 1E-14))); 

        volTensorField OmegaMod = W;
				+ 2/(0.4-2)*OmegaSS; 


        volScalarField Fw
        (
            tanh(pow4(200*(nu())/(omega_*sqr(y_))))
        );

        volScalarField x
        (      
                scalar(2)/scalar(9)*
                    (
                        scalar(1) - sqrt(magSqr(OmegaMod))/(sqrt(magSqr(S))
                        + dimensionedScalar("10e-16", dimensionSet(0, 0, -1, 0, 0), 10e-16))
                    )                
        );
        
        volScalarField etam
        (
                147.5*pow(x,5)-437.8*pow(x,4)+331.5*pow(x,3)-112*pow(x,2)+18.57*x+scalar(1)
        );

	volScalarField etaMod(max
		(
		scalar(0),
		min(etam,scalar(2.5))
		));
        
        volScalarField etaN
        (
                Fw*min(scalar(1),etaMod)+(scalar(1)-Fw)*etaMod
        );
        
        eta_ = etaN;
}

/////////////////////////////////////////////////////
    // Turbulent frequency equation
    fvScalarMatrix omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      + fvm::SuSp(-fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*fr1_
       *min
        (
            S2,
            (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*sqrt(S2))
        )
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn.relax();

    // No longer needed: matrix completes at the point of solution
    // HJ, 17/Apr/2012
//     omegaEqn.completeAssembly();

    solve(omegaEqn);
    bound(omega_, omega0_);

    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G*fr1_, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn.relax();
    solve(kEqn);
    bound(k_, k0_);

    // Re-calculate viscosity
    // Fixed sqrt(2) error.  HJ, 10/Jun/2015
    nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));

if(DhakalAndWalters_)
{
    // v2 transport equation
    fvScalarMatrix v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      + fvm::SuSp(-fvc::div(phi_), v2_)
      - fvm::laplacian(DkEff(F1), v2_)
     ==
  //   min(G*fr1_, c1_*betaStar_*v2_*omega_)
       min(v2_/k_*G*fr1_, c1_*betaStar_*v2_*omega_)
      - fvm::Sp(betaStar_*omega_, v2_)
      + fvm::SuSp(omega_/v2_*1.8*0.09*pow(eta_,2)*k_, v2_)
      - fvm::SuSp(omega_/v2_*1.8*0.09*v2_, v2_)
    );

    v2Eqn.relax();
    solve(v2Eqn);
    bound(v2_, k0_);
    nut_=sqrt(k_)*sqrt(v2_)/omega_;
}

    if(Arollabifurcation_)
    {

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volSymmTensorField S(symm(tgradU()));
    volTensorField W(skew(tgradU()));

	volSymmTensorField divS =
	    (
		fvc::ddt(S)
	       +fvc::div
		(
		    phi_, S
		)
	    );
	volSymmTensorField DSDtTMP(divS);
        volTensorField OmegaSS(((S & DSDtTMP)- (DSDtTMP & S))/(2*magSqr(S) + dimensionedScalar("small", dimensionSet(0, 0, -2, 0, 0), 1E-14)));

        volTensorField OmegaMod = W; // - OmegaSS;
        volScalarField T1(1/(betaStar_*omega_));
        volScalarField T2(scalar(6)*sqrt(nu()/(betaStar_*k_*omega_)));
        volScalarField T3(pow(pow(T1,1.625)*T2,1/2.625));
        volScalarField T = max(T1,T3); 

        eta1_ = magSqr(S)*pow(T,2.0);
        eta2_.storePrevIter();
	eta2_ = magSqr(OmegaMod)*pow(T,2.0);
	eta2_.relax();	

        volScalarField eta3(eta1_ - eta2_);      

	Cmu_.storePrevIter();
        Cmu_ = min
          (
            pow((0.04645*(mag(eta3) - (eta3)) + sqrt(1-min(0.25*(eta3), 0.99))),-1),
            2.5
          );
	Cmu_.relax();
        nut_ = Cmu_*k_/omega_;	
    }
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
