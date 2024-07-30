#include "RBSHelpers.hh"
#include "G4IonYangFluctuationModel.hh"
#include "G4IonChuFluctuationModel.hh"
#include "G4hIonEffChargeSquare.hh"
#include "G4EmCalculator.hh"

// start of functions relevant to RBS spectrum
G4double RecoilEnergy(G4double E, G4double angle, G4double M1, G4double M2)
{
	G4double k;
	G4double square = std::sqrt(std::pow(M2/M1, 2.) - std::pow(sin(angle), 2.));
	G4double M1co = cos(angle);
	G4double denominator = 1 + (M2 / M1);
	k = std::pow((M1co + square) / denominator, 2.);
	return  E * k;
}

G4double CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle, G4double Z1, G4double Z2)
{
	G4double cose = cos(angle);
	G4double sine = sin(angle);
	G4double epsilon = 55.26349406 / 1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 8 * pi * epsilon * E;
	G4double FirstTerm = std::pow((Z1 * Z2) / (FirstTermDenom), 2.);
	G4double SecondTerm = 1 / (std::pow(sine, 4.));
	G4double ThirdTerm = std::pow(M2 * cose + std::pow(std::pow(M2, 2) - std::pow(M1 * sine, 2), 0.5), 2);
	G4double ThirdTermDenom = M2 * std::pow(std::pow(M2, 2) - std::pow(M1 * sine, 2), 0.5);
	G4double fullTerm = FirstTerm * SecondTerm * ThirdTerm / ThirdTermDenom;
	return fullTerm * 10; //fm^2 conversion to milibarns
}

G4double CalcRBSYield(G4double xsec, G4double solidAngle, G4double atomDens, G4double inc_angle)
{
	G4double crossSect = xsec * 1e-27; // mbarn to cm2 = 1e−27
	G4double dist = 0.1*nm;
	G4double z = dist / cm;
	G4double atom_thickness = atomDens * z;
	G4double yield = (solidAngle * crossSect * atom_thickness / cos(inc_angle));
	return yield;
}

//energy straggling in MeV^2
G4double CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist)
{
	G4double elmcharge_squared = 1.4399764; // in units of MeV*fm
	G4double e2 = std::pow(elmcharge_squared, 2.) * 1e-26; // in units of MeV^2*cm^2
	G4double fourpi = 4 * 3.14159265358979;
	G4double fBohr = fourpi * e2;
	G4double Bohr_t = Z1 * Z1 * Z2 * fBohr * atomDens * (dist / cm);
	return Bohr_t;
}
//from http://atlas.physics.arizona.edu/~shupe/Indep_Studies_2015/Notes_Goethe_Univ/L4_Scattering_Basic.pdf
G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
{
	G4double x= 0., y=0.;
	x = (M1 / M2) * sin(angle);
	y = asin(x);
	return y + angle;
}
//energy in the CM reference frame
G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2)
{
	G4double en=0.;
	en = energy * M2 / (M1 + M2);
	return en;
}

//xsec in the CM reference system
G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2)
{
	G4double sine = sin(angleCM / 2);
	G4double epsilon = 55.26349406 / 1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 16 * pi * epsilon * E;
	G4double FirstTerm = std::pow((Z1 * Z2) / (FirstTermDenom), 2.);
	G4double SecondTerm = 1 / (std::pow(sine, 4.));
	G4double fullTerm = FirstTerm * SecondTerm;
	return fullTerm * 10;
}
//xsec in the Lab reference system from the CM reference system
G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection)
{
	G4double ratio = M1 / M2;
	G4double multiplier = std::pow(1 + (ratio * ratio) + (2 * ratio * cos(angle)), 3. / 2.) / (1 + (ratio * cos(angle)));
	G4double labXSEC = multiplier * xsection;
	return labXSEC;
}
// integrates energy loss in distance
G4double CalcTotEnLoss(G4double E, G4double distance, G4int steps, G4ParticleDefinition* fParticle, G4Material* mat)
{
	G4EmCalculator emCalculator;
	G4double stp = (distance / steps) / cm;
	for (G4int i=1; i<=steps; i++) {
		//G4double stop = emCalculator.ComputeTotalDEDX(E,fParticle,mat)/(MeV/cm);	// buvo anksciau, gal del srim?
		G4double stop = emCalculator.GetDEDX(E, fParticle,mat) / (MeV / cm); // naudoti su channeling
		E -= stop * stp;
		if (E / keV <= 10)
			break;
	}
	return E;
}

// function for Total RBS yield, combining other functions into single one
G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double solidAngle, G4double xsecmod, G4double atomDensity,G4double inc_angle)
{
	G4double CMangleX        = CalcAngleCMFrame(angle, M1, M2);
	G4double CMenergyX       = CalcEnergyCMFrame(energy, M1, M2);
	G4double CMxsecX         = CalcDiffRuthXsecCM(CMenergyX, CMangleX, Z1, Z2);
	G4double andersen_factor = CalcAndersenScreening(CMenergyX, CMangleX, Z1, Z2);
	G4double CMxsecModX      = CMxsecX * andersen_factor;
	G4double CMtoLABxsecX    = CalcDiffRuthXsecLAB(M1, M2, CMangleX, CMxsecModX);
	return CalcRBSYield((CMtoLABxsecX * xsecmod), solidAngle, atomDensity, inc_angle);
}

// yang+chu
G4double CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance)
{
	if (energy < 0)
		return 0;

	G4double straggling, nucl_strag, elec_strag;
	straggling = nucl_strag = elec_strag = 0;
	
	G4int elements      = mat->GetNumberOfElements();
	G4double *Z2        = new G4double[elements];
	G4double *M2        = new G4double[elements];
	G4double *aDensity  = new G4double[elements];

	//functions for en loss straggling evaluation
	G4hIonEffChargeSquare* eff_charge = new G4hIonEffChargeSquare(""); // gets effective charge square of ion based on energy
	G4IonYangFluctuationModel* model_yang = new G4IonYangFluctuationModel("");
	G4IonChuFluctuationModel* model_chu = new G4IonChuFluctuationModel("");
	// Particle parameters
	G4double q_squared = eff_charge->TheValue(particle, mat, energy); // effective charge squared
	G4double Z1 = particle->GetAtomicNumber();
	G4double M1 = particle->GetAtomicMass();
	G4double q_sqrt = sqrt(q_squared);
	G4double q_sqrt_z = q_sqrt / Z1;

	G4double yang_sq = model_yang->TheValue(particle, mat, energy);
	G4double chu_sq	= model_chu->TheValue(particle, mat, energy);

	G4double strag_factor = q_sqrt_z * chu_sq + yang_sq;
	const G4double *atomDensVector = mat->GetVecNbOfAtomsPerVolume();

	for (int i=0; i<elements; i++) {
		Z2[i] = mat->GetElement(i)->GetZ();
		M2[i] = mat->GetElement(i)->GetA() / (g/mole);
		aDensity[i] = atomDensVector[i] / (1 / cm3);
	}

	for (int j=0; j<elements;j++) {
		elec_strag += CalcBohrStrag(Z1, Z2[j], aDensity[j], distance) * strag_factor; // electronic straggling straggling, in MeV2
		nucl_strag += CalcNuclEnStraggling(Z1, Z2[j], M1, M2[j], aDensity[j], distance); // nuclear straggling, in MeV2

		if (elec_strag < 0)
			elec_strag = 0.;
		if (nucl_strag < 0)
			nucl_strag = 0.;
	}
	straggling = elec_strag + nucl_strag;

	delete []Z2;
	delete []M2;
	delete []aDensity;
	delete eff_charge;
	delete model_yang;
	delete model_chu;
	return straggling;
}
// energy loss in detector dead layer
G4double CalculateDeadLayerEffect(G4double energy, const G4Material* material, G4double thickness,G4ParticleDefinition* particle)
{
	G4double steps = 20.;
	G4double stp = (thickness / steps) / cm;
	G4EmCalculator emCalculator;
	for (int i=1; i<=steps; i++) {
		G4double stop = emCalculator.ComputeTotalDEDX(energy, particle, material) / (MeV / cm);
		//G4double stop = emCalculator.GetDEDX(energy,particle,material)/(MeV/cm);
		energy -= stop * stp;
		if (energy / keV <= 10)
			break;
	}	
	return energy;
}

// detector fwhm calculation
G4double CalcDetectorFWHM(G4double energy, G4double Z1)
{
	// energy must be in keV
	G4double E = energy / keV;
	G4double fwhm;
	G4double C1 = 0.0999288;
	G4double C2 = 1.1871;
	G4double C3 = 1.94699;
	G4double C4 = 0.18;
	G4double C5 = 2.70004;
	G4double C6 = 9.29965;

	G4double first_term = C1 * std::pow(Z1, C2);
	G4double second_term = std::pow(log(E), C3);
	G4double third_term = C4 * std::pow(Z1, C5);
	G4double fourth_term = std::pow(log(E), C6);

	fwhm = first_term * second_term - (third_term / fourth_term);

	if (fwhm < 10.)
		fwhm = 10.;
	return fwhm;
}

// andersen screening factor
G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2)
{
	G4double V1 = (0.04873 * Z1 * Z2 * std::sqrt(std::pow(Z1, 2. / 3.) + std::pow(Z2, 2. / 3.))) * keV; // in keV
	G4double first_term = std::pow(1 + (0.5 * (V1 / energy_cm)), 2.);
	G4double second_term = V1 / (2 * energy_cm * sin(angle_cm / 2));
	G4double second_term_sq = std::pow(second_term, 2);
	G4double denominator = std::pow(1 + (V1 / energy_cm) + second_term_sq, 2);
	return first_term / denominator;
}

// nuclear energy loss straggling
G4double CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance)
{
	return 0.26 * std::pow(Z1 * Z2, 2) * std::pow(M1 / (M1 + M2), 2) * (distance / cm) * (atdens) / (1e+24); // in MeV2
}

G4double CalcScreening_TF(G4double Z1, G4double Z2)
{
	G4double a_bohr = 52900; // in fm
	G4double zz1 = std::pow(Z1,2. / 3.);
	G4double zz2 = std::pow(Z2,2. / 3.);
	return a_bohr * 0.8853 * std::pow(zz1 + zz2, -0.5);
}

G4double CalcScreening_ZBL(G4double Z1, G4double Z2)
{
	G4double a_bohr = 52900; // in fm
	G4double zz1 = std::pow(Z1, 0.23);
	G4double zz2 = std::pow(Z2, 0.23);
	return a_bohr * 0.8853 / (zz1 + zz2);
}


G4double CalcDeltaF_K(G4double M1, G4double energy, G4double angle, G4double M2, G4double kinematic_factor)
{
	G4double delta_fk = 0.;
	G4double delta_fk1 = -2 * M1 * kinematic_factor * energy * sin(angle);
	G4double delta_fk2 = sqrt(M2 * M2 - (M1 * M1 * sin(angle) * sin(angle)));
	delta_fk = delta_fk1 / delta_fk2;
	
	return delta_fk;
}

G4double CalcScallingFactorS(G4double a, G4double b, G4double v)
{
	G4double sf = 0;
	G4double k             = b / a;
	G4double k_v           = k * (v + 1);
	G4double one_plus_k    = std::abs(1 + k);
	G4double one_plus_k_sq = std::pow(one_plus_k, v + 1);

	G4double plus_one   = (one_plus_k_sq + 1) / k_v;
	G4double minus_one  = (one_plus_k_sq - 1) / k_v;

	G4double v_plus     = v + 1;
	G4double inv_v      = 1 / v;

	if (b==0)
		sf = a;
	else if (a == 0)
		sf = b / (std::pow(v_plus, inv_v));
	else if (a * b != 0 && k < -1)
		sf = a * std::pow(std::abs(plus_one), inv_v);
	else if (a * b != 0 && k >= -1)
		sf = a * std::pow(std::abs(minus_one), inv_v);
	return sf;
}

G4double CalcScallingFactorMiu(G4double energy, G4double Z1, G4double Z2,G4double screening_radius)
{
	G4double miu = 0.;
	G4double elmcharge_squared = 1.4399764;
	G4double miu_up = (energy / MeV) * screening_radius;
	G4double miu_down = 2 * Z1 * Z2 * elmcharge_squared;
	miu = miu_up / miu_down;

	return miu;
}

// kinematic factor for multiple scattering evaluations
G4double KinematicFactor(G4double angle, G4double M1, G4double M2)
{
	G4double k = 0.;
	G4double square = std::sqrt(std::pow(M2, 2.) - std::pow(M1 * sin(angle), 2.));
	G4double M1co = M1 * cos(angle);
	G4double denominator = M1 + M2;
	k = std::pow((M1co + square) / denominator, 2.);
	return k;
}