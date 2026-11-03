// Gap-fill probe for arceciemat / PR #123 Channel S==0.
// Order: real NIST materials first (stable ElasticData), scans, then synthetic
// missing-isotope elements with SkipMissingIsotopes for stock S==0 proof.

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4HadProjectile.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <cmath>
#include <iostream>
#include <memory>

namespace
{
struct ProbeTrack
{
  std::unique_ptr<G4Step> step;
  std::unique_ptr<G4Track> track;
};

ProbeTrack MakeNeutronTrack(G4Material* mat, G4double eKin)
{
  ProbeTrack pt;
  auto* dp = new G4DynamicParticle(G4Neutron::Definition(), G4ThreeVector(0., 0., 1.), eKin);
  pt.track = std::make_unique<G4Track>(dp, 0., G4ThreeVector(0., 0., 0.));
  pt.step = std::make_unique<G4Step>();
  auto* pre = new G4StepPoint();
  pre->SetMaterial(mat);
  pt.step->SetPreStepPoint(pre);
  pt.step->SetPostStepPoint(new G4StepPoint());
  pt.track->SetStep(pt.step.get());
  return pt;
}

G4double IsotopeSum(G4ParticleHPChannel* ch, const G4HadProjectile& proj, G4Material* mat)
{
  G4ParticleHPThermalBoost boost;
  const G4double T = mat->GetTemperature();
  G4double sum = 0.;
  for (G4int i = 0; i < ch->GetNiso(); ++i) {
    if (!ch->HasAnyData(i)) continue;
    sum += ch->GetWeightedXsec(boost.GetThermalEnergy(proj, ch->GetN(i), ch->GetZ(i), T), i);
  }
  return sum;
}

void ScanMaterial(const char* label, G4Material* mat, G4ParticleHPChannel* ch,
                  G4ParticleHPElasticData* data, G4int nPts, G4double eMin_eV, G4double eMax_eV)
{
  const G4Element* el = mat->GetElement(0);
  std::cout << "\n=== Scan " << label << " T=" << mat->GetTemperature() / kelvin << " K"
            << " niso=" << ch->GetNiso() << " hasData=" << ch->HasDataInAnyFinalState() << " ===\n";

  G4int mech2 = 0, dataModel = 0, scanned = 0;
  G4double bestRatio = 0., bestE = 0., bestXs = 0., bestS = 0., bestData = 0.;
  const G4double logMin = std::log10(eMin_eV);
  const G4double logMax = std::log10(eMax_eV);

  for (G4int k = 0; k <= nPts; ++k) {
    const G4double eKin = std::pow(10., logMin + (logMax - logMin) * k / nPts) * eV;
    ProbeTrack pt = MakeNeutronTrack(mat, eKin);
    G4HadProjectile proj(*pt.track);
    G4ParticleHPThermalBoost boost;

    const G4double EthEl = boost.GetThermalEnergy(proj, el, mat->GetTemperature());
    const G4double xsEl = ch->GetXsec(EthEl);
    const G4double S = IsotopeSum(ch, proj, mat);
    // Doppler Data XS is expensive — subsample.
    G4double xsData = -1.;
    if (k % 10 == 0) {
      xsData = data->GetCrossSection(pt.track->GetDynamicParticle(), el, mat->GetTemperature());
    }
    ++scanned;

    if (xsEl > 0.) {
      const G4double ratio = (S > 0.) ? (xsEl / S) : 1e99;
      if (ratio > bestRatio) {
        bestRatio = ratio;
        bestE = eKin;
        bestXs = xsEl;
        bestS = S;
        bestData = xsData;
      }
    }
    if (xsEl > 1e-30 && S == 0.) {
      ++mech2;
      if (mech2 <= 3)
        std::cout << std::scientific << "  MECH2 Ekin=" << eKin / eV << " xsEl=" << xsEl
                  << " S=" << S << "\n";
    }
    if (xsData > 1e-30 * barn && S == 0.) {
      ++dataModel;
      if (dataModel <= 3)
        std::cout << std::scientific << "  DATA>0_S==0 Ekin=" << eKin / eV
                  << " dataXS=" << xsData / barn << " barn S=" << S << "\n";
    }
  }

  std::cout << "  scanned=" << scanned << " mech2(xsEl>0&&S==0)=" << mech2
            << " data>0&&S==0=" << dataModel << "\n";
  std::cout << std::scientific << "  strongest xsEl/S @ Ekin=" << bestE / eV << " eV xsEl=" << bestXs
            << " S=" << bestS << " ratio=" << bestRatio;
  if (bestData >= 0.) std::cout << " dataXS~" << bestData / barn << " barn";
  std::cout << "\n";
}
}  // namespace

int main()
{
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(314159);
  G4Neutron::Definition();
  G4ParticleTable::GetParticleTable()->SetReadiness();

  auto* nist = G4NistManager::Instance();
  auto* cdBase = nist->FindOrBuildMaterial("G4_Cd");
  auto* uBase = nist->FindOrBuildMaterial("G4_U");
  auto* gdBase = nist->FindOrBuildMaterial("G4_Gd");

  auto* matCd300 = new G4Material("Cd300", cdBase->GetDensity(), cdBase, kStateSolid, 300. * kelvin);
  auto* matCd1200 =
    new G4Material("Cd1200", cdBase->GetDensity(), cdBase, kStateSolid, 1200. * kelvin);
  auto* matU300 = new G4Material("U300", uBase->GetDensity(), uBase, kStateSolid, 300. * kelvin);
  auto* matU3000 = new G4Material("U3000", uBase->GetDensity(), uBase, kStateSolid, 3000. * kelvin);
  auto* matGd1200 =
    new G4Material("Gd1200", gdBase->GetDensity(), gdBase, kStateSolid, 1200. * kelvin);

  std::cout << "Phase1 elements: " << G4Element::GetNumberOfElements() << "\n";

  G4ParticleHPElastic elasticModel;
  elasticModel.BuildPhysicsTable(*G4Neutron::Definition());
  G4ParticleHPElasticData elasticData;
  elasticData.BuildPhysicsTable(*G4Neutron::Definition());

  auto* channels = G4ParticleHPManager::GetInstance()->GetElasticFinalStates();
  if (!channels) {
    std::cerr << "FATAL: no channels\n";
    return 2;
  }

  auto* chCd = (*channels)[matCd300->GetElement(0)->GetIndex()];
  auto* chU = (*channels)[matU300->GetElement(0)->GetIndex()];
  auto* chGd = (*channels)[matGd1200->GetElement(0)->GetIndex()];

  ScanMaterial("Cd", matCd300, chCd, &elasticData, 600, 1e-5, 20.);
  ScanMaterial("Cd", matCd1200, chCd, &elasticData, 600, 1e-5, 20.);
  ScanMaterial("U", matU300, chU, &elasticData, 600, 1e-5, 100.);
  ScanMaterial("U", matU3000, chU, &elasticData, 600, 1e-5, 100.);
  ScanMaterial("Gd", matGd1200, chGd, &elasticData, 600, 1e-5, 20.);

  // ---- Phase 2: synthetic missing isotopes (stable with SkipMissingIsotopes) ----
  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes(true);
  std::cout << "\n=== Phase2: SkipMissingIsotopes=true, synthetic Fe isotopes ===\n";

  auto* isoA = new G4Isotope("Fe99", 26, 99, 99.0 * g / mole);
  auto* isoB = new G4Isotope("Fe101", 26, 101, 101.0 * g / mole);
  auto* elFake = new G4Element("FakeFe", "fFe", 2);
  elFake->AddIsotope(isoA, 60. * perCent);
  elFake->AddIsotope(isoB, 40. * perCent);
  auto* matFake = new G4Material("FakeFeMat", 7.8 * g / cm3, 1, kStateSolid, 293.15 * kelvin);
  matFake->AddElement(elFake, 1.0);

  auto* isoC = new G4Isotope("Fe103", 26, 103, 103.0 * g / mole);
  auto* elFake2 = new G4Element("FakeFe2", "fFe2", 1);
  elFake2->AddIsotope(isoC, 100. * perCent);
  auto* matTwoEmpty =
    new G4Material("TwoEmpty", 7.8 * g / cm3, 2, kStateSolid, 293.15 * kelvin);
  matTwoEmpty->AddElement(elFake, 0.5);
  matTwoEmpty->AddElement(elFake2, 0.5);

  auto* iso56 = new G4Isotope("Fe56p", 26, 56, 55.934937 * g / mole);
  auto* iso99 = new G4Isotope("Fe99p", 26, 99, 99.0 * g / mole);
  auto* elPartial = new G4Element("PartialFe", "pFe", 2);
  elPartial->AddIsotope(iso56, 90. * perCent);
  elPartial->AddIsotope(iso99, 10. * perCent);
  auto* matPartial =
    new G4Material("PartialFeMat", 7.8 * g / cm3, 1, kStateSolid, 293.15 * kelvin);
  matPartial->AddElement(elPartial, 1.0);

  elasticModel.BuildPhysicsTable(*G4Neutron::Definition());
  channels = G4ParticleHPManager::GetInstance()->GetElasticFinalStates();

  auto* chFake = (*channels)[elFake->GetIndex()];
  auto* chPartial = (*channels)[elPartial->GetIndex()];

  std::cout << "\n=== STOCK S==0 (empty-data element; stock isotope-sum logic) ===\n";
  {
    ProbeTrack pt = MakeNeutronTrack(matFake, 0.0253 * eV);
    G4HadProjectile proj(*pt.track);
    const G4double S = IsotopeSum(chFake, proj, matFake);
    std::cout << "  HasData=" << chFake->HasDataInAnyFinalState() << " niso=" << chFake->GetNiso()
              << " S=" << S << "\n";
    if (S == 0.) {
      std::cout << "  STOCK_S0: Channel would enter sum==0 / G4lrint fallback.\n";
      std::cout << "  TRACKING_NOTE: with SkipMissingIsotopes, process ElasticData for this "
                   "element is expected ~0, so full tracking typically never calls the model. "
                   "S==0 is still the Channel branch if ApplyYourself(-1) is entered "
                   "(e.g. element sum==0 multi-element fallback, or a caller that bypasses Data).\n";
    }

    ProbeTrack pt2 = MakeNeutronTrack(matTwoEmpty, 1.0 * eV);
    G4HadProjectile proj2(*pt2.track);
    G4ParticleHPThermalBoost boost;
    const G4double* nav = matTwoEmpty->GetVecNbOfAtomsPerVolume();
    G4double elemSum = 0.;
    for (G4int i = 0; i < (G4int)matTwoEmpty->GetNumberOfElements(); ++i) {
      const G4Element* el = matTwoEmpty->GetElement(i);
      const G4double Eth = boost.GetThermalEnergy(proj2, el, matTwoEmpty->GetTemperature());
      elemSum += (*channels)[el->GetIndex()]->GetXsec(Eth) * nav[i];
    }
    std::cout << "  element_pick_sum(two empty)=" << elemSum
              << " → Elastic `if (sum==0 || …)` still selects an element → Channel S==0.\n";
  }

  std::cout << "\n=== Partial coverage Fe-56 + Fe-99 ===\n";
  {
    for (G4int i = 0; i < chPartial->GetNiso(); ++i) {
      std::cout << "  iso[" << i << "] A=" << chPartial->GetN(i)
                << " HasAnyData=" << chPartial->HasAnyData(i) << "\n";
    }
    ProbeTrack pt = MakeNeutronTrack(matPartial, 0.0253 * eV);
    G4HadProjectile proj(*pt.track);
    G4ParticleHPThermalBoost boost;
    const G4double Eth = boost.GetThermalEnergy(proj, elPartial, matPartial->GetTemperature());
    const G4double xsEl = chPartial->GetXsec(Eth);
    const G4double S = IsotopeSum(chPartial, proj, matPartial);
    std::cout << std::scientific << "  @thermal xsEl=" << xsEl << " S=" << S << "\n";
    std::cout << "  PARTIAL_NOTE: missing iso → x_i=0; remaining iso with data keep S≈xsEl>0. "
                 "Partial file gaps do not alone produce S==0 after a successful positive-XS pick.\n";
  }

  std::cout << "\nDone.\n";
  return 0;
}
