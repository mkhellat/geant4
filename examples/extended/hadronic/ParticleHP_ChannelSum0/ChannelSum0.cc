// ParticleHP Channel sum==0 probe (see README.md).
//  (A) process ElasticData XS > 0 while Channel isotope S == 0
//  (B) Data>0 ⇒ model ApplyYourself enters the isotope sum==0 branch
//
// Observed in G4 11.5.0 / G4NDL4.7.1:
//  - Fe-99: Channel HasData=0; Names prints Elastic CrossSection "will use" Ge;
//    but the built PhysicsVector for that element was unsafe to read (SEGV).
//  - Fe-61: CrossSection substitutes Fe-60 *and* Channel HasData=1 → S>0 (no A).
//
// Realises the Names-documented asymmetry safely:
//  1) Fake Fe-99 channel with SkipMissingIsotopes → S==0 (no FS).
//  2) Real Ge elastic Data PhysicsVector installed at FakeFe's element index.
//  3) ElasticData::GetCrossSection > 0 while IsotopeSum == 0 → OPTIONAL_A_HIT.
//  4) model.ApplyYourself → Channel sum==0 branch → OPTIONAL_B.

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4GenericIon.hh"
#include "G4HadProjectile.hh"
#include "G4IonTable.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4Nucleus.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include <iostream>
#include <memory>
#include <vector>

namespace
{
struct ProbeTrack
{
  std::unique_ptr<G4Step> step;
  std::unique_ptr<G4Track> track;
};

ProbeTrack MakeTrack(G4Material* mat, G4double eKin)
{
  ProbeTrack pt;
  auto* dp = new G4DynamicParticle(G4Neutron::Definition(), G4ThreeVector(0, 0, 1), eKin);
  pt.track = std::make_unique<G4Track>(dp, 0., G4ThreeVector());
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
  G4double sum = 0.;
  for (G4int i = 0; i < ch->GetNiso(); ++i) {
    if (!ch->HasAnyData(i)) continue;
    sum += ch->GetWeightedXsec(
      boost.GetThermalEnergy(proj, ch->GetN(i), ch->GetZ(i), mat->GetTemperature()), i);
  }
  return sum;
}

// Deep-copy a G4PhysicsVector into a new free vector (owned by caller via table).
G4PhysicsVector* ClonePhysicsVector(const G4PhysicsVector* src)
{
  if (src == nullptr || src->GetVectorLength() == 0) return nullptr;
  const auto n = (G4int)src->GetVectorLength();
  std::vector<G4double> e(n), x(n);
  for (G4int i = 0; i < n; ++i) {
    e[i] = src->Energy(i);
    x[i] = (*src)[i];
  }
  return new G4PhysicsFreeVector(e, x);
}
}  // namespace

int main()
{
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(42);

  G4Neutron::Definition();
  G4GenericIon::Definition();
  G4ParticleTable::GetParticleTable()->SetReadiness();
  G4IonTable::GetIonTable()->CreateAllIon();
  G4IonTable::GetIonTable()->CreateAllIsomer();

  // Phase 1: real Ge only — stable Data + Model tables.
  auto* geMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");
  std::cout << "Phase1: Ge element index=" << geMat->GetElement(0)->GetIndex() << "\n";

  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes(true);
  G4ParticleHPManager::GetInstance()->SetNeglectDoppler(true);
  G4ParticleHPManager::GetInstance()->SetVerboseLevel(0);

  G4ParticleHPElastic model;
  model.BuildPhysicsTable(*G4Neutron::Definition());
  G4ParticleHPElasticData data;
  data.BuildPhysicsTable(*G4Neutron::Definition());

  auto* table = G4ParticleHPManager::GetInstance()->GetElasticCrossSections();
  auto* geVec = (*table)(geMat->GetElement(0)->GetIndex());
  std::cout << "Ge Data vector length=" << (geVec ? geVec->GetVectorLength() : 0) << "\n";
  G4bool oor = false;
  const G4double geXs = geVec ? geVec->GetValue(0.0253 * eV, oor) : 0.;
  std::cout << std::scientific << "Ge table XS @ 0.0253 eV = " << geXs / barn << " barn\n";

  // Phase 2: add FakeFe-99 with SkipMissing → Channel S==0; Data slot initially empty/zero.
  auto* iso = new G4Isotope("Fe99", 26, 99, 99.0 * g / mole);
  auto* elFake = new G4Element("FakeFe", "fFe", 1);
  elFake->AddIsotope(iso, 100. * perCent);
  auto* matFake = new G4Material("FakeFeMat", 7.8 * g / cm3, 1, kStateSolid, 293.15 * kelvin);
  matFake->AddElement(elFake, 1.0);

  model.BuildPhysicsTable(*G4Neutron::Definition());  // extend channels for new element
  // Rebuild Data so FakeFe gets a (zero) slot, then overwrite with Ge clone.
  data.BuildPhysicsTable(*G4Neutron::Definition());
  table = G4ParticleHPManager::GetInstance()->GetElasticCrossSections();

  auto* channels = G4ParticleHPManager::GetInstance()->GetElasticFinalStates();
  auto* chFake = (*channels)[elFake->GetIndex()];
  std::cout << "FakeFe index=" << elFake->GetIndex()
            << " HasData=" << chFake->HasDataInAnyFinalState()
            << " HasAnyData[0]=" << chFake->HasAnyData(0) << "\n";

  // Install Ge CrossSection at FakeFe index — safe stand-in for Names "will use Germanium".
  G4PhysicsVector* geClone = ClonePhysicsVector((*table)(geMat->GetElement(0)->GetIndex()));
  if (geClone == nullptr) {
    std::cerr << "FATAL: could not clone Ge vector\n";
    return 2;
  }
  // Replace table entry (G4PhysicsTable owns vectors).
  G4PhysicsVector* old = (*table)(elFake->GetIndex());
  (*table)(elFake->GetIndex()) = geClone;
  delete old;
  G4ParticleHPManager::GetInstance()->RegisterElasticCrossSections(table);

  const G4double eKin = 0.0253 * eV;
  ProbeTrack pt = MakeTrack(matFake, eKin);
  G4HadProjectile proj(*pt.track);

  const G4double xsData =
    data.GetCrossSection(pt.track->GetDynamicParticle(), elFake, matFake->GetTemperature());
  const G4double S = IsotopeSum(chFake, proj, matFake);
  G4ParticleHPThermalBoost boost;
  const G4double xsEl =
    chFake->GetXsec(boost.GetThermalEnergy(proj, elFake, matFake->GetTemperature()));

  std::cout << std::scientific << "Ekin=" << eKin / eV << " eV\n"
            << "  ElasticData XS (FakeFe←Ge) = " << xsData / barn << " barn\n"
            << "  Channel xs_element = " << xsEl << "\n"
            << "  Channel isotope S  = " << S << "\n";

  if (xsData > 0. && S == 0.) {
    std::cout << "OPTIONAL_A_HIT: Data XS > 0 while model S == 0\n";
  }
  else {
    std::cout << "OPTIONAL_A_MISS: Data=" << xsData / barn << " S=" << S << "\n";
    return 3;
  }

  std::cout << "OPTIONAL_B: Data>0 ⇒ elastic model ApplyYourself (post process-selection)\n";
  G4Nucleus nuc;
  G4HadFinalState* fs = model.ApplyYourself(proj, nuc);
  std::cout << "  ApplyYourself returned " << (fs ? "non-null" : "null") << "\n";
  std::cout << "Done.\n";
  return 0;
}
