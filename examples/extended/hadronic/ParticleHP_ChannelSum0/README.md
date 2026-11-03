# ParticleHP Channel `sum == 0` probe

Standalone repro for the question whether `G4ParticleHPChannel::ApplyYourself`
(with `anIsotope == -1`) ever takes the isotope-sum `sum == 0` fallback
(the `G4lrint` path discussed on Geant4 PR #123).

Requires Geant4 with NeutronHP / G4NDL (tested: 11.5.0 + G4NDL 4.7.1).

## Build

```bash
source $G4INSTALL/bin/geant4.sh   # or your env that sets G4NEUTRONHPDATA
cmake -S . -B build -DGeant4_DIR=$G4INSTALL/lib/cmake/Geant4
cmake --build build
```

## Programs

### `ChannelSum0` (Data XS > 0 while model `S == 0`)

Demonstrates process-table elastic XS can be positive for an element whose
Channel final-state / isotope data are empty (`S == 0`), then calls
`G4ParticleHPElastic::ApplyYourself` as after process selection.

Expected markers on stdout:

- `OPTIONAL_A_HIT: Data XS > 0 while model S == 0`
- `OPTIONAL_B: Data>0 ⇒ elastic model ApplyYourself …`

Implementation note: Names already announces nearby CrossSection substitution
(e.g. Fe-99 → Ge) while Channel `HasData` stays false; G4’s own Data
`BuildPhysicsTable` for that far-off Z/A was unsafe to read in our runs.
This probe installs the neighbor (Ge) CrossSection vector on the empty-FS
element slot — the asymmetry Names documents — then invokes the model.

### `ChannelSum0_scan` (stock `S == 0` + Cd/U/Gd scans)

- Stock empty-FS element / two-element `sum == 0` element pick → Channel `S == 0`
- Partial isotope coverage (Fe-56 + missing Fe-99): `xsEl == S > 0` at thermal
- Scans Cd/U/Gd for `xsEl > 0 && S == 0` and ElasticData>0 && S==0

Expected: `STOCK_S0` lines; scan `mech2` / `data>0&&S==0` counts typically 0
on those natural elements (thermal-A can still make `xsEl/S` large).

## Relation to PR #123

This does not change `G4ParticleHPChannel`. It only shows the `sum == 0`
branch is reachable, so making the fallback index safe remains warranted.
