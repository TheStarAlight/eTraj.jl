using .MolecularCalculators

"Sample provider which yields electron samples through WFAT formula, matching `IonRateMethod=:WFAT`"
struct WFATSampleProvider <: ElectronSampleProvider end
