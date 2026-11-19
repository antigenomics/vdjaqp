// Native metadata regression only: shared/overlapping software strings and amounts.
#include <OpenMS/SIMULATION/DigestSimulation.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>
using Input = std::tuple<OpenMS::String, OpenMS::String, double>;
OpenMS::SimTypes::FeatureMapSim digest(const std::vector<Input>& inputs) {
  OpenMS::SimTypes::FeatureMapSim features;
  OpenMS::ProteinIdentification proteins;
  for (const auto& input : inputs) {
    OpenMS::ProteinHit protein;
    protein.setAccession(std::get<0>(input));
    protein.setSequence(std::get<1>(input));
    protein.setMetaValue("intensity", std::get<2>(input));
    proteins.insertHit(protein);
  }
  features.getProteinIdentifications().push_back(proteins);
  OpenMS::DigestSimulation digestion;
  auto params = digestion.getParameters();
  params.setValue("model_naive:missed_cleavages", 0);
  digestion.setParameters(params);
  digestion.digest(features);
  return features;
}
int main(int argc, char** argv) {
  assert(argc == 2);
  const auto identical = digest({{"fixture_a","PEPTIDER",1000},
                                 {"fixture_b","PEPTIDER",1000},
                                 {"fixture_c","PEPTIDER",1000}});
  assert(identical.size() == 1);
  const auto& hit = identical[0].getPeptideIdentifications()[0].getHits()[0];
  assert(hit.getSequence().toString() == "PEPTIDER");
  assert(identical[0].getIntensity() == 3000.0);
  assert(hit.extractProteinAccessionsSet() == std::set<OpenMS::String>({"fixture_a", "fixture_b", "fixture_c"}));
  assert(hit.getPeptideEvidences().size() == std::stoul(argv[1]));
  std::cout << "PEPTIDER intensity=3000 unique_accessions=3 evidence_rows=" << hit.getPeptideEvidences().size() << "\n";
  std::vector<Input> inputs = {{"fixture_a","PEPTIDERAAAAAK",1000},
                               {"fixture_b","PEPTIDER",2000},
                               {"fixture_c","AAAAAK",3000}};
  for (int reverse = 0; reverse < 2; ++reverse) {
    const auto features = digest(inputs);
    assert(features.size() == 2);
    for (size_t i = 0; i < features.size(); ++i) {
      const auto& feature = features[i];
      const auto& peptide = feature.getPeptideIdentifications()[0].getHits()[0];
      const OpenMS::String sequence = i == 0 ? "AAAAAK" : "PEPTIDER";
      const double intensity = i == 0 ? 4000 : 3000;
      const std::set<OpenMS::String> accessions = i == 0
        ? std::set<OpenMS::String>{"fixture_a","fixture_c"}
        : std::set<OpenMS::String>{"fixture_a","fixture_b"};
      assert(peptide.getSequence().toString() == sequence);
      assert(feature.getIntensity() == intensity);
      assert(peptide.extractProteinAccessionsSet() == accessions);
      if (std::stoul(argv[1]) == 3) {
        assert(peptide.getPeptideEvidences().size() == accessions.size());
        size_t index = 0;
        for (const auto& accession : accessions)
          assert(peptide.getPeptideEvidences()[index++].getProteinAccession() == accession);
      }
    }
    std::reverse(inputs.begin(), inputs.end());
  }
  std::cout << "overlap_order_checks=passed peptides=AAAAAK,PEPTIDER intensities=4000,3000 unique_accessions=2,2\n";
}
