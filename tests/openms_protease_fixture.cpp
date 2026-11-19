// Native API qualification: cleavage exceptions, terminal fragments and interval semantics.
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

void check(const char* enzyme, const char* sequence,
           const std::vector<std::string>& expected) {
  OpenMS::ProteaseDigestion digestion;
  digestion.setEnzyme(enzyme);
  digestion.setMissedCleavages(0);
  const auto protein = OpenMS::AASequence::fromString(sequence);
  std::vector<OpenMS::AASequence> peptides;
  std::vector<std::pair<size_t, size_t>> intervals;
  digestion.digest(protein, peptides);
  digestion.digest(protein, intervals);
  assert(peptides.size() == expected.size() && intervals.size() == expected.size());
  size_t cursor = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    assert(peptides[i].toUnmodifiedString() == expected[i]);
    assert(intervals[i].first == cursor);
    assert(intervals[i].second == cursor + expected[i].size());
    cursor = intervals[i].second;
  }
  assert(cursor == std::string(sequence).size());
  digestion.setMissedCleavages(3);
  peptides.clear();
  digestion.digest(protein, peptides);
  intervals.clear();
  digestion.digest(protein, intervals);
  assert(intervals.size() == peptides.size());
  for (size_t i = 0; i < intervals.size(); ++i)
    assert(std::string(sequence).substr(intervals[i].first, intervals[i].second - intervals[i].first)
           == peptides[i].toUnmodifiedString());
  // OpenMS 3.0 can repeat the same interval for N-terminal Asp-N cleavage.
  // Remove duplicate coordinates, never duplicate strings at distinct positions.
  const auto native_count = intervals.size();
  std::sort(intervals.begin(), intervals.end());
  intervals.erase(std::unique(intervals.begin(), intervals.end()), intervals.end());
  size_t count = 0;
  for (size_t span = 1; span <= 4 && span <= expected.size(); ++span)
    count += expected.size() - span + 1;
  if (intervals.size() != count) {
    std::cerr << enzyme << " sequence=" << sequence << " expected=" << count
              << " actual=" << peptides.size() << " products:";
    for (const auto& peptide : peptides) std::cerr << " [" << peptide.toUnmodifiedString() << "]";
    std::cerr << std::endl;
  }
  assert(intervals.size() == count);
  std::cout << enzyme << " zero_missed=" << expected.size()
            << " up_to_three_missed=" << count << " native_rows=" << native_count
            << " intervals=zero_based_half_open\n";
}

int main() {
  check("Trypsin", "AKPRAK", {"AKPR", "AK"});
  check("Trypsin", "AK", {"AK"});
  check("Chymotrypsin", "AFPLWYA", {"AFPL", "W", "Y", "A"});
  check("Asp-N", "DADD", {"DA", "D", "D"});
}
