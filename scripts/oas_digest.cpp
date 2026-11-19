// Purpose: aggregate observed OAS proteins with the installed OpenMS digestion API.
// Input is a validated TSV prepared by digest_oas_cohort.py; no protease rules are duplicated.
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

using Key = std::array<size_t, 4>;
using Counts = std::array<uint64_t, 4>; // available internal bonds, then three proteases
using Interval = std::pair<size_t, size_t>;
const std::array<std::string, 3> enzymes{{"Trypsin", "Chymotrypsin", "Asp-N"}};
const std::array<std::string, 7> regions{{"fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"}};
struct Record {
  size_t cohort, gene;
  std::string protein;
  std::array<Interval, 7> bounds;
};
struct Summary {
  std::map<Key, Counts> native, aligned;
  std::map<Key, uint64_t> lengths; // cohort, enzyme, exact missed cleavages, length
  std::map<size_t, uint64_t> records;
};
size_t integer(const std::string& s) {
  if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
    throw std::runtime_error("Expected nonnegative integer");
  return std::stoull(s);
}
Record parse(const std::string& line) {
  std::istringstream in(line); std::vector<std::string> fields; std::string s;
  while (std::getline(in, s, '\t')) fields.push_back(s);
  if (fields.size() != 17) throw std::runtime_error("Expected 17 input columns");
  Record r{integer(fields[0]), integer(fields[1]), fields[2], {}};
  if (r.protein.empty() || r.protein.find_first_not_of("ACDEFGHIKLMNPQRSTVWY") != std::string::npos)
    throw std::runtime_error("Expected canonical observed protein");
  size_t previous_end = 0;
  for (size_t i = 0; i < 7; ++i) {
    const size_t start = integer(fields[3+2*i]), end = integer(fields[4+2*i]);
    // 0,0 denotes unresolved; all other coordinates are zero-based half-open.
    if (end && (start >= end || end > r.protein.size() || start < previous_end))
      throw std::runtime_error("Invalid or overlapping inferred region");
    if (!end && start) throw std::runtime_error("Incomplete inferred region");
    r.bounds[i] = {start, end};
    if (end) previous_end = end;
  }
  return r;
}
std::vector<Interval> digest(OpenMS::ProteaseDigestion& enzyme, const OpenMS::AASequence& protein, size_t missed) {
  enzyme.setMissedCleavages(missed);
  std::vector<Interval> spans;
  enzyme.digest(protein, spans);
  std::sort(spans.begin(), spans.end());
  // OpenMS 3.0 repeats intervals for N-terminal Asp-N; positional identity is authoritative.
  spans.erase(std::unique(spans.begin(), spans.end()), spans.end());
  return spans;
}
void analyze(const Record& r, std::array<OpenMS::ProteaseDigestion, 3>& proteases, Summary& out) {
  ++out.records[r.cohort];
  const size_t n = r.protein.size();
  // Native display: 10-residue protein-length bins, 5-residue internal-bond bins.
  // Region display: gene-call rows, region-local bond offsets; no stretching or imputation.
  std::vector<int> region(n, -1);
  for (size_t k = 0; k < 7; ++k)
    for (size_t i = r.bounds[k].first; i < r.bounds[k].second; ++i) region[i] = int(k);
  for (size_t c = 1; c < n; ++c) {
    ++out.native[{r.cohort, n/10, 0, (c-1)/5}][0];
    if (region[c-1] >= 0) {
      size_t k = size_t(region[c-1]);
      ++out.aligned[{r.cohort, r.gene, k, c-r.bounds[k].first}][0];
    }
  }
  const auto protein = OpenMS::AASequence::fromString(r.protein);
  for (size_t e = 0; e < 3; ++e) {
    const auto zero = digest(proteases[e], protein, 0);
    std::vector<size_t> cuts; size_t cursor = 0;
    for (auto span : zero) {
      if (span.first != cursor || span.second <= cursor || span.second > n)
        throw std::runtime_error("Native zero-missed fragments do not partition observed protein");
      cursor = span.second;
      if (cursor < n) cuts.push_back(cursor);
    }
    if (cursor != n) throw std::runtime_error("Native digest omitted terminal sequence");
    for (size_t c : cuts) {
      ++out.native[{r.cohort, n/10, 0, (c-1)/5}][e+1];
      if (region[c-1] >= 0) {
        size_t k = size_t(region[c-1]);
        ++out.aligned[{r.cohort, r.gene, k, c-r.bounds[k].first}][e+1];
      }
    }
    const auto spans = digest(proteases[e], protein, 3);
    size_t expected = 0;
    for (size_t width = 1; width <= 4 && width <= zero.size(); ++width) expected += zero.size()-width+1;
    if (spans.size() != expected) throw std::runtime_error("Native missed-cleavage interval count differs");
    for (auto span : spans) {
      if (span.first >= span.second || span.second > n) throw std::runtime_error("Invalid native interval");
      const auto lo = std::upper_bound(cuts.begin(), cuts.end(), span.first);
      const auto hi = std::lower_bound(cuts.begin(), cuts.end(), span.second);
      const size_t missed = size_t(hi-lo);
      if (missed > 3) throw std::runtime_error("Native missed-cleavage maximum exceeded");
      ++out.lengths[{r.cohort, e, missed, span.second-span.first}];
    }
  }
}
void merge(Summary& target, const Summary& source) {
  for (const auto& x : source.records) target.records[x.first] += x.second;
  for (const auto& x : source.lengths) target.lengths[x.first] += x.second;
  for (auto pair : {std::make_pair(&target.native, &source.native), std::make_pair(&target.aligned, &source.aligned)})
    for (const auto& x : *pair.second)
      for (size_t i = 0; i < 4; ++i) (*pair.first)[x.first][i] += x.second[i];
}
int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("Usage: oas_digest INPUT_TSV NEW_OUTPUT_PREFIX CPU_COUNT");
    const size_t requested = integer(argv[3]);
    if (!requested) throw std::runtime_error("CPU count must be positive");
    const size_t threads = std::min(requested, size_t(omp_get_num_procs()));
    omp_set_dynamic(0); omp_set_num_threads(int(threads)); omp_set_max_active_levels(1);
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot read input");
    const std::string prefix = argv[2];
    std::array<std::ofstream, 4> files;
    const std::array<std::string, 4> suffix{{"native.tsv", "aligned.tsv", "lengths.tsv", "records.tsv"}};
    for (size_t i = 0; i < 4; ++i) {
      const auto name = prefix+suffix[i];
      if (std::ifstream(name)) throw std::runtime_error("Output already exists");
      files[i].open(name);
      if (!files[i]) throw std::runtime_error("Cannot create output");
    }
    std::vector<Summary> summaries(threads);
    std::vector<std::array<OpenMS::ProteaseDigestion, 3>> proteases(threads);
    for (auto& p : proteases) for (size_t e = 0; e < 3; ++e) p[e].setEnzyme(enzymes[e]);
    std::string line; size_t records = 0;
    // At most 4 MiB of source text per batch; each parallel worker receives a contiguous range.
    while (input.peek() != EOF) {
      std::vector<Record> batch; size_t bytes = 0;
      while (bytes < 4*1024*1024 && std::getline(input, line)) {
        bytes += line.size(); batch.push_back(parse(line));
      }
      std::exception_ptr error;
      #pragma omp parallel for schedule(static)
      for (size_t i = 0; i < batch.size(); ++i) {
        try { auto t = size_t(omp_get_thread_num()); analyze(batch[i], proteases[t], summaries[t]); }
        catch (...) {
          #pragma omp critical
          { if (!error) error = std::current_exception(); }
        }
      }
      if (error) std::rethrow_exception(error);
      records += batch.size();
    }
    if (!records) throw std::runtime_error("Empty input");
    Summary result;
    for (const auto& summary : summaries) merge(result, summary);
    files[0] << "cohort_id\tlength_bin10\tbond_bin5\tavailable_bonds\tTrypsin\tChymotrypsin\tAsp-N\n";
    files[1] << "cohort_id\tgene_id\tregion\tbond_offset\tavailable_bonds\tTrypsin\tChymotrypsin\tAsp-N\n";
    for (size_t type = 0; type < 2; ++type) for (const auto& x : type ? result.aligned : result.native) {
      const auto& k = x.first;
      files[type] << k[0] << '\t' << k[1] << '\t';
      if (type) files[type] << regions[k[2]] << '\t';
      files[type] << k[3];
      for (auto c : x.second) files[type] << '\t' << c;
      files[type] << '\n';
    }
    files[2] << "cohort_id\tenzyme\tmissed_cleavages\tlength_aa\tpeptide_occurrences\n";
    for (const auto& x : result.lengths) {
      const auto& k = x.first;
      files[2] << k[0] << '\t' << enzymes[k[1]] << '\t' << k[2] << '\t' << k[3] << '\t' << x.second << '\n';
    }
    files[3] << "cohort_id\tprotein_records\n";
    for (const auto& x : result.records) files[3] << x.first << '\t' << x.second << '\n';
    for (auto& f : files) { f.flush(); if (!f) throw std::runtime_error("Output write failed"); }
    std::cout << "records=" << records << " threads=" << threads << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
