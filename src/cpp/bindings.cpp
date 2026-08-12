#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <sstream>
#include <unordered_map>

#include "fastphylo/core/DistanceMatrix.hpp"
#include "fastphylo/core/SequenceTree.hpp"
#include "fastphylo/dna/DNA_b128_String.hpp"
#include "fastphylo/dna/Sequences2DistanceMatrix.hpp"
#include "fastphylo/dna/NeighborJoining.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ProtDistCalc.hpp"

namespace py = pybind11;

// ────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ────────────────────────────────────────────────────────────────────────────

static DNA_b128_String::base_frequences
aggregate_freqs(const std::vector<DNA_b128_String> &seqs) {
    DNA_b128_String::base_frequences f = seqs[0].getBaseFrequences();
    for (size_t i = 1; i < seqs.size(); ++i) {
        auto fi = seqs[i].getBaseFrequences();
        f.num_As_          += fi.num_As_;
        f.num_Cs_          += fi.num_Cs_;
        f.num_Gs_          += fi.num_Gs_;
        f.num_Ts_          += fi.num_Ts_;
        f.num_unknowns_    += fi.num_unknowns_;
        f.num_ambiguities_ += fi.num_ambiguities_;
    }
    return f;
}

// Default translation model: treat ambiguities as gaps, no ts/tv ratio override.
static sequence_translation_model default_tm() {
    sequence_translation_model tm{};
    tm.no_ambiguities      = true;
    tm.no_ambig_resolve    = true;
    tm.no_transition_probs = true;
    tm.use_base_freqs      = false;
    tm.no_tstvratio        = true;
    tm.tstvratio           = 2.0f;
    tm.pyrtvratio          = 2.0f;
    return tm;
}

// Build DNA_b128_String objects, mapping U→T so RNA works transparently.
static std::vector<DNA_b128_String>
make_b128(const std::vector<std::string> &raw) {
    std::vector<DNA_b128_String> out;
    out.reserve(raw.size());
    for (const auto &s : raw) {
        std::string mapped;
        mapped.reserve(s.size());
        for (char c : s)
            mapped += (c == 'U' || c == 'u') ? 'T' : c;
        out.emplace_back((int)mapped.size(), mapped);
    }
    return out;
}

// Build an N×N StrDblMatrix with row names set.
static StrDblMatrix make_matrix(const std::vector<std::string> &names) {
    StrDblMatrix dm(names.size());
    for (size_t i = 0; i < names.size(); ++i)
        dm.setIdentifier((int)i, names[i]);
    return dm;
}

// ────────────────────────────────────────────────────────────────────────────
// Tree extraction
// Returns (edges, leaf_names) where:
//   edges      = list[(u:int, v:int, branch_len:float)]
//   leaf_names = list[str]  — leaf_names[i] = accession of leaf with vertex i
//
// Leaf vertex IDs 0…N-1 match the original distance-matrix row order.
// Internal node IDs start at N.
// Branch lengths are -1.0 for FNJ (not computed by FastPhylo). NJ and
// BioNJ compute real branch lengths.
// ────────────────────────────────────────────────────────────────────────────

static py::tuple
extract_tree(SequenceTree &tree, const std::vector<std::string> &orig_names) {
    size_t n = orig_names.size();

    // name → leaf-ID map (preserves original matrix order)
    std::unordered_map<std::string, int> name_to_id;
    name_to_id.reserve(n);
    for (size_t i = 0; i < n; ++i)
        name_to_id[orig_names[i]] = (int)i;

    // Collect all nodes in prefix order
    SequenceTree::NodeVector nodes;
    tree.addNodesInPrefixOrder(nodes);

    // Assign integer vertex IDs
    std::unordered_map<const SequenceTree::Node *, int> node_id;
    node_id.reserve(nodes.size());
    int next_internal = (int)n;
    for (auto *nd : nodes) {
        if (nd->isLeaf()) {
            auto it = name_to_id.find(nodeName(nd));
            node_id[nd] = (it != name_to_id.end()) ? it->second : next_internal++;
        } else {
            node_id[nd] = next_internal++;
        }
    }

    // Build edge list (one edge per non-root node)
    py::list edges;
    for (auto *nd : nodes) {
        const auto *par = nd->getParent();
        if (!par) continue;
        edges.append(py::make_tuple(node_id.at(nd), node_id.at(par), nodeEdge(nd)));
    }

    return py::make_tuple(edges, orig_names);
}

static py::tuple run_nj(StrDblMatrix dm, NJ_method method) {
    size_t n = dm.getSize();
    std::vector<std::string> names(n);
    for (size_t i = 0; i < n; ++i)
        names[i] = dm.getIdentifier((int)i);

    SequenceTree tree;
    computeNJTree(dm, tree, method);
    return extract_tree(tree, names);
}

// ────────────────────────────────────────────────────────────────────────────
// Protein ML distance
//
// fastphylo_py_integration_plan.md, Phase 3: every model fastphylo-py
// supports is now native in FastPhylo's C++ (ModelMatrix.cpp, Phase 0),
// so this routes through the same build_ml_decomposition()/
// calculate_ml_dists() pair fastprot/main.cpp itself uses - FastPhylo's
// safeguarded Newton-Raphson solver, not the from-scratch Brent's-method
// optimizer this replaces (deleted here; kept, differently, in
// protein.py's pure-Python fallback path - see Phase 4).
//
// Model-name strings match protein.py's RateMatrix subclass __name__s
// exactly (the names distance_matrix()'s model= argument already
// passes through unchanged) - not FastPhylo's own fastprot -D flag
// spelling (e.g. "JTT-DCMUT"), which differs in a few places
// (underscore vs hyphen, case). Keeping fastphylo-py's own spelling
// here is what makes this a zero-Python-API-change swap.
// ────────────────────────────────────────────────────────────────────────────

static const std::unordered_map<std::string, model_type> PROTEIN_MODEL_NAMES = {
    {"WAG", wag},       {"JTT", jtt},         {"JTT_DCMut", jtt_dcmut},
    {"LG", lg},         {"VT", vt},           {"HIVb", hivb},
    {"HIVw", hivw},     {"cpREV", cprev},     {"BLOSUM62", blosum62},
    {"Dayhoff", day},   {"DCMUT", dcmut},     {"MtREV", mtrev},
    {"RtREV", rtrev},   {"PMB", pmb},
};

static StrDblMatrix compute_protein_distances_cpp(
    const std::vector<std::string> &names,
    const std::vector<std::string> &seqs,
    const std::string &model_name
) {
    if (names.size() != seqs.size())
        throw std::invalid_argument("names and seqs must have the same length");
    auto it = PROTEIN_MODEL_NAMES.find(model_name);
    if (it == PROTEIN_MODEL_NAMES.end())
        // Wording matches RateMatrix.instantiate()'s own ValueError
        // (protein.py) - the compiled path no longer routes through
        // that Python lookup at all (every model is native in C++ now,
        // Phase 0), but callers/tests should see one consistent error
        // message regardless of which path actually raised it.
        throw std::invalid_argument("Unknown protein model '" + model_name + "'");

    SeqVec sv;
    sv.reserve(names.size());
    for (size_t i = 0; i < names.size(); ++i)
        sv.emplace_back(names[i], seqs[i]);

    MatrixExpm decomp = build_ml_decomposition(it->second);
    StrDblMatrix dm;
    calculate_ml_dists(sv, dm, decomp);
    dm.setIdentifiers(names);
    return dm;
}

// ────────────────────────────────────────────────────────────────────────────
// Module
// ────────────────────────────────────────────────────────────────────────────

PYBIND11_MODULE(_fastphylo, m) {
    m.doc() = "fastphylo C++ extension — DNA/protein distances and NJ tree reconstruction, backed by libfastphylo";

    // ------------------------------------------------------------------
    // DistMatrix — thin Python wrapper around StrDblMatrix
    // ------------------------------------------------------------------
    py::class_<StrDblMatrix>(m, "DistMatrix")
        .def(py::init<size_t>(), py::arg("n"),
             "Create an N×N distance matrix initialised to zero.")
        .def("size",     &StrDblMatrix::getSize,
             "Number of taxa (rows/columns).")
        .def("get", [](const StrDblMatrix &dm, size_t i, size_t j) {
                return dm.getDistance((int)i, (int)j);
            }, py::arg("i"), py::arg("j"))
        .def("set", [](StrDblMatrix &dm, size_t i, size_t j, double v) {
                dm.setDistance((int)i, (int)j, v);
            }, py::arg("i"), py::arg("j"), py::arg("v"))
        .def("name", [](const StrDblMatrix &dm, size_t i) -> std::string {
                return dm.getIdentifier((int)i);
            }, py::arg("i"))
        .def("set_name", [](StrDblMatrix &dm, size_t i, const std::string &s) {
                dm.setIdentifier((int)i, s);
            }, py::arg("i"), py::arg("name"))
        .def("names", [](const StrDblMatrix &dm) {
                std::vector<std::string> out;
                out.reserve(dm.getSize());
                for (size_t i = 0; i < dm.getSize(); ++i)
                    out.push_back(dm.getIdentifier((int)i));
                return out;
            }, "Return all row/column names as a list.");

    // ------------------------------------------------------------------
    // DNA distance computation
    // compute_dna_distances(names, seqs, model="k2p") -> DistMatrix
    // ------------------------------------------------------------------
    m.def("compute_dna_distances",
        [](const std::vector<std::string> &names,
           const std::vector<std::string> &seqs,
           const std::string &model) -> StrDblMatrix
        {
            if (names.size() != seqs.size())
                throw std::invalid_argument(
                    "names and seqs must have the same length");
            if (names.empty())
                throw std::invalid_argument("sequence list is empty");

            auto b128 = make_b128(seqs);
            auto dm   = make_matrix(names);
            auto tm   = default_tm();

            if (model == "hamming") {
                fillMatrix_Hamming(dm, b128, tm);
            } else if (model == "jc") {
                fillMatrix_JC(dm, b128, tm);
            } else if (model == "k2p") {
                fillMatrix_K2P(dm, b128, tm);
            } else if (model == "tn93") {
                fillMatrix_TN93(dm, b128, aggregate_freqs(b128), tm);
            } else {
                throw std::invalid_argument(
                    "unknown model '" + model +
                    "'; expected: hamming, jc, k2p, tn93");
            }
            return dm;
        },
        py::arg("names"), py::arg("seqs"), py::arg("model") = "k2p",
        "Compute pairwise DNA distances.\n\n"
        "Returns a DistMatrix. model: 'hamming' | 'jc' | 'k2p' | 'tn93'.");

    // ------------------------------------------------------------------
    // Tree reconstruction
    // Each function takes a DistMatrix and returns
    //   (edges: list[(u,v,w)], leaf_names: list[str])
    // The input matrix is copied; the original is not modified.
    // NJ and BioNJ: branch lengths are computed.
    // FNJ: branch lengths are -1 (not computed by FastPhylo).
    // ------------------------------------------------------------------
    m.def("nj_tree",
        [](const StrDblMatrix &dm) { return run_nj(dm, NJ); },
        py::arg("dm"),
        "NJ tree. Returns (edges, leaf_names) with real branch lengths.");

    m.def("fnj_tree",
        [](const StrDblMatrix &dm) { return run_nj(dm, FNJ); },
        py::arg("dm"),
        "Fast NJ tree. Returns (edges, leaf_names). Branch lengths are -1.");

    m.def("bionj_tree",
        [](const StrDblMatrix &dm) { return run_nj(dm, BIONJ); },
        py::arg("dm"),
        "BioNJ tree. Returns (edges, leaf_names) with real branch lengths.");

    // ------------------------------------------------------------------
    // Protein ML distance computation
    // compute_protein_distances_cpp(names, seqs, model) -> DistMatrix
    // ------------------------------------------------------------------
    m.def("compute_protein_distances_cpp",
        &compute_protein_distances_cpp,
        py::arg("names"), py::arg("seqs"), py::arg("model"),
        "Compute pairwise protein ML distances via libfastphylo's safeguarded "
        "Newton-Raphson solver.\n\n"
        "model is one of the RateMatrix subclass names in protein.py "
        "(e.g. 'WAG', 'JTT_DCMut'). Returns a DistMatrix.");

    // ------------------------------------------------------------------
    // Newick string convenience function
    // ------------------------------------------------------------------
    m.def("newick",
        [](const StrDblMatrix &dm, const std::string &method) -> std::string {
            StrDblMatrix copy = dm;
            NJ_method m;
            if      (method == "nj")    m = NJ;
            else if (method == "fnj")   m = FNJ;
            else if (method == "bionj") m = BIONJ;
            else throw std::invalid_argument(
                "unknown method '" + method + "'; expected: nj, fnj, bionj");
            SequenceTree tree;
            computeNJTree(copy, tree, m);
            std::ostringstream ss;
            ss << tree;
            return ss.str();
        },
        py::arg("dm"), py::arg("method") = "fnj",
        "Return the Newick string for the NJ tree of dm.");
}
