#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <array>
#include <cmath>
#include <functional>
#include <sstream>
#include <unordered_map>

#include "DistanceMatrix.hpp"
#include "DNA_b128_String.hpp"
#include "Sequences2DistanceMatrix.hpp"
#include "NeighborJoining.hpp"
#include "SequenceTree.hpp"

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
// Branch lengths are -1.0 for NJ and FNJ (not computed by FastPhylo).
// BioNJ computes real branch lengths.
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
            auto it = name_to_id.find(NAME(nd));
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
        edges.append(py::make_tuple(node_id.at(nd), node_id.at(par), EDGE(nd)));
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
// Protein ML distance  (Brent bounded minimizer, no scipy)
// ────────────────────────────────────────────────────────────────────────────

// Map amino acid characters to indices 0-19 (ARNDCQEGHILKMFPSTWYV convention).
// Returns -1 for gaps or unknown characters.
static const int AA_IDX[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 0-15
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 16-31
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 32-47 (space-/)
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, // 48-63 (0-?)
    -1, 0,-1, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1, // 64-79  (A-O)
    14, 5, 1,15,16,-1,19,17,-1,18,-1,-1,-1,-1,-1,-1, // 80-95  (P-_)
    -1, 0,-1, 4, 3, 6,13, 7, 8, 9,-1,11,10,12, 2,-1, // 96-111 (a-o)
    14, 5, 1,15,16,-1,19,17,-1,18,-1,-1,-1,-1,-1,-1, // 112-127 (p-DEL)
    // 128-255: all -1 (default-initialised)
};

// 20×20 amino acid replacement count matrix from two aligned sequences.
static std::array<std::array<double,20>,20>
protein_count_matrix(const std::string &s1, const std::string &s2) {
    std::array<std::array<double,20>,20> N{};
    size_t len = std::min(s1.size(), s2.size());
    for (size_t k = 0; k < len; ++k) {
        int i = AA_IDX[(unsigned char)s1[k]];
        int j = AA_IDX[(unsigned char)s2[k]];
        if (i >= 0 && j >= 0) N[i][j] += 1.0;
    }
    return N;
}

// Negative log-likelihood of observed counts N under P(t).
// right_e[i][k], left_e[k][j], evals[k] define the model eigen decomposition.
// P(t)[i][j] = sum_k right_e[i][k] * exp(evals[k]*t) * left_e[k][j]
static double neg_log_likelihood(
    const double *right_e,  // 20×20 row-major
    const double *left_e,   // 20×20 row-major
    const double *evals,    // 20
    const std::array<std::array<double,20>,20> &N,
    double t
) {
    // exp(evals * t)
    double exp_ev[20];
    for (int k = 0; k < 20; ++k) exp_ev[k] = std::exp(evals[k] * t);

    double ll = 0.0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            double nij = N[i][j];
            if (nij == 0.0) continue;
            // P[i][j] = sum_k right_e[i*20+k] * exp_ev[k] * left_e[k*20+j]
            double p = 0.0;
            const double *ri = right_e + i * 20;
            for (int k = 0; k < 20; ++k)
                p += ri[k] * exp_ev[k] * left_e[k * 20 + j];
            if (p < 1e-300) p = 1e-300;
            ll += nij * std::log(p);
        }
    }
    return -ll;
}

// Brent's bounded minimizer (mirrors scipy.optimize.minimize_scalar method="bounded").
static double brent_minimize(
    std::function<double(double)> f,
    double xa, double xb,
    double xtol = 0.02,
    int maxiter = 20
) {
    const double mintol = 1.0e-11;
    const double cg = 0.3819660112501051;  // (3 - sqrt(5)) / 2

    double fulc = xa + cg * (xb - xa);
    double nfc = fulc, xf = fulc;
    double rat = 0.0, e = 0.0;
    double x = xf, fx = f(x);
    double ffulc = fx, fnfc = fx;
    double xm = 0.5 * (xa + xb);
    double tol1 = xtol * std::abs(x) + mintol;
    double tol2 = 2.0 * tol1;

    for (int iter = 0; iter < maxiter && std::abs(x - xm) > tol2 - 0.5 * (xb - xa); ++iter) {
        bool golden = true;
        double p = 0.0, q = 0.0, r_val = 0.0;

        if (std::abs(e) > tol1) {
            r_val = (x - nfc) * (fx - ffulc);
            q     = (x - fulc) * (fx - fnfc);
            p     = (x - fulc) * q - (x - nfc) * r_val;
            q     = 2.0 * (q - r_val);
            if (q > 0.0) p = -p; else q = -q;
            r_val = e;
            e     = rat;
            if (std::abs(p) < std::abs(0.5 * q * r_val) &&
                p > q * (xa - x) && p < q * (xb - x)) {
                rat = p / q;
                double xnew = x + rat;
                if ((xnew - xa) < tol2 || (xb - xnew) < tol2)
                    rat = (xm > x) ? tol1 : -tol1;
                golden = false;
            }
        }
        if (golden) {
            e   = (x >= xm) ? xa - x : xb - x;
            rat = cg * e;
        }

        double u = (std::abs(rat) >= tol1) ? x + rat : x + (rat > 0 ? tol1 : -tol1);
        double fu = f(u);

        if (fu <= fx) {
            if (u < x) xb = x; else xa = x;
            fulc = nfc; ffulc = fnfc;
            nfc  = x;   fnfc  = fx;
            x    = u;   fx    = fu;
        } else {
            if (u < x) xa = u; else xb = u;
            if (fu <= fnfc || nfc == x) {
                fulc = nfc; ffulc = fnfc;
                nfc  = u;   fnfc  = fu;
            } else if (fu <= ffulc || fulc == x || fulc == nfc) {
                fulc = u; ffulc = fu;
            }
        }
        xm   = 0.5 * (xa + xb);
        tol1 = xtol * std::abs(x) + mintol;
        tol2 = 2.0 * tol1;
    }
    return x;
}

// Compute all pairwise protein ML distances.
// right_e, left_e: (20,20) C-contiguous double arrays.
// evals: (20,) C-contiguous double array.
static StrDblMatrix compute_protein_distances_cpp(
    const std::vector<std::string> &names,
    const std::vector<std::string> &seqs,
    py::array_t<double, py::array::c_style | py::array::forcecast> right_e,
    py::array_t<double, py::array::c_style | py::array::forcecast> left_e,
    py::array_t<double, py::array::c_style | py::array::forcecast> evals
) {
    const double *re  = right_e.data();
    const double *le  = left_e.data();
    const double *ev  = evals.data();
    const double DELTA    = 0.0001;
    const double MAX_DIST = 3.0;
    const double XTOL     = 0.02;
    const int    MAXITER  = 20;

    size_t n = names.size();
    auto dm = make_matrix(names);

    for (size_t i = 0; i < n; ++i) {
        dm.setDistance((int)i, (int)i, 0.0);
        for (size_t j = i + 1; j < n; ++j) {
            auto N = protein_count_matrix(seqs[i], seqs[j]);
            double d = brent_minimize(
                [&](double t){ return neg_log_likelihood(re, le, ev, N, t); },
                DELTA, MAX_DIST, XTOL, MAXITER
            );
            dm.setDistance((int)i, (int)j, d);
            dm.setDistance((int)j, (int)i, d);
        }
    }
    return dm;
}

// ────────────────────────────────────────────────────────────────────────────
// Module
// ────────────────────────────────────────────────────────────────────────────

PYBIND11_MODULE(_fastphylo, m) {
    m.doc() = "fastphylo C++ extension — DNA distances and NJ tree reconstruction";

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
    // NJ and FNJ: branch lengths are -1 (not computed by FastPhylo).
    // BioNJ: branch lengths are computed.
    // ------------------------------------------------------------------
    m.def("nj_tree",
        [](const StrDblMatrix &dm) { return run_nj(dm, NJ); },
        py::arg("dm"),
        "NJ tree. Returns (edges, leaf_names). Branch lengths are -1.");

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
    // compute_protein_distances_cpp(names, seqs, right_e, left_e, evals) -> DistMatrix
    // ------------------------------------------------------------------
    m.def("compute_protein_distances_cpp",
        &compute_protein_distances_cpp,
        py::arg("names"), py::arg("seqs"),
        py::arg("right_e"), py::arg("left_e"), py::arg("evals"),
        "Compute pairwise protein ML distances using a pre-computed eigen decomposition.\n\n"
        "right_e and left_e are (20,20) float64 arrays; evals is a (20,) float64 array.\n"
        "Returns a DistMatrix.");

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
            tree.printOn(ss);
            return ss.str();
        },
        py::arg("dm"), py::arg("method") = "fnj",
        "Return the Newick string for the NJ tree of dm.");
}
