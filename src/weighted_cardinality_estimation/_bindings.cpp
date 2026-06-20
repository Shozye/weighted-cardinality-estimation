#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "memory_flag.hpp"
#include "exp_sketch.hpp"
#include "exp_sketch_float32.hpp"
#include "fast_exp_sketch.hpp"
#include "fast_exp_sketch_t.hpp"
#include "fastgm_exp_sketch.hpp"
#include "q_sketch_dyn.hpp"
#include "q_sketch.hpp"
#include "fast_k_q_sketch.hpp"
#include "fast_k_q_sketch_rounding.hpp"
#include "k_q_sketch_rounded_dyn.hpp"
#include "k_q_sketch_shifted.hpp"
#include "weighted_min_hash.hpp"
#include "fast_exp_sketch_custom_float.hpp"
#include "log_exp_sketch_slow_no_shifted.hpp"
#include "log_exp_sketch_slow_shifted.hpp"
#include "log_exp_sketch_fast_no_shifted.hpp"
#include "log_exp_sketch_fast_shifted.hpp"
#include "min_hash.hpp"
#include "martingale_min_hash.hpp"
#include "hyper_log_log.hpp"
#include "weighted_hyper_log_log.hpp"
#include "weighted_hyper_log_log_custom_float.hpp"
#include "rng_engine_type.hpp"

namespace py = pybind11;

// ─── Method binders ──────────────────────────────────────────────────────────

// Unweighted sketches: only get_registers (add/add_many/estimate/memory come from CardinalitySketch)
template <typename PyClass>
PyClass& bind_unweighted_base(PyClass& cls) {
    return cls.def("get_registers", &PyClass::type::get_registers);
}

// Weighted sketches: add(x, weight) + add_many(elems, weights) overloads, plus get_registers
template <typename PyClass>
PyClass& bind_sketch_base(PyClass& cls) {
    using Cls = typename PyClass::type;
    return cls
        .def("add",      static_cast<void (Cls::*)(const std::string&, double)>(&Cls::add),
             py::arg("x"), py::arg("weight"))
        .def("add",      static_cast<void (Sketch::*)(const std::string&)>(&Sketch::add),
             py::arg("x"))
        .def("add_many", static_cast<void (Cls::*)(const std::vector<std::string>&, const std::vector<double>&)>(&Cls::add_many),
             py::arg("elems"), py::arg("weights"))
        .def("add_many", static_cast<void (Sketch::*)(const std::vector<std::string>&)>(&Sketch::add_many),
             py::arg("elems"))
        .def("get_registers", &Cls::get_registers);
}

// ─── Pickle helpers ──────────────────────────────────────────────────────────

// Shape: (m, master_seed, registers)
template <typename Cls, typename RegT, typename... Bases>
void bind_pickle_regs(py::class_<Cls, Bases...>& cls) {
    cls.def(py::pickle(
        [](const Cls& p) {
            return py::make_tuple(p.get_sketch_size(), p.get_master_seed(), p.get_registers());
        },
        [](const py::tuple& t) {
            if (t.size() != 3) throw std::runtime_error("Invalid pickle state!");
            return Cls(t[0].cast<std::size_t>(),
                       t[1].cast<std::uint64_t>(),
                       t[2].cast<std::vector<RegT>>());
        }
    ));
}

// Shape: (m, master_seed, amount_bits, registers)
template <typename Cls>
void bind_pickle_q(py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin>& cls) {
    cls.def(py::pickle(
        [](const Cls& p) {
            return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                  p.get_amount_bits(), p.get_registers());
        },
        [](const py::tuple& t) {
            if (t.size() != 4) throw std::runtime_error("Invalid pickle state!");
            return Cls(t[0].cast<std::size_t>(),
                       t[1].cast<std::uint64_t>(),
                       t[2].cast<std::uint8_t>(),
                       t[3].cast<std::vector<int>>());
        }
    ));
}

// Shape: (m, master_seed, amount_bits, registers, log_base) — note: ctor takes log_base before registers
template <typename Cls, typename... Bases>
void bind_pickle_log_exp(py::class_<Cls, Bases...>& cls) {
    cls.def(py::pickle(
        [](const Cls& p) {
            return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                  p.get_amount_bits(), p.get_registers(), p.get_logarithm_base());
        },
        [](const py::tuple& t) {
            if (t.size() != 5) throw std::runtime_error("Invalid pickle state!");
            return Cls(t[0].cast<std::size_t>(),
                       t[1].cast<std::uint64_t>(),
                       t[2].cast<std::uint8_t>(),
                       t[4].cast<float>(),
                       t[3].cast<std::vector<int>>());
        }
    ));
}

// ─── Family binders ──────────────────────────────────────────────────────────

// ExpSketch family: (m, seed) ctor, Jaccard + merge + regs pickle
template <typename Cls, typename RegT>
void bind_jaccard_sketch(py::module_& m, const char* name) {
    auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, name)
        .def(py::init<std::size_t, std::uint64_t>(), py::arg("m"), py::arg("seed"));
    bind_sketch_base(cls)
        .def("jaccard_struct", &Cls::jaccard_struct)
        .def("merge", &Cls::merge, py::arg("other"));
    bind_pickle_regs<Cls, RegT>(cls);
}

// WeightedMinHash: (m, seed) ctor, merge only, regs pickle
template <typename Cls, typename RegT>
void bind_mergeable_sketch(py::module_& m, const char* name) {
    auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin>(m, name)
        .def(py::init<std::size_t, std::uint64_t>(), py::arg("m"), py::arg("seed"));
    bind_sketch_base(cls)
        .def("merge", &Cls::merge, py::arg("other"));
    bind_pickle_regs<Cls, RegT>(cls);
}

// ─── Module definition ───────────────────────────────────────────────────────

PYBIND11_MODULE(_core, m) {

    // ── Base marker classes ──────────────────────────────────────────────────
    py::class_<CardinalitySketch>(m, "CardinalitySketch")
        .def("add",      [](CardinalitySketch& self, const std::string& x) {
            static_cast<SketchBase&>(self).add(x);
        }, py::arg("x"))
        .def("add_many", [](CardinalitySketch& self, const std::vector<std::string>& elems) {
            static_cast<SketchBase&>(self).add_many(elems);
        }, py::arg("elems"))
        .def("estimate",              [](const CardinalitySketch& self) { return static_cast<const SketchBase&>(self).estimate(); })
        .def("memory_usage", [](const CardinalitySketch& self, uint64_t flags) {
            return static_cast<const SketchBase&>(self).memory_usage(flags);
        }, py::arg("flags"));
    py::class_<WeightedMixin>(m, "WeightedMixin");
    py::class_<MergeableMixin>(m, "MergeableMixin");
    py::class_<JaccardMixin>(m, "JaccardMixin");
    py::class_<NewtonMixin>(m, "NewtonMixin");

    auto mf = m.def_submodule("MemoryFlag");
    mf.attr("NOTHING")                    = MemoryFlag::NOTHING;
    mf.attr("REGISTERS")                  = MemoryFlag::REGISTERS;
    mf.attr("ALL_WRITE_NO_REGISTERS")     = MemoryFlag::ALL_WRITE_NO_REGISTERS;
    mf.attr("FISHER_YATES_PERM_INIT")     = MemoryFlag::FISHER_YATES_PERM_INIT;
    mf.attr("FISHER_YATES_NON_PERM_INIT") = MemoryFlag::FISHER_YATES_NON_PERM_INIT;
    mf.attr("SEEDS")                      = MemoryFlag::SEEDS;
    mf.attr("TOTAL")                      = MemoryFlag::TOTAL;
    mf.attr("ALL_WRITE")                  = MemoryFlag::ALL_WRITE;
    mf.attr("FISHER_YATES")               = MemoryFlag::FISHER_YATES;

    // ── RngEngine enum ───────────────────────────────────────────────────────
    py::enum_<RngEngine>(m, "RngEngine")
        .value("PCG64",        RngEngine::PCG64)
        .value("MT19937",      RngEngine::MT19937)
        .value("XOSHIRO128PP", RngEngine::XOSHIRO128PP)
        .value("XOSHIRO256PP", RngEngine::XOSHIRO256PP);

    // ── ExpSketch family ─────────────────────────────────────────────────────
    bind_jaccard_sketch<ExpSketchT<double>, double>(m, "ExpSketch");
    bind_jaccard_sketch<ExpSketchT<float>,  float >(m, "ExpSketchFloat32");
    {
        using Cls = FastExpSketchT<float>;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "FastExpSketchFloat32")
            .def(py::init<std::size_t, std::uint64_t, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        bind_pickle_regs<Cls, float>(cls);
    }
    {
        using Cls = FastExpSketchT<double>;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "FastExpSketch")
            .def(py::init<std::size_t, std::uint64_t, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        bind_pickle_regs<Cls, double>(cls);
    }
    {
        using Cls = FastGMExpSketch;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "FastGMExpSketch")
            .def(py::init<std::size_t, std::uint64_t, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        bind_pickle_regs<Cls, double>(cls);
    }
    bind_mergeable_sketch<WeightedMinHash, double>(m, "WeightedMinHash");
    bind_mergeable_sketch<WeightedHyperLogLog, double>(m, "WeightedHyperLogLog");
    bind_mergeable_sketch<WeightedHyperLogLogFloat32, float>(m, "WeightedHyperLogLogFloat32");

    {
        using Cls = MinHash;
        auto cls = py::class_<Cls, CardinalitySketch, MergeableMixin, JaccardMixin>(m, "MinHash")
            .def(py::init<std::size_t, std::uint64_t>(), py::arg("m"), py::arg("seed"));
        bind_unweighted_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        bind_pickle_regs<Cls, double>(cls);
    }

    {
        using Cls = MartingaleMinHash;
        auto cls = py::class_<Cls, CardinalitySketch>(m, "MartingaleMinHash")
            .def(py::init<std::size_t, std::uint64_t>(), py::arg("m"), py::arg("seed"));
        bind_unweighted_base(cls);
        cls.def("get_E", &Cls::get_E);
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_registers(), p.get_E());
            },
            [](const py::tuple& t) {
                if (t.size() != 4) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::vector<double>>(),
                           t[3].cast<double>());
            }
        ));
    }

    {
        using Cls = HyperLogLog;
        auto cls = py::class_<Cls, CardinalitySketch, MergeableMixin>(m, "HyperLogLog")
            .def(py::init<std::size_t, std::uint64_t>(), py::arg("m"), py::arg("seed"));
        bind_unweighted_base(cls)
            .def("merge", &Cls::merge, py::arg("other"));
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(), p.get_registers());
            },
            [](const py::tuple& t) {
                if (t.size() != 3) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::vector<uint8_t>>());
            }
        ));
    }

    // ── QSketch family ───────────────────────────────────────────────────────
    {
        using Cls = QSketch;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin>(m, "QSketch")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("merge", &Cls::merge, py::arg("other"));
        bind_pickle_q(cls);
    }

    {
        using Cls = QSketchDyn;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin>(m, "QSketchDyn")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, std::uint32_t>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("g_seed") = 42);
        bind_sketch_base(cls)
            .def(py::pickle(
                [](const Cls& p) {
                    return py::make_tuple(
                        p.get_sketch_size(), p.get_amount_bits(), p.get_g_seed(),
                        p.get_master_seed(), p.get_registers(),
                        p.get_t_histogram(), p.get_cardinality());
                },
                [](const py::tuple& t) {
                    if (t.size() != 7) throw std::runtime_error("Invalid pickle state!");
                    return Cls(t[0].cast<std::size_t>(),
                               t[1].cast<std::uint8_t>(),
                               t[2].cast<std::uint32_t>(),
                               t[3].cast<std::uint64_t>(),
                               t[4].cast<std::vector<int>>(),
                               t[5].cast<std::vector<std::uint32_t>>(),
                               t[6].cast<double>());
                }
            ));
    }

    // ── kQSketch family ─────────────────────────────────────────────────────
    {
        using Cls = kQSketch;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, NewtonMixin>(m, "kQSketch")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, float, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("logarithm_base"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("merge", &Cls::merge, py::arg("other"))
            .def("estimate_direct", &Cls::estimate_direct)
            .def("estimate_newton_cold", &Cls::estimate_newton_cold)
            .def("estimate_newton_warm", &Cls::estimate_newton_warm)
            .def("estimate_newton_cold_iterations", &Cls::estimate_newton_cold_iterations)
            .def("estimate_newton_warm_iterations", &Cls::estimate_newton_warm_iterations);
        bind_pickle_log_exp(cls);
    }

    {
        using Cls = kQSketchRounding;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, NewtonMixin>(m, "kQSketchRounding")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, float, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("logarithm_base"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("merge", &Cls::merge, py::arg("other"))
            .def("estimate_corrected", &Cls::estimate_corrected)
            .def("estimate_direct", &Cls::estimate_direct)
            .def("estimate_newton_cold", &Cls::estimate_newton_cold)
            .def("estimate_newton_warm", &Cls::estimate_newton_warm)
            .def("estimate_newton_cold_iterations", &Cls::estimate_newton_cold_iterations)
            .def("estimate_newton_warm_iterations", &Cls::estimate_newton_warm_iterations);
        bind_pickle_log_exp(cls);
    }

    // ── kQSketchRoundedDyn (martingale) ─────────────────────────────────────
    {
        using Cls = kQSketchRoundedDyn;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, NewtonMixin>(m, "kQSketchRoundedDyn")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, float, std::uint32_t>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("logarithm_base"), py::arg("g_seed") = 42);
        bind_sketch_base(cls)
            .def("estimate_direct", &Cls::estimate_direct)
            .def("estimate_newton_cold", &Cls::estimate_newton_cold)
            .def("estimate_newton_warm", &Cls::estimate_newton_warm)
            .def("estimate_newton_cold_iterations", &Cls::estimate_newton_cold_iterations)
            .def("estimate_newton_warm_iterations", &Cls::estimate_newton_warm_iterations)
            .def(py::pickle(
                [](const Cls& p) {
                    return py::make_tuple(
                        p.get_sketch_size(), p.get_amount_bits(), p.get_logarithm_base(),
                        p.get_g_seed(), p.get_master_seed(), p.get_registers(),
                        p.get_t_histogram(), p.get_cardinality());
                },
                [](const py::tuple& t) {
                    if (t.size() != 8) throw std::runtime_error("Invalid pickle state!");
                    return Cls(t[0].cast<std::size_t>(),
                               t[1].cast<std::uint8_t>(),
                               t[2].cast<float>(),
                               t[3].cast<std::uint32_t>(),
                               t[4].cast<std::uint64_t>(),
                               t[5].cast<std::vector<int>>(),
                               t[6].cast<std::vector<std::uint32_t>>(),
                               t[7].cast<double>());
                }
            ));
    }


    // ── kQSketchShifted (fast) ───────────────────────────────────────────────
    {
        using Cls = kQSketchShifted;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, NewtonMixin>(m, "kQSketchShifted")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, float, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("logarithm_base"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("merge", &Cls::merge, py::arg("other"))
            .def("get_offset", &Cls::get_offset)
            .def("estimate_direct", &Cls::estimate_direct)
            .def("estimate_newton_cold", &Cls::estimate_newton_cold)
            .def("estimate_newton_warm", &Cls::estimate_newton_warm)
            .def("estimate_newton_cold_iterations", &Cls::estimate_newton_cold_iterations)
            .def("estimate_newton_warm_iterations", &Cls::estimate_newton_warm_iterations);
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_amount_bits(), p.get_registers(),
                                      p.get_logarithm_base(), p.get_offset());
            },
            [](const py::tuple& t) {
                if (t.size() != 6) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::uint8_t>(),
                           t[4].cast<float>(),
                           t[3].cast<std::vector<int>>(),
                           t[5].cast<int>());
            }
        ));
    }

    // ── QuantizationMode enum ────────────────────────────────────────────────
    py::enum_<QuantizationMode>(m, "QuantizationMode")
        .value("ALL_NORMAL",      QuantizationMode::ALL_NORMAL)
        .value("WITH_SUBNORMALS", QuantizationMode::WITH_SUBNORMALS)
        .value("LINEAR",          QuantizationMode::LINEAR)
        .value("LOGARITHMIC",     QuantizationMode::LOGARITHMIC);

    // ── LogExpSketchSlowNoShifted ─────────────────────────────────────────────────────────
    {
        using Cls = LogExpSketchSlowNoShifted;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "LogExpSketchSlowNoShifted")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, double>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("v_max"));
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_amount_bits(), p.get_v_max(),
                                      p.get_registers());
            },
            [](const py::tuple& t) {
                if (t.size() != 5) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::uint8_t>(),
                           t[3].cast<double>(),
                           t[4].cast<std::vector<int>>());
            }
        ));
    }

    // ── LogExpSketchSlowShifted ──────────────────────────────────────────────────────────
    {
        using Cls = LogExpSketchSlowShifted;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "LogExpSketchSlowShifted")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, double>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("v_max"));
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"))
            .def("get_offset", &Cls::get_offset);
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_amount_bits(), p.get_v_max(),
                                      p.get_registers(), p.get_offset());
            },
            [](const py::tuple& t) {
                if (t.size() != 6) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::uint8_t>(),
                           t[3].cast<double>(),
                           t[4].cast<std::vector<int>>(),
                           t[5].cast<int>());
            }
        ));
    }

    // ── WeightedHyperLogLogCustomFloat ───────────────────────────────────────
    {
        using Cls = WeightedHyperLogLogCustomFloat;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin>(m, "WeightedHyperLogLogCustomFloat")
            .def(py::init<std::size_t, std::uint64_t, int, int>(),
                 py::arg("m"), py::arg("seed"),
                 py::arg("exp_bits"), py::arg("mant_bits"));
        bind_sketch_base(cls)
            .def("merge", &Cls::merge, py::arg("other"))
            .def_property_readonly("exp_bits",  &Cls::get_exp_bits)
            .def_property_readonly("mant_bits", &Cls::get_mant_bits)
            .def(py::pickle(
                [](const Cls& p) {
                    return py::make_tuple(
                        p.get_sketch_size(), p.get_master_seed(),
                        p.get_exp_bits(), p.get_mant_bits(),
                        p.get_registers());
                },
                [](const py::tuple& t) {
                    if (t.size() != 5) throw std::runtime_error("Invalid pickle state!");
                    return Cls(t[0].cast<std::size_t>(),
                               t[1].cast<std::uint64_t>(),
                               t[2].cast<int>(),
                               t[3].cast<int>(),
                               t[4].cast<std::vector<double>>());
                }
            ));
    }

    // ── FastExpSketchCustomFloat ─────────────────────────────────────────────
    {
        using Cls = FastExpSketchCustomFloat;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "FastExpSketchCustomFloat")
            .def(py::init<std::size_t, std::uint64_t, int, int, RngEngine>(),
                 py::arg("m"), py::arg("seed"),
                 py::arg("exp_bits"), py::arg("mant_bits"),
                 py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge",      &Cls::merge, py::arg("other"))
            .def("clone_with", &Cls::clone_with, py::arg("exp_bits"), py::arg("mant_bits"))
            .def_property_readonly("exp_bits",  &Cls::get_exp_bits)
            .def_property_readonly("mant_bits", &Cls::get_mant_bits)
            .def(py::pickle(
                [](const Cls& p) {
                    return py::make_tuple(
                        p.get_sketch_size(), p.get_master_seed(),
                        p.get_exp_bits(), p.get_mant_bits(),
                        p.get_registers());
                },
                [](const py::tuple& t) {
                    if (t.size() != 5) throw std::runtime_error("Invalid pickle state!");
                    return Cls(t[0].cast<std::size_t>(),
                               t[1].cast<std::uint64_t>(),
                               t[2].cast<int>(),
                               t[3].cast<int>(),
                               t[4].cast<std::vector<double>>());
                }
            ));
    }

    // ── LogExpSketchFastNoShifted ─────────────────────────────────────────────────────
    {
        using Cls = LogExpSketchFastNoShifted;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "LogExpSketchFastNoShifted")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, double, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("v_max"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"));
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_amount_bits(), p.get_v_max(),
                                      p.get_registers());
            },
            [](const py::tuple& t) {
                if (t.size() != 5) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::uint8_t>(),
                           t[3].cast<double>(),
                           t[4].cast<std::vector<int>>());
            }
        ));
    }

    // ── LogExpSketchFastShifted ───────────────────────────────────────────────────────
    {
        using Cls = LogExpSketchFastShifted;
        auto cls = py::class_<Cls, CardinalitySketch, WeightedMixin, MergeableMixin, JaccardMixin>(m, "LogExpSketchFastShifted")
            .def(py::init<std::size_t, std::uint64_t, std::uint8_t, double, RngEngine>(),
                 py::arg("m"), py::arg("seed"), py::arg("amount_bits"), py::arg("v_max"), py::arg("rng_engine") = kDefaultRngEngine);
        bind_sketch_base(cls)
            .def("jaccard_struct", &Cls::jaccard_struct)
            .def("merge", &Cls::merge, py::arg("other"))
            .def("get_offset", &Cls::get_offset);
        cls.def(py::pickle(
            [](const Cls& p) {
                return py::make_tuple(p.get_sketch_size(), p.get_master_seed(),
                                      p.get_amount_bits(), p.get_v_max(),
                                      p.get_registers(), p.get_offset());
            },
            [](const py::tuple& t) {
                if (t.size() != 6) throw std::runtime_error("Invalid pickle state!");
                return Cls(t[0].cast<std::size_t>(),
                           t[1].cast<std::uint64_t>(),
                           t[2].cast<std::uint8_t>(),
                           t[3].cast<double>(),
                           t[4].cast<std::vector<int>>(),
                           t[5].cast<int>());
            }
        ));
    }
}
