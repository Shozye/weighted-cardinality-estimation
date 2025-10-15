#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "base_log_exp_sketch.hpp"
#include "base_q_sketch.hpp"
#include "exp_sketch.hpp"
#include "fast_exp_sketch.hpp"
#include "fast_q_sketch.hpp"
#include "fastgm_exp_sketch.hpp"
#include "q_sketch_dyn.hpp"
#include "q_sketch.hpp"
#include "fast_log_exp_sketch.hpp"
#include "base_shifted_log_exp_sketch.hpp"
#include "fast_shifted_log_exp_sketch.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    py::class_<ExpSketch>(m, "ExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&>(),
             py::arg("m"), py::arg("seeds"))
        .def("add", &ExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("estimate", &ExpSketch::estimate)
        .def("add_many", &ExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("jaccard_struct", &ExpSketch::jaccard_struct)
        .def("memory_usage_total", &ExpSketch::memory_usage_total)
        .def("memory_usage_write", &ExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &ExpSketch::memory_usage_estimate)
        .def(py::pickle(
            [](const ExpSketch &p) {return py::make_tuple(p.get_sketch_size(), p.get_seeds(), p.get_registers());},
            [](const py::tuple& t) {
                if (t.size() != 3) {
                    throw std::runtime_error("Invalid state for ExpSketch pickle!");
                }
                return ExpSketch( 
                    t[0].cast<std::size_t>(),
                    t[1].cast<std::vector<std::uint32_t>>(),
                    t[2].cast<std::vector<double>>()
                );
            }
        ));
 
    py::class_<FastExpSketch>(m, "FastExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&>())
        .def("add", &FastExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("estimate", &FastExpSketch::estimate)
        .def("add_many", &FastExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("jaccard_struct", &FastExpSketch::jaccard_struct)
        .def("memory_usage_total", &FastExpSketch::memory_usage_total)
        .def("memory_usage_write", &FastExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &FastExpSketch::memory_usage_estimate)
        .def(py::pickle(
    [](const FastExpSketch &p) {
        return py::make_tuple(
            p.get_sketch_size(),
            p.get_seeds(),
            p.get_registers()
        );
    },
    [](const py::tuple& t) {
        if (t.size() != 3) {
            throw std::runtime_error("Invalid state for FastExpSketch pickle!");
        }
        return FastExpSketch(
            t[0].cast<std::size_t>(),
            t[1].cast<std::vector<std::uint32_t>>(),
            t[2].cast<std::vector<double>>()
        );
    }
    ));;

    py::class_<FastGMExpSketch>(m, "FastGMExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&>())
        .def("add", &FastGMExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("estimate", &FastGMExpSketch::estimate)
        .def("add_many", &FastGMExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("jaccard_struct", &FastGMExpSketch::jaccard_struct)
        .def("memory_usage_total", &FastGMExpSketch::memory_usage_total)
        .def("memory_usage_write", &FastGMExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &FastGMExpSketch::memory_usage_estimate)
        .def(py::pickle(
    [](const FastGMExpSketch &p) {
        return py::make_tuple(
            p.get_sketch_size(),
            p.get_seeds(),
            p.get_registers()
        );
    },
    [](const py::tuple& t) {
        if (t.size() != 3) {
            throw std::runtime_error("Invalid state for FastGMExpSketch pickle!");
        }
        return FastGMExpSketch(
            t[0].cast<std::size_t>(),
            t[1].cast<std::vector<std::uint32_t>>(),
            t[2].cast<std::vector<double>>()
        );
    }
    ));;

    py::class_<FastQSketch>(m, "FastQSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"))
        .def("add", &FastQSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &FastQSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &FastQSketch::estimate)
        .def("memory_usage_total", &FastQSketch::memory_usage_total)
        .def("memory_usage_write", &FastQSketch::memory_usage_write)
        .def("memory_usage_estimate", &FastQSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const FastQSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 4) {
                throw std::runtime_error("Invalid state for FastQSketch pickle!");
            }
            return FastQSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[3].cast<std::vector<int>>()
            );
        }
    ));;

    py::class_<QSketchDyn>(m, "QSketchDyn")
    .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t, std::uint32_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"), py::arg("g_seed") = 42)
    .def("add", &QSketchDyn::add, py::arg("elem"), py::arg("weight") = 1.0)
    .def("add_many", &QSketchDyn::add_many, py::arg("elems"), py::arg("weights"))
    .def("estimate", &QSketchDyn::estimate)
    .def("memory_usage_total", &QSketchDyn::memory_usage_total)
    .def("memory_usage_write", &QSketchDyn::memory_usage_write)
    .def("memory_usage_estimate", &QSketchDyn::memory_usage_estimate)
    .def(py::pickle(
        [](const QSketchDyn &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_amount_bits(),
                p.get_g_seed(),
                p.get_seeds(),
                p.get_registers(),
                p.get_t_histogram(),
                p.get_cardinality()
            );
        },
        [](const py::tuple &t) {
            if (t.size() != 7) { throw std::runtime_error("Invalid state for QSketchDyn pickle!"); }
            return QSketchDyn(
                t[0].cast<std::size_t>(),
                t[1].cast<std::uint8_t>(),
                t[2].cast<std::uint32_t>(),
                t[3].cast<std::vector<std::uint32_t>>(),
                t[4].cast<std::vector<int>>(),
                t[5].cast<std::vector<int>>(),
                t[6].cast<double>()
            );
        }
    ));

    py::class_<BaseQSketch>(m, "BaseQSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"))
        .def("add", &BaseQSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &BaseQSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &BaseQSketch::estimate)
        .def("memory_usage_total", &BaseQSketch::memory_usage_total)
        .def("memory_usage_write", &BaseQSketch::memory_usage_write)
        .def("memory_usage_estimate", &BaseQSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const BaseQSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 4) {
                throw std::runtime_error("Invalid state for BaseQSketch pickle!");
            }
            return BaseQSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[3].cast<std::vector<int>>()
            );
        }
    ));

    py::class_<QSketch>(m, "QSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"))
        .def("add", &QSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &QSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &QSketch::estimate)
        .def("memory_usage_total", &QSketch::memory_usage_total)
        .def("memory_usage_write", &QSketch::memory_usage_write)
        .def("memory_usage_estimate", &QSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const QSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 4) {
                throw std::runtime_error("Invalid state for QSketch pickle!");
            }
            return QSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[3].cast<std::vector<int>>()
            );
        }
    ));

    py::class_<FastLogExpSketch>(m, "FastLogExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t, float>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"), py::arg("logarithm_base"))
        .def("add", &FastLogExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &FastLogExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &FastLogExpSketch::estimate)
        .def("memory_usage_total", &FastLogExpSketch::memory_usage_total)
        .def("memory_usage_write", &FastLogExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &FastLogExpSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const FastLogExpSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers(),
                p.get_logarithm_base()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 5) {
                throw std::runtime_error("Invalid state for FastLogExpSketch pickle!");
            }
            return FastLogExpSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[4].cast<float>(),
                t[3].cast<std::vector<int>>()
            );
        }
    ));;

    py::class_<BaseLogExpSketch>(m, "BaseLogExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t, float>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"), py::arg("logarithm_base"))
        .def("add", &BaseLogExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &BaseLogExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &BaseLogExpSketch::estimate)
        .def("memory_usage_total", &BaseLogExpSketch::memory_usage_total)
        .def("memory_usage_write", &BaseLogExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &BaseLogExpSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const BaseLogExpSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers(),
                p.get_logarithm_base()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 5) {
                throw std::runtime_error("Invalid state for BaseLogExpSketch pickle!");
            }
            return BaseLogExpSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[4].cast<float>(),
                t[3].cast<std::vector<int>>()
            );
        }
    ));;
    py::class_<BaseShiftedLogExpSketch>(m, "BaseShiftedLogExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t, float>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"), py::arg("logarithm_base"))
        .def("add", &BaseShiftedLogExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &BaseShiftedLogExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &BaseShiftedLogExpSketch::estimate)
        .def("memory_usage_total", &BaseShiftedLogExpSketch::memory_usage_total)
        .def("memory_usage_write", &BaseShiftedLogExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &BaseShiftedLogExpSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const BaseShiftedLogExpSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_logarithm_base(),
                p.get_registers(),
                p.get_offset()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 6) {
                throw std::runtime_error("Invalid state for BaseShiftedLogExpSketch pickle!");
            }
            return BaseShiftedLogExpSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[3].cast<float>(),
                t[4].cast<std::vector<uint32_t>>(),
                t[5].cast<int>()
            );
        }
    ));;

    py::class_<FastShiftedLogExpSketch>(m, "FastShiftedLogExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t, float>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"), py::arg("logarithm_base"))
        .def("add", &FastShiftedLogExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &FastShiftedLogExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &FastShiftedLogExpSketch::estimate)
        .def("memory_usage_total", &FastShiftedLogExpSketch::memory_usage_total)
        .def("memory_usage_write", &FastShiftedLogExpSketch::memory_usage_write)
        .def("memory_usage_estimate", &FastShiftedLogExpSketch::memory_usage_estimate)
        .def(py::pickle(
        [](const FastShiftedLogExpSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_logarithm_base(),
                p.get_registers(),
                p.get_offset()
            );
        },
        [](const py::tuple& t) {
            if (t.size() != 6) {
                throw std::runtime_error("Invalid state for FastShiftedLogExpSketch pickle!");
            }
            return FastShiftedLogExpSketch(
                t[0].cast<std::size_t>(),
                t[1].cast<std::vector<std::uint32_t>>(),
                t[2].cast<std::uint8_t>(),
                t[3].cast<float>(),
                t[4].cast<std::vector<uint32_t>>(),
                t[5].cast<int>()
            );
        }
    ));;
}
