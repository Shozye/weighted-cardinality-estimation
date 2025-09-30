#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "exp_sketch.hpp"
#include "fast_exp_sketch.hpp"
#include "fast_q_sketch.hpp"


namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    py::class_<ExpSketch>(m, "ExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&>(),
             py::arg("m"), py::arg("seeds"))
        .def("add", &ExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("estimate", &ExpSketch::estimate)
        .def("add_many", &ExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("jaccard_struct", &ExpSketch::jaccard_struct)
        .def(py::pickle(
            [](const ExpSketch &p) {return py::make_tuple(p.get_m(), p.get_seeds(), p.get_registers());},
            [](py::tuple t) {
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
        .def(py::pickle(
    // __getstate__: Zwróć tylko 3 fundamentalne pola
    [](const FastExpSketch &p) {
        return py::make_tuple(
            p.get_m(),
            p.get_seeds(),
            p.get_registers()
        );
    },
    // __setstate__: Wywołaj nowy, inteligentny konstruktor
    [](py::tuple t) {
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

    py::class_<FastQSketch>(m, "FastQSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_bits"))
        .def("add", &FastQSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &FastQSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &FastQSketch::estimate)
        .def(py::pickle(
        // __getstate__: Zwróć tylko 3 fundamentalne pola
        [](const FastQSketch &p) {
            return py::make_tuple(
                p.get_sketch_size(),
                p.get_seeds(),
                p.get_amount_bits(),
                p.get_registers()
            );
        },
        // __setstate__: Wywołaj nowy, inteligentny konstruktor
        [](py::tuple t) {
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
}