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
        .def("jaccard_struct", &ExpSketch::jaccard_struct);
 
    py::class_<FastExpSketch>(m, "FastExpSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&>())
        .def("add", &FastExpSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("estimate", &FastExpSketch::estimate)
        .def("add_many", &FastExpSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("jaccard_struct", &FastExpSketch::jaccard_struct);

    py::class_<FastQSketch>(m, "FastQSketch")
        .def(py::init<std::size_t, const std::vector<std::uint32_t>&, std::uint8_t>(),
            py::arg("m"), py::arg("seeds"), py::arg("amount_seeds"))
        .def("add", &FastQSketch::add, py::arg("x"), py::arg("weight") = 1.0)
        .def("add_many", &FastQSketch::add_many, py::arg("elems"), py::arg("weights"))
        .def("estimate", &FastQSketch::estimate);
}