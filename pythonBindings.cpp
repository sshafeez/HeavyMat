//g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` pythonBindings.cpp -o HeavyMat`python3-config --extension-suffix`
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include "matrix.h"
#include <string>
using namespace std;
using namespace pybind11::literals;
namespace py = pybind11;

PYBIND11_MODULE(HeavyMat, m) {
    py::class_<matrix>(m, "matrix")
        //Constructors
        .def(py::init< int, int, bool >(), "rows"_a, "cols"_a, "initRandom"_a = false)
        .def(py::init< const matrix& >())
        .def(py::init< vector<vector<double>>& >())

        //inner data access
        .def_readwrite("grid", &matrix::grid)
        .def("at", &matrix::at, py::return_value_policy::reference)
        .def("set", &matrix::set)
        .def("print", &matrix::print)
        .def("size",&matrix::size)
        .def("__repr__", 
            [](matrix& mat){
                string result="";
                result += to_string(mat.grid.size()) + " : " + to_string(mat.grid[0].size()) + "\n";
                for(vector<double>& row : mat.grid){
                    for(double val : row){
                        result += to_string(val) + " ";		
                    }
                result += "\n";
                }
                return result;
            }
        )

        //data modifications
        .def("append", &matrix::append)
        .def("invert", &matrix::invert)
        .def("transpose", &matrix::transpose)
    ;

    //multiplication
    m.def("multiply",&multiply);
    m.def("fastMultiply",&fastMultiply);

    //Heavy Mat Class
    py::class_<heavy_matrix,matrix>(m, "heavy_matrix")
        //Constructors
        .def(py::init< int, int, bool >(), "rows"_a, "cols"_a, "initRandom"_a = false)
        .def(py::init< const heavy_matrix& >())
        .def(py::init< vector<vector<double>>& >())

        //compression
        .def("cache", &heavy_matrix::cache, "search for linear dependencies", py::arg("error")=epsilon)
        .def("writeback", &heavy_matrix::writeback)

    ;

}