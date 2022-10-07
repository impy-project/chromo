#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Pythia8/Event.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/ParticleData.h>
#include <array>

namespace py = pybind11;
using namespace Pythia8;
using namespace pybind11::literals;

struct Hepevt
{
    void fill(const Event &event)
    {
        ++nevhep;
        if (nhep > py::len(idhep))
        {
            int maxsize = nhep * 1.6;
            idhep = py::array_t<int>({maxsize});
            isthep = py::array_t<int>({maxsize});
            phep = py::array_t<double>({5, maxsize});
            vhep = py::array_t<double>({4, maxsize});
            jmohep = py::array_t<int>({2, maxsize});
        }
        nhep = event.size();

        // TODO copy stuff over
    }

    int nevhep = 0;
    int nhep;
    py::array_t<int> idhep;
    py::array_t<int> isthep;
    py::array_t<double> phep;
    py::array_t<double> vhep;
    py::array_t<int> jmohep;
} hepevt;

PYBIND11_MODULE(_pythia8, m)
{
    py::class_<Hepevt>(m, "Hepevt")
        // .def(py::init<const Event &>())
        .def_readwrite("nevhep", &Hepevt::nevhep)
        .def_readwrite("nhep", &Hepevt::nhep)
        .def_readwrite("idhep", &Hepevt::idhep)
        .def_readwrite("isthep", &Hepevt::isthep)
        .def_readwrite("phep", &Hepevt::phep)
        .def_readwrite("vhep", &Hepevt::vhep)
        .def_readwrite("jmohep", &Hepevt::jmohep)

        ;

    py::class_<Pythia>(m, "Pythia")
        .def(py::init<string, bool>())
        .def("init", &Pythia::init)
        .def("next", [](Pythia &self)
             {
            bool ok = self.next();
            if (ok)
                hepevt.fill(self.event);
            return ok; })
        .def("readString", &Pythia::readString, "setting"_a, "warn"_a = true)
        // .def_readonly("info", &Pythia::info)
        .def_readwrite("particleData", &Pythia::particleData)
        .def_property_readonly("hepevt", [](Pythia &self)
                               { return hepevt; })

        ;

    py::class_<ParticleData>(m, "ParticleData")
        .def("mayDecay", py::overload_cast<int, bool>(&ParticleData::mayDecay))

        ;
}