#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Pythia8/Event.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/ParticleData.h>
#include <array>
#include <cassert>

namespace py = pybind11;
using namespace Pythia8;
using namespace pybind11::literals;

struct Hepevt
{
    Hepevt() { reset(10000); }

    void reset(int size)
    {
        idhep = py::array_t<int>({size});
        isthep = py::array_t<int>({size});
        phep = py::array_t<double>({5, size});
        vhep = py::array_t<double>({4, size});
        jmohep = py::array_t<int>({2, size});
    }

    void fill(const Event &event)
    {
        ++nevhep;
        nhep = event.size();
        if (nhep > py::len(idhep))
            reset(nhep);

        auto pid = idhep.mutable_unchecked<1>();
        auto sta = isthep.mutable_unchecked<1>();
        auto p = phep.mutable_unchecked<2>();
        auto v = vhep.mutable_unchecked<2>();
        auto par = jmohep.mutable_unchecked<2>();
        int k = 0;
        for (const auto &particle : event)
        {
            // record starts with internal particle 90,
            // skip this to start with beam particles
            if (k == 0 && particle.id() == 90)
            {
                --nhep;
                continue;
            }
            pid(k) = particle.id();
            sta(k) = particle.statusHepMC();
            p(0, k) = particle.px();
            p(1, k) = particle.py();
            p(2, k) = particle.pz();
            p(3, k) = particle.e();
            p(4, k) = particle.m();
            v(0, k) = particle.xProd();
            v(1, k) = particle.yProd();
            v(2, k) = particle.zProd();
            v(3, k) = particle.tProd();
            par(0, k) = particle.mother1();
            par(1, k) = particle.mother2();
            ++k;
        }
    }

    int nevhep = -1;
    int nhep = 0;
    py::array_t<int, py::array::f_style> idhep;
    py::array_t<int, py::array::f_style> isthep;
    py::array_t<double, py::array::f_style> phep;
    py::array_t<double, py::array::f_style> vhep;
    py::array_t<int, py::array::f_style> jmohep;
};

int charge_from_pid(const ParticleData &pd, int pid)
{
    auto pptr = pd.findParticle(pid);
    assert(pptr); // never be a nullptr if charge_from_pid is used on particles produced by Pythia
    return pptr->charge();
}

PYBIND11_MODULE(_pythia8, m)
{
    py::class_<Hepevt>(m, "Hepevt")
        .def(py::init<>())
        .def("fill", &Hepevt::fill)
        .def_readwrite("nevhep", &Hepevt::nevhep)
        .def_readwrite("nhep", &Hepevt::nhep)
        .def_readwrite("idhep", &Hepevt::idhep)
        .def_readwrite("isthep", &Hepevt::isthep)
        .def_readwrite("phep", &Hepevt::phep)
        .def_readwrite("vhep", &Hepevt::vhep)
        .def_readwrite("jmohep", &Hepevt::jmohep)

        ;

    py::class_<ParticleData>(m, "ParticleData")
        .def("mayDecay", py::overload_cast<int, bool>(&ParticleData::mayDecay))
        .def("addParticle", py::overload_cast<int, string, string, int, int, int, double, double, double, double, double, bool>(&ParticleData::addParticle), "pdgid"_a, "name"_a, "antiname"_a, "spin"_a = 0, "charge"_a = 0, "color"_a = 0, "mass"_a = 0, "width"_a = 0, "min"_a = 0, "max"_a = 0, "tau"_a = 0, "var_width"_a = false)
        .def("isParticle", &ParticleData::isParticle)
        .def("findParticle", py::overload_cast<int>(&ParticleData::findParticle))

        ;

    py::class_<ParticleDataEntry, ParticleDataEntryPtr>(m, "ParticleDataEntry")
        .def_property_readonly("id", &ParticleDataEntry::id)
        .def_property_readonly("antiId", &ParticleDataEntry::antiId)
        .def_property_readonly("hasAnti", &ParticleDataEntry::hasAnti)
        .def_property_readonly("name", [](const ParticleDataEntry &self)
                               { return self.name(); })
        .def_property_readonly("charge", [](const ParticleDataEntry &self)
                               { return self.charge(); })
        .def_property_readonly("spinType", &ParticleDataEntry::spinType)
        .def_property_readonly("chargeType", [](const ParticleDataEntry &self)
                               { return self.chargeType(); })
        .def_property_readonly("colType", [](const ParticleDataEntry &self)
                               { return self.colType(); })
        .def_property_readonly("m0", &ParticleDataEntry::m0)
        .def_property_readonly("mWidth", &ParticleDataEntry::mWidth)
        .def_property_readonly("mMin", &ParticleDataEntry::mMin)
        .def_property_readonly("mMax", &ParticleDataEntry::mMax)
        .def_property_readonly("m0Min", &ParticleDataEntry::m0Min)
        .def_property_readonly("m0Max", &ParticleDataEntry::m0Max)
        .def_property_readonly("tau0", &ParticleDataEntry::tau0)
        .def_property_readonly("isResonance", &ParticleDataEntry::isResonance)
        .def_property_readonly("varWidth", &ParticleDataEntry::varWidth)
        .def_property_readonly("mayDecay", &ParticleDataEntry::mayDecay)
        .def_property_readonly("tauCalc", &ParticleDataEntry::tauCalc)
        .def_property_readonly("doExternalDecay", &ParticleDataEntry::doExternalDecay)
        .def_property_readonly("isVisible", &ParticleDataEntry::isVisible)
        .def_property_readonly("doForceWidth", &ParticleDataEntry::doForceWidth)
        .def_property_readonly("hasChanged", &ParticleDataEntry::hasChanged)
        .def_property_readonly("hasChangedMMin", &ParticleDataEntry::hasChangedMMin)
        .def_property_readonly("hasChangedMMax", &ParticleDataEntry::hasChangedMMax)
        .def_property_readonly("initBWmass", &ParticleDataEntry::initBWmass)
        .def_property_readonly("constituentMass", &ParticleDataEntry::constituentMass)
        .def_property_readonly("mSel", &ParticleDataEntry::mSel)
        .def("mRun", &ParticleDataEntry::mRun)
        .def_property_readonly("useBreitWigner", &ParticleDataEntry::useBreitWigner)
        .def_property_readonly("canDecay", &ParticleDataEntry::canDecay)
        .def_property_readonly("isLepton", &ParticleDataEntry::isLepton)
        .def_property_readonly("isQuark", &ParticleDataEntry::isQuark)
        .def_property_readonly("isGluon", &ParticleDataEntry::isGluon)
        .def_property_readonly("isDiquark", &ParticleDataEntry::isDiquark)
        .def_property_readonly("isParton", &ParticleDataEntry::isParton)
        .def_property_readonly("isHadron", &ParticleDataEntry::isHadron)
        .def_property_readonly("isMeson", &ParticleDataEntry::isMeson)
        .def_property_readonly("isBaryon", &ParticleDataEntry::isBaryon)
        .def_property_readonly("isOnium", &ParticleDataEntry::isOnium)
        .def_property_readonly("isOctetHadron", &ParticleDataEntry::isOctetHadron)
        .def_property_readonly("heaviestQuark", [](const ParticleDataEntry &self)
                               { return self.heaviestQuark(); })
        .def_property_readonly("baryonNumberType", [](const ParticleDataEntry &self)
                               { return self.baryonNumberType(); })
        .def("nQuarksInCode", &ParticleDataEntry::nQuarksInCode)

        ;

    py::class_<Event>(m, "Event");

    py::class_<Pythia>(m, "Pythia")
        .def(py::init<string, bool>())
        .def("init", &Pythia::init)
        .def("next", py::overload_cast<>(&Pythia::next))
        .def("readString", &Pythia::readString, "setting"_a, "warn"_a = true)
        .def_readwrite("particleData", &Pythia::particleData)
        // .def_readonly("info", &Pythia::info)
        .def_readwrite("event", &Pythia::event)

        ;

    // only vectorize over pid, not over ParticleData
    m.def("charge_from_pid",
          [](const ParticleData &pd, py::array_t<int> pid)
          { return py::vectorize([&pd](int pid)
                                 { return charge_from_pid(pd, pid); })(std::move(pid)); });
}