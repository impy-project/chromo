#include <Pythia8/Event.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/Info.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/PythiaStdlib.h>
#include <Pythia8/Basics.h>
#include <Pythia8Plugins/PythiaCascade.h>
#include <array>
#include <cassert>
#include <limits>
#include <private_access.hpp>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Pythia8;
using namespace pybind11::literals;

float charge_from_pid(const ParticleData &pd, int pid)
{
    auto pptr = pd.findParticle(pid);

    // return NaN if unknown pid is met
    if (pptr == nullptr)
        return std::numeric_limits<float>::quiet_NaN();

    // ParticleData returns partice even if anti-particle pid is used
    return pid == pptr->id() ? pptr->charge() : -pptr->charge();
}

PRIVATE_ACCESS_MEMBER(Particle, idSave, int)
PRIVATE_ACCESS_MEMBER(Particle, statusSave, int)
PRIVATE_ACCESS_MEMBER(Particle, pSave, Vec4)
PRIVATE_ACCESS_MEMBER(Particle, mSave, double)
PRIVATE_ACCESS_MEMBER(Particle, vProdSave, Vec4)
PRIVATE_ACCESS_MEMBER(Particle, mother1Save, int)
PRIVATE_ACCESS_MEMBER(Particle, mother2Save, int)
PRIVATE_ACCESS_MEMBER(Particle, daughter1Save, int)
PRIVATE_ACCESS_MEMBER(Particle, daughter2Save, int)
PRIVATE_ACCESS_MEMBER(Vec4, xx, double)
PRIVATE_ACCESS_MEMBER(Vec4, yy, double)
PRIVATE_ACCESS_MEMBER(Vec4, zz, double)
PRIVATE_ACCESS_MEMBER(Vec4, tt, double)
PRIVATE_ACCESS_MEMBER(PythiaCascade, pythiaMain, Pythia)
PRIVATE_ACCESS_MEMBER(PythiaCascade, pythiaColl, Pythia)

template <class Accessor>
auto event_array(Event &event)
{
    // skip first pseudoparticle
    auto *ptr = &private_access::member<Accessor>(event[1]);
    int number = event.size() - 1;
    int shape[1] = {number};
    int strides[1] = {sizeof(Particle)};
    return py::array_t<std::decay_t<decltype(*ptr)>>(shape, strides, ptr);
}

auto event_status(Event &event)
{
    // skip first pseudoparticle
    py::array_t<int> result(event.size() - 1);
    int *ptr = result.mutable_data();
    for (auto pit = event.begin() + 1; pit != event.end(); ++pit)
        *ptr++ = pit->statusHepMC();
    return result;
}

template <class Accessor>
py::array_t<double> event_array_p(Event &event)
{
    // skip first pseudoparticle
    auto &pref = private_access::member<Particle_pSave>(event[1]);
    double *ptr = &private_access::member<Accessor>(pref);
    int number = event.size() - 1;
    int shape[1] = {number};
    int strides[1] = {sizeof(Particle)};
    return py::array_t<double>(shape, strides, ptr);
}

template <class Accessor>
py::array_t<double> event_array_v(Event &event)
{
    // skip first pseudoparticle
    auto &pref = private_access::member<Particle_vProdSave>(event[1]);
    double *ptr = &private_access::member<Accessor>(pref);
    int number = event.size() - 1;
    int shape[1] = {number};
    int strides[1] = {sizeof(Particle)};
    return py::array_t<double>(shape, strides, ptr);
}

py::array_t<int> event_array_mothers(Event &event)
{
    // skip first pseudoparticle
    auto ptr1 = &private_access::member<Particle_mother1Save>(event[1]);
    auto ptr2 = &private_access::member<Particle_mother2Save>(event[1]);
    int number = event.size() - 1;
    int shape[2] = {number, 2};
    int strides[2] = {sizeof(Particle), sizeof(int)};
    return py::array_t<int>(shape, strides, ptr1);
}

py::array_t<int> event_array_daughters(Event &event)
{
    // skip first pseudoparticle
    auto ptr1 = &private_access::member<Particle_daughter1Save>(event[1]);
    auto ptr2 = &private_access::member<Particle_daughter2Save>(event[1]);
    int number = event.size() - 1;
    int shape[2] = {number, 2};
    int strides[2] = {sizeof(Particle), sizeof(int)};
    return py::array_t<int>(shape, strides, ptr1);
}

// refills "event" stack with particles
void fill(Event &event,
          py::array_t<int> &pid,
          py::array_t<int> &status,
          py::array_t<double> &px,
          py::array_t<double> &py,
          py::array_t<double> &pz,
          py::array_t<double> &energy,
          py::array_t<double> &mass)
{
    // Get a raw reference to numpy array
    auto pid_ = pid.unchecked<1>();
    auto status_ = status.unchecked<1>();
    auto px_ = px.unchecked<1>();
    auto py_ = py.unchecked<1>();
    auto pz_ = pz.unchecked<1>();
    auto energy_ = energy.unchecked<1>();
    auto mass_ = mass.unchecked<1>();

    event.reset();
    for (int i = 0; i != pid.size(); ++i)
    {
        event.append(pid_[i], status_[i], 0, 0,
                     px_[i], py_[i], pz_[i], energy_[i], mass_[i]);
    }
}

// High-level wrapper around PythiaCascade for use from chromo.
// Owns the PythiaCascade internally, manages the sigmaSetuphN→nextColl
// sequence, and returns particle data as pre-copied numpy arrays so
// there is no Event& lifetime issue on the Python side.
class PythiaCascadeForChromo
{
public:
    bool init(double eKinMin, double enhanceSDtarget, string initFile,
              bool rapidDecays, double smallTau0, bool slowDecays,
              bool listFinalOnly)
    {
        return _cascade.init(eKinMin, enhanceSDtarget, initFile,
                             rapidDecays, smallTau0, slowDecays, listFinalOnly);
    }

    // Attempt one hadron+nucleus inelastic collision.
    // Returns py::none() when the projectile cannot interact (low energy,
    // unknown species, or empty event generated).  Otherwise returns a
    // tuple (pid, status, px, py, pz, en, m, vx, vy, vz, vt, mothers,
    // daughters) of numpy arrays.  These are views into the internal Event
    // buffer and will be overwritten by the next call; callers must copy
    // the arrays if they need to persist across calls.
    py::object next_coll(int id, double px, double py, double pz, double e,
                         double m, int Z, int A)
    {
        Vec4 pVec(px, py, pz, e);
        if (!_cascade.sigmaSetuphN(id, pVec, m))
            return py::none();
        Event &ev = _cascade.nextColl(Z, A);
        if (ev.size() == 0)
            return py::none();
        return _extract(ev);
    }

    // Returns the inelastic h+A cross section (mb) for the given projectile
    // and lab-frame four-momentum.  Returns 0 if not available.
    double sigma_hA(int id, double px, double py, double pz, double e,
                    double m, int A)
    {
        Vec4 pVec(px, py, pz, e);
        if (!_cascade.sigmaSetuphN(id, pVec, m))
            return 0.0;
        return _cascade.sigmahA(A);
    }

    int  n_collisions()        { return _cascade.nCollisions(); }
    int  first_collision_code(){ return _cascade.firstCollisionCode(); }
    void stat()                { _cascade.stat(); }

    ParticleData &particle_data() { return _cascade.particleData(); }
    
    // Set mayDecay flag for a particle in both internal Pythia instances.
    void set_may_decay(int pdgid, bool may_decay) {
        private_access::member<PythiaCascade_pythiaMain>(_cascade).particleData.mayDecay(pdgid, may_decay);
        private_access::member<PythiaCascade_pythiaColl>(_cascade).particleData.mayDecay(pdgid, may_decay);
    }

    // Get RNG state from both internal Pythia instances.
    py::dict getRndmState() {
        auto &main = private_access::member<PythiaCascade_pythiaMain>(_cascade);
        auto &coll = private_access::member<PythiaCascade_pythiaColl>(_cascade);
        py::dict result;
        result["main"] = _getRndmStateFrom(main);
        result["coll"] = _getRndmStateFrom(coll);
        return result;
    }

    // Set RNG state on both internal Pythia instances.
    void setRndmState(py::dict state_dict) {
        auto &main = private_access::member<PythiaCascade_pythiaMain>(_cascade);
        auto &coll = private_access::member<PythiaCascade_pythiaColl>(_cascade);
        _setRndmStateOn(main, state_dict["main"].cast<py::dict>());
        _setRndmStateOn(coll, state_dict["coll"].cast<py::dict>());
    }

    py::array_t<float> charge(py::object result)
    {
        // result is the tuple returned by next_coll; pid is element 0
        auto pid = result.attr("__getitem__")(0).cast<py::array_t<int>>();
        auto &pd = _cascade.particleData();
        int size = pid.size();
        py::array_t<float> out(size);
        float *ptr = out.mutable_data();
        auto pid_ = pid.unchecked<1>();
        for (int i = 0; i < size; ++i)
            ptr[i] = charge_from_pid(pd, pid_[i]);
        return out;
    }

      
private:
    PythiaCascade _cascade;

    static py::dict _getRndmStateFrom(Pythia &p) {
        RndmState state = p.rndm.getState();
        py::dict result;
        result["i97"] = state.i97;
        result["j97"] = state.j97;
        result["seed"] = state.seed;
        result["sequence"] = state.sequence;
        result["c"] = state.c;
        result["cd"] = state.cd;
        result["cm"] = state.cm;
        py::list u_array;
        for (int i = 0; i < 97; i++)
            u_array.append(state.u[i]);
        result["u"] = u_array;
        return result;
    }

    static void _setRndmStateOn(Pythia &p, py::dict state_dict) {
        RndmState state;
        state.i97 = state_dict["i97"].cast<int>();
        state.j97 = state_dict["j97"].cast<int>();
        state.seed = state_dict["seed"].cast<int>();
        state.sequence = state_dict["sequence"].cast<long>();
        state.c = state_dict["c"].cast<double>();
        state.cd = state_dict["cd"].cast<double>();
        state.cm = state_dict["cm"].cast<double>();
        py::list u_list = state_dict["u"].cast<py::list>();
        for (int i = 0; i < 97; i++)
            state.u[i] = u_list[i].cast<double>();
        p.rndm.setState(state);
    }

    // Reuses the event_array<> helpers already defined in this file.
    static py::tuple _extract(Event &ev)
    {
        return py::make_tuple(
            event_array<Particle_idSave>(ev),
            event_status(ev),
            event_array_p<Vec4_xx>(ev),
            event_array_p<Vec4_yy>(ev),
            event_array_p<Vec4_zz>(ev),
            event_array_p<Vec4_tt>(ev),
            event_array<Particle_mSave>(ev),
            event_array_v<Vec4_xx>(ev),
            event_array_v<Vec4_yy>(ev),
            event_array_v<Vec4_zz>(ev),
            event_array_v<Vec4_tt>(ev),
            event_array_mothers(ev),
            event_array_daughters(ev));
    }
};

PYBIND11_MODULE(_pythia8, m)
{
    py::class_<ParticleData>(m, "ParticleData")
        .def("mayDecay", py::overload_cast<int, bool>(&ParticleData::mayDecay))
        .def("addParticle", py::overload_cast<int, string, string, int, int, int, double, double, double, double, double, bool>(&ParticleData::addParticle), "pdgid"_a, "name"_a, "antiname"_a, "spinType"_a = 0, "chargeType"_a = 0, "colType"_a = 0, "m0"_a = 0, "mWidth"_a = 0, "mMin"_a = 0, "mMax"_a = 0, "tau0"_a = 0, "varWidth"_a = false)
        .def("isParticle", &ParticleData::isParticle)
        .def("findParticle", py::overload_cast<int>(&ParticleData::findParticle))
        // TODO better style is .def("__iter__", ...) and return make_value_iterator
        // here, but this leads to segfault which I cannot fix right now
        .def("all",
             [](ParticleData &self)
             {
                 py::list pl;
                 for (auto p : self)
                     pl.append(p.second);
                 return pl;
             });

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

    py::class_<HIInfo>(m, "HIInfo")
        .def_property_readonly("nPartProj", &HIInfo::nPartProj)
        .def_property_readonly("nPartTarg", &HIInfo::nPartTarg)
        .def_property_readonly("b", &HIInfo::b)
        .def("glauberTot",   &HIInfo::glauberTot)
        .def("glauberINEL",  &HIInfo::glauberINEL)
        .def("glauberEL",    &HIInfo::glauberEL)
        .def("glauberND",    &HIInfo::glauberND)
        .def("glauberReset", &HIInfo::glauberReset)

        ;

    py::class_<SigmaTotal>(m, "SigmaTotal")
        .def("calc", &SigmaTotal::calc)
        .def_property_readonly("sigmaTot", &SigmaTotal::sigmaTot)
        .def_property_readonly("sigmaEl", &SigmaTotal::sigmaEl)
        .def_property_readonly("sigmaXB", &SigmaTotal::sigmaXB)
        .def_property_readonly("sigmaAX", &SigmaTotal::sigmaAX)
        .def_property_readonly("sigmaXX", &SigmaTotal::sigmaXX)
        .def_property_readonly("sigmaAXB", &SigmaTotal::sigmaAXB)
        .def_property_readonly("sigmaND", &SigmaTotal::sigmaND)
        .def_property_readonly("rho", &SigmaTotal::rho)

        ;

    py::class_<Info>(m, "Info")
        .def_property_readonly("hiInfo", [](Info &self)
                               { return self.hiInfo; })
        .def_property_readonly("sigmaTot", [](Info &self)
                               { return self.sigmaTotPtr; })
        .def_property_readonly("isDiffractiveA", &Info::isDiffractiveA)
        .def_property_readonly("isDiffractiveB", &Info::isDiffractiveB)
        .def_property_readonly("isDiffractiveC", &Info::isDiffractiveC)
        .def_property_readonly("isNonDiffractive", &Info::isNonDiffractive)

        ;

    py::class_<Settings>(m, "Settings")
        .def("resetAll", &Settings::resetAll)
        .def("listChanged", [](Settings &self) {
            py::scoped_ostream_redirect stream(
                std::cout,
                py::module_::import("sys").attr("stdout")
            );
            self.listChanged();
        })
        .def("listAll", [](Settings &self) {
            py::scoped_ostream_redirect stream(
                std::cout,
                py::module_::import("sys").attr("stdout")
            );
            self.listAll();
        })
        .def("writeFile", py::overload_cast<string, bool>(&Settings::writeFile),
             "toFile"_a, "writeAll"_a = false)

        ;

    py::class_<Pythia>(m, "Pythia")
        .def(py::init<string, bool>(), py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
        .def("init", &Pythia::init, py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
        .def("next", py::overload_cast<>(&Pythia::next))
        .def("readString", &Pythia::readString, "line"_a, "warn"_a = true, "subrun"_a = SUBRUNDEFAULT)
        .def("forceHadronLevel", &Pythia::forceHadronLevel, "find_junctions"_a = true)
        .def_readwrite("particleData", &Pythia::particleData)
        .def_readwrite("settings", &Pythia::settings)
        .def_readwrite("event", &Pythia::event)
        .def_property_readonly("info", [](Pythia &self)
                               { return self.info; })
        .def("getRndmState", [](Pythia &self) -> py::dict {
            // Get the internal RNG state from Pythia8
            RndmState state = self.rndm.getState();
            py::dict result;
            result["i97"] = state.i97;
            result["j97"] = state.j97;
            result["seed"] = state.seed;
            result["sequence"] = state.sequence;
            result["c"] = state.c;
            result["cd"] = state.cd;
            result["cm"] = state.cm;
            
            // Convert the u array to a Python list
            py::list u_array;
            for (int i = 0; i < 97; i++) {
                u_array.append(state.u[i]);
            }
            result["u"] = u_array;
            
            return result;
        })
        .def("setRndmState", [](Pythia &self, py::dict state_dict) {
            // Set the internal RNG state in Pythia8
            RndmState state;
            state.i97 = state_dict["i97"].cast<int>();
            state.j97 = state_dict["j97"].cast<int>();
            state.seed = state_dict["seed"].cast<int>();
            state.sequence = state_dict["sequence"].cast<long>();
            state.c = state_dict["c"].cast<double>();
            state.cd = state_dict["cd"].cast<double>();
            state.cm = state_dict["cm"].cast<double>();
            
            // Convert Python list back to array
            py::list u_list = state_dict["u"].cast<py::list>();
            for (int i = 0; i < 97; i++) {
                state.u[i] = u_list[i].cast<double>();
            }
            
            self.rndm.setState(state);
        })
        .def("charge",
             [](Pythia &self)
             {
                 // skip first pseudoparticle
                 int size = self.event.size() - 1;
                 py::array_t<float> result(size);
                 float *ptr = result.mutable_data();
                 for (auto pit = self.event.begin() + 1; pit != self.event.end(); ++pit)
                     *ptr++ = charge_from_pid(self.particleData, pit->id());
                 return result;
             })
        .def("readFile",
             py::overload_cast<string, bool, int>(&Pythia::readFile),
             "fileName"_a, "warn"_a = true, "subrun"_a = SUBRUNDEFAULT,
             py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
        .def("setBeamIDs", &Pythia::setBeamIDs, "idA"_a, "idB"_a = 0)
        .def("setKinematics", py::overload_cast<double>(&Pythia::setKinematics), "eCM"_a);

    py::class_<Event>(m, "Event")
        .def_property_readonly("size", [](Event &self)
                               { return self.size() - 1; })
        .def("pid", event_array<Particle_idSave>)
        .def("status", event_status)
        .def("m", event_array<Particle_mSave>)
        .def("px", event_array_p<Vec4_xx>)
        .def("py", event_array_p<Vec4_yy>)
        .def("pz", event_array_p<Vec4_zz>)
        .def("en", event_array_p<Vec4_tt>)
        .def("vx", event_array_v<Vec4_xx>)
        .def("vy", event_array_v<Vec4_yy>)
        .def("vz", event_array_v<Vec4_zz>)
        .def("vt", event_array_v<Vec4_tt>)
        .def("mothers", event_array_mothers)
        .def("daughters", event_array_daughters)
        .def("reset", &Event::reset)
        .def("list", py::overload_cast<bool, bool, int>(&Event::list, py::const_), "showScaleAndVertex"_a = false, "showMothersAndDaughters"_a = false, "precision"_a = 3)
        .def("append", py::overload_cast<int, int, int, int, double, double, double, double, double, double, double>(&Event::append), "pdgid"_a, "status"_a, "col"_a, "acol"_a, "px"_a, "py"_a, "pz"_a, "e"_a, "m"_a = 0, "scale"_a = 0, "pol"_a = 9.)
        .def("fill", [](Event &self, py::array_t<int> pid, py::array_t<int> status, py::array_t<double> px, py::array_t<double> py, py::array_t<double> pz, py::array_t<double> energy, py::array_t<double> mass)
             { fill(self, pid, status, px, py, pz, energy, mass); }, "pid"_a, "status"_a, "px"_a, "py"_a, "pz"_a, "energy"_a, "mass"_a);

    py::class_<PythiaCascadeForChromo>(m, "PythiaCascadeForChromo")
        .def(py::init<>())
        .def("init", &PythiaCascadeForChromo::init,
             py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
             "eKinMin"_a = 0.3, "enhanceSDtarget"_a = 0.5,
             "initFile"_a = "../share/Pythia8/setups/InitDefaultMPI.cmnd",
             "rapidDecays"_a = false, "smallTau0"_a = 1e-10,
             "slowDecays"_a = true, "listFinalOnly"_a = false)
        .def("next_coll", &PythiaCascadeForChromo::next_coll,
             "id"_a, "px"_a, "py"_a, "pz"_a, "e"_a, "m"_a, "Z"_a, "A"_a)
        .def("sigma_hA", &PythiaCascadeForChromo::sigma_hA,
             "id"_a, "px"_a, "py"_a, "pz"_a, "e"_a, "m"_a, "A"_a)
        .def("n_collisions",        &PythiaCascadeForChromo::n_collisions)
        .def("first_collision_code",&PythiaCascadeForChromo::first_collision_code)
        .def("stat",                &PythiaCascadeForChromo::stat)
        .def("particle_data", &PythiaCascadeForChromo::particle_data,
             py::return_value_policy::reference_internal)
        .def("set_may_decay", &PythiaCascadeForChromo::set_may_decay,
             "pdgid"_a, "may_decay"_a)
        .def("charge", &PythiaCascadeForChromo::charge, "result"_a)
        .def("getRndmState", &PythiaCascadeForChromo::getRndmState)
        .def("setRndmState", &PythiaCascadeForChromo::setRndmState);
}