#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include "HepMC/FourVector.h"
#include "HepMC/Data/SmartPointer.h"
#include "HepMC/GenHeavyIon.h"
#include "HepMC/GenPdfInfo.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/GenRunInfo.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Units.h"
#include "HepMC/Print.h"
#include "HepMC/WriterAscii.h"
#include "HepMC/ReaderAscii.h"
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

PYBIND11_DECLARE_HOLDER_TYPE(T, HepMC::SmartPointer<T>);

// need to customize getter for SmartPointer
namespace pybind11 { namespace detail {
    template <typename T>
    struct holder_helper<HepMC::SmartPointer<T>> { // <-- specialization
        static const T *get(const HepMC::SmartPointer<T> &p) { return &(*p); }
    };
}}

namespace HepMC {

using GenRunInfoPtr = std::shared_ptr<GenRunInfo>;

// equality comparions used by unit tests
bool operator==(const GenRunInfo::ToolInfo& a,
                const GenRunInfo::ToolInfo& b) {
    return a.name == b.name &&
        a.version == b.version &&
        a.description == b.description;
}

bool is_close(const FourVector& a, const FourVector& b,
              double rel_eps = 1e-7) {
    auto is_close = [rel_eps](double a, double b) {
            return std::abs(a-b) < rel_eps;
        };
    return is_close(a.x(), b.x()) &&
        is_close(a.y(), b.y()) &&
        is_close(a.z(), b.z()) &&
        is_close(a.t(), b.t());
}

bool operator==(const GenParticle& a, const GenParticle& b) {
    return a.id() == b.id() && a.pid() == b.pid() &&
        a.status() == b.status() && is_close(a.momentum(), b.momentum());
}

bool operator==(const GenVertex& a, const GenVertex& b) {
    auto equal_id = [](const GenParticlePtr& a,
                       const GenParticlePtr& b) { return a->id() == b->id(); };
    return a.id() == b.id() && a.status() == b.status() &&
        is_close(a.position(), b.position()) &&
        std::equal(a.particles_in().begin(), a.particles_in().end(),
                   b.particles_in().begin(), b.particles_in().end(),
                   equal_id) &&
        std::equal(a.particles_out().begin(), a.particles_out().end(),
                   b.particles_out().begin(), b.particles_out().end(),
                   equal_id);
}

bool operator==(const std::vector<GenVertexPtr>& a,
                const std::vector<GenVertexPtr>& b)
{
    if (a.size() != b.size())
        return false;
    int i = 0;
    for (; i < a.size() && *a[i] == *b[i]; ++i);
    return i == a.size();
}

bool operator==(const std::vector<GenParticlePtr>& a,
                const std::vector<GenParticlePtr>& b)
{
    if (a.size() != b.size())
        return false;
    int i = 0;
    for (; i < a.size() && *a[i] == *b[i]; ++i);
    return i == a.size();
}

bool operator==(const GenEvent& a, const GenEvent& b) {
    // incomplete:
    // missing comparison of GenHeavyIon, GenPdfInfo, GenCrossSection

    if (a.event_number() != b.event_number() ||
        a.momentum_unit() != b.momentum_unit() ||
        a.length_unit() != b.length_unit())
        return false;

    auto a_tools = std::vector<GenRunInfo::ToolInfo>();
    auto b_tools = std::vector<GenRunInfo::ToolInfo>();
    if (a.run_info())
        a_tools = a.run_info()->tools();
    if (b.run_info())
        b_tools = b.run_info()->tools();
    if (!std::equal(a_tools.begin(), a_tools.end(),
                    b_tools.begin(), b_tools.end()))
        return false;

    return a.vertices() == b.vertices() && a.particles() == b.particles();
}

} // namespace HepMC

namespace py = pybind11;
using namespace py::literals;

#define FUNC(name) m.def(#name, name)
#define PROP_RO(name, cls) .def_property_readonly(#name, &cls::name)
#define PROP(name, cls) .def_property(#name, &cls::name, &cls::set_##name)
#define METH(name, cls) .def(#name, &cls::name)
#define METH_OL(name, cls, args) .def(#name, py::overload_cast<args>(&cls::name))

PYBIND11_MODULE(pyhepmc, m) {
    using namespace HepMC;

    m.doc() = R"pbdoc(
        pyhepmc plugin
        --------------
        .. currentmodule:: pyhepmc
        .. autosummary::
           :toctree: _generate
           Units
           FourVector
           GenEvent
           GenParticle
           GenVertex
           print_content
           print_listing
    )pbdoc";

    // m.def("dummy", []() { return 0; }, R"pbdoc(
    //     Return zero
    //     Some other explanation about the dummy function.
    // )pbdoc");

    py::class_<Units> clsUnits(m, "Units");

    py::enum_<Units::MomentumUnit>(clsUnits, "MomentumUnit")
        .value("MEV", Units::MomentumUnit::MEV)
        .value("GEV", Units::MomentumUnit::GEV)
        .export_values();

    py::enum_<Units::LengthUnit>(clsUnits, "LengthUnit")
        .value("CM", Units::LengthUnit::CM)
        .value("MM", Units::LengthUnit::MM)
        .export_values();

    py::class_<FourVector>(m, "FourVector")
        .def(py::init<>())
        .def(py::init<double, double, double, double>(),
             "x"_a, "y"_a, "z"_a, "t"_a)
        .def(py::init([](py::sequence seq) {
                const double x = py::cast<double>(seq[0]);
                const double y = py::cast<double>(seq[1]);
                const double z = py::cast<double>(seq[2]);
                const double t = py::cast<double>(seq[3]);
                return new FourVector(x, y, z, t);
            }))
        .def_property("x", &FourVector::x, &FourVector::setX)
        .def_property("y", &FourVector::y, &FourVector::setY)
        .def_property("z", &FourVector::z, &FourVector::setZ)
        .def_property("t", &FourVector::t, &FourVector::setT)
        .def_property("px", &FourVector::px, &FourVector::setPx)
        .def_property("py", &FourVector::py, &FourVector::setPy)
        .def_property("pz", &FourVector::pz, &FourVector::setPz)
        .def_property("e", &FourVector::e, &FourVector::setE)
        METH(length2, FourVector)
        METH(length, FourVector)
        METH(perp2, FourVector)
        METH(perp, FourVector)
        METH(interval, FourVector)
        METH(pt, FourVector)
        METH(m2, FourVector)
        METH(m, FourVector)
        METH(phi, FourVector)
        METH(theta, FourVector)
        METH(eta, FourVector)
        METH(rap, FourVector)
        METH(abs_eta, FourVector)
        METH(abs_rap, FourVector)
        METH(is_zero, FourVector)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())
        .def(py::self /= double())
        ;

    py::implicitly_convertible<py::sequence, FourVector>();

    FUNC(delta_phi);
    FUNC(delta_eta);
    FUNC(delta_rap);
    FUNC(delta_r2_eta);
    FUNC(delta_r_eta);
    FUNC(delta_r2_rap);
    FUNC(delta_r_rap);

    py::class_<GenRunInfo, GenRunInfoPtr> clsGenRunInfo(m, "GenRunInfo");
    clsGenRunInfo
        .def(py::init<>())
        .def_property("tools",
                (std::vector<GenRunInfo::ToolInfo>&(GenRunInfo::*)()) &GenRunInfo::tools,
                [](GenRunInfo& self, py::sequence seq) {
                    self.tools() = py::cast<std::vector<GenRunInfo::ToolInfo>>(seq);
                }
            )
        ;

    py::class_<GenRunInfo::ToolInfo>(clsGenRunInfo, "ToolInfo")
        .def(py::init<std::string, std::string, std::string>(),
             "name"_a, "version"_a, "description"_a)
        .def(py::init([](py::sequence seq) {
                return new GenRunInfo::ToolInfo({
                    py::cast<std::string>(seq[0]),
                    py::cast<std::string>(seq[1]),
                    py::cast<std::string>(seq[2])
                });
            }))
        .def(py::self == py::self)
        ;

    py::implicitly_convertible<py::sequence, GenRunInfo::ToolInfo>();

    py::class_<GenHeavyIon, GenHeavyIonPtr>(m, "GenHeavyIon")
        .def(py::init([](int nh, int np, int nt, int nc, int ns, int nsp,
                         int nnw=0, int nwn=0, int nwnw=0,
                         float im=0., float pl=0., float ec=0., float s=0., float cent=0.) {
                auto x = GenHeavyIonPtr(new GenHeavyIon());
                x->set(nh, np, nt, nc, ns, nsp, nnw, nwn, nwnw, im, pl, ec, s, cent);
                return x;
            }),
            "n_coll_hard"_a,
            "n_part_proj"_a,
            "n_part_targ"_a,
            "n_coll"_a,
            "n_spec_neut"_a,
            "n_spec_prot"_a,
            "nnw"_a = 0, "nwn"_a = 0, "nwnw"_a = 0,
            "impact_parameter"_a = 0.f, "event_plane_angle"_a = 0.f,
            "eccentricity"_a = 0.f, "sigma_inel_NN"_a = 0.f,
            "centrality"_a = 0.f)
        ;

    py::class_<GenEvent>(m, "GenEvent")
        .def(py::init<std::shared_ptr<GenRunInfo>, Units::MomentumUnit, Units::LengthUnit>(),
             "run"_a, "momentum_unit"_a = Units::GEV, "length_unit"_a = Units::MM)
        .def(py::init<Units::MomentumUnit, Units::LengthUnit>(),
             "momentum_unit"_a = Units::GEV, "length_unit"_a = Units::MM)
        .def_property("particles",
                (std::vector<GenParticlePtr>&(GenEvent::*)()) &GenEvent::particles,
                [](GenEvent& self, py::sequence seq) {
                    self.particles() = py::cast<std::vector<GenParticlePtr>>(seq);;
                }
            )
        .def_property("vertices",
                (std::vector<GenVertexPtr>&(GenEvent::*)()) &GenEvent::vertices,
                [](GenEvent& self, py::sequence seq) {
                    self.vertices() = py::cast<std::vector<GenVertexPtr>>(seq);
                }
            )
        PROP(event_number, GenEvent)
        PROP_RO(momentum_unit, GenEvent)
        PROP_RO(length_unit, GenEvent)
        METH(set_units, GenEvent)
        .def_property("heavy_ion", &GenEvent::heavy_ion, (void (GenEvent::*)(const GenHeavyIonPtr&))&GenEvent::set_heavy_ion)
        .def_property("pdf_info", &GenEvent::pdf_info, (void (GenEvent::*)(const GenPdfInfoPtr&))&GenEvent::set_pdf_info)
        .def_property("cross_section", &GenEvent::cross_section, (void (GenEvent::*)(const GenCrossSectionPtr&))&GenEvent::set_cross_section)
        METH(event_pos, GenEvent)
        METH(beams, GenEvent)
        METH_OL(add_vertex, GenEvent, GenVertexPtr)
        METH_OL(add_particle, GenEvent, GenParticlePtr)
        METH_OL(remove_vertex, GenEvent, GenVertexPtr)
        METH_OL(remove_particle, GenEvent, GenParticlePtr)
        .def("reserve", &GenEvent::reserve, "particles"_a, "vertices"_a = 0)
        METH(clear, GenEvent)
        .def(py::self == py::self)
        ;

    py::class_<GenParticle, GenParticlePtr>(m, "GenParticle")
        .def(py::init<const FourVector&, int, int>(),
             "momentum"_a = FourVector::ZERO_VECTOR(), "pid"_a = 0, "status"_a = 0)
        PROP_RO(in_event, GenParticle)
        PROP_RO(parent_event, GenParticle)
        PROP_RO(id, GenParticle)
        // PROP_RO(data, GenParticle)
        .def_property_readonly("production_vertex", (GenVertexPtr (GenParticle::*)()) &GenParticle::production_vertex)
        .def_property_readonly("end_vertex", (GenVertexPtr (GenParticle::*)()) &GenParticle::end_vertex)
        PROP_RO(parents, GenParticle)
        PROP_RO(children, GenParticle)
        PROP_RO(ancestors, GenParticle)
        PROP_RO(descendants, GenParticle)
        PROP(pid, GenParticle)
        PROP(status, GenParticle)
        PROP(momentum, GenParticle)
        METH(is_generated_mass_set, GenParticle)
        METH(generated_mass, GenParticle)
        METH(set_generated_mass, GenParticle)
        METH(unset_generated_mass, GenParticle)
        .def(py::self == py::self)
        ;

    py::class_<GenVertex, GenVertexPtr>(m, "GenVertex")
        .def(py::init<const FourVector&>(),
             "position"_a = FourVector::ZERO_VECTOR())
        PROP_RO(parent_event, GenVertex)
        PROP_RO(in_event, GenVertex)
        PROP_RO(id, GenVertex)
        PROP(status, GenVertex)
        // PROP_RO(data, GenVertex)
        METH_OL(add_particle_in, GenVertex, GenParticlePtr)
        METH_OL(add_particle_out, GenVertex, GenParticlePtr)
        METH_OL(remove_particle_in, GenVertex, GenParticlePtr)
        METH_OL(remove_particle_out, GenVertex, GenParticlePtr)
        METH(particles, GenVertex)
        PROP_RO(particles_in, GenVertex)
        PROP_RO(particles_out, GenVertex)
        PROP(position, GenVertex)
        METH(has_set_position, GenVertex)
        .def(py::self == py::self)
        ;

    // py::class_<GenParticleData>(m, "GenParticleData");
    // py::class_<GenVertexData>(m, "GenVertexData");

    // py::class_<std::ostream>("ostream");

    py::class_<std::ostringstream>(m, "ostringstream")
        .def(py::init<>())
        .def("str", (std::string (std::ostringstream::*)() const) &std::ostringstream::str)
        METH(flush, std::ostringstream)
        ;

    py::class_<WriterAscii>(m, "WriterAscii")
        .def(py::init<const std::string&, GenRunInfoPtr>(),
             "filename"_a, "run"_a = nullptr)
        .def(py::init<std::ostringstream&, GenRunInfoPtr>(),
             "ostringstream"_a, "run"_a = nullptr,
             py::keep_alive<1, 2>())
        METH(write_event, WriterAscii)
        METH(write_run_info, WriterAscii)
        METH(failed, WriterAscii)
        METH(close, WriterAscii)
        METH(set_precision, WriterAscii)
        ;

    py::class_<ReaderAscii>(m, "ReaderAscii")
        .def(py::init<const std::string>(), "filename"_a)
        METH(read_event, ReaderAscii)
        METH(failed, ReaderAscii)
        METH(close, ReaderAscii)
        ;

    m.def("print_content", [](const GenEvent& event) { Print::content(event); });
    m.def("print_listing", [](const GenEvent& event) { Print::listing(event); });

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}