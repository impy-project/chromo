#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

PYBIND11_MODULE(_util, m)
{
  // def _select_mothers(mask, mothers):
  //     # attach parentless particles to beam particles,
  //     # unless those are also removed
  //     fallback = (-1, -1)
  //     if mask[0] and mask[1]:
  //         fallback = (0, 1)

  //     n = len(mothers)
  //     indices = np.arange(n)[mask]
  //     result = mothers[mask]
  //     mapping = {old: i for i, old in enumerate(indices)}

  //     n = len(result)
  //     for i in range(n):
  //         a = result[i, 0]
  //         if a == -1:
  //             continue
  //         p = mapping.get(a, -1)
  //         if p == -1:
  //             a, b = fallback
  //             result[i, 0] = a
  //             result[i, 1] = b
  //         elif p != a:
  //             q = -1
  //             b = result[i, 1]
  //             if b > -1:
  //                 q = mapping.get(b, -1)
  //             result[i, 0] = p
  //             result[i, 1] = q
  //     return result

  m.def("select_mothers", )
}