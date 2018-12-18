from ._solver_cpp import _SolverCPP
from ..objects import Fiber


class Solver(_SolverCPP):

    def get_fibers(self):
        fiber_bundles = _SolverCPP.get_fibers(self)
        for fb in fiber_bundles:
            for f in fb:
                f = Fiber(f.points, f.radii)

        return fiber_bundles
