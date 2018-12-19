from ._solver_cpp import _SolverCPP
from ..objects import Fiber


class Solver(_SolverCPP):

    def get_fiber_bundles(self):
        # casting FiberCPP to Fiber
        fiber_bundles = _SolverCPP.get_fiber_bundles(self)
        for i in range(len(fiber_bundles)):
            for j in range(len(fiber_bundles[i])):
                fiber_bundles[i][j] = Fiber(
                    fiber_bundles[i][j].points,
                    fiber_bundles[i][j].radii)

        return fiber_bundles
