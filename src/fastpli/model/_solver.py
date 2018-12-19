from ._solver_cpp import _SolverCPP
from ..objects import Fiber


class Solver(_SolverCPP):

    def get_fiber_bundles(self):
        # casting FiberCPP to Fiber
        fiber_bundles = super().get_fiber_bundles()
        for i in range(len(fiber_bundles)):
            for j in range(len(fiber_bundles[i])):
                fiber_bundles[i][j] = Fiber(
                    fiber_bundles[i][j].points,
                    fiber_bundles[i][j].radii)
        return fiber_bundles

    @property
    def fiber_bundles(self):
        return self.get_fiber_bundles()

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        super().set_fiber_bundles(fbs)

    @property
    def parameters(self):
        return super().get_parameters()

    @parameters.setter
    def parameters(self, tuple):
        super().set_parameters(tuple[0], tuple[1], tuple[2])
