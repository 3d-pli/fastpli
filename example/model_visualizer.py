import fastpli

fiber_bundles = [[fastpli.objects.Fiber(
    [0, 0, 10, 1, 1, 10, 1, 2, 10], [1, 1, 0.5])]]

vis = fastpli.model.Vis()
vis.set_data(fiber_bundles)
vis.run()
