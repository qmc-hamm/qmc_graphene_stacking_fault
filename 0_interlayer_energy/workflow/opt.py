'''
Define wave function optimization strategies we may wish to use.
'''
import nexus 

def basic(samples=100000, max_loop=12):
    return [
        nexus.loop(
            max = max_loop,
            qmc = nexus.linear(
                minmethod = 'OneShiftOnly',
                minwalkers = 0.3,
                samples = samples,
                timestep = 0.8,
                substeps = 10,
                warmupsteps = 0,
                blocks = 50
                )
            )
        ]

def split(samples=20000, max_loop=12):
    return [
        nexus.loop(
            max = 6,
            qmc = nexus.linear(
                minmethod = 'OneShiftOnly',
                minwalkers = 1e-4,
                samples = 10000,
                timestep = 0.8,
                substeps = 4,
                warmupsteps = 5,
                blocks = 25,
                )
            ),
        nexus.loop(
            max = max_loop,
            qmc = nexus.linear(
                minmethod = 'OneShiftOnly',
                minwalkers = 0.5,
                samples = samples,
                timestep = 0.8,
                substeps = 4,
                warmupsteps = 2,
                blocks = 50)
            )
        ]

def strats(opt_method, samples, max_loop):
    if opt_method == 'basic':
        return basic(samples=samples, max_loop=max_loop)
    elif opt_method == 'split':
        return split(samples=samples, max_loop=max_loop)

