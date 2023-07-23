'''
Define the computational resources we may wish to use
'''
import nexus

def apply_settings(generate_only=False):
    nexus.settings(
        pseudo_dir = 'pplib',
        sleep = 10,
        generate_only = generate_only,
        status_only = 0,
        machine = 'summit',
        account = 'mat221',
        command_line = False
        )

def pw_job(pw_nodes, minutes):
    return nexus.job(
        minutes = minutes,
        nodes = pw_nodes,
        app = '/gpfs/alpine/mat221/proj-shared/yeluo_qmcpack/qmcpack/external_codes/quantum_espresso/qe-6.8/build/bin/pw.x',
        presub = 'module load gcc/9.3.0 spectrum-mpi cuda/11.0.3 essl netlib-lapack netlib-scalapack hdf5/1.10.7 fftw; module use /gpfs/alpine/mat221/proj-shared/compile/llvm_recipes/olcf_summit/modules; module load llvm/main-20210811',
        app_options = '-nk 6',
        alloc_flags = 'smt1'
        )

def p2q_job(p2q_nodes, minutes):
    return nexus.job(
        minutes = minutes,
        nodes = p2q_nodes,
        app = '/gpfs/alpine/mat221/proj-shared/yeluo_qmcpack/qmcpack/external_codes/quantum_espresso/qe-6.8/build/bin/pw2qmcpack.x',
        presub = 'module load gcc/9.3.0 spectrum-mpi cuda/11.0.3 essl netlib-lapack netlib-scalapack hdf5/1.10.7 fftw; module use /gpfs/alpine/mat221/proj-shared/compile/llvm_recipes/olcf_summit/modules; module load llvm/main-20210811',
        alloc_flags = 'smt1'
        )

def qmc_job(qmc_nodes, minutes):
    return nexus.job(
        minutes = minutes,
        nodes = qmc_nodes,
        app = '/gpfs/alpine/mat221/proj-shared/mick_graphene/bin/qmcpack-3.9.2_newsystem',
        presub = 'module load gcc/9.3.0 spectrum-mpi cuda/11.0.3 essl netlib-lapack hdf5/1.10.7 fftw; module use /gpfs/alpine/mat221/proj-shared/compile/llvm_recipes/olcf_summit/modules; module load llvm/main-20210811',
        run_options = f'-b rs -c 7 -g 1 -d packed -n {qmc_nodes*6} -r 6 -a 1 --smpiargs="-disable_gpu_hooks"',
        threads = 7,
        alloc_flags = 'smt1'
        )
