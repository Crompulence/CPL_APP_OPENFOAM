# Example uses of pytest

These two folders give an example of the two possible ways to test coupled simulations, 
using Python to drive the testing of OpenFOAM with [pytest][https://docs.pytest.org/en/latest/].
The two methods are:

    1) coupled_to_pytest - Running an instance of pytest AND an instance of OpenFOAM with information exchange.
    2) pytest_runs_subprocess - Running a range of OpenFOAM cases each linked to different dummy scripts.

## Why Two Methods?

The use of MPI requires jobs to be started with mpiexec, specifying the number of processes to be used (e.g. 20) as follows,

```bash
    mpiexec -n  20 ./some_job_name
```
Whereas, pytest typically works by finding python scripts starting with "test" and running then,
```bash
pytest test_case.py
```
However, a pytest cannot become part of a parallel run as it was not started with mpiexec.

Instead, a parallel run has to be started with pytest as part of the run,
```bash
    mpiexec -n  20 ./some_job_name : -n 1 pytest test_case.py
```
which is option 1), or we need to start an mpiexec instance using subprocess inside a pytest,
```python
    subprocess.Popen("mpiexec -n  20 ./some_job_name : -n 1 python coupled_run.py", shell=True)
```
which is option 2), run with something like
```bash
pytest test_subprocess.py
```
and this creates a range of runs.

## Which One to Use?
Which one to use depends on the nature of the development you'd like to test. They both have advantages and disadvantages.

By using a direct coupling with pytest and OpenFOAM, as in option 1), the effect of injecting controlled information into OpenFOAM may be used to test fundamental behaviour.
As an example, a given force field or boundary condition could be sent from the dummy script using CPL_send. 
The outcome from CPL_recv is obtained directly in the pytest code and can then be tested against a range of assert statements.
This allows any variable in OpenFOAM (provided it can be packaged as a 4D array) to be sent directly to the pytest code and checked.
The downside of this approach, is that only information that can be specified during a single run can be tested, with no changes to initial conditions.
In addition, the OpenFOAM run has persistent information, that is, after a send/recv statement, the OpenFOAM solver will potentially have changed in some way because of the sent information.
The next test must take this into account.
 The example provided in "coupled_to_pytest" is the extreme example of using persistence of state, where the test is simply checking each timestep that thesolver evolves in the expected way (vs. analytical solution) when no information is sent.

In option 2), a far more flexible testing method is put forward, using multiple instances of OpenFOAM and a dummy MD script to be started by pytest as subprocesses. 
This allows the dummy MD script, the input for OpenFOAM and the processor topology to be adjusted and a wide range of cases explored.
The parameter space is navigated using [SimWrapLib][https://github.com/edwardsmith999/SimWrapPy] to change script and OpenFOAM input values.

The disadvantage of this approach is the coupled run is a subprocess so the actual driving pytest does not have any access to information generated during the simulation. 
Writing to disk and reading at the end (or during) the run is the only way to assert the coupled run is working as expected.

