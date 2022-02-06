# Using in FORTRAN (2003)

It starts with

```FORTRAN
use, intrinsic :: iso_c_binding     ! To translate C API into Fortran (possible in 2003)
include 'fftw3.f03'                 ! The definitions for translation of FFTW modules
```

An Example starting script

```FORTRAN
type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), dimension(1024,1000) :: in, out
plan = fftw_plan_dft_2d(1000,1024, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
...
call fftw_execute_dft(plan, in, out)
...
call fftw_destroy_plan(plan)
```

A few <mark>important things</mark> to keep in mind are:

- <mark>FFTW plans are type(C_PTR)</mark>. Other C types are mapped in the obvious way via the iso_c_binding standard: <mark>int turns into integer(C_INT), fftw_complex turns
into complex(C_DOUBLE_COMPLEX), double turns into real(C_DOUBLE)</mark>, and so on. See Section 7.3 [FFTW Fortran type reference], page 80.

- Functions in C become functions in Fortran if they have a return value, and subroutines
in Fortran otherwise.

- The ordering of the Fortran array dimensions must be reversed when they are passed
to the FFTW plan creation, thanks to differences in array indexing conventions (see Section 3.2 [Multi-dimensional Array Format], page 15). This is unlike the legacy
Fortran interface (see Section 8.1 [Fortran-interface routines], page 87), which reversed
the dimensions for you. See Section 7.2 [Reversing array dimensions], page 78.

- Using ordinary Fortran array declarations like this works, but may yield suboptimal
performance because the data may not be not aligned to exploit SIMD instructions on
modern proessors (see Section 3.1 [SIMD alignment and fftw malloc], page 15). Better
performance will often be obtained by allocating with ‘fftw_alloc’. See Section 7.5
[Allocating aligned memory in Fortran], page 82.

- Similar to the legacy Fortran interface (see Section 8.3 [FFTW Execution in Fortran], page 88), we currently recommend not using fftw_execute but rather using the more
specialized functions like fftw_execute_dft (see Section 4.6 [New-array Execute Func-
tions], page 38). However, you should execute the plan on the same arrays as the ones
for which you created the plan, unless you are especially careful. See Section 7.4 [Plan execution in Fortran], page 81. To prevent you from using fftw_execute by mistake, the fftw3.f03 file does not provide an fftw_execute interface declaration.

- Multiple planner flags are combined with ior (equivalent to ‘|’ in C). e.g.
FFTW_MEASURE | FFTW_DESTROY_INPUT becomes ior(FFTW_MEASURE, FFTW_DESTROY_
INPUT). (You can also use ‘+’ as long as you don’t try to include a given flag more
than once.)

## Sources to be added

- [Fortran Complex Array](https://www.physicsforums.com/threads/fortran-complex-array-assignment.592750/)