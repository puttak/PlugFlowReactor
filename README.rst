Welcome to PlugFlowReactor
==========================

This is a plug-flow reactor implemented with Cantera. Once this project reaches
its maturity (mostly debugged and documented, with surface chemistry enabled
and possibility to connect to reactor networks), it will be deprecated and I
intend to make a contribution to Cantera.

Currently only supported for Linux. Check my project CanteraPFR_ that gave
origin to this simpler and concise package if you wish to try building Cantera
under Cygwin (Windows) to be able to build the project (check specially the
script `get_cygwin.bat`). Please do not ask me about Visual Studio, I personally
abhorror Windows and already did great effort automatizing Cygwin install with
Cantera compilation.

.. _CanteraPFR: https://github.com/waltermateriais/CanteraPFR/

Build instructions
-------------------

Under Linux (or Cygwin) you will need most classical build tools such as a
C++ compiler with C++14 compatibility (any recent version of `g++`), `make`
build manager, and a `python` distribution above Python 3.6.5.

**Note:** if building the shared library with `make`, consider modifying the
values of `ROOT_CANTERA` and `ROOT_SUNDIALS` (and other variables you think
need to be changed) to the ones corresponding to your machine.

I already provide the library `libPlugFlowReactor_shared.so` with the package
as compiled for Linux-64 systems (you need Cantera 2.4.0 and Sundials 2.7 to
be in your library path or declared under `LD_LIBRARY_PATH`). If you do not
trust me (you should avoid binaries from unknown sources, anyways), you will
need to run `make` to build the project. Otherwise, all you need to do is
install the package with `python setup.py install` or just manually place it
in your `PYTHONPATH`.

Once everything is done and if I am able to make a pull request to Cantera, you
will finally have the peace of community support for other platforms than Linux.
Anyways, if you are doing chemical kinetics and/or fluid dynamics with open
source libraries, you should consider moving to Linux.

;)

Walter
