# MIproject
A collection of programs and scripts to solve and analyze the Mott-Hubbard transition in a variety of dynamical mean-field settings.

## Requirements
The project relies on several inter-dependent[^1] external libraries:

[^1]: **Please build them in the exact order reported here.** More specifically: DMFT-tools depends on SciFortran, DMFT-ED depends on SciFortran, CDMFT-ED depends on both SciFortran and DMFT-tools. DMFT-LAB depends on everything.

- [SciFortran](dependencies/scifor)
- [DMFT-tools](dependencies/dmft-tools)
- [DMFT-ED](dependencies/dmft-ed)
- [CDMFT-ED](dependencies/cdmft-ed)
- [DMFT-LAB](dependencies/dmft-lab)

These dependencies are handled through git-submodule embedding. To clone the project shipping also all the correct versions of such libraries, just run:

```
git clone --recursive https://github.com/bellomia/MIproject.git MIproject
```

To pull upstream changes of one of the submodules you can run:

```
git submodule update --remote --merge <submodule-name>
```

To update instead all the submodules, use the `foreach` command:

```
git submodule foreach `git pull origin`
git submodule update
```

Complete documentation of the `git submodule` tools can be found [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

**Warning:** The build and installation of the libraries is not automatized (yet?), thus you will need to enter each submodule directory and follow the provided instructions. All the upstream requirements have to be met, even the "optional" ones (e.g. MPI related stuff).
