## Master Script to Test All Builds and Run 3D Live Coding Demo

To start everything using a single command:

```bash
bash scripts/e2e-test-builds.bash
```
This will git clone the repo again in a temporary folder, builds, and run the 3D live coding demo on a local server. (only uses git-pushed changes). (Fetches the files first time, but will reuse them).

Then, open your browser at the location shown, for example:
`http://11.111.11.111:8000/mp5_json_code.html`

This script is refined to run on a remote dev machine (ssh via terminal). But originally designed to run on MacOS too.

If you use a remote ssh-linux setup: Also, to access the local web-server in this remote setup, the you will need to open the firewall `sudo ufw allow ...` to allow incoming http on port 8000.

## Prerequisites:
(May require manual steps if not installed)
* Docker
* python3
* git, bash, Linux or MacOS.
* +
    * Only E2E live demo: Linux/ssh remote setup:
       * ufw (if remote setup on Linux/ssh )
    * Only python route, only on a MacOS local machine: Not needed `e2e-test-builds.bash` route.
       * `brew` (if MacOS, for some python packages).

## Python Development Setup:

For Python development (not the above route): The following packages need a manual instalation ebfore the automated pip installtions; Their installation is non-trivial (not just pip): May need brew on MacOS. See `sandbox/sympy-experiment/run_script.bash`
* mayavi
* qt5
* vtk
Other packages will be automatically installed automatically via pip on venv.

## Other routes

* There exists a NodeJS target too (for notejs tests, npm package publish, and cli-style 3D live-coding in nodejs/npm/yarn setup).

* The Browser-based 3D designer: https://github.com/sohale/mp5-private

* The standalone slicer: https://github.com/sohale/mp5slicer

* The full-platform ( aka WeDesign.Live ), including the live-collaborateive 3D: See https://github.com/sohale/mp5-private
