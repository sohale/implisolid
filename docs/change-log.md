## Notes as it progresses

#### Fri 29 Apr 2022 22:10:47 BST

Attmmpting to replace the ace-build with submodule with a specific tag

Example deleted file: `js_iteration_1/ace-builds/src-min/ace.js`
Corresponds to `https://github.com/ajaxorg/ace-builds/blob/v1.2.6/src-min/ace.js`

Branch.tag info: v1.2.6 2268d21 2268d21c5893320330c9daf10b70b1973ca45eba

Cannot use a specific tag.

Failed attempt: on `~/cs/implisolid/js_iteration_1/`: `git submodulegit submodule add -b v1.2.6 git@github.com:ajaxorg/ace-builds.git ./ace-builds`



#### Sat Apr 19 09:32:39 BST 2025
  Migrating to tagged "emsdk", new URLs for dependencies, and installation on a new place (yet another configuration) + improving docs for better clarity on how to build/deploy.

  * June 12, 2022
         * Last meaningful commit where scrips were tested working and being developed

  * April 18 2025
     * explicitly specifying docker tag
           emscripten/emsdk 3.1.14  `2022-06-20T15:58:58.450669Z`

     * Technical debt from the "latest" `emsdk` docker:
          * Latest no longer support: `DEMANGLE_SUPPORT`
          * Latest no longer support: `EXTRA_EXPORTED_RUNTIME_METHODS`
          ```txt
            `-sEXTRA_EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap']`: No longer supported, use EXPORTED_RUNTIME_METHODS
          ```

     * Dependencies:
        * Updated the URL to download Boost (older one lo longer available)
        * submodules are cloned shallow (ace, Eigen)

     * Fixed: The `assert_env_nonempty` needs `"`

     * documented env variables need to run E2E
     ```bash
      export DEPLOY_LOCATION=/dataneura/implisolid/build/b2/build2/demo1
      export LIB_FOLDER=/dataneura/implisolid/build/lib
      export BUILD_LOCATION=/dataneura/implisolid/build
      export CACHE_TEMP=/dataneura/implisolid/cache-temp
      export IMPLISOLID=/dataneura/implisolid
      export IMPLISOLID_REPO=/dataneura/implisolid
      ```
    * Tested on `/dataneura/implisolid`
    * Scripts asummary:

      * `build-clonepull.sh`
          * clones (if not cloned): has two uses: for installing cloned, and, for re-cloning (a speparate)
          * I avoid cMake. (Currenlty the build scripts are for JS target (and soon Wasm). As a C++ library, ImpliSolid yet to be packaged as a C++ library (for use in native C++ projects) -- althought it is a native C++ library, but the build scripts are fine-tuned for use as JS in HTML.)
      * `build-emscripten.sh`
          * The main buid one
      * `e2e-test-builds.bash`
          * re-clones (a sparate copy), and runs, to make sure it is deployable, E2E
          * Uses folder `IMPLISOLID_REPO` rather than `IMPLISOLID`

      * Not needed: ( internal use: bash-utils.sh, deprecated: build_configuration.sh )
