How is the github pages for implisolid [sohale.github.io/implisolid/](https://sohale.github.io/implisolid/) published?


* The main page is a [readme](https://github.com/sohale/implisolid/blob/revival-sohale/docs/readme.md) file that is copied (how) from this folder: [./docs/readme.md](https://github.com/sohale/implisolid/blob/revival-sohale/docs/readme.md)
* Submodules are used: The repo [implisolid-build](https://github.com/sohale/implisolid-build/tree/master) is used as a submodule in the implisolid as well as sohale.github.io 's folder [demos](https://github.com/sohale/sohale.github.io/tree/master/demos). Which one leads to publishing the readme?
* The submodule is extracted using the following githubworkflow:
* [.github/ workflows/jekyll-with-submodules.yml](https://github.com/sohale/sohale.github.io/blob/master/.github/workflows/jekyll-with-submodules.yml)
* Then? The readme is not there is the [demos](https://github.com/sohale/sohale.github.io/tree/master/demos) or the  [implisolid-build](https://github.com/sohale/implisolid-build/tree/master) instantiated (as submodule) in it.
* 
