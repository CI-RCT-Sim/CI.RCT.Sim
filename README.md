# CI.RCT.Sim

This repository contains code for the simulation study investigating causial
inference (CI) methods used to account for intercurrent events (ICE) in
randomized controlled trials (RCT).

## Setup Development Environment

* install git
* clone this repository
* open the rstudio project
* install dependencies with `renv::restore()`

### Project Setup Troubleshooting

If you're not using Rstudio: start R with the root folder of the project as
working directory. renv should be loaded automatically and you should see a
message like `Project '...' loaded. [renv 1.1.5]` and you should be able to use
`renv::restore()`. If renv is not automatically loaded use
`source("renv/activate.R")` to load it.

## Project Structure

Reusable functions, methods, etc. will be available as an r-package (built from
this repository), that can be tested, built, installed, and distributed. Bespoke
code for the simulation study will be included in the directory `scripts`.

## Contribute Code

To contribute code:

* switch to the main branch `git checkout main`
* make sure you are on the latest version of the main branch `git pull main`
* create a new branch from there `git checkout -b <new_branch_name>` (replace <new_branch_name> with the name of your branch)
* push your code to github with `git push -u origin <new_branch_name>` for the first time and just `git push` afterwards
* to merge code into the main branch, open a pull request on github. Navigate to your branch use the "Compare & pull request" button. 
