# Project setup and git workflow

## Setup Development Environment

- install git
- clone this repository
- open the rstudio project
- install dependencies with
  [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)

### Project Setup Troubleshooting

If you’re not using Rstudio: start R with the root folder of the project
as working directory. renv should be loaded automatically and you should
see a message like `Project '...' loaded. [renv 1.1.5]` and you should
be able to use
[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html).
If renv is not automatically loaded use `source("renv/activate.R")` to
load it.

## Project Structure

Reusable functions, methods, etc. will be available as an r-package
(built from this repository), that can be tested, built, installed, and
distributed. Bespoke code for the simulation study will be included in
the directory `scripts`.

## Contribute Code / git workflow

To contribute code:

1.  switch to the main branch `git checkout main`
2.  make sure you are on the latest version of the main branch
    `git pull main`
3.  create a new branch from there `git checkout -b <new_branch_name>`
    (replace with the name of your branch)
4.  write your code
5.  Run
    [`renv::status()`](https://rstudio.github.io/renv/reference/status.html)
    to check if any pacakges need to be installed, updated or added to
    renv. If necessary use
    [`renv::snapshot`](https://rstudio.github.io/renv/reference/snapshot.html)
    to add packages you installed to renv or
    [`renv::restore`](https://rstudio.github.io/renv/reference/restore.html)
    if you are missing pacakges others have added.
6.  run `devtools::document()` to generate documentation and update
    `DESCRIPTION` and `NAMESPACE`
7.  run `devtools::test()` and `devtools::check()`, if there are any
    errors, fix them and iterate from step 5. if all errors are taken
    care of continue to 8.
8.  commit your code to git through the Rstudio userinterface or with
    `git add .` followed by `git commit`
9.  push your code to github with `git push -u origin <new_branch_name>`
    for the first time and just `git push` afterwards
10. to merge code into the main branch, open a pull request on github.
    Navigate to your branch use the “Compare & pull request” button.

# Coding Style

Please approximately keep to the [tidyverse style
guide](https://style.tidyverse.org/) to keep the code easily readable.
