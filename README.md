Efficient Cumulative Incidence Estimation in Biobank Studies Using All Prevalent and Incident Events

This repository contains an R package illdthCIF1 which implements the methods in the following paper:

Zucker, D. M., and Gorfine, M. (2026). Efficient cumulative incidence estimation in biobank studies using all prevalent and incident events. arXiv 2606.19041.

The package contains a single function called "cifcmp.full".

Preferred package installation method:

install.packages("remotes")

remotes::install_github("david-zucker/illdthCIF1")

If that fails, follow the steps below:

1. Click the green button "Code" in the upper right of the main package page and then in the dropdown select "Download ZIP".

2. In R, set your working directory to the directory into which you dowloaded the zip file.

3. In R, enter the following command:

remotes::install_local("illdthCIF1-main.zip")

The companion branch "docs" contains the arxiv paper and the package manual.
