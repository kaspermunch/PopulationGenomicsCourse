# course-template repo

This repository is a webpage template for courses.

## Requirements

In order to create a webpage from this repo, you will need to have installed [conda or pip](https://docs.conda.io/projects/conda/en/latest/index.html). Then install the following packages using the terminal:

```bash
pip install mkdocs
pip install mkdocs-material
pip install mkdocs-video
```

## Usage

I would recommend that all the webpage source materials should be saved in an independent branch of the repo (maybe call it webapage). The template should look like this:

```
course-template/
    |
    |-> develop/ # This folder contains the source material
    |   |-> images/ # This folder contains all the images for the webpage as well as the favicon and the logo
    |   |
    |   |-> stylesheets/
    |   |  |-> extra.css # This file contains formatting of the webpage: colors and fonts
    |   |   
    |   |-> *.md # These will be individual webpages that should be inside the mkdocs.yml navigation tree
    |
    |-> docs/ # This folder contains the built webpage
    |
    |-> mkdocs.yml # This is the configuration file of the webpage
```

You should develop the materials for the webpage in markdown format (.md) in the `develop` folder. Mkdocs should also handle Jupypter notebooks (.ipynb) files as well. Ideally, all the images of the .md files will link to the `develop/images` folder.

Once you have created your materials inside the `develop` folder, you can test the website using:

```bash
mkdocs serve
```
**NOTE**: remember to set the working directory to {course-name} general folder, otherwise mkdocs will not be able to find the mkdocs.yml file.

This will prompt a local link were you can check how the webpage will look (usually it is http://127.0.0.1:8000/). 
Then, if you are happy with how it looks, you can build the website:

```bash
mkdocs build
```

This will create the webpage in the `docs` folder. It is necessary that the webpage is created in this folder so that Github Pages can fetch it properly. To set up the Github Pages, go to the **`Settings`** of the repo in [github.com](https://github.com/). Then click on **`Pages`** and set up the **`Build and deployment section`** as below:

![](./develop/images/git_pages.png)

Where:
-   **`Source`** is *`Deploy from a branch`*
-   In **`Branch`**:
    -   **`Select branch`** is *`webpage`* or the name of the branch were this material is
    -   **`Select folder`** is *`/docs`*

When you click on **`Save`**, the webpage should be activated shortly afterwards in the displayed at the top.
