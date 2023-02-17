# Splashscreen for jupyterlab

Show an intro widget in the laucher that displays information about your own content and code. The code is taken from the [ENKI-intro splashscreen reposiroty](https://gitlab.com/ENKI-portal/jupyterlab_enkiintro/-/tree/master/).


## Prerequisites

* JupyterLab > 3.0

## How to use it

- Clone this repository in place into your project. 

- Go into the folder `style` to change the `index.css` file with the splashscreen design, which uses also images saved in the same folder. 

- Go into the folder `src` to modify the file `index.ts`, containing (apart from various settings that can be left untouched), the text appearing into the splash screen. This is modified using the HTML language.

## Installation

Once you have modified the files above, try out your splashscreen. Install it locally by going into the main folder of the splashscreen, and write the `bash` command
```bash
jupyter labextension install .
```
Now you should be able to launch jupyterlab and see the splashscreen.

### Uninstall

To remove the splashscreen:
```bash
jupyter labextension uninstall jupyterlab_test
```
You do not need this to reinstall the splashscreen, you can just install it again to test changes.