import {
  JupyterFrontEnd,
  JupyterFrontEndPlugin
} from '@jupyterlab/application';

/**
 * Initialization data for the @enki-portal/enkiintro extension.
 */
const extension: JupyterFrontEndPlugin<void> = {
  id: '@enki-portal/enkiintro',
  autoStart: true,
  activate: (app: JupyterFrontEnd) => {
    console.log('JupyterLab extension @enki-portal/enkiintro is activated!');
  }
};

export default extension;
