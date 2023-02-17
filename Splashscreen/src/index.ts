import {
  JupyterFrontEnd, JupyterFrontEndPlugin, ILayoutRestorer
} from '@jupyterlab/application';

import {
  Dialog, ICommandPalette, WidgetTracker, ISplashScreen
} from '@jupyterlab/apputils';

import {
	ILauncher
} from '@jupyterlab/launcher'

import {
	JSONExt
} from '@lumino/coreutils'

import {
  Widget
} from '@lumino/widgets';

import {
  DisposableDelegate, IDisposable
} from '@lumino/disposable';

import '../style/index.css';

const SPLASH_RECOVER_TIMEOUT = 120000;

class ENKIWidget extends Widget {
	constructor() {
		super();
		this.id = 'enkiinfo-jupyterlab';
  		this.title.label = 'Health Data Science Sandbox';
  		this.title.closable = true;
  		this.addClass('jp-ENKIWidget');

    	this.img = Private.createSplash();
    	this.node.appendChild(this.img);
	}

	readonly img: HTMLElement;
};

let splasher: IDisposable;

function activate(app: JupyterFrontEnd, palette: ICommandPalette, launcher: ILauncher, restorer: ILayoutRestorer, splash: ISplashScreen) {
	console.log('JupyterLab extension Splashscreen_sandbox is activated!');

	let widget: ENKIWidget;

	const command: string = 'enkiinfo:open';
	app.commands.addCommand(command, {
		label: args => (args['isPalette'] ? 'ENKI Information' : 'Guidelines'),
    caption: 'Course Information',
    iconClass: args => (args['isPalette'] ? '' : 'jp-ENKIicon'),
		execute: args => {
			if (!widget) {
				widget = new ENKIWidget();
				widget.update();
			}
			if (!tracker.has(widget)) {
				tracker.add(widget);
			}
			if (!widget.isAttached) {
				app.shell.add(widget);
			} else {
				widget.update();
			}
			app.shell.activateById(widget.id);
			return widget;
		}
	});

	launcher.add({
    command: command,
		category: 'Other',
		rank: 0
	});

	palette.addItem({command, category: 'Help', args: { isPalette: true } });

	let tracker = new WidgetTracker<Widget>({ namespace: 'Guidelines'});
	restorer.restore(tracker, {
		command,
		args: () => JSONExt.emptyObject,
		name: () => 'Guidelines'
	});

	splasher = splash.show();

};

const palette: JupyterFrontEndPlugin<void> = {
  id: 'jupyterlab_enkiintro:palette',
  autoStart: true,
  requires: [ICommandPalette, ILauncher, ILayoutRestorer, ISplashScreen],
  activate: activate
};

const splash: JupyterFrontEndPlugin<ISplashScreen> = {
  id: 'jupyterlab_enkiintro:splash',
  autoStart: true,
  provides: ISplashScreen,
  activate: app => {
    return {
      show: () => {
        const { restored } = app;
        const recovery = () => { splasher.dispose(); };

        return Private.showSplash(restored, recovery);
      }
    };
  }
};

const plugins: JupyterFrontEndPlugin<any>[] = [
  palette, splash
];
export default plugins;



let mainContent =
'<h1>Welcome to the material for the Population Genetics Course</h1>'
+ '<p>This page contains some basic information and resources to help you Navigating the material effectively.</p>'
+ '<br>'
+ '<p>For more info about the material visit, see <a href="https://hds-sandbox.github.io/Popgen_course_aarhus" target="_blank">the course website</a>.</p>'
+ '<br>'
+ '<h2>Find the interactive notebooks</h2>'
+ '<p>You can find the code for the course into jupyter notebooks. Those are located'
+ 'in the folder <b>Course_Material/Notebooks</b> located on the left-side browser of your window </p>'
+ '<br>'
+ '<h2>Saving the datasets and notebooks on uCloud</h2>'
+ '<p>If you are using this material on the computing cluster uCLoud, the data and notebooks will not be automatically saved'
+ 'into your account. If you launch the app again, you will start from scratch running the code, and any change you made will not be available.'
+ 'To see how you can reuse your material and datasets, please follow the simple instructions '
+ 'contained in the <a href="https://docs.cloud.sdu.dk/Apps/genomics.html#copying-folders-from-your-session" target="_blank">Documentation page .</p>';

/*
let mainContent = 
'<h1>Welcome to the ENKI server</h1>'
+ '<p>This page displays some basic information and links to resources that will help '
+ 'you use the server effectively.  The information on this page can always be redisplayed '
+ 'by executing the <strong>ENKI Information</strong> command on the <strong>Commands</strong> '
+ 'palette or by clicking on the <strong>ENKI info</strong> launcher button.</p>'
+ '<p>The ENKI server is built on top of a Jupyter Lab computing environment. If you '
+ 'have never used Jupyter Lab, please consult the excellent '
+ '<a href="http://jupyterlab.readthedocs.io/en/latest/" target="_blank">Jupyter Lab '
+ 'User Guide</a> found at Read The Docs.</p>'
+ '<p>There are a series of videos on the '
+ '<a href="https://www.youtube.com/channel/UCJSeYpnbGcxv8WLrDEN383A" target="_blank">'
+ 'Youtube (ENKI-portal)</a> that describe how to use software provided on the server. </p>' 
+ '<p>Jupyter notebooks that illustrate how to perform thermodynamic calculations '
+ 'using the ENKI infrastructure are accessible from the ENKI tab (visible on the '
+ 'left edge of the browser window once the splash screen is closed).</p>'
+ '<p>The principal <strong>software repository</strong> for ENKI may be found at ' 
+ '<a href="https://gitlab.com/ENKI-portal" target="_blank">Gitlab (ENKI-portal)</a>. '
+ 'The code base is open source and in development. Until '
+ 'ENKI is officially released, permission to access the software requires a login and '
+ 'approval.  Please request access from the '
+ '<a href="mailto:ghiorso@ofm-research.org" target="_blank">ENKI PIs</a>. If you have '
+ 'GitLab credentials, you may interact with GitLab repositories via the '
+ 'GitLab tab.  Consult the Youtube How-to video for proper configuration of ENKI related '
+ 'GitLab settings.</p>'
+ '<p>The link <a href="http://enki-portal.org" target="_blank">ENKI Portal</a> points to '
+ 'the project website.</p>'
+ '<p>There are two main software and data repositories that support the ENKI '
+ 'infrastructure. </p>'
+ '<dl>'
+ '<dt><strong>ThermoEngine</strong> '
+ 'A Python package for calculating thermodynamic properties from various databases '
+ 'and for performing equilibrium calculations (includes interfaces to Berman, Holland and '
+ 'Powell, Stixrude, MELTS and DEW)</dt>'
+ '<dd>- <a href="https://enki-portal.gitlab.io/ThermoEngine" target="_blank">Documentation</a></dd>'
+ '<dd>- <a href="https://gitlab.com/ENKI-portal/ThermoEngine" target="_blank">Code '
+ 'respository</a> (requires GitLab login)</dd>'
+ '</dl>'
+ '<dl>'
+ '<dt><strong>Geothermodat</strong> '
+ 'A Python package and data store for accumulating, manipulating and querying '
+ 'phase equilibrium data used in calibrating thermodynamic models of minerals and melts</dt>'
+ '<dd>- <a href="https://enki-portal.gitlab.io/geothermodat" target="_blank">Documentation</a></dd>'
+ '<dd>- <a href="https://gitlab.com/ENKI-portal/geothermodat" target="_blank">Code '
+ 'respository</a> (requires GitLab login)</dd>'
+ '</dl>';
*/

namespace Private {
  export
  function createSplash(): HTMLElement {
      const splash = document.createElement('div');
      splash.id = 'enki-splash';
      const container = document.createElement('div');
      container.id = 'container';

      const header = document.createElement('header');
      header.innerHTML = '<h1>Guidelines</h1>';
      container.appendChild(header);

      const logo = document.createElement('nav');
      //logo.id = 'main-logo';
      container.appendChild(logo);

      const mainText = document.createElement('article');
      mainText.innerHTML = mainContent;
      container.appendChild(mainText);

      const footer = document.createElement('footer');
      footer.innerHTML = 'Copyright &copy; HDS Sandbox';
      container.appendChild(footer);

      splash.appendChild(container);

      return splash;
  }

  let debouncer = 0;
  let dialog: Dialog<any>;

  function recover(fn: () => void): void {
    if (dialog) {
      return;
    }

    dialog = new Dialog({
      title: 'Loading...',
      body: `The loading screen is taking a long time.
        Would you like to clear the workspace or keep waiting?`,
      buttons: [
        Dialog.cancelButton({ label: 'Keep Waiting' }),
        Dialog.warnButton({ label: 'Clear Workspace' })
      ]
    });

    dialog.launch().then(result => {
      if (result.button.accept) {
        return fn();
      }

      dialog.dispose();
      dialog = null;

      debouncer = window.setTimeout(() => {
        recover(fn);
      }, SPLASH_RECOVER_TIMEOUT);
    });
  }

  const splash = createSplash();

  let splashCount = 0;
  let splashButton = 0;

  /**
   * Show the splash element.
   *
   * @param ready - A promise that must be resolved before splash disappears.
   *
   * @param recovery - A function that recovers from a hanging splash.
   */
  export
  function showSplash(ready: Promise<any>, recovery: () => void): IDisposable {
    splash.classList.remove('splash-fade');
    splashCount++;

    if (debouncer) {
      window.clearTimeout(debouncer);
    }
    debouncer = window.setTimeout(() => {
      recover(recovery);
    }, SPLASH_RECOVER_TIMEOUT);

    if (splashButton == 0) {
    	let buttonCont = document.createElement('div');
    	buttonCont.setAttribute('style', 'position: relative;');
    	let button = document.createElement("button");
    	button.id = 'enki-button';
    	button.textContent = 'Close this screen';
    	button.onclick = function() {
      		console.log('Button clicked.')
      		splasher.dispose();
    	};
    	buttonCont.appendChild(button);
    	splash.appendChild(buttonCont);
    	splashButton = 1;
	}

	document.body.appendChild(splash);

    return new DisposableDelegate(() => {
      ready.then(() => {
        if (--splashCount === 0) {
          if (debouncer) {
            window.clearTimeout(debouncer);
            debouncer = 0;
          }

          if (dialog) {
            dialog.dispose();
            dialog = null;
          }

          splash.classList.add('splash-fade');
          window.setTimeout(() => { document.body.removeChild(splash); }, 500);
        }
      });
    });
  }

  export
  function noOp() { /* no-op */ }
}
