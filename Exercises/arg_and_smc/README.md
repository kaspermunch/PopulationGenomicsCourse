# Ancestral Recombination Graphs

## Install the ARG dashboard

### Plan A: Conda

Create a conda environment on your own computer:

```bash
conda create -n arg-dashboard -c conda-forge -c plotly -c kaspermunch popgen-dashboards=1.1.4
```

Acitvate the environment:

```bash
conda activate arg-dashboard
```

Run the dashboard app:

    arg-dashboard

if it does not show up in your browser, you can right-click (open in new tab) or paste this address into your browser : [http://127.0.0.1:8050](http://127.0.0.1:8050)

### Plan B: Docker

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) on your own machine. 

Start the Docker Desktop application (you may be prompted to create a login and password for DockerHub)

Run the dashboard app in your in the Terminal on Mac or in Anaconda PowerShell Prompt on Windows:

```bash
docker run --rm -i -t -p 8050:8050 kaspermunch/arg-dashboard-linux-amd64:1.1.5
```

To view the dashboard, right-click (open in new tab) or paste this address into your browser : [http://127.0.0.1:8050](http://127.0.0.1:8050)


> If you have a small screen, you may need to zoom out a bit to see en entire dashboard. On Chrome, you click the top right three dots and select a zoom level of 80%.

## Exercise 1: Simulate some ARGs

1. The *Main* panel shows a simulated graph. Click `New` to generate a new graph. Dropdown menus control the type of graph, the number of samples and the sequence length. Choose "ARG" for the simulation, a sample size of "5" and a sequence length of "2kb". 
2. Simulate a lot of ARGS. Just keep clicking `New` to see the variation of ARGS.
3. Pick and ARG with 2-3 recombincombination events.
4. The *Coalesce and recombination events* panel controls the number of events shown in the graph. The *Recombination points* panel shows the points of recombination in the sequence. Use the *Coalesce and recombination events* panel to reconstruct the graph step by step by moving the slider from left to right. Make sure you understand what happens at each coalescence and recombination event.
5. Use the two sliders in the *Recombination points* to retrict the view to a smaller part of the sequence. Watch what happens to the arg as you include more or fewer recombination points on the sequence.

## Exercise 2: Ancestral sequences and marginal trees

1. Click `New` to find an ARG with 2-3 recombincombination events where some of the nodes are not red. Ancestral node colors reflect the proportion of the sequence that is represented in the ancestor. Ancestral material is the part of the sequence that has a descendants among the sampled sequences. Mouse-over a node to see the type of event and the proportion of ancestal proportion sequence at each node. 
2. Notice how mouse-over also activates the *Ancestral sequences* and *Marginal trees* panels. The *Ancestral sequences* panel shows how ancestral material is merged at coalescence events and divided at recombination events. Sequence segments separated by recombination events are shown with separate colors. The genalogy for each segment is shown in the *Marginal tree* panel with a color matching the segment. Non-ancestral segments are shown in white. Mouse-over the root node of the ARG to see all marginal trees (genealogies) for the ARG. Mouse-over other nodes to see the marginal trees below that node.
3. Use the slider in the  *Recombination points* panel to show the graph for only a subset of the sequence. Notice how this affects the shown ARG and the marginal trees when you mouse-over a node. Sequence outside the range specified in the *Recombination points* panel is shown as gray in the *Ancestral sequences* panel.
4. The marginal trees look a lot like each other. Try to understand how each one if different from its neighbor. You can make the range in the *Recombination points* slider really small so it shows only one marginal tree at a time. Now move it along the sequence to see how the ARG changes. Can you see that only one branch detaches and is moved somewhere else?

## Exercise 3: "Captured" non-ancestral sequence

1. See if you can find any nodes where non-ancestral sequence is "captured" between two ancestral (colored) segments. Figure out how the segment you found got trapped between two ancestral sequence segments.
2. See if you can find it i the SMC. If you cannot, why do you think that is?

## Exercise 4: "Diamond" recombinations

1. Simulating from the ARG, you may have noticed that some recombinations events produce two lineages that coalesce with each other rather than with other lineages. The following coalescence nullifies the preceeding recombination event. This produces producing an "eye" or "diamond" in the graph. Such recombination events cannot be inferred from data (if you mouse over the nodes you can see why), but they are part of the evolutionary history of a sample none the less.
2. Try the SMC. Do you find any "diamonds"?. Why?
3. Now try the SMC'. Do you find any "diamonds"?. Why?
4. What do you think assuming diamonds do not exist does to the estimation of the graph?

## Exercise 5: Imagine complexity

1. Imagine and ARG for a human chromosome - 100,000 times larger than the 2kb we are considering here.


<!-- Log into [UCloud](https://cloud.sdu.dk/app/dashboard) and complete this part of the exercise there. -->


<!-- Set this up on your own machine

```
conda create --name popgen-dashboards -c conda-forge -c plotly -c kaspermunch popgen-dashboards
```

First, clone the following github: 

git clone https://github.com/kaspermunch/popgen-dashboards/

Then download the notebook by right-clicking <a href="https://raw.githubusercontent.com/kaspermunch/PopulationGenomicsCourse/master/Notebooks/arg-dashboard.ipynb" download="arg-dashboard.ipynb">
this link
</a> and "choose save link as". Place it in the popgen_dashboards folder, and run it using jupyter notebook -e popgen-dashboards -->
