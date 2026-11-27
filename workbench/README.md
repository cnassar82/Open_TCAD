# Open TCAD Workbench

Launch the desktop app from the repository root.

```sh
python3 workbench/app.py
```

The editor uses Python's Tkinter. Install NumPy and Matplotlib in the same Python environment for plotting and the existing structure viewer. The app uses the checkout's validated SUPREM4 launcher and local Genius runtime. Set `TCAD_SUPREM4` or `TCAD_GENIUS` to select a different executable. Set `TCAD_GENIUS_DIR` to the runtime directory that contains the Genius `lib` folder.

GDS import uses [gdstk](https://heitzmann.github.io/gdstk/library/gdstk.read_gds.html). For a local environment with the same plotting libraries as your current Python, run these commands from the repository root.

```sh
python3 -m venv --system-site-packages tmp/workbench-venv
tmp/workbench-venv/bin/python -m pip install 'gdstk>=0.9.61,<2'
tmp/workbench-venv/bin/python workbench/app.py
```

## Edit and run

- Open a SUPREM4 or Genius deck from the library or use **Open…** to browse. **New** loads a diode starter for the selected simulator. **Insert command** adds an editable card.
- Click **Run deck** or press F5 to run the editor contents, including unsaved changes. Use **Run inputs…** to stage relative inputs such as included decks, masks, or imported device meshes. Save the deck first so the workbench can preserve their paths relative to the deck's directory.
- Use relative output names so every run writes to its own case folder. Deck commands that refer to absolute paths or paths outside the case folder still resolve on the local filesystem.

## Build a process step

Choose **Process wizards → Set up starting structure** to start a deck from a rectangular 2D substrate. Set its width, depth, mesh spacing, material, crystal orientation, background doping, and surface and backside boundaries.

Choose **Process wizards → Import GDS mask edges** and select a GDS file, cell, and layer/datatype. Enter start and end points for a straight cut in layout coordinates (µm). The cut maps into deck X as distance from the start plus your selected deck origin. Choose whether drawn polygons are openings or protected areas. **Compute cut** previews the opening intervals. **Mesh lines** merges their edges into the starting X mesh. **Mask opening etches** inserts polygon etch cards at the caret. Place those cards after depositing the mask material and set the etch top and bottom to span that film. **Copy intervals** copies the deck X coordinates as CSV.

Choose **SUPREM4 → Process wizards** to open **Deposit**, **Diffuse**, **Oxidize**, or **Etch**. The forms read the active runtime's compiled command catalog, including available materials, ambients, model selectors, and parameter descriptions. Advanced fields accept simulator expressions where the card allows them. Blank optional fields leave the simulator's defaults in effect.

- **Deposit** sets the film material, thickness, mesh divisions, and optional doping or contour file.
- **Implant** sets species, dose, energy, and the distribution, angle, or damage options provided by the SUPREM4 command catalog.
- **Diffuse** sets the furnace temperature, duration, and anneal or dopant ambient.
- **Oxidize** generates a `diffuse` card with dry or wet oxygen.
- **Etch** supports whole-material removal, directional depth, either side of a line, material-dependent rates, a polygon, or a contour file. Enter polygon vertices in µm, one X/Y pair per line.

Use **Models and coefficients** to choose a diffusion or oxidation model on the `method` card, or configure an `oxide`, species, or material card. **Queue model card** lets you configure another instance of the same card. The preview shows the model cards before the process step, and **Insert into deck** adds those editable lines at the editor caret. Stage any contour files with **Run inputs…** before running.

## Save and resume a process checkpoint

1. Put the caret after the step you want to save. Choose **Insert command → Save checkpoint…** to add a SUPREM4 `structure outfile` card. Use **Save checkpoint and stop…** to end the run there. Place the save between completed process steps.
2. Run the deck. In **Run results → Artifacts**, select the resulting `.str` file and click **New deck from checkpoint**. A new deck opens with a `structure infile` card and a place to add the next process step. To load an existing `.str` into the current deck, use **Insert command → Load checkpoint…**.
3. Reapply any custom `method`, model, and process option cards needed by the resumed steps. A structure file carries the mesh, materials, and saved solution fields. A new run starts with the simulator's runtime settings. The generated load card refers to the original checkpoint path, so keep the source file available.

## Inspect process steps and compare checkpoints

- Open **Process wizards → Inspect process** to see the current editor deck in execution order. Select a complete process card, such as the last polygon `etch done` that opens a resist mask, and click **Run through step**. The run executes the deck prefix, saves a checkpoint, then stops. **Insert capture** instead adds a `structure outfile` card to the editor so a later full run can save the same step. Use **Run with step log** when running the full deck to record command echoes.
- In **Run results**, select a SUPREM4 case and click **Inspect process**. Select a checkpoint row to see its cross-section. A step without a saved structure can show its preceding checkpoint. The **Step log** tab shows matched command output when the run recorded echoes and the recent case log otherwise.
- Ctrl-click two saved checkpoint rows and click **Compare checkpoints** for common X and Y limits and material colors. **Branch from checkpoint** opens a continuation deck. **Measure checkpoints** applies the same recipe to the selected structures.
- **Add views to report** inserts a structure or side-by-side comparison in the open deck. Use **Report → Build editor report with selected case** to render it. To inspect older outputs, use **Process wizards → Inspect saved process folder** and select the saved deck alongside its output files.

## Sweep a parameter

1. Select a number such as the final `60` in `diffuse time=60 temp=1000`. Click **Mark factor…** and name it `anneal`. The editor inserts `{{anneal=60}}`. A regular run substitutes the default value.
2. Click **Design DOE…** and enter levels such as `5, 10, 20`. With multiple factors, the workbench runs every combination. Cases in a batch run sequentially, up to 128 per batch.
3. Select a case in **Run results** to see its log and files. Open `.str` or `.vtu` files in the existing mesh and field viewer. Open `.csv`, `.tsv`, or Genius `.dat` tables to choose and plot numeric columns. Select a table and click **Compare DOE curves** to overlay the same file across completed cases.

For an electrical example, choose **Genius**, click **New**, and sweep `stop_voltage` over `0.1, 0.2`. Compare `forward_iv.dat` in the resulting batch.

The workbench stores histories, rendered decks, logs, and outputs under `tmp/workbench/`. **Open result…** also opens files from earlier runs outside the workbench.

## Run a process-to-device DOE

1. Open or write a SUPREM4 process deck. Stage any relative process inputs with **Run inputs**. Use **Assign contacts** on a representative final structure and save its contact map.
2. Choose **Process-to-device DOE → New DOE project**. Select the contact map. The process can initialize its own substrate or load a starting checkpoint at an `@checkpoint@` line. The runner writes the final structure as `study-device.str`.
3. Mark numeric factors as `{{name=default}}` in process cards, Genius options, or numerical electrical test settings. A factor name used in several places shares one setting and must have the same default everywhere.
4. Add electrical tests for an ID–VDS family, ID–VGS transfer with optional constant-current VT, or BVDSS. Each test uses a distinct output prefix and runs in a separate Genius job using the case's exported mesh. Use **Measurements → Add saved recipe** to include an optional scalar response from any matching case-relative artifact.
5. In **DOE design**, use **Find factors** to populate the level list. Full factorial runs every combination of levels. Latin hypercube uses one fixed value or two range endpoints per factor and a reproducible seed. Check the case count, save the project, and run it.
6. Open **Run results → DOE results → Open response table** to compare factors, test status, VT, BVDSS, and optional measurement responses. Double-click a case to inspect its artifacts. **Export response CSV** saves the table. **Rerun failed or canceled cases** starts a separate batch using the old run's staged inputs and preserves prior results.

Batches run locally and sequentially, with up to 128 cases. Identical rendered process decks in a batch reuse the first successful process output. Retries can also reuse a previously completed process stage. Each response shows its status. A missing criterion crossing leaves the voltage blank, and a recipe that returns several values records the individual rows in the case's `doe-measurements.json`. Existing implant-dose projects can be reopened in the new editor, and old results retain their dose plot and CSV.

## Save and reuse a measurement

- In **Run results → Artifacts**, select a structure or an electrical `.measurement.json` and open **Measurements → Reusable measurements**. Choose a doping point, doping cutline, junction depth, material thickness, VT, or BVDSS. **Open files** selects multiple independent checkpoints, and **Across DOE cases** applies the recipe to the matching artifact from each completed case.
- **Save recipe** stores the settings in JSON. **Load recipe** restores them later. **Evaluate** fills the result table, **Plot cutlines** overlays the profiles, and **Export CSV** saves all evaluated rows.
- Junction depth is estimated from signed net-doping crossings within the specified crystalline-silicon depth interval. Film thickness measures each matching layer at the given X. VT and BVDSS reevaluate the existing curve with the recipe's current criterion.
- **Add to report** embeds the recipe in the deck as an `@measure` annotation. A report generated from a matching case uses the same settings and case-relative artifact. A DOE project can also import the saved recipe for scalar responses.

## Write a report alongside a deck

Use the **Report** menu to insert Markdown or a data block at the editor caret. Each report line starts with `#|`, which keeps it a simulator comment. Ordinary deck lines run as usual. **Preview editor layout** renders the writing and shows placeholders for data that will come from a case.

```text
#| # Implant study
#| This section records the final structure and a channel measurement.
structure outfile=final.str
#| @structure final.str title="Final materials"
#| @profile final.str 1.5 0.005 1.5 1 field=NetDoping scale=log
#| @point final.str 1.5 0.1 field=NetDoping
#| @slide
#| ## Electrical results
#| @value electrical/case-001/transfer.measurement.json result.VT_V title="Threshold voltage" unit=V
#| @factor dose title="Implant dose" unit="cm⁻²"
```

- Use case-relative file paths for structures and measurement JSON. A profile takes two X/Y endpoint pairs in µm. A point reading takes one X/Y pair. If a field name or file path has spaces, use quotes.
- Use `@image picture.png title="Caption"` or Markdown image syntax for a local PNG or JPEG beside the deck or in the case. Figures and images are embedded in the HTML.
- Select a finished case and choose **Build report from selected case** to render the run's saved deck. **Build editor report with selected case** applies the current notes to the same output and writes a separate file. Both reports appear under **Artifacts**. Missing files or unavailable fields are shown in the report.
- Choose **Slides** in the browser to present, with Previous/Next or the arrow keys. A Markdown H1 or `@slide` starts a new slide. Use **Report** for a continuous view and **Print** to print or save the current view to PDF.

The report renderer uses Mistune 3. If it is absent, install `mistune>=3,<4` in the Python environment you use to launch the app. Use the same option for the checkout-local virtual environment when running with GDS support.

## Measure device results

For a process checkpoint, select a SUPREM4 `.str` artifact and click **Assign contacts**. To browse an earlier structure file, use **Assign contacts to checkpoint** from either wizard menu. Choose **MOSFET** for Source, Drain, Gate, and Body, **Diode** for Anode and Cathode, **BJT** for Emitter, Base, and Collector, or **Custom**.

Select a metal or polysilicon region, or use **Edge** or **Box edges** to select exposed silicon for an ohmic contact. Click **Add from selection** and choose the terminal name. Set the gate work function when assigning a gate. Connected metal and polysilicon regions form one ideal conductor in the exported mesh. The exporter checks for duplicate terminal names, overlapping assignments, and unassigned conductive islands.

**Generate Genius deck** saves a paired `.inp` and `.tif` and opens the deck with its mesh staged. Click **Run deck** to solve equilibrium and produce a structure for inspection before adding biased electrical tests.

The starter deck uses `METHOD Type=DDML1 NS=Basic LS=MUMPS MaxIt=30`. It sets the contact types and gate work function, and Genius supplies its default device physics. Add region-specific `MODEL` and `PMI` cards before the solve when you want to control the physical models or material parameters explicitly. The terminal preset only supplies contact names and types.

For a Genius MOSFET deck, open **Electrical wizards**. Insert the cards after device, contact, model, and solver setup and before `END`. Set the source and body at 0 V. Use the electrode labels already present on your device and a distinct output prefix for each measurement.

- **ID–VDS family** sweeps VDS at each requested fixed VGS. **ID–VGS transfer** sweeps VGS at a fixed VDS. Fill in the optional **VT current criterion** to extract constant-current VT automatically.
- **BVDSS** sweeps VDS with VGS at zero. Enter the total drain-current criterion in amperes and choose an impact-ionization option for the semiconductor region, or use the model already set in the deck. The report interpolates between neighboring converged samples that bracket the criterion. It reports the last sampled voltage and leaves the breakdown value blank when the criterion is unreached.
- Open the resulting `.measurement.json` from **Run results → Artifacts** for the grouped curves or BVDSS result. Individual curve tables are saved alongside it.

To extract from an existing DAT, CSV, or TSV curve, open it and choose **Extract VT / BVDSS**. Set the plot's X column to VGS for VT or VDS for BVDSS. Set Y to drain current and enter the absolute current criterion in amperes. Use a transfer curve at fixed nonzero VDS for VT, or a sweep with VGS = 0 V for BVDSS. Choose a criterion appropriate to the device width and test bias. The extraction linearly interpolates between adjacent samples, retains the voltage sign, and leaves the voltage blank when the data do not resolve a crossing. You can export one result per plotted curve, including DOE overlays.

Select a saved `.str` or `.vtu` structure in **Run results → Artifacts** and click **Measure doping**. Enter one X/Y pair per line in µm. Choose signed net doping, donor, acceptor, or a supported saved species. Sample one structure or the matching output across completed DOE cases and export the result as CSV. At a material interface, the table can report a row for each material. SUPREM4 totals use the validated runtime's default dopant types and prefer active concentrations where saved. Select individual species when using custom donor/acceptor assignments.
