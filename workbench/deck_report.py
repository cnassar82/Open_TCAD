"""Render Markdown annotations and case artifacts into a portable report."""

from dataclasses import dataclass
from base64 import b64encode
from html import escape
from io import BytesIO
import json
from math import isfinite
from pathlib import Path
import re
import shlex
from string import Template
from tempfile import NamedTemporaryFile
from urllib.parse import unquote

NOTE = "#|"
FILENAME = "workbench-report.html"
EDITOR_FILENAME = "workbench-editor-report.html"
REPORT_FILES = {FILENAME, EDITOR_FILENAME}
DOPING_FIELDS = {"NetDoping", "Donor", "Acceptor"}


@dataclass(frozen=True)
class Block:
    kind: str
    content: str
    line: int


def markdown_lines(text):
    return "\n".join(f"{NOTE} {line}" if line else NOTE for line in text.splitlines())


def annotated_blocks(deck):
    blocks = []
    contents = []
    start = 0
    fence = None

    def flush():
        if contents and any(line.strip() for line in contents):
            blocks.append(Block("markdown", "\n".join(contents), start))
        contents.clear()

    for i, line in enumerate(deck.splitlines(), 1):
        visible = line.lstrip()
        if not visible.startswith(NOTE):
            if contents and fence is None and contents[-1]:
                contents.append("")
            continue
        annotation = visible[len(NOTE):].removeprefix(" ")
        match = re.match(r"^ {0,3}(`{3,}|~{3,})(.*)$", annotation)
        if match:
            token = match[1]
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence) and not match[2].strip():
                fence = None
        elif fence is None and annotation.startswith("@"):
            flush()
            blocks.append(Block("directive", annotation, i))
            continue
        elif fence is None and re.match(r"^#\s+\S", annotation):
            flush()
            blocks.append(Block("slide", "", i))
        if not contents:
            start = i
        contents.append(annotation)
    flush()
    return blocks


class ReportRenderer:
    def __init__(self, case_dir=None, deck_dir=None, factors=None):
        self.case_dir = Path(case_dir).resolve() if case_dir is not None else None
        self.deck_dir = Path(deck_dir).resolve() if deck_dir is not None else None
        self.factors = factors or {}
        self.structures = {}
        self.warnings = []

    def file(self, name, *, authored=False):
        path = Path(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Use a file path relative to the case or deck folder.")
        roots = [self.deck_dir, self.case_dir] if authored else [self.case_dir]
        for root in roots:
            if root is None:
                continue
            target = (root / path).resolve()
            if target.is_relative_to(root) and target.is_file():
                return target
        raise FileNotFoundError(f"Could not find {name}.")

    def structure(self, name):
        from doping import load_structure

        path = self.file(name)
        if path.suffix.lower() not in {".str", ".vtu"}:
            raise ValueError("Choose a .str or .vtu structure.")
        if path not in self.structures:
            self.structures[path] = load_structure(path)
        return self.structures[path]

    @staticmethod
    def figure(canvas, caption):
        data = BytesIO()
        canvas.savefig(data, format="png", dpi=155, bbox_inches="tight", pad_inches=.15)
        return image_tag(data.getvalue(), "image/png", caption)

    @staticmethod
    def coordinates(values):
        parsed = [float(value) for value in values]
        if not all(isfinite(value) for value in parsed):
            raise ValueError("Coordinates must be finite.")
        return parsed

    def structure_figure(self, args, opts):
        if len(args) != 1:
            raise ValueError("Use @structure FILE with optional field= and scale=.")
        structure = self.structure(args[0])
        import numpy as np
        from matplotlib import colormaps
        from matplotlib.figure import Figure
        from matplotlib.tri import Triangulation
        from suprem_viewer import add_material_legend, material_entries, scalar_values, transform_values

        field = viewer_field(opts.get("field", "Material"))
        if field not in structure.fields:
            raise ValueError(f"{field} is unavailable in {structure.path.name}.")
        figure = Figure(figsize=(8, 4.1), layout="constrained")
        ax = figure.subplots()
        coords = structure.points.copy()
        # Sparse SUPREM coordinate IDs may leave unreferenced NaN rows.
        coords[~np.isfinite(coords)] = 0.0
        mesh = Triangulation(coords[:, 0], coords[:, 1], structure.triangles)
        values, label = transform_values(scalar_values(structure, field), field, scale_log(opts))
        if not np.isfinite(values).any():
            raise ValueError(f"{field} has no finite values in {structure.path.name}.")
        if field == "Material":
            cmap = colormaps["tab20"].resampled(max(1, len(material_entries(structure))))
            ax.tripcolor(mesh, facecolors=values, cmap=cmap, shading="flat", vmin=-.5,
                         vmax=max(.5, len(material_entries(structure)) - .5))
            add_material_legend(ax, structure, cmap)
        else:
            artist = ax.tripcolor(mesh, facecolors=values, cmap="viridis", shading="flat")
            figure.colorbar(artist, ax=ax, label=label)
        ax.set_xlabel("X (µm)")
        ax.set_ylabel("Y (µm)")
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        return self.figure(figure, opts.get("title", f"{structure.path.name} · {label}"))

    def compare_checkpoints(self, args, opts):
        if len(args) != 2:
            raise ValueError("Use @compare FIRST.str SECOND.str.")
        figure = checkpoint_figure([self.structure(name) for name in args])
        return self.figure(figure, opts.get("title", "Checkpoint comparison"))

    def profile_figure(self, args, opts):
        if len(args) != 5:
            raise ValueError("Use @profile FILE X0 Y0 X1 Y1.")
        structure = self.structure(args[0])
        import numpy as np
        from matplotlib.figure import Figure
        from suprem_viewer import sample_line_profile

        coords = self.coordinates(args[1:])
        if coords[:2] == coords[2:]:
            raise ValueError("Choose two different profile endpoints.")
        field = viewer_field(opts.get("field", "NetDoping"))
        if field == "Material" or field not in structure.fields:
            raise ValueError("Choose a saved scalar field for the profile.")
        samples = int(opts.get("samples", "200"))
        if not 2 <= samples <= 2000:
            raise ValueError("Use between 2 and 2000 profile samples.")
        dist, values, _, _, label = sample_line_profile(structure, field, scale_log(opts), *coords, samples)
        if not np.isfinite(values).any():
            raise ValueError("The line has no finite field samples. Check the cut coordinates and field.")
        figure = Figure(figsize=(8, 3.75), layout="constrained")
        ax = figure.subplots()
        ax.plot(dist, values, color="#2356a8", lw=1.7)
        ax.set_xlabel("Distance along cut (µm)")
        ax.set_ylabel(label)
        ax.grid(alpha=.27)
        return self.figure(figure, opts.get("title", f"{structure.path.name} · ({coords[0]:g}, {coords[1]:g}) "
                                             f"to ({coords[2]:g}, {coords[3]:g}) µm"))

    def point(self, args, opts):
        if len(args) != 3:
            raise ValueError("Use @point FILE X Y.")
        structure = self.structure(args[0])
        from doping import ACCEPTORS, DONORS, sample_doping

        field = viewer_field(opts.get("field", "NetDoping"))
        if field == "Material" or field not in structure.fields:
            raise ValueError("Choose a saved scalar field for the probe.")
        x, y = self.coordinates(args[1:])
        values = sample_doping(structure, [(x, y)], "Net doping" if field == "NetDoping" else field)
        unit = " (cm⁻³)" if field in DOPING_FIELDS or field.removesuffix("Active") in DONORS + ACCEPTORS else ""
        rows = "".join(
            "<tr>" + "".join(f"<td>{escape(item)}</td>" for item in
                             (material or "—", f"{value:.6g}" if value is not None else "—", status))
            + "</tr>" for _, _, material, value, status in values
        )
        title = opts.get("title", f"{field} at ({x:g}, {y:g}) µm · {structure.path.name}")
        return (f"<table><caption>{escape(title)}</caption><thead><tr><th>Material</th>"
                f"<th>{escape(field + unit)}</th><th>Status</th></tr></thead><tbody>{rows}</tbody></table>")

    def saved_measurement(self, args, opts):
        from measurements import Measurement, cell, evaluate

        if len(args) != 2:
            raise ValueError("Insert a measurement recipe with Add to report.")
        name, recipe = args
        spec = Measurement.from_dict(json.loads(recipe))
        caption = opts.get("title", f"{spec.kind} · {name}")
        if spec.kind == "Doping cutline":
            a, b, c, d, field, samples, scale = spec.settings
            return self.profile_figure([name, a, b, c, d],
                                       dict(field=field, samples=samples, scale=scale, title=caption))
        path = self.file(name)
        columns, rows = evaluate(spec, path, structure=self.structure(name) if spec.kind not in {"VT", "BVDSS"} else None)
        headings = "".join(f"<th>{escape(heading)}</th>" for heading in columns)
        content = "".join("<tr>" + "".join(f"<td>{escape(cell(value))}</td>" for value in row) + "</tr>"
                          for row in rows)
        return (f"<table><caption>{escape(caption)}</caption><thead><tr>{headings}</tr></thead>"
                f"<tbody>{content}</tbody></table>")

    def scalar_value(self, args, opts):
        if len(args) != 2:
            raise ValueError("Use @value FILE field.path.")
        path = self.file(args[0])
        if path.suffix.lower() != ".json":
            raise ValueError("Choose a JSON result for @value.")
        value = json.loads(path.read_text(encoding="utf-8"))
        for part in args[1].split("."):
            if not isinstance(value, dict) or part not in value:
                raise ValueError(f"{args[1]} is unavailable in {path.name}.")
            value = value[part]
        if value is None or isinstance(value, (dict, list)):
            raise ValueError(f"{args[1]} does not contain a scalar value.")
        if isinstance(value, float):
            value = f"{value:.7g}"
        return metric(opts.get("title", args[1]), str(value), opts.get("unit", ""))

    def local_image(self, args, opts):
        return self.read_image(args, opts, inline=False)

    def read_image(self, args, opts, *, inline):
        if len(args) != 1:
            raise ValueError("Use @image relative/path.png.")
        path = self.file(args[0], authored=True)
        from PIL import Image

        if path.stat().st_size > 15_000_000:
            raise ValueError("Use an image smaller than 15 MB.")
        try:
            with Image.open(path) as image:
                formats = {"PNG": "image/png", "JPEG": "image/jpeg"}
                if image.format not in formats:
                    raise ValueError("Use a PNG or JPEG image.")
                mime = formats[image.format]
                image.verify()
        except Image.DecompressionBombError as error:
            raise ValueError("Use an image with fewer pixels.") from error
        return image_tag(path.read_bytes(), mime, opts.get("title", path.name), inline=inline)

    def markdown_image(self, url, caption):
        decoded = unquote(url)
        if ":" in decoded:
            return f'<span class="inline-error">Image: {escape(caption or url)} (use a relative path)</span>'
        try:
            return self.read_image([decoded], {"title": caption or Path(decoded).name}, inline=True)
        except (OSError, ValueError, ImportError) as error:
            self.warnings.append(str(error))
            return f'<span class="inline-error">Image: {escape(caption or decoded)} · {escape(str(error))}</span>'

    @staticmethod
    def placeholder(name, label, detail):
        return f'<p class="note"><strong>{escape(name)} · {escape(label)}</strong><br>{escape(detail)}</p>'

    def directive(self, block):
        try:
            first, *tokens = shlex.split(block.content)
            name = first[1:]
            args, opts = [], {}
            for token in tokens:
                key, separator, value = token.partition("=")
                if separator and key in {"field", "scale", "samples", "title", "unit"}:
                    if key in opts or not value:
                        raise ValueError(f"Enter one value for {key}.")
                    opts[key] = value
                else:
                    args.append(token)
            handlers = {"structure": self.structure_figure, "profile": self.profile_figure,
                        "point": self.point, "value": self.scalar_value, "image": self.local_image,
                        "compare": self.compare_checkpoints, "measure": self.saved_measurement}
            if name == "factor":
                if len(args) != 1:
                    raise ValueError("Use @factor NAME.")
                if args[0] in self.factors:
                    return metric(opts.get("title", args[0]), str(self.factors[args[0]]), opts.get("unit", ""))
                if self.case_dir is None:
                    return self.placeholder("Factor", args[0], "Filled from the selected case.")
                raise ValueError(f"The selected case has no factor named {args[0]}.")
            if name not in handlers:
                raise ValueError(f"Unknown report directive @{name}.")
            if self.case_dir is None and name != "image":
                return self.placeholder(name.capitalize(), " ".join(args), "Filled from the selected case.")
            return handlers[name](args, opts)
        except (OSError, ValueError, TypeError, RuntimeError, KeyError, IndexError, ImportError) as error:
            self.warnings.append(f"Line {block.line}: {error}")
            return self.placeholder(f"Deck line {block.line}", block.content, str(error))


def checkpoint_figure(structures):
    import numpy as np
    from matplotlib import colormaps
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch
    from matplotlib.tri import Triangulation

    labels = [sorted({region.material for region in structure.regions.values()}) for structure in structures]
    materials = sorted(set().union(*labels))
    colors = colormaps["tab20"].resampled(len(materials))
    palette = {name: colors(i) for i, name in enumerate(materials)}
    figure = Figure(figsize=(9, 4.2), layout="constrained")
    axes = figure.subplots(1, len(structures), sharex=True, sharey=True, squeeze=False)[0]
    bounds = []
    for ax, structure in zip(axes, structures):
        coords = np.nan_to_num(structure.points, nan=0.0, posinf=0.0, neginf=0.0)
        mesh = Triangulation(coords[:, 0], coords[:, 1], structure.triangles)
        values = [materials.index(structure.regions[int(region)].material) for region in structure.triangle_regions]
        ax.tripcolor(mesh, facecolors=np.array(values), shading="flat",
                     cmap=ListedColormap(list(palette.values())), vmin=-.5, vmax=len(materials) - .5)
        bounds.append(coords[structure.triangles].reshape(-1, 2))
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(structure.path.name)
        ax.set_xlabel("X (µm)")
    low = np.min(np.concatenate(bounds), axis=0)
    high = np.max(np.concatenate(bounds), axis=0)
    padding = np.maximum((high - low) * .035, 1e-6)
    axes[0].set(xlim=(low[0] - padding[0], high[0] + padding[0]),
                ylim=(high[1] + padding[1], low[1] - padding[1]), ylabel="Y (µm)")
    figure.legend(handles=[Patch(color=color, label=name) for name, color in palette.items()],
                  loc="outside lower center", ncols=min(5, len(materials)), fontsize=8)
    return figure


def viewer_field(field):
    return "NetDoping" if field.casefold().replace(" ", "") == "netdoping" else field


def scale_log(opts):
    scale = opts.get("scale", "linear")
    if scale not in {"linear", "log"}:
        raise ValueError("Use scale=linear or scale=log.")
    return scale == "log"


def image_tag(data, mime, caption, *, inline=False):
    image = b64encode(data).decode("ascii")
    if inline:
        return f'<img class="inline-image" src="data:{mime};base64,{image}" alt="{escape(caption)}">'
    return (f'<figure><img src="data:{mime};base64,{image}" alt="{escape(caption)}">'
            f"<figcaption>{escape(caption)}</figcaption></figure>")


def metric(title, value, unit=""):
    label = escape(f"{value} {unit}".strip())
    return f'<div class="metric">{escape(title)}<strong>{label}</strong></div>'


def render_report(deck, *, case_dir=None, deck_dir=None, factors=None, title="Deck report", metadata=()):
    import mistune

    renderer = ReportRenderer(case_dir, deck_dir, factors)

    class MarkdownRenderer(mistune.HTMLRenderer):
        def image(self, text, url, title=None):
            return renderer.markdown_image(url, title or text)

    md = mistune.create_markdown(renderer=MarkdownRenderer(escape=True), plugins=["table"])
    pages = [[]]
    for block in annotated_blocks(deck):
        if block.kind == "slide" or (block.kind == "directive" and block.content.strip() == "@slide"):
            if pages[-1]:
                pages.append([])
        elif block.kind == "markdown":
            pages[-1].append(md(block.content))
        else:
            pages[-1].append(renderer.directive(block))
    if not any(pages):
        raise ValueError("Add report notes with the Report menu or prefix Markdown lines with #|.")
    sections = "\n".join(f'<section class="page">{"".join(page)}</section>' for page in pages if page)
    header = "".join(f"<span>{escape(str(item))}</span>" for item in metadata if item)
    template = Template(Path(__file__).with_name("report_template.html").read_text(encoding="utf-8"))
    html = template.substitute(title=escape(title), pages=sections,
                               metadata=f'<p class="metadata">{header}</p>' if header else "")
    return html, renderer.warnings


def write_report(path, deck, **options):
    html, warnings = render_report(deck, **options)
    with NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False, mode="w",
                            encoding="utf-8") as output:
        temporary = Path(output.name)
        try:
            output.write(html)
            output.close()
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
    return warnings
