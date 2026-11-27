"""Generate Genius electrical measurements and collect their result artifacts."""

from math import copysign, isfinite
import json
from pathlib import Path
import re

from curves import read_table

MEASUREMENTS = ("ID–VDS family", "ID–VGS transfer", "BVDSS")
MARKER = "# workbench-measurement "
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_.:-]*")
PREFIX = re.compile(r"[A-Za-z_][A-Za-z0-9_-]{0,63}")


def measurement_cards(kind, values):
    def number(name):
        try:
            result = float(values[name])
        except ValueError:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} requires a number.") from None
        if not isfinite(result):
            raise ValueError("Measurement voltages, steps, and currents must be finite.")
        return result

    def magnitude(name):
        value = number(name)
        if value <= 0:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be positive.")
        return value

    def fmt(value):
        return format(value, ".12g")

    def sweep(contact, start, stop, step, output):
        if start == stop:
            return
        increment = copysign(min(step, abs(stop - start)), stop - start)
        cards.append(f"SOLVE Type=DCSWEEP VScan={contact} VStart={fmt(start)} "
                     f"VStep={fmt(increment)} VStop={fmt(stop)} out.prefix={output}")

    drain, gate, source = (values[name] for name in ("drain", "gate", "source"))
    if any(not IDENTIFIER.fullmatch(name) for name in (drain, gate, source)):
        raise ValueError("Enter the existing Genius drain, gate, and source electrode labels.")
    if len(set(name.casefold() for name in (drain, gate, source))) != 3:
        raise ValueError("The drain, gate, and source labels must be distinct.")
    prefix = values["prefix"]
    if not PREFIX.fullmatch(prefix):
        raise ValueError("Use a filename prefix starting with a letter or underscore.")
    sweep_contact = gate if kind == MEASUREMENTS[1] else drain
    meta = {"kind": kind, "prefix": prefix, "drain": drain, "gate": gate, "source": source,
            "x": f"Vapp({sweep_contact}) [V]", "y": f"I({drain}) [A]", "curves": []}
    cards = []
    stop = number("stop")
    step = magnitude("sweep_step")
    if kind == "BVDSS" and values["impact"] != "Existing model":
        region = values["region"]
        if not IDENTIFIER.fullmatch(region):
            raise ValueError("Enter the Genius semiconductor region label for the impact model.")
        cards.append(f"MODEL Region={region} II={values['impact']} II.Force=EVector")
    cards.append("SOLVE Type=Equilibrium")
    if kind == MEASUREMENTS[0]:
        if stop == 0:
            raise ValueError("VDS stop must differ from zero.")
        try:
            levels = list(dict.fromkeys(float(token) for token in re.split(r"[ ,\t]+", values["levels"].strip())))
        except ValueError:
            raise ValueError("Enter gate biases separated by commas or spaces.") from None
        if not 1 <= len(levels) <= 32 or not all(map(isfinite, levels)):
            raise ValueError("Enter between 1 and 32 finite gate biases.")
        ramp = magnitude("ramp_step")
        previous = 0.0
        for index, bias in enumerate(levels, start=1):
            name = f"{prefix}_vg{index:02}"
            sweep(gate, previous, bias, ramp, name + "_gate")
            sweep(drain, 0, stop, step, name)
            if index < len(levels):
                sweep(drain, stop, 0, step, name + "_return")
            meta["curves"].append({"name": f"VGS = {fmt(bias)} V", "file": name + ".dat"})
            previous = bias
    elif kind == MEASUREMENTS[1]:
        start, bias = number("start"), number("drain_bias")
        if start == stop:
            raise ValueError("The VGS start and stop must differ.")
        ramp = magnitude("ramp_step")
        sweep(drain, 0, bias, ramp, prefix + "_prebias")
        sweep(gate, 0, start, ramp, prefix + "_gate")
        sweep(gate, start, stop, step, prefix)
        if values.get("criterion"):
            meta["criterion_A"] = magnitude("criterion")
            if bias == 0:
                raise ValueError("Use a nonzero fixed VDS to extract VT.")
        meta["curves"].append({"name": f"VDS = {fmt(bias)} V", "file": prefix + ".dat"})
    else:
        if stop == 0:
            raise ValueError("BVDSS requires a nonzero drain sweep.")
        meta["criterion_A"] = magnitude("criterion")
        meta["target_V"] = stop
        meta["curves"].append({"name": "VGS = 0 V", "file": prefix + ".dat"})
        sweep(drain, 0, stop, step, prefix)
    return "\n".join([MARKER + json.dumps(meta, separators=(",", ":")), *cards])


def current_crossing(volts, currents, criterion, *, low_drive_first=False):
    """Interpolate the first rising crossing of an absolute-current criterion."""
    if criterion <= 0 or not isfinite(criterion):
        raise ValueError("Enter a finite, positive current criterion in amperes.")
    points = [(v, abs(current)) if isfinite(v) and isfinite(current) else None
              for v, current in zip(volts, currents)]
    valid = [point for point in points if point is not None]
    result = {"criterion_A": criterion}
    if not valid:
        return dict(result, status="no_converged_points")
    if low_drive_first and abs(valid[0][0]) > abs(valid[-1][0]):
        points.reverse()
        valid.reverse()
    result.update(last_voltage_V=valid[-1][0], max_abs_Id_A=max(i for _, i in valid))
    first_v, first_i = valid[0]
    if first_i > criterion:
        return dict(result, status="criterion_already_exceeded")
    if first_i == criterion:
        return dict(result, status="criterion_reached", voltage_V=first_v)
    for previous, current in zip(points, points[1:]):
        if previous is None or current is None:
            continue
        (v0, i0), (v1, i1) = previous, current
        if i0 < criterion <= i1 and v0 != v1:
            voltage = v0 + (v1 - v0) * (criterion - i0) / (i1 - i0)
            return dict(result, status="criterion_reached", voltage_V=voltage)
    status = "crossing_unresolved" if result["max_abs_Id_A"] >= criterion else "criterion_not_reached"
    return dict(result, status=status)


def column_index(columns, name):
    names = [column.casefold() for column in columns]
    key = name.casefold()
    if key in names:
        return names.index(key)
    # Imported Genius contacts can be reported as Electrode[boundary].
    terminal, separator, unit = key.rpartition(")")
    matches = [i for i, column in enumerate(names)
               if separator and column.startswith(terminal + "[") and column.endswith("])" + unit)]
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError(f"Multiple boundary columns match {name}. Choose the intended curve columns.")
    raise ValueError(f"Missing column {name}.")


def table_column(columns, data, name):
    return data[column_index(columns, name)]


def terminal_voltages(columns, data, measurement):
    return tuple(table_column(columns, data, f"Vapp({measurement[name]}) [V]")
                 for name in ("drain", "gate", "source"))


def breakdown(columns, data, measurement):
    drain, gate, source = terminal_voltages(columns, data, measurement)
    if any(not (isfinite(vg) and isfinite(vs)) or abs(vg - vs) > 1e-8 or abs(vs) > 1e-8
           for vg, vs in zip(gate, source)):
        return {"status": "gate_or_source_not_at_zero", "criterion_A": measurement["criterion_A"]}
    result = current_crossing((vd - vs for vd, vs in zip(drain, source)),
                              table_column(columns, data, measurement["y"]), measurement["criterion_A"],
                              low_drive_first=True)
    if "voltage_V" in result:
        result["BVDSS_V"] = result.pop("voltage_V")
    if "last_voltage_V" in result:
        result["last_Vds_V"] = result.pop("last_voltage_V")
    if "target_V" in measurement:
        result["target_V"] = measurement["target_V"]
    return result


def threshold(columns, data, measurement):
    drain, gate, source = terminal_voltages(columns, data, measurement)
    bias = [vd - vs for vd, vs in zip(drain, source)]
    if not bias or abs(bias[0]) <= 1e-8 or any(not isfinite(v) or abs(v - bias[0]) > 1e-8 for v in bias):
        return {"status": "vds_not_fixed_and_nonzero", "criterion_A": measurement["criterion_A"]}
    result = current_crossing((vg - vs for vg, vs in zip(gate, source)),
                              table_column(columns, data, measurement["y"]), measurement["criterion_A"],
                              low_drive_first=True)
    result["Vds_V"] = bias[0]
    if "voltage_V" in result:
        result["VT_V"] = result.pop("voltage_V")
    return result


def collect_measurements(directory, deck, simulation_status):
    warnings = []
    with deck.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if not line.startswith(MARKER):
                continue
            try:
                measurement = json.loads(line[len(MARKER):])
                prefix = measurement["prefix"]
                if not PREFIX.fullmatch(prefix) or measurement["kind"] not in MEASUREMENTS:
                    raise ValueError("Invalid measurement specification.")
                curves = [curve for curve in measurement["curves"]
                          if Path(curve["file"]).name == curve["file"] and
                          (directory / curve["file"]).is_file()]
                measurement.update(simulator_status=simulation_status, expected_curves=len(measurement["curves"]),
                                   curves=curves)
                if "criterion_A" in measurement and measurement["kind"] in MEASUREMENTS[1:]:
                    extract = breakdown if measurement["kind"] == "BVDSS" else threshold
                    if curves:
                        try:
                            measurement["result"] = extract(*read_table(directory / curves[0]["file"]), measurement)
                        except (OSError, ValueError) as error:
                            measurement["result"] = {"status": "unavailable", "detail": str(error)}
                    else:
                        measurement["result"] = {"status": "no_curve"}
                (directory / f"{prefix}.measurement.json").write_text(
                    json.dumps(measurement, indent=2) + "\n", encoding="utf-8",
                )
            except (OSError, ValueError, KeyError, TypeError) as error:
                warnings.append(str(error))
    return warnings
