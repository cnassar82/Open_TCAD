# Withdrawn cavity-removal experiment

This is diagnostic evidence, not an accepted engine implementation.
The active source was restored to the pre-experiment build snapshot
`suprem4/build/native64_b43fm60y/src/`, preserving earlier workspace changes.
No launcher or built binary was promoted or rewritten.

The directed-star / exact-ear triangulation experiment eventually passed
independent geometry checks on the captured PMOS activation step, but the
engine rejected its exposed-point bookkeeping:

```text
exposed point 3890 has no gas node
panic after mesh self test
```

See `decks/subcmos2014/results/pmos_exact_ear_replay/result.json`,
`panic.str`, and `gas_backtrace.txt`. The failing candidate was
`native64_44j683ve`; passing geometry does not make that candidate acceptable.

`src/` preserves experimental files, `patches/` records their differences
from the restored snapshot, and `manifest.json` records hashes. `tests/`
preserves experimental helper tests; they are not part of the active suite
and refer to helpers absent from the restored engine.

Restoring the old coarsener is a withdrawal of an unvalidated change, not a
claim that the old algorithm is reliable. Follow
`suprem4/doc/mesh_patch_replacement_contract.md` before replacing it.
