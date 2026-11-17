# Review of the standalone mesh proposal reference

Reviewed 2026-11-15 against `tests/mesh_patch_oracle.py`, its five tests,
the reduced RIT failure, and the mesh replacement contract.

Later implementation/evidence: [material domains and coupled-solver reliability](material_domain_reliability.md).
That report distinguishes the original removal failure, artificial oxide-domain
constraints, cancellation in finite-volume weights, and the separate ILD offset
boundary failure. This review remains the historical reference review, not the
current implementation status.

## Verdict

Keep this as a useful reference for a simple interior cavity in one region.
The geometric approach is reasonable within that domain. It is insufficient
for the RIT repair and should not be ported directly into the removal loop as
though it implements the replacement contract.

The captured failure is an exposed boundary point selected by edge collapse.
The reference returns UNSUPPORTED for every open boundary fan, including that
case. The passing RIT test proves refusal, not successful processing of the
activation step. An engine integration still needs to avoid an invalid removal
while making progress on short edges and evolving interfaces.

## What is sound

Directed triangle incidence determines the cavity boundary. Exact rational
predicates avoid angle-sorting and tolerance ambiguity for the given coordinates.
The algorithm prepares output without mutating its input. It checks positive
triangles, input/output overlap, boundary shape and output edge multiplicity.
Ear clipping is a reasonable initial triangulator for the supported simple
polygon; I did not find an accepted geometry defect in that domain during this
review.

I reran the five reference tests and three fixture tests: all pass. I also used
an independent exact polygon-clipping calculation to compare the old and new
unions in both directions. It passed 72 generated interior stars with 4–12 rays,
ID/order permutations, 432 output triangles and unchanged inputs. This is useful
additional evidence, not a proof or a production acceptance result.

Reproduce that review from the repository root:

```sh
python3 suprem4/experiments/mesh_proposal_review_20261115/coverage_review.py
```

## Required changes before production integration

1. **Define the actual operation.** `rem_edg` currently selects a surface point
   for deletion and falls through to cavity removal. This reference only handles
   an interior point. Build a constrained edge-collapse or local-remeshing
   proposal for the complete incident patch of both endpoints. The constraints
   must decide which points can move or disappear before velocity/quality
   preferences choose a direction.

2. **Carry region and constraint identities.** The function receives vertex IDs
   and coordinates, with no region labels, material nodes, constrained segments,
   GAS ownership or external attachment map. Its single-region assumption is a
   caller precondition, not an enforced guarantee. A caller can feed triangles
   from different regions without the function detecting the semantic error.
   Preserve the geometric kernel as a component, but put an explicitly validated
   multi-material patch around it.

3. **Add quality acceptance.** The first available ear is selected without a
   quality objective. The existing 64-ray circle test passes while minimum angle
   falls from about 5.625 degrees to 2.8125 degrees. This does not itself prove
   solver failure, but demonstrates that OK means geometrically valid, not an
   improved numerical mesh. Use application-specific quality/error limits and
   consider valid diagonal flips or other candidate triangulations. When a good
   replacement cannot be formed, retain/refine the patch and report the reason.

4. **Implement field transfer and atomic commit.** Area preservation is not
   dopant conservation: removing a node changes the finite-element field.
   Specify transfer rules per field and material, verify the relevant integrals,
   and prepare point/node ownership and connectivity before retiring entities.
   Inject preparation/commit failures and verify unchanged state on rejection.
   This reference implements none of those operations.

5. **Make invalid-input behavior explicit.** NaN, infinity and wrong-dimensional
   coordinates currently raise ValueError/OverflowError rather than return an
   INVALID proposal. Validate the input schema and finiteness before arithmetic.
   In production, prevent malformed inputs and allocation failures from turning
   a supposedly nonmutating proposal into an uncontrolled engine failure.

## Test and reasoning improvements

The old-star area equals the directed polygon area algebraically when interior
rays cancel. That equality is a useful consistency check, but cannot independently
prove domain coverage. Keep the separate intersection/boundary checks and use
independent clipping-based union tests such as the review experiment.

The production C port must be checked against this rational reference on actual
binary32 inputs. Current generated examples also include binary64 circle
coordinates. Rounding, construction of new points, overflow handling and exact
predicate compilation conditions need their own tests. Rational Python should
remain a bounded verification oracle; its pairwise checks and Fraction arithmetic
are not evidence of production scalability.

Add constrained interfaces, same-material/different-region boundaries,
point-only contacts, rejection-state comparisons, field-transfer tests, external
attachment tests, and failure injection when their corresponding implementation
exists. Do not label tests that merely expect UNSUPPORTED as repaired process
cases.

## Concrete next experiment

In the captured case, point 3796 is exposed and point 3567 is interior. Before
building another general triangulator, capture the complete patch of both
endpoints and evaluate the opposite contraction direction: retain surface point
3796 and remove/merge interior point 3567. This is a hypothesis to test, not an
approved repair. Verify the topological link condition, positive triangles,
unchanged material boundary, attachments, ownership and field transfer. If that
direction is invalid too, expand the local remeshing patch while preserving the
surface constraint.

This focuses the next implementation on the actual failing operation. It also
avoids assuming that safe refusal alone will let oxidation progress.

## Documentation corrections to keep in mind

The early claim of overlapping live oxide triangles was based on a panic export
that also serialized retired triangles. The subsequent live-triangle audit
passed at the failing operation. Use `rit_first_invalid_removal.md` as the current
evidence; the older overlap assertion still present in the replacement contract
must not be used as proof of an earlier live-mesh corruption.

The temporary loss of GAS occurs inside an incomplete point removal. It explains
why partial mutation is unsafe when cavity construction subsequently fails; it
does not prove that this same deletion caused the different missing-GAS failure
in the withdrawn exact-ear candidate. Keep those two causal claims separate.
