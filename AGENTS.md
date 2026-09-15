# Codex Working Agreement

## Purpose

These instructions define how Codex should plan, implement, verify, and report code changes in this repository. The goal is to keep changes explicit, reviewable, and easy to inspect using GitHub Desktop or another external Git client.

## Repository Location

Before making changes:

1. Identify the absolute path of the repository being edited.
2. State whether it is:
   - the user's normal local checkout;
   - a Codex worktree; or
   - another checkout.
3. When the user requests changes in their normal local checkout, do not make the changes only in a Codex worktree.
4. Do not duplicate the same uncommitted change across multiple checkouts.
5. If the requested checkout is unclear, ask before editing.

## Planning and Scope

Routine documentation, formatting, and small mechanical changes may proceed
when directly requested; reserve a separate approval step for substantive
design or scientific behavior changes.

For nontrivial changes, or whenever the user requests a plan first:

1. Do not edit files during the planning phase.
2. Create a checklist with stable, descriptive change identifiers, such as:
   - `AUC-SECOND-ORDER`
   - `RANK-NORMALIZATION`
   - `TIE-BREAKING-NOISE`
   - `PARTITION-ALIGNMENT`
3. For each change identifier, state:
   - intended behavior;
   - exact files expected to change;
   - functions or sections expected to change;
   - validation or tests to perform;
   - important assumptions or risks.
4. Wait for explicit approval before implementing the proposed changes.
5. Do not renumber or silently reinterpret approved change identifiers later in the conversation.

Before editing, restate the approved scope briefly:

```text
Approved for implementation:
- [CHANGE-ID]: description

Not approved or still pending:
- [CHANGE-ID]: description
```

## Local Repository and Current Branch Only

1. Work only in the user's normal local checkout on its current branch unless
   the user explicitly requests otherwise.
2. Do not create or switch branches, create worktrees, or move work to another
   checkout unless explicitly requested by the user.
3. Do not use online GitHub to inspect or modify this repository. Do not use
   GitHub websites, APIs, or remote repository tools unless explicitly requested.
4. Do not fetch, pull, push, or create pull requests unless explicitly requested.
   Local Git inspection, such as status, diff, and log, is allowed.

## Code Style

1. Prefer compact, readable code that makes use of a standard line of roughly
   80 characters. Use `R/tdhia_stat_tests.R` as a reference for compact layout.
2. Keep multiple arguments on each line in function calls and definitions, and
   multiple elements on each line when defining lists, when they fit comfortably.
   Do not default to one argument or list element per line. Wrap near 80
   characters at sensible boundaries without sacrificing readability.
3. Aim for roughly one code comment line per four lines of code. Treat this as
   a general guideline, using judgment where more explanation is needed.

## Argument Validation

Keep argument validation concise and appropriate for a research-grade R
package. Prefer relying on downstream functions to validate their inputs.

Add explicit checks only when downstream calls would not catch the problem,
or would produce a misleading error, silently incorrect result, or unintended
side effect. Focus on assumptions specific to the function, such as matching
sample identifiers or compatible dimensions.

Avoid duplicating downstream type, length, range, or missing-value checks.
Do not add extensive defensive validation for hypothetical misuse.

## Research Package Priorities

This package supports research and student contributors. Prefer clear,
practical implementations over elaborate frameworks or unnecessary
standardization. Follow established repository conventions unless there
is a concrete reason to change them.

Keep changes focused. Avoid unrelated cleanup, new dependencies, or broad
refactoring during a targeted fix.

## Scientific Behavior

Treat changes to statistical methods, model formulas, normalization,
filtering, missing-value handling, and multiple-testing correction as
scientific behavior changes. Explain the old and proposed behavior and
obtain explicit approval before changing it.

Preserve sample and feature alignment. Where relevant, document matrix
orientation, identifier columns or row names, methylation scale, and
whether filtering occurs before or after statistical adjustment.

Do not silently change scientific behavior while refactoring or updating
documentation. If the implementation appears incorrect, report it and
propose a separate fix.

## Public Function Compatibility

Preserve exported function names, argument names, defaults, return fields,
and output column names unless a change is explicitly approved.

Before changing a public interface, check its callers, examples, and
documentation. Explain any migration needed by existing users.

Keep internal helpers unexported unless they are intentionally useful as
public functions.

## Documentation

When changing a public function, update its roxygen documentation in the
same change. Describe inputs, returned objects, missing-value handling,
and relevant side effects accurately.

Regenerate affected Rd files with roxygen2. Do not edit generated Rd files
by hand. Report if regeneration cannot be completed.

Prefer small examples using synthetic data. Avoid examples that require
private datasets, lengthy analyses, downloads, or credentials.

## Verification

Match verification effort to the change. For scientific behavior changes
and bug fixes, prefer a small reproducible example that distinguishes the
correct result from the previous behavior.

Check identifiers and ordering as well as numeric results when modifying
data processing. Use a fixed random seed when reproducibility matters.

Do not run expensive analyses, download datasets, or regenerate large
outputs merely to check a small change. State what was checked and what
remains unverified.

## Debugging and Research Data

Default db_flag to FALSE. When a wrapper exposes db_flag, propagate that
choice to downstream calls rather than forcing debugging on.

Document debug-file locations and contents. Avoid overwriting snapshots
from parallel workers.

Do not commit participant data, debug workspaces, credentials, or local
analysis outputs. Use synthetic or explicitly approved public data for
examples and tests.

## Supporting Student Contributors

Explain the purpose of a change and any important R or statistical
concept needed to maintain it. Keep explanations concise and concrete.

Use comments to explain assumptions and non-obvious decisions. Prefer
readable, explicit code when a shorter expression would be harder for
a student to understand.

When reporting a problem, distinguish confirmed behavior from suspected
issues and include a small reproduction when practical.
