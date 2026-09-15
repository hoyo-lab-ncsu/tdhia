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
