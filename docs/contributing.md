---
title: "Contributing to the Docs"
nav_order: 89
layout: default
render_with_liquid: false
---

# Contributing to the Docs

The site is built with [Jekyll](https://jekyllrb.com) and the [Just the Docs](https://just-the-docs.com) theme. Every page is a Markdown file in `docs/`. Pushing to `master` triggers an automatic GitHub Actions build and deploy.

## Access

The repo is open to collaborators with write access. If you don't have access yet, write to **Luca Greco** to be added.

---

## Installing prerequisites

You need **Ruby ≥ 3.1** and **Bundler**. Git is assumed to already be installed.

**macOS**

The system Ruby is too old. Install a current version via [Homebrew](https://brew.sh):

```bash
brew install ruby
```

Then add Ruby to your PATH as Homebrew instructs (something like `export PATH="/opt/homebrew/opt/ruby/bin:$PATH"`), and install Bundler:

```bash
gem install bundler
```

**Linux (Debian/Ubuntu)**

```bash
sudo apt install ruby-full build-essential zlib1g-dev
gem install bundler
```

**Windows**

Use [RubyInstaller](https://rubyinstaller.org/) (select the version with Devkit). Bundler is included.

---

## Local preview

```bash
git clone git@github.com:<username>/<fork-name>.git
cd fork-name
bundle install        # only needed once, or after Gemfile changes
bundle exec jekyll serve --baseurl /fork-name
```

Open `http://localhost:4000/fork-name/`. Changes to `.md` files are picked up on reload without restarting the server.

---

## Making a pull request

1. Clone the fork (once) and keep it up to date (see [Developments](/code_development/development_workflow.md) ):
   ```bash
   git clone git@github.com:<username>/<fork-name>.git
   git pull origin develop
   ```
2. Create a branch from `develop`:
   ```bash
   git checkout -b my-topic
   ```
3. Add or edit files, commit, and push the branch:
   ```bash
   git add docs/my_page.md
   git commit -m "Add documentation for X"
   git push origin my-topic
   ```
4. Open a pull request against `develop` of the main repository on GitHub. Once merged, the site redeploys automatically — no further action needed.

---

## Front matter reference

Every page starts with a YAML block delimited by `---`. YAML is whitespace-sensitive and uses `key: value` pairs. A few rules worth knowing:

- **Strings** generally do not need quotes. Use double quotes if the value contains a colon, `#`, or leading/trailing spaces: `title: "Solver: Details"`.
- **Booleans** are bare `true` / `false` — no quotes.
- **Integers** are bare numbers — no quotes.
- Indentation must be spaces, not tabs.

The keys used in this site are listed below. For the full reference see the [Just the Docs navigation docs](https://just-the-docs.com/docs/navigation/).

| Key | Type | Required | Description |
|---|---|---|---|
| `title` | string | yes | Page title shown in the browser tab and sidebar. |
| `layout` | string | yes | Always `default` for content pages. |
| `render_with_liquid` | bool | yes | Set to `false` to prevent Jekyll from interpreting `{{` and `{%` as template syntax — important for pages containing code examples. |
| `nav_order` | integer | no | Position within the parent group in the sidebar. Lower = higher up. Defaults to the page title alphabetically if omitted. |
| `parent` | string | no | Title of the parent page, exactly as written in that page's `title` field. Nests this page under it in the sidebar. |
| `grand_parent` | string | no | Title of the grandparent page. Required when `parent` is itself a child (third-level pages). |
| `has_children` | bool | no | Marks the page as a section that can contain children. |
| `nav_fold` | bool | no | Collapses the section in the sidebar by default. Only meaningful when `has_children: true`. |
| `nav_exclude` | bool | no | Hides the page from the sidebar. The page is still built and accessible via its URL. |

---

## Creating a page

Every `.md` file under `docs/` is a page. It needs a YAML front matter block:

```yaml
---
title: "My Page"
layout: default
render_with_liquid: false
---

Content goes here.
```

That's the minimum. The page will be built and accessible via its URL but will not appear in the sidebar unless you add navigation fields (see below).

To **exclude a page from the sidebar** while still making it linkable:

```yaml
nav_exclude: true
```

---

## Adding a sidebar section

A sidebar section is just a page with `has_children: true` and `nav_fold: true`. By convention, section landing pages are named `cat_<name>.md` and live in their own subdirectory:

```
docs/
  mysection/
    cat_mysection.md   ← section landing page
    page_a.md
    page_b.md
```

`cat_mysection.md` front matter:

```yaml
---
title: "My Section"
nav_order: 8
has_children: true
nav_fold: true
layout: default
render_with_liquid: false
---
```

Child pages reference the section by its `title`:

```yaml
---
title: "Page A"
nav_order: 1
parent: "My Section"
layout: default
render_with_liquid: false
---
```

---

## Managing sidebar elements

| Front matter key | Effect |
|---|---|
| `nav_order` | Position within the parent group (lower = higher up). Applies at every level. |
| `parent` | Nests the page under a section. Value must exactly match the section's `title`. |
| `grand_parent` | Required for third-level pages. Value must match the section title (the grandparent). |
| `has_children: true` | Marks a page as a section that can have children. |
| `nav_fold: true` | Starts the section collapsed in the sidebar. |
| `nav_exclude: true` | Hides the page from the sidebar entirely. |

The navigation supports multiple levels of nesting. The `grand_parent` key always refers to the direct grandparent (one step above `parent`), not necessarily the top-level section.

To **rename a section**, change the `title` in its `cat_*.md` file and update the `parent` field of every child page to match.

---

## Linking between pages

Links use standard Markdown syntax with `.md` extensions — the build converts them to `.html` automatically:

```markdown
See [Normalization](normalization.md) for details.
```

For pages in a different directory, use a path relative to the current file:

```markdown
[Compiling](compiling/compilation_and_libraries/compiling.md)
```

Absolute paths from the repo root also work:

```markdown
[Compiling](compiling/compilation_and_libraries/compiling.md)
```

Do not hard-code `.html` extensions — they will break if the site baseurl changes.

---

## Images

Keep images next to the page that uses them. Create an `assets/` folder in the same directory as the `.md` file, with a subfolder named after the file:

```
docs/howto/physics_options/physical_terms/
  diamag.md
  assets/
    diamag/
      my_figure.png
```

Reference them with a path relative to the `.md` file:

```markdown
![Alt text](assets/diamag/my_figure.png)
```

To control size, use an HTML `<img>` tag:

```html
<img src="assets/diamag/my_figure.png" alt="Alt text" width="600">
```

Prefer PNG or SVG for diagrams and screenshots; keep file sizes reasonable — large binaries slow down the repository for everyone.

---

## Maths

Inline LaTeX: `$...$`. Display math: `$$...$$`. MathJax renders both. Standard macros work; avoid bare `||` inside inline math — use `\Vert` instead.
