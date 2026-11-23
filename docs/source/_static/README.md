# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
html_static_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images

## Logo

The logo font is "Avenir Next" bold.

See below for logo colors.

The PDF image `imdclient-logo.pdf` is editable. 
The PNG image `imdclient-logo.png` is rendered at 
300dpi on a solid white background.

## Theme colors

| Role | Hex | RGB |
|---|---:|---:|
| Deep navy | `#062B4A` | `rgb(6, 43, 74)` |
| Stream blue | `#0788D1` | `rgb(7, 136, 209)` |
| Neutral gray | `#8C99A5` | `rgb(140, 153, 165)` |
| White | `#FFFFFF` | `rgb(255, 255, 255)` |
| Background | Transparent | `rgba(0, 0, 0, 0)` |

### `mdanalysis-sphinx-theme` mapping

These map onto variables in the theme’s
[`site.sass`](https://github.com/MDAnalysis/mdanalysis-sphinx-theme/blob/main/mdanalysis_sphinx_theme/sass/site.sass):

| `site.sass` | IMDClient color |
|---|---|
| `$primary` | Stream blue (`#0788D1`) |
| `$blue` | Deep navy (`#062B4A`) |
| `$accent` (`color_accent`) | Deep navy (`#062B4A`) |

`color_accent` and `mobile_navbar_background` are set in `conf.py` via
`html_theme_options`. `$primary` / `$blue` are hardcoded in the theme’s
compiled CSS, so brand overrides live in `custom.css` and are loaded with:

```python
html_css_files = ["custom.css"]
```

