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

## Theme colors

| Role | Hex | RGB |
|---|---:|---:|
| Deep navy | `#062B4A` | `rgb(6, 43, 74)` |
| Stream blue | `#0788D1` | `rgb(7, 136, 209)` |
| Neutral gray | `#8C99A5` | `rgb(140, 153, 165)` |
| White | `#FFFFFF` | `rgb(255, 255, 255)` |
| Background | Transparent | `rgba(0, 0, 0, 0)` |

