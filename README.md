# sourmash_plugin_abundhist

## Installation

```
pip install sourmash_plugin_abundhist
```

## Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
