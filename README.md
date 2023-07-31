# sourmash_plugin_abundhist

## Installation

```
pip install sourmash_plugin_abundhist
```

## Usage

### Example


#### Basic command-line usage:
```
% sourmash scripts abundhist example/SRR606249-abund-100k.sig.zip

== This is sourmash version 4.8.3.dev0. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

loaded 1 total that matched ksize & molecule type

36   [3487]  ****************************************
72   [ 485]  ******
107  [ 171]  **
143  [  38]  *
178  [   5]
214  [   3]
249  [   7]  *
285  [   0]
320  [   2]
356  [   2]
```

#### Create a nice histogram figure:

```
% sourmash scripts abundhist example/SRR606249-abund-100k.sig.zip --figure hist.png
```
will create this figure:

![](example/hist.png)

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`abundhist` is developed at https://github.com/ctb/sourmash_plugin_abundhist.

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
