# Third-party source and redistribution notices

METIS and GKlib are separate bundled dependencies. Their original source-file
copyright headers, license files and license inventories are retained.

| Component | Source revision | Upstream notices |
| --- | --- | --- |
| METIS 5.2.1 | `f5ae915a84d3bbf1508b529af90292dd7085b9ec` | [Copyright and license](METIS/src/LICENSE) |
| GKlib | `3b7d61b9f885063c89901f3901fb4426f9cfb58f` | [Copyright and license](GKlib/LICENSE.txt), [per-file inventory](GKlib/LICENSES.md) |

The default license of both libraries is Apache-2.0. GKlib also contains
LGPL-2.1-or-later GNU C Library derivatives and a BSD-3-Clause Mersenne Twister
implementation. See the upstream inventory for exact files and build conditions;
do not describe the entire dependency set as Apache-only.

Full license texts are supplied in [licenses/Apache-2.0.txt](licenses/Apache-2.0.txt),
[GKlib/LICENSES/LGPL-2.1.txt](GKlib/LICENSES/LGPL-2.1.txt), and
[GKlib/LICENSES/BSD-3-Clause-MT.txt](GKlib/LICENSES/BSD-3-Clause-MT.txt).

nticLBM integrates the libraries as separate shared libraries. GKlib's LGPL
sort templates are also instantiated in METIS, so preserve the LGPL notices
and corresponding-source availability for both libraries when packaging them.
Do not restrict replacement/debugging of the LGPL-covered library components
contrary to their license terms.

For redistribution, include these notices, all applicable license texts,
upstream copyright notices and the exact corresponding third-party sources
and build recipes (including the generated-header settings). The install
rules include the bundled source and integration recipes under
`share/nticLBM/third-party-sources/`; a binary package must retain those materials and
the shared libraries rather than copying only the application executable.
If upstream files are changed later, document the changes prominently in
those files and update `SOURCE_VERSION.txt`; retain existing notices.

The vendored upstream source files are currently unmodified. Integration
recipes under `cmake/3rdparties/` and version manifests are nticLBM additions.
Shared libraries under `3rdParties/lib/` and generated headers under `build/`
contain build products and are not source revisions.
