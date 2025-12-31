# Sandlersteam

> Digitized steam tables from Sandler's 5th ed.

Sandlersteam implements a python interface to the steam tables found in Appendix III of _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.

The interface operates similarly to the IAPWS steam tables.

## Installation 

Sandlersteam is available via `pip`:

```sh
pip install sandlersteam
```

## Usage

The steam tables provided in _Chemical, Biochemical, and Engineering Thermodynamics_ use temperature in Celsius, pressure in MPa (or kPa for very low pressures in the saturated steam tables), and report specific properties per kg in kJ (enthalpy and internal energy), kJ/K (entropy), or m<sup>3</sup> (volume).

### Command-line interface

Querying properties of steam and water:

```bash
$ sandlersteam state -TC 800 -P 40
THERMODYNAMIC STATE OF UNSATURATED STEAM/WATER:
T =  1073.2 K =  800.0 C
P =  40.00 MPa =  400.00 bar
u =  3517.8 kJ/kg =  63374.2 J/mol
v =  0.011523 m3/kg =  0.00020759 m3/mol
s =  6.6662 kJ/kg-K =  120.093 J/mol-K
h =  3978.7 kJ/kg =  71677.4 J/mol
$ sandlersteam state -TC 200 -x 0.5 
THERMODYNAMIC STATE OF SATURATED STEAM/WATER:
TC =  200.0 C =  473.1 K
P  =  1.55 MPa =  15.54 bar
u  =  1722.98 kJ/kg =  31039.9 J/mol
v  =  0.0642585 m3/kg =  0.00115763 m3/mol
s  =  4.3816 kJ/kg-K =  78.9358 J/mol-K
h  =  1822.82 kJ/kg =  32838.7 J/mol
x  =  0.50 kg vapor/kg
uL =  850.65 kJ/kg =  15324.7 J/mol
uV =  2595.3 kJ/kg =  46755.1 J/mol
vL =  0.001157 m3/kg =  2.08437e-05 m3/mol
vV =  0.12736 m3/kg =  0.00229443 m3/mol
sL =  2.3309 kJ/kg-K =  41.9918 J/mol-K
sV =  6.4323 kJ/kg-K =  115.88 J/mol-K
hL =  852.45 kJ/kg =  15357.1 J/mol
hV =  2793.2 kJ/kg =  50320.3 J/mol
```
<!-- 200     1.5538  0.001157 0.12736   850.65 1744.7 2595.3  852.45 1940.7 2793.2 2.3309 4.1014 6.4323 -->
Generating LaTeX source for a table or a block:

```bash
$ sandlersteam latex --suphP 4.0 --output tmp.tex
Request completed: wrote output to tmp.tex
$ cat tmp.tex
```
```latex
\clearpage
\noindent THERMODYNAMIC PROPERTIES OF STEAM (Selected)\\*[1cm]
Superheated steam:\\*[0mm]
\noindent\begin{minipage}{0.6\textwidth}
\footnotesize\vspace{5mm}
\begin{center}
$P$ = 4.0 MPa\\*[1ex]
\begin{tabular}{>{\raggedleft}p{8mm}@{}p{5mm}>{\raggedleft}p{4mm}@{}p{10mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft\arraybackslash}p{3mm}@{}p{8mm}}
\toprule
\multicolumn{2}{c}{$T$~($^\circ$C)} & \multicolumn{2}{c}{$\hat{V}$} & \multicolumn{2}{c}{$\hat{U}$} & \multicolumn{2}{c}{$\hat{H}$} & \multicolumn{2}{c}{$\hat{S}$}\\
\toprule
\midrule
250 & .4 & 0 & .04978 & 2602 & .3 & 2801 & .4 & 6 & .0701 \\
275 &  & 0 & .05457 & 2667 & .9 & 2886 & .2 & 6 & .2285 \\
300 &  & 0 & .05884 & 2725 & .3 & 2960 & .7 & 6 & .3615 \\
350 &  & 0 & .06645 & 2826 & .7 & 3092 & .5 & 6 & .5821 \\
400 &  & 0 & .07341 & 2919 & .9 & 3213 & .6 & 6 & .7690 \\
450 &  & 0 & .08002 & 3010 & .2 & 3330 & .3 & 6 & .9363 \\
500 &  & 0 & .08643 & 3099 & .5 & 3445 & .3 & 7 & .0901 \\
600 &  & 0 & .09885 & 3279 & .1 & 3674 & .4 & 7 & .3688 \\
700 &  & 0 & .11095 & 3462 & .1 & 3905 & .9 & 7 & .6198 \\
800 &  & 0 & .12287 & 3650 & .0 & 4141 & .5 & 7 & .8502 \\
900 &  & 0 & .13469 & 3843 & .6 & 4382 & .3 & 8 & .0647 \\
1000 &  & 0 & .14645 & 4042 & .9 & 4628 & .7 & 8 & .2662 \\
1100 &  & 0 & .15817 & 4248 & .0 & 4880 & .6 & 8 & .4567 \\
1200 &  & 0 & .16987 & 4458 & .6 & 5138 & .1 & 8 & .6376 \\
1300 &  & 0 & .18156 & 4674 & .3 & 5400 & .5 & 8 & .8100 \\
\bottomrule
\end{tabular}
\end{center}
\end{minipage}

\noindent $\hat{V}\ [=]\ \mbox{m$^3$/kg}$; $\hat{U}\ [=]\ \mbox{kJ/kg}$; $\hat{H}\ [=]\ \mbox{kJ/kg}$; $\hat{S}\ [=]\ \mbox{kJ/kg-K}$\\*[1cm]
```

### API

Below we create a `State` object to define a thermodynamic state for steam at 100 deg. C and 0.1 MPa:

```python
>>> from sandlersteam.state import State
>>> state1 = State(TC=100.0, P=0.1)
>>> state1.h  # enthalpy in kJ/kg
2676.2
>>> state1.u  # internal energy in kJ/kg
2506.7
>>> state1.v  # volume in m3/kg
1.6958
>>> state1.s  # entropy in kJ/kg-K
7.3614
```

Specifying a state requires values for two independent state variables.  The state variables recognized by sandlersteam are:

* `TC` temperature in C (alternatively `TK` for temperature in K)
* `P` pressure in MPa
* `u` specific internal energy in kJ/kg
* `v` specific volume in m<sup>3</sup>/kg
* `h` specific enthalpy in kJ/kg
* `s` specific entropy in kJ/kg
* `x` quality; mass fraction of vapor in a saturated vapor/liquid system (between 0 and 1)

Initializing a `State` instance with any two of these values set by keyword parameters results in the other
properties receiving values by (bi)linear interpolation.

When specifying quality, a `State` objects acquires `Liquid` and `Vapor` attributes that each hold intensive, saturated single-phase property values.  The property value attributes owned directly by the `State` object reflect the quality-weighted sum of the respective single-phase values:

```python
>>> s = State(TC=100,x=0.5)
>>> s.P
0.10135
>>> s.v
0.836972
>>> s.Vapor.v
1.6729
>>> s.Liquid.v
0.001044
>>> 0.5*(s.Vapor.v+s.Liquid.v)
0.836972
```
One can also import the `SteamTables` dictionary from the state `state` module and then generate LaTeX-compatible versions of either blocks in the superheated/subcooled steam tables or entire saturated steam stables, listed by temperature or pressure. For example:

```python
>>> from sandlersteam.state import SteamTables as st
>>> print(st['suph'].to_latex(P=1.0))  # generates latex for the 1 MPa block of the superheated steam table
```

```latex
\begin{minipage}{0.6\textwidth}
\footnotesize\vspace{5mm}
\begin{center}
$P$ = 1.0 MPa\\*[1ex]
\begin{tabular}{>{\raggedleft}p{8mm}@{}p{5mm}>{\raggedleft}p{4mm}@{}p{10mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft\arraybackslash}p{3mm}@{}p{8mm}}
\toprule
\multicolumn{2}{c}{$T$~($^\circ$C)} & \multicolumn{2}{c}{$\hat{V}$} & \multicolumn{2}{c}{$\hat{U}$} & \multicolumn{2}{c}{$\hat{H}$} & \multicolumn{2}{c}{$\hat{S}$}\\
\toprule
\midrule
179 & .91 & 0 & .19444 & 2583 & .6 & 2778 & .1 & 6 & .5865 \\
200 &  & 0 & .2060 & 2621 & .9 & 2827 & .9 & 6 & .6940 \\
250 &  & 0 & .2327 & 2709 & .9 & 2942 & .6 & 6 & .9247 \\
300 &  & 0 & .2579 & 2793 & .2 & 3051 & .2 & 7 & .1229 \\
350 &  & 0 & .2825 & 2875 & .2 & 3157 & .7 & 7 & .3011 \\
400 &  & 0 & .3066 & 2957 & .3 & 3263 & .9 & 7 & .4651 \\
500 &  & 0 & .3541 & 3124 & .4 & 3478 & .5 & 7 & .7622 \\
600 &  & 0 & .4011 & 3296 & .8 & 3697 & .9 & 8 & .0290 \\
700 &  & 0 & .4478 & 3475 & .3 & 3923 & .1 & 8 & .2731 \\
800 &  & 0 & .4943 & 3660 & .4 & 4154 & .7 & 8 & .4996 \\
900 &  & 0 & .5407 & 3852 & .2 & 4392 & .9 & 8 & .7118 \\
1000 &  & 0 & .5871 & 4050 & .5 & 4637 & .6 & 8 & .9119 \\
1100 &  & 0 & .6335 & 4255 & .1 & 4888 & .6 & 9 & .1017 \\
1200 &  & 0 & .6798 & 4465 & .6 & 5145 & .4 & 9 & .2822 \\
1300 &  & 0 & .7261 & 4681 & .3 & 5407 & .4 & 9 & .4543 \\
\bottomrule
\end{tabular}
\end{center}
\end{minipage}
```

This renders as

![a steam table block](https://github.com/cameronabrams/Sandlersteam/raw/main/stimage.png)
<!-- ![a steam table block](stimage.png "1 MPa superheated steam table block") -->

## Release History

* 0.7.0
    * Temperature must be specified via `TC` or `TK`; `T` is not longer recognized.
* 0.6.1
    * fixed `README`
* 0.6.0
    * uses `StateReporter`
    * reports in kg and molar units
    * default input temperature in Kelvins; T in C allowed via `-TC` argument
* 0.5.1:
    * units bugfix after cli introduction
* 0.5.0:
    * command-line interface
* 0.4.2:
    * `RandomState` class introduced
* 0.3.3:
    * Included subcooled liquid data in the `Request` capability
* 0.2.1
    * bugfix:  set `left` and `right` parameters in all `numpy.interp()` to `numpy.nan` to override the pinning default for extrapolated values
* 0.2.0
    * Added the `Request` class for dynamically selecting and outputting (as LaTeX) steam table blocks requested by (for example) exam problems
* 0.1.9
    * bugfix: allows for specification of `x` as 0
* 0.1.8
    * Update readme
* 0.1.7
    * Updated pandas indexing for saturated table LaTeX printing
* 0.1.5
    * Updated interpolators
* 0.1.1
    * Updated pyproject.toml and README.md
* 0.1.0
    * Initial version

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/sandlersteam/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
