# BeamStructures

`BeamStructures` ist eine Julia-Bibliothek zur Simulation ebener Balkenstrukturen.
Sie unterstützt gerade und gekrümmte Balken, verschiedene Randbedingungen sowie
automatische Differentiation und nichtlineare Solver.

## Installation

Im Projektverzeichnis:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Danach kann das Paket geladen werden:

```julia
using BeamStructures
```

## Kurzes Beispiel

```julia
using BeamStructures

beam = Beam(30.0, 1.6, 10.0, 0.0; E = 1.4e5)
fixed = Clamp(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
free = Branch(30.0, 0.0, 0.0, 0.0, 0.0, 0.0)

beams, nodes = BeamStructures.prepare(fixed, beam, free)
structure = Structure([0 1; 1 0])
```

Für eine vollständige nichtlineare Berechnung und Visualisierung siehe
[`examples/Cross_axis_Flexible_Pivot.jl`](examples/Cross_axis_Flexible_Pivot.jl).

## Beispiel: flexibler Kreuzachsen-Pivot

Das Beispiel zeigt einen flexiblen Pivot aus geraden (`Beam`) und gekrümmten
(`CurvedBeam`) Balken. Eine wirkende Kraft wird gradientenbasiert optimiert;
die Berechnung erfolgt mit `NonlinearSolve` und wird mit `GLMakie` visualisiert.

**Ergebnis:**

[Video `optimisation_curved_zero.mp4` ansehen](examples/optimisation_curved_zero.mp4)

<video controls width="720">
  <source src="examples/optimisation_curved_zero.mp4" type="video/mp4">
</video>

## Tests

```julia
using Pkg
Pkg.test()
```

## Lizenz

Lizenziert unter der [GNU Affero General Public License v3.0](LICENSE.txt).
