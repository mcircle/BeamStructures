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

## Tangentensteifigkeit für Clamp und Branch

`BeamStructures.admittance_matrix(states, adjacency, structure, values(beams))`
gibt ein `BeamStructures.TangentSystem` zurück. Es enthält gekoppelte
Knotenverschiebungen und innere Balkenkräfte; der Rückgabewert ist **keine**
gewöhnliche `3n × 3n`-Matrix mehr. Dieser API-Wechsel betrifft bestehende
Aufrufer, die den Rückgabewert direkt indizieren oder als Array verarbeiten.

```julia
system = BeamStructures.admittance_matrix(states, adjacency, structure, values(beams))
steifigkeit = BeamStructures.effective_stiffness(system, feste_dofs => belastete_dofs)
nachgiebigkeit = BeamStructures.effective_compliance(system, feste_dofs => belastete_dofs)
```

Die Knotenfreiheitsgrade sind `(ϕ, x, y)` je Knoten; feste und belastete DOFs
beziehen sich nur darauf. Die internen Balkenkraft-DOFs werden automatisch
berücksichtigt. Eine SVD prüft Lösbarkeit und Eindeutigkeit der belasteten
Verschiebung; bei unvereinbaren oder nicht eindeutigen Fällen wird ein Fehler
gemeldet. Die absolute Maschinepsilon-Regularisierung im Kehrwert sowie die
SVD-Rangschwelle sind numerische Konventionen, keine Stabilitätsvalidierung.
Exakt verschwindende Kanten werden aus dem Gleichungssystem entfernt; dort ist
die Topologieauswahl nicht differenzierbar.

Die aktuelle Balken-ODE enthält keinen aktiven axialen Dehnungsterm `F/(EA)`.
Geprüft wurden gerade Balken mit `Clamp` und `Branch`; bei nicht erfülltem
Gleichgewicht ist die berechnete Linearisierung nicht automatisch eine
Gleichgewichtstangente. Der Pullback der Assemblierung liefert Ableitungen nach
der Adjazenz bei festem Zustand und festen Balken. Zustands- und
Geometrieableitungen der Tangentenassemblierung fehlen weiterhin, sodass der
Steifigkeitsverlust darüber keinen vollständigen Netzgradienten liefert.

Die Kraftresiduen für `Branch` verwenden nun das eintreffende Balkenende und
den ausgehenden Balkenanfang. Damit sind ältere Ergebnisse mit dieser
Vorwärtsrechnung nicht ohne erneute Prüfung vergleichbar.

## Tests

```julia
using Pkg
Pkg.test()
```

## Lizenz

Lizenziert unter der [GNU Affero General Public License v3.0](LICENSE.txt).
