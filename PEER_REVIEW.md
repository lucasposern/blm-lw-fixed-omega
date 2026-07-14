# Peer Review

**Manuskript:** *Black–Litterman with Ledoit–Wolf Shrinkage: Fixed-Ω Leakage and the Critical Correlation Threshold* (fixed_omega.tex, Stand Juli 2026)
**Review-Modus:** Vollprüfung nach Journal-Maßstab (Quantitative Finance / Journal of Asset Management-Niveau), bewusst streng.
**Verifikationsbasis:** Ich habe nicht nur den Text geprüft, sondern `code/paper_numbers.py` gegen die eingefrorenen Daten in `data/` ausgeführt. **Alle im Paper zitierten Zahlen reproduzieren exakt** (α, μ_s, δ, κ, Δπ, ρ*-Tabelle, ρ**, alle 28 Paar-Shifts, 2×2-Beispiel, Backtest-Tabelle, Memmel-Test, absolute Views +12.9/−0.5, Basket +5.8, c* = 0.446, Floor). Zusätzlich habe ich jede Herleitung von Hand nachgerechnet und die Linearitätsbehauptung numerisch getestet. Das Ergebnis der Rechnung steht in Abschnitt 1.

---

## 1 Wissenschaftliche Korrektheit

### 1.1 Verifiziert korrekt

Die folgenden Kernbestandteile habe ich Schritt für Schritt nachgerechnet; sie sind fehlerfrei.

- **Gl. (1), BL-Posterior:** Standardform, konsistent mit He–Litterman-Masterformel (der Code nutzt die Woodbury-Form, äquivalent).
- **Idzorek-Closed-Form** ω = (1−c)/c · τ · pᵀΣp und die exakte Kürzung c_eff = c: korrekt; τ und v kürzen sich.
- **Herleitung von (⋆) und Gl. (5), Single-View-Identität:** jede Zeile korrekt, inkl. des Prefaktor-Schritts τ(1−c_eff)/ω = c_eff/v.
- **Gl. (6), Gewichtsidentität** und Fußnote 1 („exactly fully invested" für relative Views): korrekt; 1ᵀp = 0 impliziert 1ᵀw = 1ᵀw_eq. Im Code-Output bestätigt (net = 1.0000 in allen drei Varianten).
- **Δc_eff-Algebra (Gl. 2 und 3):** Beide Umformungen stimmen; Substitution ω = (1−c)/c·τ·v^Std kürzt τ vollständig.
- **c\* = 1/(1+√(v^LW/v^Std)):** korrekt (Nullstelle von c²D + 2cb − b = 0 nachgerechnet); numerisch bestätigt (c\* = 0.446, Shift +10.86 pp).
- **Remark 2 (Range):** Grenzen (−c(1−c)α/(1−cα), 1−c), strikte Monotonie in v^Std, Ableitungszähler −c(1−c)αμ_s‖p‖²: alles korrekt nachgerechnet. Stetige Fortsetzung auf den PD-Rand korrekt argumentiert.
- **Gl. (4), ρ\*:** Herleitung korrekt; alle 16 Tabellenwerte unabhängig nachgerechnet (Abweichung < 0.001). Die Behauptung, die Tabelle hänge nur von μ_s/σ² und c ab, folgt aus der Formel und stimmt.
- **ρ\*\*-Existenzbedingung** ε < c(1−c)α/(1−cα) und die Werte −0.32 bzw. −25.3: korrekt.
- **Δπ-Formel** (Abschnitt 5): algebraisch korrekt, Δπ = αδ(μ_s I − Σ̂_s)w_eq komponentenweise ausgeschrieben. (Aber: die verbale Interpretation ist falsch, siehe E2.)
- **Zerlegung Return-/Precision-Channel** (Abschnitt 7): exakte Identität (add-and-subtract), korrekt.
- **Eigenwertaussage** λ_k^LW = (1−α)λ_k + αμ_s: korrekt, weil Shrinkage Richtung μ_s·I die Eigenbasis erhält. Der „Faktor drei" für 1/λ_n ist konsistent mit κ: 341.8/108.0 = 3.2.
- **Remark 1 (Mehr-View-Absorption):** A = S(S+Ω)⁻¹ korrekt hergeleitet; die Diagnostik-Zahlen (Residuen ≤ 0.535 % < 0.6 %, Absorption innerhalb 1 pp von 50 %, Bewegung < 0.3 pp) reproduzieren.
- **Verallgemeinerung auf beliebige p** via ‖p‖²: korrekt (pᵀI p = ‖p‖²); Assertion im Code bestätigt die Identität maschinengenau.

### 1.2 Gefundene Fehler

**E1 (sachlicher Fehler, mehrfach im Text): „both exactly linear in α" ist falsch.**
Betroffen: Abstract („the remaining effect on the posterior runs through the equilibrium prior and the view projection, both exactly linear in α"), Abschnitt 5 („Both effects are exactly linear in α"), Conclusion („prior and projection shifts that are exactly linear in α"), und implizit Abschnitt 7 („Δπ = O(α) is the only remaining term").
Sachverhalt: π(α) = δ[(1−α)Σ̂_s + αμ_s I]w_eq ist exakt linear in α. Die normierte Projektion Σ(α)p/v(α) ist es nicht: Zähler und Nenner sind affin in α, der Quotient ist eine rationale Funktion und nur zu erster Ordnung linear (exakt linear nur, wenn p Eigenvektor von Σ̂_s ist). Ich habe das numerisch bestätigt: für μ_BL(α) unter Idzorek-Ω gilt max|μ(0.5) − ½(μ(0)+μ(1))| ≈ 4.1·10⁻³, während der reine Prior maschinengenau linear ist (10⁻¹⁸).
Korrektur: „the prior shift is exactly linear in α; the view projection Σ(α)p/v(α) is a smooth rational function of α, linear to first order." Alle vier Stellen anpassen. Zusätzlich widerspricht die Formulierung in Abschnitt 7 („the only remaining term") der eigenen Diagnostik in Abschnitt 8, wonach vom maximalen Posterior-Shift von 0.42 % nur höchstens 0.21 % aus dem Prior stammen, der Rest aus der Projektion. Der Projektionsterm ist also gleich groß wie der Priorterm und darf nicht wegdefiniert werden.

**E2 (sachlicher Fehler, von den eigenen Daten widerlegt): Interpretation von Δπ.**
Zitat (Abschnitt 5): „Assets with above-average variance (σ_i² > μ_s) receive a lower prior, assets with below-average variance receive a higher one."
Diese Aussage beschreibt nur den Diagonalterm (μ_s − σ_i²)w_eq,i und ignoriert den Kreuzterm −Σ_{j≠i} Σ̂_s,ij w_eq,j, der in equity-lastigen Universen dominiert. Auf den eigenen Daten: IWF hat σ² = 0.0241 < μ_s = 0.0282 (unterdurchschnittlich), aber Δπ_IWF = **−0.113 %** < 0; ebenso IWD (σ² = 0.0244 < μ_s, Δπ = −0.114 %). Zwei der acht Assets widerlegen den Satz.
Korrektur: Δπ_i = αδ(μ_s w_eq,i − Cov(r_i, r_p)) mit r_p = w_eqᵀr. Shrinkage ersetzt die Kovarianz jedes Assets mit dem Gleichgewichtsportfolio anteilig durch μ_s w_eq,i; Assets mit hoher Kovarianz zum Marktportfolio verlieren Prior-Rendite, Assets mit niedriger (Bonds) gewinnen. Diese Lesart erklärt das Vorzeichenmuster der Daten (AGG +0.042, BWX +0.009, alle Equities negativ) vollständig.

**E3 (falsche Genauigkeitsangabe): „Both predictions match the numerical Δμ and Δw to within a few hundredths of a percentage point."**
Für Δμ stimmt es (Abweichung ≈ 0.01 pp). Für Δw nicht: Vorhersage ±30.0 pp, gemessen (+30.2, −29.8) pp, Abweichung 0.2 pp, also Zehntel, nicht Hundertstel. Korrekt wäre: „to within ≈ 0.01 pp for Δμ and ≈ 0.2 pp for Δw".

**E4 (inkonsistente Rundung): Deflations-Floor.**
Remark 2 behauptet „no view in any universe can lose more than 1.1 percentage points" bei α = 0.041, c = 0.50. Der Floor ist c(1−c)α/(1−cα) = 1.05 pp (mit dem tatsächlichen α = 0.0406: 1.04 pp). Als obere Schranke ist 1.1 nicht falsch, aber die natürliche Rundung wäre 1.0; im selben Satz wird −0.4 pp auf eine Dezimale gerundet. Einheitlich runden oder 1.05 schreiben.

**E5 (methodische Unschärfe in 6.6): „Differencing … at fixed π, Σ and ω".**
Step 1 leitet Δμ_BL = Δc_eff(Q − pᵀπ)Σp/v „bei festem Σ" ab, wendet das Ergebnis aber auf das Experiment an, in dem Σ von Σ̂_s auf Σ̂_LW wechselt und damit die Richtung Σp/v selbst. Exakt gilt über den Estimatorwechsel hinweg nur die Projektion pᵀΔμ = Δc_eff(Q − pᵀπ) (weil pᵀ(Σp/v) = 1 für jedes Σ). Die Vektorgleichung ist eine Approximation mit Fehler von der Ordnung α mal Volatilitäts-Heterogenität. Der Text vermischt zwei Näherungen (fixe Richtung; Σp/v ≈ p/‖p‖²) und deklariert in Step 2 die gemessene Größe als „exact partial derivative at fixed Σ", obwohl die Messung den Richtungseffekt enthält (sichtbar an der Asymmetrie +30.2/−29.8). Empfehlung: erst die exakte Along-p-Identität festhalten, dann die Vektorapproximation mit explizitem Fehlerterm, dann die numerische Bestätigung.

**E6 (interner Widerspruch): „The answer depends entirely on a bookkeeping question … what happens to Ω?" (Einleitung).**
Abschnitt 7 stellt selbst fest, dass der größte Effekt des Estimatorwechsels im Precision-Channel Σ⁻¹ liegt, der von Ω unabhängig ist. Die Antwort hängt also gerade nicht „entirely" an Ω, sondern nur im View-Channel. Formulierung präzisieren: „Within the view channel, the answer depends entirely on …".

### 1.3 Präzisierungsbedarf (nicht falsch, aber unter Journal-Standard)

- **Das probabilistische Modell fehlt.** Gl. (1) fällt vom Himmel. Ein Journal erwartet die Bayes-Struktur: Prior μ ~ N(π, τΣ), Views Q = Pμ + ε, ε ~ N(0, Ω), ε ⊥ μ, Σ als bekannt behandelt. Gerade Letzteres ist hier zentral: Das Paper untersucht Kovarianz-*Schätzung* in einem Modell, dessen Bayes-Herleitung Σ als bekannt voraussetzt. Diese Spannung (Plug-in-Bayes) muss einmal explizit benannt werden.
- **„reliable only if T ≫ n; the relative estimation error grows with n/T"** ist Folklore-Niveau. Entweder ein präzises Resultat zitieren (z. B. LW 2004, Lemma/Theorem zur erwarteten Frobenius-Verlustrate; oder Marchenko–Pastur für die Eigenwertspreizung) oder als Heuristik kennzeichnen.
- **Oracle-α:** „minimises E‖Σ̂ − Σ‖_F²" ist verkürzt; LW 2004 minimieren über der Konvexkombination mit Zielgrößen aus Populationsmomenten, und der analytische Schätzer ersetzt vier Momente konsistent. Die Schätzformel für α̂ fehlt komplett (nur „scikit-learn"); für Selbstständigkeit des Papers gehört sie in einen Appendix, inkl. sklearn-Versionspinning.
- **Idzorek-Äquivalenz:** „the closed form is its exact single-view solution" wird behauptet, nicht gezeigt. Der Beweis ist zweizeilig (Tilt(ω) = c_eff · Tilt_100 via Gl. (6); Matching c·Tilt_100 ⟺ ω = (1−c)/c·τv) und gehört in einen Appendix, zumal die eigene Codebasis (views.py) die Kalibrierung ausdrücklich als „vereinfacht, nicht das iterative Tilt-Matching" kennzeichnet.
- **α = 0.041 vs. 0.0406:** Der Code liefert 0.0406. Einmal exakt angeben („α̂ = 0.0406, im Folgenden 0.041"), sonst wirken Ableitungen wie 1.05 pp vs. 1.04 pp inkonsistent.
- **σ_i² := Σ̂_s,ii** wird beiläufig in einer Klammer definiert; in die Notationsfestlegung von Abschnitt 2 ziehen.

---

## 2 Inhaltliche Logik

Der Gesamtaufbau ist schlüssig: Modell → Problem → Estimator → inerter Fall → Leck-Fall → Propagation → Kanal-Trennung → Backtest. Die Trennung View-Channel vs. Precision-Channel ist die stärkste konzeptionelle Leistung des Papers. Konkrete Brüche:

1. **Der Begriff „leakage" wird informell eingeführt** („We call the resulting silent shift of the effective inputs leakage", Einleitung), aber nie formal definiert. Der natürliche Kandidat ist Δc_eff ≠ 0 bei konstantem c. Eine Definition (ein Satz, ggf. als Definition-Umgebung) würde Abstract, Einleitung und Abschnitt 6 verankern.
2. **Konzeptionelle Lücke, die ein Referee sicher anspricht: ω-Primitiv vs. c-Primitiv.** Das Paper unterstellt durchgehend, das Primitiv des Analysten sei die Konfidenz c, und ein eingefrorenes ω „verrate" c. Wenn aber ω selbst das Primitiv ist (eine absolute Varianz des View-Rauschens, unabhängig vom Kovarianzmodell), dann ist das Einfrieren von ω Bayesianisch kohärent, und die Änderung von c_eff ist die korrekte Antwort auf eine bessere Kovarianzschätzung; „recompute Ω" wäre dann der Fehler. Die normative Empfehlung in 6.9 gilt nur unter c-Primat. Das Paper muss diese Alternative benennen und begründen, warum es c als Primitiv setzt (Idzoreks Verfahren und die Risk-Governance-Formulierung in 6.9 liefern die Begründung, sie muss nur explizit werden). Ohne diese Diskussion liest sich die Kernbotschaft angreifbar.
3. **Motivationszirkel:** Einleitung: „Two workflows coexist in practice." Abschnitt 6: „We are not aware of survey evidence on how widespread this workflow is." Der zweite Satz ist ehrlich, der erste damit unbelegt. Entweder Belege liefern (naheliegend: Software-Schnittstellen, in denen Ω als nutzerdefinierte Konstante übergeben wird, während Σ separat geschätzt wird, z. B. gängige Open-Source-BL-Implementierungen) oder die Behauptung konditional formulieren.
4. **Spannung zwischen Abschnitt 5 und 8** (siehe E1/E6): Abschnitt 5 verkauft den Idzorek-Fall als „inert bis auf linearen Prior", die Diagnostik in 8 zeigt, dass der Projektionsterm gleich groß ist wie der Priorterm. Kein Widerspruch in den Zahlen, aber in der Rhetorik.
5. **Parameterwechsel mitten in 6.5:** Tabelle 1 nutzt Illustrationsparameter (μ_s = 4·10⁻³, σ² = 3.5·10⁻³), Tabelle 2 die Datenwerte. Der Wechsel ist begründet, aber der Leser stolpert; ein Satz, warum nicht die Datenwerte in Tabelle 1 stehen (Antwort: nur die Ratio zählt, und 8/7 ist repräsentativ), steht zwar da, dennoch wäre eine Zeile mit den Daten-kalibrierten ρ*-Werten die überzeugendere Wahl.
6. **Gedankensprung in 6.2:** „As v_k grows, so does the absorption c_eff,k." Das gilt bei festem ω, was genau der Punkt ist, steht aber erst im Satz danach implizit. Ein „at fixed ω_k" im Satz schließt die Lücke.

---

## 3 Vollständigkeit

Fehlende Elemente, jeweils mit Begründung, Umfang und Ort:

1. **Formales Modell + Annahmenblock** (siehe 1.3). Warum: ohne Bayes-Struktur sind Ω, τ und „confidence" nicht wohldefiniert; die Plug-in-Spannung ist der konzeptionelle Kern des Papers. Umfang: halbe Seite. Ort: Abschnitt 2, vor Gl. (1).
2. **Proposition/Beweis-Struktur.** Die Theorem-Umgebungen sind im Präambel definiert und werden für die Hauptresultate nie benutzt; Gl. (3), (4), (5) und Remark 2 sind die Sätze des Papers und sollten als Propositionen mit expliziten Voraussetzungen (Σ ≻ 0, ω > 0, K = 1, p ≠ 0) formuliert werden. Warum: Prüfbarkeit und Zitierbarkeit. Umfang: Umformatierung, kein neuer Inhalt. Ort: Abschnitte 5–6.
3. **α̂-Formel + Versionspinning** (Appendix, halbe Seite). Warum: Reproduzierbarkeit ohne sklearn-Quellcode-Lektüre; die (T−1)/T-Konvention ist bereits vorbildlich dokumentiert, die Schätzformel fehlt.
4. **Beweis der Idzorek-Äquivalenz** (Appendix, fünf Zeilen). Warum: tragende Behauptung von Abschnitt 5.
5. **Long-only-Ergebnisse.** Abschnitt 8 behauptet: „Imposing a long-only constraint instead compresses both variants toward each other and toward far tamer risk figures" — ohne eine einzige Zahl. Der Code berechnet die Long-only-Zwillinge (cvxpy), das Paper berichtet sie nicht. Entweder Tabelle ergänzen (eine Zeile pro Variante genügt) oder Satz streichen. Warum: unbelegte empirische Aussage ist ein klassischer Ablehnungsgrund.
6. **Robustheit/Sensitivität.** Fehlt vollständig: (a) Rolling-Origin-Backtest statt einem Split; (b) Sensitivität der Backtest-Ergebnisse in α (z. B. α ∈ {0.02, 0.041, 0.10, 0.30}), c und Q; (c) alternative View-Sets (systematisches Gitter oder zufällige Paare), um die Hindsight-Abhängigkeit der drei gewählten Views zu neutralisieren. Umfang: 1–2 Seiten + 1 Abbildung. Ort: nach 8.
7. **Schätzunsicherheit der Leckage selbst.** Δc_eff ist eine deterministische Funktion von Σ̂ und α̂, aber beide sind Schätzer; ein Bootstrap (Block-Bootstrap über Monate) liefert Konfidenzbänder für Δc_eff pro Paar und beantwortet die Frage, ob +10.7 pp von +8.8 pp unterscheidbar ist. Umfang: halbe Seite + Tabellenspalte. Ort: 6.5.
8. **w_eq-Werte.** Die Idzorek-Marktgewichte werden benutzt, aber nirgends abgedruckt. Tabelle in Abschnitt 2 oder Appendix. Warum: Selbstständigkeit; π und alle Gewichte hängen daran.
9. **Datenkonstruktion:** Wie wurde die 13-Wochen-T-Bill-Rendite in eine Monatsrate umgerechnet (Diskontsatz-Konvention, /12 oder (1+y)^(1/12)−1)? Eine Fußnote. Warum: Excess-Return-Definition betrifft jede Zahl.
10. **Abbildungen.** Null Abbildungen im gesamten Paper. Mindestens: (i) Δc_eff als Funktion von v^Std/μ_s‖p‖² mit den 28 Paaren als Punkte (macht Sign-Kriterium, Asymmetrie und Bond-Paar-Pointe auf einen Blick sichtbar); (ii) kumulierte Out-of-sample-Pfade der vier Portfolios inkl. Drawdown-Schattierung. Warum: Journal-Standard; die Bond-Paar-Botschaft (AGG/BWX bei ρ = 0.57) trägt visuell weiter als jede Tabelle.
11. **Diskussion Live-Rebalancing.** Das Paper framet die Leckage als einmaligen Estimatorwechsel. In der Praxis wird Σ̂ rollierend neu geschätzt; ohne Rekalibrierung von Ω tritt die Leckage dann an jedem Rebalancing-Termin auf. Ein Absatz in der Conclusion erweitert die Relevanz erheblich.
12. **Related-Work-Absatz zur BL-Kovarianz-Kombination:** Es gibt empirische Literatur, die BL mit Shrinkage kombiniert (siehe Abschnitt 6 dieses Reviews); die Aussage „receives comparatively little attention" braucht diese Abgrenzung, sonst wirkt der Neuheitsanspruch unrecherchiert.

---

## 4 Quantitative Finance

- **Black–Litterman korrekt dargestellt?** Ja, in der „canonical reference model"-Variante. Die Wahl w = (δΣ)⁻¹μ_BL statt Σ_post = Σ + M ist offengelegt (Fußnote 1) und konsistent durchgehalten, inkl. der Fully-invested-Konsequenz. Gut. Es fehlt die Einordnung, dass diese Wahl die Ergebnisse quantitativ beeinflusst (mit Σ + M wäre c_eff um den Faktor ≈ 1/(1+τ) verschoben); ein Satz genügt.
- **Herleitung vollständig?** Die Single-View-Kette ja (und sauber). Die Bayes-Herleitung von Gl. (1) fehlt (siehe 3.1). Die Mehr-View-Welt existiert nur als Remark; für die Kernaussage akzeptabel, weil Scope klar begrenzt wird.
- **Bayesianische Interpretation sauber?** Funktional ja (c_eff als Posterior-Gewicht entlang der View-Richtung ist korrekt und schön interpretiert), formal nein (Modell nie aufgeschrieben; „confidence" nie als Modellgröße definiert, nur operativ über Idzorek).
- **Ledoit–Wolf korrekt beschrieben?** Ja; die ddof-Konvention und die (T−1)/T-Falle von sklearn sind vorbildlich dokumentiert (Abschnitt 4 und Reproducibility). Warum Shrinkage funktioniert (Bias-Varianz-Abtausch, Eigenwertspreizung) wird zitiert, aber nicht erklärt; zwei Sätze mehr wären angemessen.
- **Wann wird Shrinkage problematisch?** Indirekt behandelt (Deflationsregime, PD-Rand, Remark 2: „LW regularises exactly the near-degenerate view directions"), das ist eine der besten Beobachtungen des Papers. Nicht behandelt: Fälle, in denen das Identitäts-Target selbst schlecht ist (stark heterogene Varianzen, Faktorstruktur), obwohl genau dann μ_s‖p‖² als Schwelle ökonomisch fragwürdig wird. Ein Absatz.
- **Mathematik vs. Ökonomie getrennt?** Überwiegend ja; Verstöße sind markierte Overclaims (E1, E2, E6) und die unbelegte Long-only-Aussage.
- **Alle Modell-Einschränkungen genannt?** Fehlend: Normalität, Σ-bekannt-Spannung, Ein-Perioden-Rahmen, Proxy-Inkonsistenz (δ aus S&P 500, w_eq aus Idzorek-Assetklassen, Universum aus ETFs: drei verschiedene „Märkte" in einem Gleichgewichtsargument; als bewusste Praxis-Approximation kennzeichnen).

---

## 5 Empirie

- **Datensatz:** Ein Vendor (Yahoo, adjusted closes), keine Qualitätsprüfung berichtet (Splits/Dividenden-Adjustierung, Stale Prices). Frozen CSVs im Repo sind gut (Reproduzierbarkeit), ersetzen aber keine Validierung gegen eine zweite Quelle (CRSP, Bloomberg, oder mindestens Vergleich der ETF-NAV-Returns vom Anbieter). Umfang: eine Fußnote plus Robustheits-Satz.
- **Stichprobe:** T = 132, n = 8. Der Selbstwiderspruch ist adressiert (α klein, „the analytic intensity is modest"), aber die Konsequenz fehlt: Bei n/T = 0.06 ist die Motivationserzählung aus Abschnitt 3 (Konditionierungsproblem) schwach; die Leckage ist hier groß, *obwohl* α klein ist, weil v^Std ≪ 2μ_s. Diese Pointe (Leckage-Treiber ist die Varianzstruktur, nicht die Dimensionalität) steht implizit in 6.5 und sollte prominenter werden, weil sie die Wahl des kleinen Universums rechtfertigt.
- **Look-ahead/Survivorship:** w_eq (2004-Vintage) ist sauber look-ahead-frei, gut argumentiert. ETF-Universum hat kein Survivorship-Problem im engeren Sinn (alle acht existieren durchgehend). Aber: **die View-Auswahl ist ex post.** Das Paper wurde 2026 geschrieben; „deliberately not tuned to the test window" ist eine Absichtserklärung, kein Design. Zwei von drei Views stimmen, und der dominante IWF>IWD-View liegt mit +9.1 % realisiert weit über den behaupteten +3 %. Für die Kernaussage (Vergleich frozen vs. recomputed bei identischen Views) ist das unkritisch, für die Niveaus (Sharpe 0.72) nicht. Lösung: View-Gitter oder Zufallsviews als Robustheitscheck, und die Niveaus explizit als nicht interpretierbar kennzeichnen (steht teilweise da, kann schärfer).
- **Parameterwahl:** τ = 0.05, c = 0.50 begründet und (für das Kernresultat) irrelevant, gut. δ = 2.44 aus S&P statt aus dem Universum: dokumentieren, warum (Konvention), und die Sensitivität ist trivial (δ skaliert alle Gewichte), ein Satz genügt.
- **Benchmarks:** Market-Equilibrium-Spalte ja; es fehlen 1/N und Minimum-Variance als Kontext (DeMiguel et al. werden zitiert, aber nicht als Benchmark umgesetzt). Optional, aber ein Referee fragt danach.
- **Metriken:** Arithmetische Annualisierung + MDD auf kumuliertem Excess-Pfad: offengelegt, aber MDD auf Excess-Basis ist unüblich (Drawdown wird ökonomisch auf Total-Return-Pfaden erlebt); mindestens Fußnote, besser zusätzlich Total-Return-MDD.
- **Statistische Aussagekraft:** Ein Memmel-Test für ein Spaltenpaar. Fehlend: (a) Test frozen vs. recomputed (die eigentliche Vergleichsgröße des Papers!); (b) HAC-robuste Inferenz nach Ledoit–Wolf (2008) statt Memmel, da 2019–2024 den COVID-Crash enthält und Monatsrenditen weder iid noch normal sind (LW 2008 wird zitiert, aber nicht benutzt); (c) Konfidenzintervalle für die Kernzahlen (Bootstrap, siehe 3.7).
- **Data Snooping:** Ein Split, ein Universum, ein Target, ein View-Set, keine Multiple-Testing-Fragen; das Design ist zu klein für formales Snooping, aber die View-Hindsight-Frage (oben) bleibt.
- **Gesamturteil Empirie:** illustrativ solide und vollständig reproduzierbar (großes Plus, selten in Einreichungen), aber als empirische Evidenz dünn: eine Beobachtung pro Aussage.

---

## 6 Literatur

**Formales:** Alle 18 Einträge habe ich geprüft; Bände, Hefte und Seitenzahlen stimmen (inkl. Memmel, Finance Letters 1(1), und der a/b-Disambiguierung LW 2004). Attributionen im Text sind korrekt (Walters für die Closed Form, He–Litterman für die Gewichtszerlegung, Jagannathan–Ma S. 1652 passt zum Constraints-Argument).

**Fehlende Standardreferenzen**, mit Begründung und Ort:

| Referenz | Warum nötig | Ort |
|---|---|---|
| Jorion (1986), *Bayes–Stein estimation for portfolio analysis*, JFQA 21(3) | Das Paper handelt von Shrinkage im Mittel (BL) × Shrinkage in der Kovarianz (LW); Bayes–Stein ist der kanonische Vorläufer der Mean-Shrinkage-Seite | Einleitung, Satz zur Input-Sensitivität |
| Frost & Savarino (1986), *An empirical Bayes approach to efficient portfolio selection*, JFQA 21(3) | dito, Prior auf beiden Inputs | Einleitung |
| Kan & Zhou (2007), *Optimal portfolio choice with parameter uncertainty*, JFQA 42(3) | Formalisiert genau die Plug-in-Problematik (Σ̂ statt Σ), die das Paper implizit behandelt | Abschnitt 2 oder 3 |
| Meucci (2010), *The Black–Litterman approach: original model and extensions*, in Encyclopedia of Quantitative Finance | Standard-Survey; ohne ihn wirkt der BL-Literaturüberblick lückenhaft | Einleitung |
| Kolm & Ritter (2017), *On the Bayesian interpretation of Black–Litterman*, EJOR 258(2) | Das Paper argumentiert Bayesianisch, ohne die maßgebliche Klärung der Bayes-Struktur zu zitieren; direkt relevant für die ω-Primitiv-Frage | Abschnitt 2 |
| Michaud, Esch & Michaud (2013), *Deconstructing Black–Litterman*, Journal of Investment Management 11(1) | Prominente Kritik; Evenhandedness des Literaturbilds | Einleitung |
| Bevan & Winkelmann (1998), *Using the Black–Litterman global asset allocation model* (Goldman Sachs) | Praxis-Workflow-Beleg; stützt die Motivationslücke „two workflows coexist" | Einleitung/Abschnitt 6 |
| Herold (2003), *Portfolio construction with qualitative forecasts*, JPM 30(1) | Alternative Ω-Spezifikationen; zeigt, dass Ω-Wahl ein eigenes Literaturfeld ist | Abschnitt 5 |
| Bessler, Opfer & Wolff (2017), *Multi-asset portfolio optimization and out-of-sample performance…*, European Journal of Finance 23(1) | Empirische BL-Studie inkl. Kovarianz-Varianten; nötig für die Behauptung, die Kombination sei wenig beachtet | Einleitung |
| Ledoit & Wolf (2020), *Analytical nonlinear shrinkage*, Annals of Statistics 48(5) oder LW (2022), *The power of (non-)linear shrinking*, Journal of Financial Econometrics 20(1) | Aktueller Stand der Shrinkage-Literatur; das Paper zitiert 2017 als letzten Stand | Einleitung/Conclusion |
| DeMiguel, Martín-Utrera & Nogales (2013), *Size matters: Optimal calibration of shrinkage estimators*, JBF 37(8) | Kalibrierung von Shrinkage speziell für Portfoliozwecke; nahe am Thema Intensitätswahl | Abschnitt 4 |

**Einbindung:** Die vorhandene Literatur ist gut verwoben (kein Zitat-Dumping). Mit den Ergänzungen sollte die Einleitung einen kurzen expliziten Related-Work-Absatz bekommen (5–8 Sätze), der drei Stränge trennt: BL-Spezifikationssensitivität (τ, Ω), Kovarianz-Shrinkage im Portfoliokontext, Parameterunsicherheit/Plug-in. Das Paper positioniert sich dann als Schnittpunkt.

---

## 7 Sprache

Gesamturteil vorab: Das Englisch ist auf muttersprachlichem Fachniveau, die Sätze tragen Information, es gibt fast keine Grammatikfehler. Das Problem ist nicht Kompetenz, sondern Register: Der Text ist stellenweise Essay, nicht Journal-Artikel. Konkrete Befunde (Prüfung nach stop-slop-Kriterien):

**7.1 Formelhafte Kontrastkonstruktionen („not X but Y" / „X rather than Y")** treten gehäuft auf und wirken ab dem dritten Auftreten wie ein Tick:
- „The practical lesson is procedural rather than statistical."
- „The larger effect … lies not in the return estimator but in the inverse"
- „a structural sketch rather than the general statement"
- „the list is informed conjecture rather than a tested claim"
- „is a property of frozen scalars, not of a particular calibration formula"
- „The point is not academic"
Jeweils zwei oder drei stehen lassen, den Rest direkt formulieren. Beispiel: „The point is not academic: the analyst who freezes Ω has…" → „An analyst who freezes Ω has therewith raised the confidence of every affected view…".

**7.2 Rhetorische Frage:** „Where does LW earn its keep?" → „Three configurations plausibly favour LW." (Der Folgesatz existiert bereits; die Frage ist verzichtbar.)

**7.3 Pull-Quote-Sätze** (zu feuilletonistisch für ein Quant-Journal, je einzeln vertretbar, in Summe zu viel):
- „The model believes the view more strongly than the analyst stated."
- „…and the leakage rides exactly this lever."
- „The analyst who freezes Ω forfeits close to half of the de-risking they think they bought."
- „This is the dimension where LW carries its own weight"
Vorschläge: „The effective confidence exceeds the stated confidence." / „The leakage operates through the same factor 1/(δv_k)." / „Freezing Ω forfeits roughly 40 % of the exposure and drawdown reduction." (nebenbei: „close to half" ≠ 40 %, siehe T8) / „The drawdown reduction is the dimension in which LW adds value."

**7.4 Kolloquialismen:**
- „A different workflow is easy to fall into." → „A second workflow arises when …"
- „We can keep this case brief." → streichen; der Abschnitt ist ohnehin kurz.
- „made the crossing from journals to practice" → „is widely used in practice"
- „The absolute position sizes are themselves a comment on unconstrained Black–Litterman" → „The absolute position sizes illustrate a known property of unconstrained Black–Litterman with highly correlated pairs"
- „prices the frozen-Ω workflow directly" / „prices the leakage" (zweimal) → „quantifies"
- „on the input form" → einmal in Ordnung als Bild, nicht zweimal.

**7.5 Wiederholungen:** „Two specialisations are worth recording" (6.2) und „Three consequences are worth recording" (6.3) direkt hintereinander; eine der beiden Formulierungen ersetzen („Two special cases follow immediately").

**7.6 Em-Dash-Dichte:** Der Text setzt Gedankenstriche als Standard-Interpunktion ein („substantially---on our data…", „five percentage points---among them…", „by +10.7---the model operates…"). In dieser Dichte ermüdet es; die Hälfte durch Punkt, Komma oder Klammer ersetzen.

**7.7 Passiv:** Der Data-Absatz reiht Passiva („are taken", „are converted", „is reserved", „are the … weights published"). Für Methodenprosa akzeptabel; zwei aktive Umformulierungen lockern („We take monthly total returns from Yahoo Finance and convert them …").

**7.8 Adverbien:** „exactly" ist als Fachwort (exakt vs. approximativ) gerechtfertigt und bleibt. Streichbar: „squarely" („sits squarely in the homogeneous special case"), „deliberately", „silently" (zweimal; einmal reicht), „routinely" (Abstract, erster Satz), „essentially" („stays essentially put").

**7.9 Grammatik/Mikro:**
- „the workflow behind Table 2" — „behind" ist salopp; „the design underlying Table 2".
- „such as the unusually wide eigenvalue spreads of 2008--2009) that flows into the sample covariance and propagates into the test window" — die Anomalie propagiert nicht ins Testfenster, sondern in die Gewichte; umformulieren („…that enters the sample covariance and distorts the weights held during the test window").
- Einheitlich „percentage points"/„pp": Abstract nutzt einmal „points" („by +12.9 points").

---

## 8 Wissenschaftlicher Stil

**Unbelegte oder übertriebene Behauptungen** (jeweils abschwächen oder belegen):
1. „one of the few pieces of portfolio theory that made the crossing from journals to practice" — anekdotisch; entweder Beleg (Umfragen zur Praxisverbreitung, z. B. in Meucci 2010 referenziert) oder neutral: „is widely implemented in practice".
2. „The standard remedy is the shrinkage estimator of Ledoit and Wolf (2004)" — *ein* Standard, nicht *der*; Faktor-Modelle, EWMA und Sparse-Methoden sind gleichrangige Alternativen. „A standard remedy" genügt als Fix.
3. „Two workflows coexist in practice" — siehe 2.3.
4. „depends entirely on a bookkeeping question" — siehe E6.
5. „a correlation level routine for style pairs within one equity market" — stimmt, aber ein Zitat oder der Verweis auf die eigene Tabelle 2 (ρ = 0.92 IWF/IWD) macht es selbsttragend.
6. „a drawdown above 50 % is not investable for many institutional mandates" — plausibel, unbelegt; „regardless of the realised return" verschärft zusätzlich. Entweder Referenz (z. B. Grossman–Zhou-artige Drawdown-Constraints-Literatur) oder als Einschätzung kennzeichnen.
7. „The intensity α itself grows with n/T" — als Faktum formuliert, ist aber eine Heuristik über den analytischen Schätzer; „tends to grow" plus Verweis auf die Asymptotik in LW (2004).

**Positiv hervorzuheben** (bitte beibehalten): die konsequente Kennzeichnung von Unsicherheit („We are not aware of survey evidence", „informed conjecture rather than a tested claim", „dated but look-ahead-free"), die Offenlegung aller schmeichelnden Annahmen (Kosten, Leverage, Financing) und das Reproducibility-Statement. Das ist überdurchschnittliche wissenschaftliche Hygiene.

**Objektivität/Neutralität:** insgesamt gewahrt; die normative Sprache („the clean fix", „self-healing") ist unter c-Primat vertretbar, braucht aber die Diskussion aus 2.2.

---

## 9 Format

- **Abbildungen:** keine. Für ein Quant-Journal ein Ausschlusskriterium; siehe 3.10.
- **Struktur:** 9 Abschnitte, klare Gliederung, Roadmap-Absatz vorhanden. Abschnitte 3 und 4 sind je ~10 Zeilen; zu „Preliminaries" zusammenlegen.
- **Theoremumgebungen:** definiert, aber nur für zwei Remarks genutzt; Hauptresultate als Proposition 1 (Idzorek-Invarianz), Proposition 2 (Δc_eff, Sign, τ-Freiheit), Corollary (ρ*), Proposition 3 (Single-View-Identität) setzen.
- **Notation:**
  - **k-Kollision:** k ist View-Index (ω_k, c_k, v_k), in Tabelle 2 aber Volatilitätsratio „k = σ_i/σ_j". Ratio umbenennen (z. B. κ_vol wäre wiederum mit κ(Σ) kollidiert; nimm r oder h).
  - **Superskript-Inkonsistenz:** Matrizen tragen Subskripte (Σ̂_s, Σ̂_LW), Skalare Superskripte (v^Std, v^LW), und „Std" vs. „s" bezeichnen dasselbe. Vereinheitlichen (v_s, v_LW).
  - ρ** wird im Fließtext eingeführt; bei Formalisierung als Teil der Proposition ausweisen.
  - (⋆)-Tag: bei Nummerierung aller referenzierten Gleichungen auflösen.
  - Symbolverzeichnis: bei dieser Dichte (α, τ, δ, c, c_eff, μ_s, v_k, ω_k, κ, ρ*, ρ**) eine halbe Seite im Appendix wert; von QF-Journals nicht verlangt, vom Leser dankbar angenommen.
- **Tabellen:** technisch sauber (booktabs, Captions vollständig, Konventionen in den Captions wiederholt — gut). Tabelle 1: Spaltenkopf „α = 0.166" ist ein unerklärter Wert; entweder begründen (woher stammt 0.166? vermutlich ein Thesis-Artefakt) oder auf runde Werte wechseln.
- **Pseudocode:** ein Python-Listing im Fließtext ist für JPM/JAM in Ordnung, für QF/JBF unüblich; in den Appendix oder aufs Repo verweisen.
- **Zitierstil:** natbib round, konsistent; Bibliographie handgepflegt mit korrektem Hinweis auf Umstellung bei Submission. \date{\today} vor Submission fixieren.
- **Abstract:** ~230 Wörter; QF erlaubt das, viele Journals kappen bei 150–200. Die Mehr-View-Aussage und die +12.9-Zahl können raus, ohne dass die Kernbotschaft leidet.
- **Titel/Keywords/JEL:** präzise und passend (G11, C13, C58 korrekt gewählt).

---

## 10 Originalität

**Beitrag:** Die Identifikation, dass die Ω-Behandlung beim Estimatorwechsel eine exakt quantifizierbare, τ-freie Konfidenzverschiebung erzeugt, samt Sign-Kriterium (v vs. μ_s‖p‖²), Range-Asymmetrie und Korrelationsschwelle, ist meines Wissens in dieser Isolation neu. Ich kenne keine Arbeit, die diesen Mechanismus formalisiert; die nächstliegenden (Walters 2014 zur (τ,Ω)-Sensitivität, Satchell–Scowcroft) beschreiben Sensitivität, keine Closed Form für den Workflow-Bruch. *Unsicherheitsvermerk:* Eine systematische Prüfung gegen die graue Literatur (SSRN, Praktiker-Notizen zu PyPortfolioOpt/Bloomberg-PORT-Workflows) habe ich nicht durchgeführt; der Autor sollte mit den Suchbegriffen „Black-Litterman omega recalibration / stale omega / view uncertainty shrinkage" gegenprüfen.

**Grenzen des Beitrags:** Die Mathematik ist elementar (Bruchrechnung über gemeinsamem Nenner, eine Sherman-Morrison-freie Single-View-Identität, die in der Substanz bei He–Litterman/Walters vorliegt). Der Wert liegt in der Identifikation, Exaktheit und Quantifizierung, nicht in der Technik. Das ist ehrlich einzuordnen: Es trägt ein prägnantes Practitioner-Paper oder eine Research Note, kein Full Article in einem Top-Field-Journal. Die interessanteste Einzelbeobachtung ist die Bond-Paar-Pointe (AGG/BWX: großer Shift bei moderater Korrelation) und der Absolute-View-Fall (+12.9 pp auf AGG), weil beide die naive Intuition „nur hohe Korrelation ist gefährlich" widerlegen.

**Relevanz/Zielgruppe:** Portfoliokonstrukteure, BL-Software-Autoren (jede API, die Ω als Nutzer-Input entgegennimmt, während Σ intern geschätzt wird, hat genau diese Bruchstelle), Risk-Governance. Die Zielgruppe ist real; die Motivationslücke (Verbreitung des Workflows) bleibt der wunde Punkt.

**Verbesserungen mit Hebel:** (a) ω-Primitiv-Diskussion (2.2) macht aus einer Rezeptur ein Argument; (b) Multi-View-Closed-Form via Sherman–Morrison (statt sie der Zukunft zu überlassen) würde den Beitrag auf Full-Paper-Niveau heben; (c) Software-Beleg für den Frozen-Workflow schließt die Motivationslücke.

---

## 11 Journal-Tauglichkeit (1–10)

| Kategorie | Score | Begründung (eine Zeile) |
|---|---|---|
| Einleitung | 7 | Klar motiviert, gute Kanal-Vorschau; Overclaims (E6, „coexist in practice") und fehlender Related-Work-Absatz |
| Literaturüberblick | 5 | Vorhandenes korrekt und gut verwoben; 8–11 Standardreferenzen fehlen (Jorion, Meucci, Kolm–Ritter, Bessler et al.) |
| Methodik | 7 | Saubere Kanaltrennung, exakte Closed Forms, vorbildliche Konventions-Pinnung; ω-Primitiv-Lücke, kein formales Modell |
| Mathematik | 8 | Alle Herleitungen korrekt nachgerechnet; ein falscher Linearitäts-Claim (E1), Interpretationsfehler (E2), Präsentationsrigor fehlt |
| Ergebnisse | 6 | Vollständig reproduzierbar (selten!); aber ein Split, ein Universum, keine Unsicherheitsquantifizierung |
| Diskussion | 6 | Limitationen benannt und ehrlich; Long-only unbelegt, Live-Rebalancing und ω-Primitiv fehlen |
| Schlussfolgerung | 7 | Präzise Zusammenfassung, gute Extensions-Liste; erbt den Linearitätsfehler |

---

## 12 Referee Report (English, as submitted to the editor)

### Summary

The paper studies what happens inside the Black–Litterman (BL) model when the covariance estimator is switched from the sample estimator to Ledoit–Wolf (LW) linear shrinkage. The authors distinguish two workflows for the view-uncertainty matrix Ω. When Ω is recalibrated from Σ via Idzorek's confidence mapping, the switch leaves the view absorption of a single view exactly unchanged. When Ω is calibrated once on the sample covariance and frozen, the effective confidence shifts by an exact, τ-independent closed form; its sign is governed by the view variance relative to μ_s‖p‖², and in the symmetric case a closed-form critical correlation ρ* results. An eight-asset ETF application quantifies the shifts (up to +10.7 pp across pairs, +12.9 pp for an absolute view on the lowest-variance asset) and a 2019–2024 backtest attributes the performance differences to the precision channel, showing that the frozen-Ω workflow forfeits roughly 40 % of the de-risking the estimator switch was meant to deliver.

The mechanism is, to my knowledge, new in this isolated form, the algebra is correct (I verified every derivation and reproduced every reported number from the authors' code and frozen data, which deserves explicit credit), and the paper is honest about its limitations. The contribution is nevertheless narrow, the presentation is informal for a research journal, and the empirical section rests on a single split with ex-post-selected views.

### Major Comments

**M1. A repeated mathematical claim is false.** The abstract, Section 5, and the conclusion state that under Idzorek calibration the remaining posterior impact runs through prior and view projection, "both exactly linear in α". The prior is exactly linear; the normalised projection Σ(α)p/v(α) is a ratio of affine functions of α and hence not linear (it is linear only to first order, or exactly when p is an eigenvector of Σ̂_s). Section 7's "Δπ = O(α) is the only remaining term" contradicts the authors' own diagnostics in Section 8 (max posterior shift 0.42 % of which the prior contributes at most 0.21 %). Please correct all four locations.

**M2. The economic interpretation of the prior shift is contradicted by the authors' own data.** Section 5 claims below-average-variance assets receive a higher prior under LW. IWF and IWD have below-average variance yet negative Δπ (−0.113 %, −0.114 %), because the cross-covariance term dominates. Restate via Δπ_i = αδ(μ_s w_i − Cov(r_i, r_p)).

**M3. The normative message presumes confidence-primacy without saying so.** If an analyst's primitive is the view-noise variance ω itself, freezing ω under an estimator upgrade is Bayesianly coherent and the induced change in c_eff is the correct posterior response; "recompute Ω" would then be the error. The paper's recommendation is valid only when the stated confidence c is the primitive (as under Idzorek's procedure and typical risk-governance sign-off). This needs an explicit paragraph; without it the central recommendation is open to a one-line objection.

**M4. The prevalence of the frozen-Ω workflow is asserted, then admitted to be undocumented.** "Two workflows coexist in practice" (Introduction) vs. "We are not aware of survey evidence" (Section 6). Either provide evidence — e.g., open-source BL implementations whose interfaces accept a user-supplied Ω while Σ is estimated internally — or make the motivation explicitly hypothetical.

**M5. Unsupported empirical claim.** "Imposing a long-only constraint instead compresses both variants toward each other and toward far tamer risk figures" is stated without any reported numbers, although the accompanying code computes long-only twins. Report the table or delete the claim.

**M6. The backtest is a single, hindsight-exposed design.** One train/test split, one universe, one view set chosen by authors who know the test window; the dominant view (IWF>IWD) happens to realise +9.1 % p.a. The frozen-vs-recomputed comparison is internally valid, but the levels are not interpretable and no robustness is offered. At minimum: (i) a rolling-origin protocol or several splits; (ii) sensitivity of the key comparison in α, c, and the view set (e.g., a grid of pair views); (iii) HAC-robust inference (Ledoit–Wolf 2008, which the paper cites but does not use; the test window contains March 2020) also for the frozen-vs-recomputed pair, not only sample-vs-LW.

**M7. No probabilistic model is stated.** Equation (1) is used without the underlying Bayesian structure (normal prior on μ, view noise independent, Σ treated as known). Since the paper's subject is precisely the estimation of Σ inside a model that conditions on Σ, the plug-in tension must be acknowledged once, with references (e.g., Kolm & Ritter 2017; Kan & Zhou 2007).

**M8. Presentation is below journal form.** Main results should be stated as propositions with hypotheses (Σ ≻ 0, ω > 0, K = 1) and short proofs; the paper defines theorem environments and never uses them. There is not a single figure; the sign criterion and the pair table beg for one. The Idzorek closed-form equivalence ("its exact single-view solution") is asserted, not proven — the proof is three lines and belongs in an appendix, together with the analytic α̂ formula (currently outsourced to scikit-learn).

### Minor Comments

1. "The answer depends entirely on a bookkeeping question" overstates: the precision channel (Section 7) is Ω-independent. Restrict to the view channel.
2. Abstract states the multi-view deviation "stays an order of magnitude below the effects studied here" as a general fact; it is a finding on one three-view example. Qualify.
3. Report α̂ = 0.0406 once; "0.041" everywhere else is fine.
4. Table 1: the column α = 0.166 is unexplained; use a round value or explain its origin.
5. δ is implied by S&P 500 excess returns while w_eq comes from Idzorek's 2004 asset-class weights and the universe is ETFs; three different market proxies enter one equilibrium argument. State this as a deliberate approximation.
6. Notation: k doubles as view index and volatility ratio (Table 2); "Std"/"s" superscripts vs. subscripts inconsistent; "percentage points" vs. "points".
7. Maximum drawdown on the compounded *excess*-return path is non-standard; add total-return drawdowns or a footnote.
8. "forfeits close to half" (conclusion, Section 8) vs. the computed ≈ 40 %; keep "roughly 40 %".
9. Data construction: state how the 13-week T-bill yield was de-annualised; print the w_eq vector.
10. The deflation floor is 1.05 pp, stated as "1.1"; harmonise rounding.
11. Section 6.6, Step 1: the differencing is taken "at fixed Σ" but applied across the estimator switch; only the projection pᵀΔμ is exact. Restructure (exact along-p identity first, vector approximation with explicit error second).
12. "Both predictions match … to within a few hundredths of a percentage point": true for Δμ (0.01 pp), false for Δw (0.2 pp).
13. Tone down essayistic phrasing ("easy to fall into", "earns its keep", "rides exactly this lever", "prices the leakage", "Where does LW earn its keep?"); one rhetorical question and several pull-quote sentences exceed the register of this journal.
14. References: add Jorion (1986), Frost & Savarino (1986), Kan & Zhou (2007), Meucci (2010), Kolm & Ritter (2017), Michaud et al. (2013), Herold (2003), Bevan & Winkelmann (1998), Bessler et al. (2017), Ledoit & Wolf (2020 or 2022), DeMiguel et al. (2013).

### Technical Corrections

T1. Abstract/S5/S7/Conclusion: replace "exactly linear in α" as per M1.
T2. Section 5: replace the Δπ interpretation sentence as per M2.
T3. Section 6.7: correct the accuracy claim (0.01 pp for Δμ, 0.2 pp for Δw).
T4. Remark 2: floor is 1.05 pp (1.04 pp at α̂ = 0.0406), not 1.1.
T5. Section 6: "Σ moves in the numerator of c_eff,k … while ω_k does not move in the denominator" — v_k moves in numerator *and* denominator; rewrite ("v_k moves while ω_k stays fixed").
T6. Abstract: "+12.9 points" → "percentage points".
T7. Section 8: "propagates into the test window" → the anomaly propagates into the *weights* held during the test window.
T8. "close to half" → "roughly 40 %" (twice).
T9. Footnote 1 or Section 2: note that the Σ-vs-Σ+M choice rescales c_eff by ≈ 1/(1+τ), so results are mapping-dependent in the third decimal.
T10. Fix \date{\today}; resolve the (⋆) tag if referenced.

### Recommendation

**Major Revision.**

The paper contains a correct, novel, and practically relevant observation, and its reproducibility discipline is exemplary. It is not publishable in its current form: a repeated analytical claim is false (M1), an interpretation is contradicted by the paper's own data (M2), the central recommendation lacks its conceptual premise (M3), the motivating workflow is undocumented (M4), one empirical claim is unsupported (M5), and the empirical design carries no robustness (M6). All of these are fixable without new theory; M1–M5 are days of work, M6 is a few weeks. Scope-wise the paper fits *Journal of Asset Management*, the practitioner side of *Quantitative Finance*, or *Journal of Portfolio Management*; for a top field journal the multi-view closed form (currently future work) would be needed. I would be willing to review a revised version.

---

## 13 Satz-für-Satz-Prüfung

Format je Eintrag: **Zitat** → Problem → Warum problematisch → Bessere Version. Sätze, die hier nicht auftauchen, haben die mathematische und die stilistische Prüfung bestanden (das gilt insbesondere für die gesamte Herleitungskette in 6.6, die ich Zeile für Zeile verifiziert habe).

### Abstract

**A1.** *„The Black–Litterman model is routinely combined with shrinkage estimators of the covariance matrix…"*
Problem: „routinely" ist eine empirische Behauptung ohne Beleg. Warum: Der erste Satz eines Abstracts trägt Beweislast. Besser: „The Black–Litterman model is often combined with shrinkage estimators of the covariance matrix, most prominently the Ledoit–Wolf estimator."

**A2.** *„with several views the deviation is driven by the cross-view covariances and stays an order of magnitude below the effects studied here."*
Problem: Ein Befund aus einem Drei-View-Beispiel wird als allgemeine Tatsache formuliert. Warum: Für hinreichend korrelierte Views ist die Abweichung nicht beschränkt; der Abstract verspricht mehr, als Remark 1 liefert. Besser: „…and, in our three-view application, stays an order of magnitude below the effects studied here."

**A3.** *„the remaining effect on the posterior runs through the equilibrium prior and the view projection, both exactly linear in α"*
Problem: falsch (siehe E1). Warum: Σ(α)p/v(α) ist rational, nicht linear in α; numerisch belegt. Besser: „the remaining effect runs through the equilibrium prior, exactly linear in the shrinkage intensity α, and the view projection, linear in α to first order."

**A4.** *„an absolute view on the lowest-variance asset shifts by +12.9 points"*
Problem: „points" statt „percentage points". Warum: Einheiten-Konsistenz; im selben Abstract steht „percentage points". Besser: „…by +12.9 percentage points".

### Einleitung (Abschnitt 1)

**I1.** *„is one of the few pieces of portfolio theory that made the crossing from journals to practice."*
Problem: feuilletonistische Metapher plus unbelegtes „one of the few". Warum: Erster Satz des Papers; ein Referee liest ihn als Stilprobe. Besser: „is among the most widely implemented models of quantitative portfolio construction."

**I2.** *„The standard remedy is the shrinkage estimator of Ledoit and Wolf (2004)…"*
Problem: bestimmter Artikel. Warum: Faktor-Modelle, EWMA, Sparse- und Nonlinear-Shrinkage sind gleichrangige Alternativen; „the standard" ist nicht haltbar. Besser: „A standard remedy is the shrinkage estimator of Ledoit and Wolf (2004)…"

**I3.** *„The answer depends entirely on a bookkeeping question that, to our knowledge, has not been isolated in the literature: what happens to Ω?"*
Problem: „entirely" widerspricht Abschnitt 7 (Precision-Channel ist Ω-unabhängig und laut Paper der größere Effekt). Warum: interner Widerspruch im selben Dokument. Besser: „Within the view channel, the answer turns on a bookkeeping question that, to our knowledge, has not been isolated in the literature: the treatment of Ω."

**I4.** *„Two workflows coexist in practice."*
Problem: Behauptung, die Abschnitt 6 selbst als unbelegt ausweist. Warum: Die Motivation des Papers hängt daran. Besser: „Two workflows arise naturally." plus Beleg (Software-Schnittstellen, die Ω als Nutzerkonstante akzeptieren) oder explizit konditional.

**I5.** *„We call the resulting silent shift of the effective inputs \emph{leakage}."*
Problem: Der zentrale Begriff wird informell und vor der zugehörigen Größe definiert. Warum: „leakage" ist der Titelbegriff; er braucht ein formales Objekt. Besser: „We call the resulting shift of the effective view confidence, c_eff − c ≠ 0 at frozen ω, \emph{leakage}."

**I6.** *„…the remaining effect on the posterior runs through the equilibrium prior and the view projection, both exactly linear in α."*
Problem/Warum: wie A3 (dritte Nennung). Besser: wie A3.

**I7.** *„with the analytic shrinkage intensity α = 0.041"*
Problem: Der Schätzwert ist 0.0406. Warum: Reproduzierbarkeit auf die genannte Ziffer; die Rundung erzeugt Folge-Inkonsistenzen (Floor 1.04 vs. 1.05). Besser: einmal „α̂ = 0.0406 (0.041 im Folgenden)".

**I8.** *„which shows that the intuition a reader might bring, that only high correlations are affected, is incomplete."*
Problem: verschachtelt, drei Kommata, „a reader might bring" ist Füllmaterial. Warum: Lesbarkeit an einer Pointe-Stelle. Besser: „which shows that high correlation is not necessary for a material shift."

**I9.** *„The practical lesson is procedural rather than statistical."*
Problem: Kontrastformel als Absatzauftakt; die Substanz folgt erst im nächsten Satz. Warum: rhetorisches Muster, das im Text noch fünfmal vorkommt. Besser: Satz streichen und direkt beginnen: „Covariance estimator, confidence and view uncertainty form one unit; …"

### Abschnitt 2

**S2-1.** *„On the training window the implied risk aversion of the S&P 500 excess return series is δ = 2.44 (the ratio of annualised mean to annualised variance)."*
Problem: (a) Der Annualisierungsfaktor kürzt sich, die Klammer suggeriert eine Skalenabhängigkeit, die nicht existiert; (b) δ stammt vom S&P 500, w_eq von Idzoreks Assetklassen, das Universum sind ETFs, also drei Markt-Proxies in einem Gleichgewichtsargument, unkommentiert. Warum: Konsistenz des Equilibrium-Arguments ist ein Kernbaustein von π. Besser: „…δ = 2.44, the ratio of mean to variance of monthly S&P 500 excess returns (the annualisation factor cancels). We use the S&P 500 as the proxy for δ while the equilibrium weights come from Idzorek's asset classes; we accept this mismatch as the customary approximation."

### Abschnitt 3

**S3-1.** *„For n assets and T observations this is reliable only if T ≫ n; the relative estimation error grows with n/T."*
Problem: unquantifizierte Folklore. Warum: Ein Paper über Kovarianzschätzung muss die eigene Prämisse präzise fassen. Besser: „Its expected loss grows with n/T (Ledoit and Wolf, 2004); once n/T is non-negligible, the eigenvalues of Σ̂_s spread systematically…" (schließt direkt an den vorhandenen, korrekten Eigenwert-Satz an).

### Abschnitt 4

**S4-1.** *„where the oracle shrinkage coefficient α ∈ [0,1] minimises the expected squared Frobenius error E[‖Σ̂ − Σ‖_F²]"*
Problem: Das Minimierungsproblem ist unterspezifiziert (worüber wird minimiert, welches Target, Populations- vs. Stichprobengrößen). Warum: Genau diese Feinheiten (T vs. T−1, Populationstarget) behandelt das Paper an anderer Stelle vorbildlich; hier fällt es dahinter zurück. Besser: „where the oracle coefficient α minimises E[‖(1−α)Σ̂_s + αμI_n − Σ‖_F²] over α ∈ [0,1], with μ the population counterpart of μ_s; the analytic estimator replaces the four unobservable moments by consistent plug-ins (formula in Appendix B)."

### Abschnitt 5

**S5-1.** *„We can keep this case brief."*
Problem: Meta-Kommentar ohne Information. Warum: Der Abschnitt ist ohnehin kurz; der Satz kostet die Zeile, die er zu sparen behauptet. Besser: ersatzlos streichen.

**S5-2.** *„(Idzorek's original procedure matches weight tilts iteratively; the closed form is its exact single-view solution)."*
Problem: tragende Äquivalenzbehauptung ohne Beweis oder präzise Referenz. Warum: Die Inertness-Aussage des ganzen Abschnitts hängt daran; die eigene Codebasis nennt die Kalibrierung ausdrücklich „vereinfacht". Besser: „(Idzorek's original procedure matches weight tilts iteratively; Appendix A verifies that the closed form solves it exactly for a single view under the weight mapping of Section 2.)"

**S5-3.** *„Both effects are exactly linear in α."*
Problem: falsch (E1, vierte Nennung im Dokument mitgezählt). Besser: „The prior shift is exactly linear in α; the projection Σp/v_k is a smooth rational function of α, linear to first order."

**S5-4.** *„Assets with above-average variance (σ_i² := Σ̂_{s,ii} > μ_s) receive a lower prior, assets with below-average variance receive a higher one."*
Problem: falsch auf den eigenen Daten (E2: IWF, IWD). Warum: Die Aussage unterschlägt den dominanten Kreuzterm; ein Leser, der sie auf Tabelle im Code-Output anwendet, scheitert an 2 von 8 Assets. Besser: „The diagonal term raises the prior of below-average-variance assets, but the removal of cross-covariances dominates: Δπ_i = αδ(μ_s w_{eq,i} − Cov(r_i, r_p)) with r_p = w_eqᵀr, so assets that co-move strongly with the equilibrium portfolio lose prior return. On our data all six equity ETFs have Δπ_i < 0, including IWF and IWD despite below-average variances; the bond ETFs gain."

**S5-5.** *„Idzorek's Ω matches only diag(S), so off-diagonal elements of S break both the exact absorption c and its estimator-invariance."* (Remark 1)
Problem: „matches diag(S)" gilt wörtlich nur für c = 0.5. Warum: Präzision in einer formalen Remark. Besser: „Idzorek's Ω matches diag(S) up to the factor (1−c_k)/c_k, so off-diagonal elements of S break…"

### Abschnitt 6

**S6-1.** *„A different workflow is easy to fall into."*
Problem: umgangssprachlich, personalisiert. Warum: Registerbruch am Beginn des Hauptabschnitts. Besser: „A second workflow stores ω_k once and leaves it untouched."

**S6-2.** *„…the point of this section is that whoever follows it changes the model's inputs without seeing it."*
Problem: „without seeing it" kolloquial; „whoever" locker. Besser: „…the point of this section is that this workflow changes the model's effective inputs without any visible change at the input stage."

**S6-3.** *„The view absorption shifts implicitly, because Σ moves in the numerator of c_eff,k (defined below) while ω_k does not move in the denominator."*
Problem: sachlich schief, v_k steht in Zähler *und* Nenner. Warum: Der Satz beschreibt den Kernmechanismus; wer die Formel noch nicht kennt, lernt sie hier falsch. Besser: „The view absorption shifts implicitly: the view variance v_k moves with the estimator while ω_k stays fixed."

**S6-4.** *„As v_k grows, so does the absorption c_eff,k."*
Problem: gilt nur bei festem ω_k, was hier der Punkt ist, aber ungesagt bleibt. Besser: „At fixed ω_k, the absorption c_eff,k increases with v_k."

**S6-5.** *„Two specialisations are worth recording." / „Three consequences are worth recording."*
Problem: identische Floskel in zwei aufeinanderfolgenden Unterabschnitten. Besser: eine der beiden ersetzen, z. B. „Two special cases follow immediately."

**S6-6.** *„…no view in any universe can lose more than 1.1 percentage points of effective confidence"* (Remark 2)
Problem: Floor ist 1.05 pp (bzw. 1.04 bei α̂ = 0.0406); im selben Satz wird auf eine Dezimale gerundet. Besser: „…more than 1.05 percentage points".

**S6-7.** *„The symmetric case is a structural sketch rather than the general statement:"*
Problem: Kontrastformel; „structural sketch" unscharf. Besser: „The symmetric case serves as a rule of thumb; for σ_i ≠ σ_j we evaluate the exact criterion (3), and the pair-level results below do so."

**S6-8.** *„The parameters … are illustrative, chosen of the same order as the monthly quantities in our data"*
Problem: Grammatik („chosen of the same order"). Besser: „chosen to match the order of magnitude of the monthly quantities in our data".

**S6-9.** *„a correlation level routine for style pairs within one equity market."*
Problem: Behauptung ohne Anker. Besser: „a correlation level typical of style pairs within one equity market (ρ = 0.92 for IWF/IWD, Table 2)."

**S6-10.** *„The intensity α itself grows with n/T, so larger universes on the same window sit further to the right in the table."*
Problem: Heuristik über den analytischen Schätzer als Faktum. Besser: „The analytic intensity tends to grow with n/T (Ledoit and Wolf, 2004), so larger universes on the same window typically sit further to the right in the table."

**S6-11.** *„…a frozen Ω lets the shrinkage-induced widening of the spread variance flow one-for-one into extra view confidence."*
Problem: „one-for-one" ist quantitativ falsch; die Übersetzung läuft über Gl. (3) und ist nichtlinear. Besser: „…lets the shrinkage-induced widening of the spread variance flow directly into extra view confidence, at the rate given by (3)."

**S6-12.** *„The model believes the view more strongly than the analyst stated."*
Problem: Anthropomorphismus, Pull-Quote. Besser: „The effective confidence exceeds the stated confidence."

**S6-13.** *„Differencing \eqref{eq:bl-single-exact} at fixed π, Σ and ω:"* (Step 1)
Problem: Die Ableitung fixiert Σ, die Anwendung wechselt Σ (E5); exakt bleibt nur die Projektion auf p. Warum: Der Unterschied ist im 2×2-Beispiel messbar (Asymmetrie +30.2/−29.8). Besser: „Projecting the exact identity on p gives pᵀΔμ_BL = Δc_eff (Q − pᵀπ) for any pair of estimators, since pᵀ(Σp/v_k) = 1 for each. The full vector shift additionally moves the direction Σp/v_k; at fixed Σ it reduces to Δμ_BL = Δc_eff (Q − pᵀπ) Σp/v_k, and for σ_i ≈ σ_j …" (dann wie gehabt).

**S6-14.** *„\eqref{eq:dw-dc} is the exact partial derivative at fixed Σ, i.e. the pure view-absorption channel. This is the quantity measured numerically in the 2×2 example below…"*
Problem: Das Beispiel misst zusätzlich den Richtungsshift; „This is the quantity measured" stimmt nur approximativ. Besser: „…the pure view-absorption channel. The 2×2 example below measures it up to the direction shift Σ_LW p/v^LW − Σ_s p/v^Std, which the shared Σ̂_s⁻¹ does not remove; the residual is visible as the asymmetry of the measured (+30.2, −29.8)."

**S6-15.** *„Real Example: Clean 2×2 System with IWF and IWD"* (Überschrift 6.7)
Problem: „Real" und „Clean" sind Werbevokabeln. Besser: „A Worked 2×2 Example: IWF and IWD".

**S6-16.** *„…the pair sits squarely in the homogeneous special case, so the symmetric approximations of Section 6.6 should be accurate here."*
Problem: Adverb „squarely"; „should be" vage vor einer Stelle, die es dann exakt prüft. Besser: „…the pair lies in the homogeneous special case (k = 0.99), so the symmetric approximations of Section 6.6 apply."

**S6-17.** *„Both predictions match the numerical Δμ and Δw to within a few hundredths of a percentage point."*
Problem: falsch für Δw (0.2 pp, E3). Besser: „The Δμ prediction matches to within 0.01 pp, the Δw prediction to within 0.2 pp; both gaps stem from the volatility asymmetry that the symmetric split ignores."

**S6-18.** *„The chain runs in one direction."*
Problem: dramatisches Fragment ohne eigenständigen Inhalt. Besser: mit dem Folgesatz verschmelzen: „The chain runs from the correlation drop to the weights: the correlation falls from 0.920 to 0.882, v_k grows…"

**S6-19.** *„The absolute position sizes are themselves a comment on unconstrained Black–Litterman with a ρ = 0.92 pair: … and the leakage rides exactly this lever."*
Problem: „are a comment on", „rides this lever" kolloquial. Besser: „The absolute position sizes illustrate a known property of unconstrained Black–Litterman with highly correlated pairs: the translation factor 1/(δv_k) ≈ 105 turns any view into triple-digit tilts, and the leakage is amplified by the same factor."

**S6-20.** *„The point is not academic: the analyst who freezes Ω has, without knowing it, raised the confidence of every affected view, inconsistent with whatever risk-governance process produced the original confidences."* (6.9)
Problem: „The point is not academic" ist Kontrastfloskel; „whatever" abschätzig. Besser: „An analyst who freezes Ω has raised the confidence of every affected view without a decision to do so, in conflict with the risk-governance process that produced the original confidences."

### Abschnitt 7

**S7-1.** *„The return channel is small because Idzorek fixes view absorption at c and Δπ = O(α) is the only remaining term."*
Problem: „the only remaining term" unterschlägt den Projektionsterm, der laut eigener Diagnostik (Abschnitt 8) gleich groß ist. Besser: „The return channel is small because Idzorek fixes the view absorption at c; the remaining terms, the prior shift and the projection shift, are both O(α) and of comparable size on our data (Section 8)."

### Abschnitt 8

**S8-1.** *„The view set is in the spirit of Idzorek's examples—relative pair views with moderate spreads—and is deliberately not tuned to the test window."*
Problem: „deliberately not tuned" ist als Absichtserklärung nicht prüfbar; die Autoren kennen das Testfenster. Warum: Hindsight-Verdacht ist die naheliegendste Referee-Frage an diesen Backtest. Besser: „The view set follows Idzorek's examples (relative pair views with moderate spreads). Since we know the test window, the frozen-vs-recomputed comparison is repeated over a grid of pair views in Section X; the conclusions are unchanged." (setzt die Robustheitsergänzung voraus).

**S8-2.** *„Where does LW earn its keep? Three configurations plausibly favour it; the list is informed conjecture rather than a tested claim."*
Problem: rhetorische Frage plus Idiom plus Kontrastformel in zwei Sätzen. Besser: „Three configurations plausibly favour LW; we state them as conjectures."

**S8-3.** *„The analyst who freezes Ω forfeits close to half of the de-risking they think they bought."*
Problem: „close to half" vs. gemessene ≈ 40 %; „they think they bought" kolloquial. Besser: „Freezing Ω forfeits roughly 40 % of the exposure and drawdown reduction that the estimator switch delivers."

**S8-4.** *„This is the dimension where LW carries its own weight; a drawdown above 50 % is not investable for many institutional mandates, regardless of the realised return."*
Problem: Idiom; die Investierbarkeits-Aussage ist unbelegt und mit „regardless" verschärft. Besser: „The drawdown reduction is where LW adds value in this window; drawdowns above 50 % breach common institutional loss limits." (Referenz ergänzen oder als Einschätzung kennzeichnen.)

**S8-5.** *„Imposing a long-only constraint instead compresses both variants toward each other and toward far tamer risk figures, consistent with the observation of Jagannathan and Ma (2003, p. 1652)…"*
Problem: empirische Aussage ohne berichtete Zahlen (M5), obwohl der Code sie berechnet. Besser: Tabelle ergänzen und schreiben: „Table X reports the long-only counterparts: both variants compress to …% volatility and …% drawdown, consistent with Jagannathan and Ma (2003)."

**S8-6.** *„…(such as the unusually wide eigenvalue spreads of 2008–2009) that flows into the sample covariance and propagates into the test window."*
Problem: Die Anomalie propagiert in die Gewichte, nicht ins Testfenster. Besser: „…that enters the sample covariance and distorts the weights held during the test window."

### Conclusion

**C1.** *„Whether the Ledoit–Wolf switch is harmless inside Black–Litterman is decided by the treatment of Ω, and both regimes now have closed forms."*
Problem: erbt E6; der Precision-Channel entscheidet laut Abschnitt 7 mit. Besser: „Within the view channel, whether the Ledoit–Wolf switch is harmless is decided by the treatment of Ω; both regimes now have closed forms."

**C2.** *„…the posterior impact is confined to prior and projection shifts that are exactly linear in α."*
Problem: E1, letzte Nennung. Besser: „…confined to prior and projection shifts, the former exactly linear and the latter first-order linear in α."

**C3.** *„…the frozen-Ω workflow forfeits close to half of the exposure and drawdown reduction that the estimator switch was meant to deliver."*
Problem: wie S8-3. Besser: „…forfeits roughly 40 % of the exposure and drawdown reduction."

Alle übrigen Sätze des Manuskripts, einschließlich der vollständigen Herleitungen in 6.6, der Remark-2-Argumentation und des Reproducibility-Statements, sind sachlich korrekt und sprachlich journaltauglich.

---

## 14 Priorisierte To-do-Liste

### Kritisch vor Einreichung

1. **„Exactly linear in α" an allen vier Stellen korrigieren** (Abstract, Abschnitt 5, Abschnitt 7 „only remaining term", Conclusion). Ein nachweislich falscher mathematischer Satz im Abstract ist ein eigenständiger Ablehnungsgrund. (E1/M1)
2. **Δπ-Interpretationssatz ersetzen** durch die Kovarianz-zum-Portfolio-Lesart; die aktuelle Version wird von IWF/IWD widerlegt. (E2/M2)
3. **ω-Primitiv vs. c-Primitiv diskutieren** (ein Absatz in Abschnitt 6 oder 6.9): begründen, warum c das Primitiv ist, und einräumen, dass unter ω-Primat das Einfrieren kohärent wäre. Ohne diesen Absatz kippt die normative Kernaussage mit einem Einzeiler. (M3)
4. **Long-only: Zahlen berichten oder Satz streichen.** Der Code berechnet die Zwillinge bereits. (M5)
5. **Motivationslücke schließen:** Beleg für den Frozen-Ω-Workflow (Software-Schnittstellen mit nutzerdefiniertem Ω, Praxisreferenz Bevan–Winkelmann) oder konditionale Umformulierung von „Two workflows coexist in practice"; zugleich „depends entirely" (E6) und „is decided by" (C1) auf den View-Channel einschränken. (M4)
6. **Formales Modell + Propositionen:** Bayes-Struktur vor Gl. (1), Hauptresultate als Proposition 1–3 mit Voraussetzungen, Idzorek-Äquivalenzbeweis und α̂-Formel in den Appendix. (M7/M8)
7. **6.6/6.7 präzisieren:** exakte Along-p-Identität von der Vektorapproximation trennen; Genauigkeitsangabe korrigieren (0.01 pp / 0.2 pp); Floor 1.05 pp. (E3/E4/E5)

### Wichtig

8. **Robustheit des Backtests:** Rolling-Origin oder mehrere Splits; View-Gitter über alle 28 Paare; Sensitivität in α (0.02–0.30) und c; Frozen-vs-Recomputed zusätzlich testen; Inferenz auf Ledoit–Wolf (2008) HAC umstellen (COVID im Testfenster). (M6)
9. **Bootstrap-Konfidenzbänder für Δc_eff** (Block-Bootstrap über Monate), als Spalte in Tabelle 2; beantwortet, ob die Paar-Rangfolge statistisch belastbar ist.
10. **Zwei Abbildungen:** (i) Δc_eff gegen v^Std/(μ_s‖p‖²) mit den 28 Paaren plus den absoluten Views als Punkten (macht Sign-Kriterium, Asymmetrie und die AGG/BWX-Pointe sichtbar); (ii) kumulierte Out-of-sample-Pfade mit Drawdown-Schattierung.
11. **Literatur ergänzen** (11 Referenzen aus Abschnitt 6 dieses Reviews) und einen Related-Work-Absatz mit den drei Strängen (BL-Spezifikation, Shrinkage im Portfolio, Plug-in-Unsicherheit) einziehen.
12. **Appendix-Paket:** Idzorek-Äquivalenz (5 Zeilen), α̂-Formel + sklearn-Version, w_eq-Tabelle, T-Bill-Umrechnungskonvention.
13. **Notation bereinigen:** k-Doppelbelegung (Tabelle 2), Std/s-Vereinheitlichung, „points"→„pp", (⋆)-Tag, σ_i²-Definition nach Abschnitt 2, optional Symbolverzeichnis.
14. **Sprachpass nach Abschnitt 7 dieses Reviews:** Kontrastformeln reduzieren, rhetorische Frage raus, Pull-Quotes entschärfen, Em-Dash-Dichte halbieren, Abschnitte 3+4 zusammenlegen.
15. **Abstract straffen** (< 200 Wörter): Mehr-View-Satz qualifizieren (A2), +12.9-Detail optional streichen.

### Optional

16. **Multi-View-Closed-Form via Sherman–Morrison** ausarbeiten statt als Future Work zu deklarieren; hebt das Paper vom Note- auf Full-Article-Niveau und wäre die wirksamste Einzelinvestition.
17. **Live-Rebalancing-Absatz:** Leckage tritt bei rollierender Neuschätzung an jedem Termin auf; erweitert die Relevanz über den einmaligen Estimatorwechsel hinaus.
18. **1/N- und Minimum-Varianz-Benchmarks** in die Backtest-Tabelle (DeMiguel et al. werden ohnehin zitiert).
19. **Zweite Datenquelle** zur Validierung der Yahoo-Reihen (oder NAV-Returns der Anbieter).
20. **Long-only-Zwilling auch für die Frozen-Ω-Variante** (der Code rechnet bisher nur sample/LW-recomputed).
21. **Tabelle 1:** α = 0.166 erklären oder durch runden Wert ersetzen; optional eine Zeile mit datenkalibriertem μ_s/σ².
22. **Total-Return-Drawdowns** zusätzlich zu den Excess-Pfad-Drawdowns ausweisen.

---

## 15 Abschlussbewertung

| Kategorie | Score (1–10) | Kurzbegründung |
|---|---|---|
| Wissenschaftliche Qualität | 6.5 | Echter, sauber isolierter Mechanismus; konzeptionelle Lücke (ω-Primat) und dünne Evidenzbasis |
| Mathematische Korrektheit | 8 | Alle Herleitungen und alle 40+ Zahlen unabhängig verifiziert; ein falscher Linearitäts-Claim, zwei Präzisionsmängel |
| Quantitative Finance | 7 | BL- und LW-Darstellung korrekt und konventionsfest; fehlendes Bayes-Fundament, Proxy-Inkonsistenz unkommentiert |
| Methodik | 7 | Kanaltrennung und Konventions-Pinnung vorbildlich; kein formales Modell, Approximationslogik in 6.6 unscharf |
| Empirie | 5 | Vollständig reproduzierbar (herausragend), aber ein Split, ein Universum, Ex-post-Views, keine Unsicherheitsquantifizierung |
| Literatur | 5 | Vorhandenes korrekt, Einbindung gut; rund ein Dutzend Standardreferenzen fehlt |
| Originalität | 6 | Mechanismus in dieser Isolation m. W. neu; Mathematik elementar, Beitrag schmal, aber praktisch relevant |
| Sprache | 7.5 | Muttersprachliches Niveau, präzise; Register stellenweise essayistisch, formelhafte Muster |
| Struktur | 6 | Logischer Aufbau, klare Scope-Abgrenzung; keine Abbildungen, ungenutzte Theorem-Umgebungen, Mini-Abschnitte |
| Journal Readiness | 4.5 | Vor Behebung von E1/E2/M5 nicht einreichbar; danach solide Major-Revision-Kandidatur |

**Gesamtbewertung: 6.3 / 10** (ungewichtetes Mittel; gewichtet man mathematische Korrektheit und Originalität höher, ändert sich wenig).

**Annahmewahrscheinlichkeiten (Einschätzung):**

- **In aktueller Form bei einem guten Quant-Finance-Journal** (Quantitative Finance, Journal of Empirical Finance, JBF): **≈ 5 %.** Realistisch ist Desk-Reject (30–40 % Risiko wegen Scope und fehlender Abbildungen/Propositionen) oder Major Revision; eine Annahme ohne Revision scheidet wegen E1/E2/M5 aus.
- **Nach Umsetzung aller kritischen und wichtigen Punkte:** **≈ 35–45 %** bei Quantitative Finance; **≈ 60–70 %** bei Journal of Asset Management, Journal of Portfolio Management oder Financial Analysts Journal, deren Leserschaft die prozedurale Pointe direkt adressiert. Mit der Multi-View-Erweiterung (Punkt 16) steigt die QF-Schätzung auf **≈ 50–60 %**.

**Schlussurteil:** Das Manuskript ist deutlich besser als der Durchschnitt studentischer Einreichungen: Die Kernrechnung stimmt, jede Zahl reproduziert per Skript, Limitationen werden benannt statt versteckt, und die Bond-Paar- und Absolute-View-Befunde tragen eine echte, nicht offensichtliche Pointe. Was fehlt, ist die Härtung zum Journal-Artikel: ein korrektes Linearitäts-Statement, eine daten-konsistente Δπ-Interpretation, das konzeptionelle Fundament der Empfehlung (c-Primat), belegte statt behaupteter Motivation, berichtete statt versprochener Long-only-Zahlen, Robustheit und Form (Propositionen, Abbildungen, Literatur). Nichts davon erfordert neue Theorie. Empfehlung an den Autor: erst Punkte 1–7, dann Venue-Entscheidung; für QF zusätzlich Punkt 16 ernsthaft erwägen.

---

*Review erstellt am 14.07.2026. Verifikation: `code/paper_numbers.py` ausgeführt gegen `data/` (sklearn 1.7.2, α̂ = 0.0406); alle Paper-Zahlen reproduziert; Linearitätstest (E1) und Δπ-Gegenbeispiel (E2) numerisch bestätigt.*



