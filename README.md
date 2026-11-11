# Gain-Scheduled Youla–Kučera Anti-Windup für LPV-Systeme

[🇬🇧 English](README.md) | 🇩🇪 Deutsch

Ein kleines, in sich geschlossenes Python-Projekt, das einen **beobachterbasierten, verstärkungsgeplanten Regler in Youla–Kučera-Form (YK) bzw. in koprimer Faktorisierung mit Anti-Windup** für eine lineare parametervariante (LPV) Strecke mit Stellgrößenbegrenzung implementiert.

Das Projekt vergleicht drei Fälle an derselben Strecke und demselben Szenario:

1. die **unbeschränkte lineare** Antwort (Referenzverhalten),
2. den **gesättigten** Regelkreis **ohne** Anti-Windup (Integrator-Windup),
3. den **gesättigten** Regelkreis **mit** YK-Anti-Windup.

Die Reglerverstärkungen lassen sich auf zwei Arten entwerfen:

| Methode | Beschreibung | Stabilitätsaussage |
|---|---|---|
| `lqr` | LQR (Zustandsrückführung) + Kalman (Beobachter) auf einem θ-Gitter, linear interpoliert | nur für *eingefrorenes* θ stabil |
| `lmi` | eine gemeinsame Lyapunov-Matrix über dem Polytop (Eckpunkt-LMIs, gelöst mit CVXPY) | quadratisch stabil für *beliebig schnelles* θ(t) (B, C konstant) |

---

## Inhalt

- [Idee in Kürze](#idee-in-kürze)
- [Projektstruktur](#projektstruktur)
- [Installation](#installation)
- [Verwendung](#verwendung)
- [Beispielstrecke und Szenario](#beispielstrecke-und-szenario)
- [Erwartete Ergebnisse](#erwartete-ergebnisse)
- [Anpassung an die eigene Strecke](#anpassung-an-die-eigene-strecke)
- [Einschränkungen](#einschränkungen)
- [Fehlersuche](#fehlersuche)
- [Literatur](#literatur)
- [Lizenz](#lizenz)

---

## Idee in Kürze

**Strecke (LPV):** $\dot x = A(\theta)x + Bu,\; y = Cx$ mit $A(\theta)=A_0+\theta A_1$, $\theta\in[0,1]$.

**Führungsverhalten:** Ein Integrator $\dot\xi = r-y$ wird ergänzt; es entsteht die erweiterte Strecke $(A_a,B_a)$.

**Verstärkungen:** $F(\theta)$ mit stabilem $A_a+B_aF$ und $L(\theta)$ mit stabilem $A-LC$ (Vorzeichenkonvention aus Dings Buch).

**Koprime Faktoren** (numerisch über die Bezout-Identität an mehreren eingefrorenen θ geprüft):

$$G = NM^{-1} = \tilde M^{-1}\tilde N,\qquad XM-YN=I,\qquad K_0 = X^{-1}Y$$

**Beobachterbasierter Regler** mit Zuständen $x_k=[\hat x;\xi]$:

$$\dot x_k = A_k x_k + B_k\begin{bmatrix}y\\ r\end{bmatrix} + H\,(u_{sat}-u_c),\qquad u_c = C_k x_k,\quad u_{sat}=\mathrm{sat}(u_c)$$

**Anti-Windup-Verstärkung** $H=\begin{bmatrix}B\\ k_{aw}F_i^{-1}\end{bmatrix}$. Damit ist

$$A_k-HC_k=\begin{bmatrix}A-LC&0\\ -h_iF_x&-k_{aw}I\end{bmatrix}$$

für jedes θ Hurwitz. Der Beobachter wird mit der *tatsächlichen* (gesättigten) Stellgröße gespeist, und der Integratorzustand wird während der Sättigung mit der Rate $k_{aw}$ zurückgeführt. Ohne Sättigung ist der Regler exakt der nominale verstärkungsgeplante Regler $K_0(\theta)$.

**Warum die Planung in dieser Struktur funktioniert:** Der Schätzfehler genügt $\dot e=(A(\theta)-L(\theta)C)e$ unabhängig von $F$ (Separationsprinzip). Daher können $F$ und $L$ unabhängig voneinander geplant werden.

---

## Projektstruktur

```
lpv_yk_project/
├── plant.py        # LPV-Strecke (A0, A1, B, C), A_of(θ), Integrator-Erweiterung, Szenario θ(t) und r(t)
├── design.py       # LQR/Kalman, GainSchedule (Interpolation), design_lqr(), design_lmi()
├── controller.py   # Koprime Faktoren / Bezout-Prüfung, Reglermatrizen inkl. Anti-Windup-Verstärkung H
├── simulate.py     # RK4-Simulation des Regelkreises mit Sättigung, Kennzahlen, Plots
├── main.py         # Einstiegspunkt: Entwurf → Prüfungen → Simulation → Kennzahlen → Plot
└── plot_only.py    # zeichnet gespeicherte Logs neu, ohne erneute Simulation
```

Importkette (keine zirkulären Importe): `plant → design → controller → simulate → main`.

---

## Installation

Erfordert **Python 3.9+**.

```bash
git clone <ihre-repo-url>
cd lpv_yk_project

python3 -m venv .venv
source .venv/bin/activate          # Windows PowerShell: .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install numpy scipy matplotlib
pip install cvxpy                  # optional, nur für die Methode `lmi` nötig
```

Falls `pip install cvxpy` in ein Timeout läuft oder versucht, aus dem Quelltext zu bauen:

```bash
pip install --only-binary :all: --default-timeout=120 --retries 10 cvxpy
```

---

## Verwendung

Im Projektordner mit aktivierter virtueller Umgebung ausführen.

```bash
python main.py lqr     # LQR/Kalman-Verstärkungen auf einem 5-Punkte-θ-Gitter
python main.py lmi     # LMI-Verstärkungen mit gemeinsamer Lyapunov-Matrix (benötigt cvxpy)
```

Jeder Lauf:

1. entwirft die Verstärkungsplanung,
2. gibt die Bezout-/Faktorisierungsfehler bei θ = 0, 0,5, 1 aus (sollten ≈ 1e-15 sein),
3. gibt die Pole des stabilen Anti-Windup-Faktors bei θ = 0,5 aus,
4. simuliert die drei Fälle (einige Minuten — reines Python-RK4, je 60 000 Schritte),
5. gibt IAE und maximalen Integratorzustand aus,
6. speichert `logs_<methode>.pkl` und `results_<methode>.png`.

Neu zeichnen ohne erneute Simulation:

```bash
python plot_only.py lqr
python plot_only.py lmi
```

> **Hinweis:** `plot_only.py <methode>` funktioniert erst, nachdem `main.py <methode>` einmal gelaufen ist, da `main.py` die Log-Datei erzeugt. Laden Sie nur `.pkl`-Dateien, die Sie selbst erzeugt haben.

Für schnellere Testläufe die Schrittweite in `main.py` erhöhen, z. B. `simulate(..., dt=2e-3)`.

---

## Beispielstrecke und Szenario

$$\ddot x = a(\theta)\,x - 0{,}4\,\dot x + u,\qquad a(\theta) = -1 + 0{,}9\,\theta \in[-1,-0{,}1]$$

- Zustände $x=[x_1,\dot x_1]$, Ausgang $y=x_1$, skalare Stellgröße $u$.
- Stellgrößenbegrenzung $|u|\le 1$.
- Planung: $\theta(t)=0{,}5+0{,}5\sin(2\pi t/40)$.
- Referenz: $r=2$ für $2\le t<22$, $r=-1{,}5$ für $22\le t<42$, danach $r=0{,}5$.

Die Strecke ist bewusst **offen-schleifig stabil**: Bei beschränkter Stellgröße kann kein Anti-Windup-Verfahren eine instabile Strecke global stabilisieren.

---

## Erwartete Ergebnisse

IAE = $\int|r-y|\,dt$; maximaler Integrator = $\max|\xi|$.

**LQR-Verstärkungen** (Referenzlauf, reproduzierbar):

| Fall | IAE | Max. Integrator |
|---|---|---|
| Linear (ohne Sättigung) | 13,40 | 3,59 |
| Gesättigt, ohne Anti-Windup | 23,22 | 7,35 |
| Gesättigt, YK-Anti-Windup | 17,85 | 3,59 |

**LMI-Verstärkungen** (ungefähr; Werte können je nach Solver-Version leicht abweichen):

| Fall | IAE | Max. Integrator |
|---|---|---|
| Linear (ohne Sättigung) | ≈ 15,8 | ≈ 4,2 |
| Gesättigt, ohne Anti-Windup | ≈ 24,9 | ≈ 7,5 |
| Gesättigt, YK-Anti-Windup | ≈ 19,8 | ≈ 4,2 |

Mit Anti-Windup bleibt der Integratorzustand auf dem Niveau der linearen Antwort, und der Regelfehler liegt deutlich näher am unbeschränkten Fall. Die LMI-Verstärkungen sind etwas langsamer als LQR, weil sie auf einen garantierten Stabilitätsabstand ausgelegt sind und nicht auf ein Gütemaß.

Anti-Windup-Pole bei θ = 0,5: LQR ≈ `[-6, -1.263, -2.732]`, LMI ≈ `[-6, -2 ± 1.83j]`.

---

## Anpassung an die eigene Strecke

**`plant.py`** bearbeiten:

- `A0, A1, B, C` und die Dimensionen `n, m, p` ersetzen.
- `B` und `C` konstant und `A` affin in θ halten, wenn die LMI-Garantie gelten soll.
- `th_fun(t)` und `r_fun(t)` für das eigene Szenario anpassen.

Stellschrauben:

| Ort | Parameter | Wirkung |
|---|---|---|
| `design.py` → `design_lqr` | `Qw, Rw` | Regelgeschwindigkeit vs. Stellaufwand |
| `design.py` → `design_lqr` | `W, V` | Beobachtergeschwindigkeit vs. Rauschen |
| `design.py` → `design_lmi` | `alpha_F, alpha_L` | garantierte Abklingraten (zu groß → unlösbar) |
| `controller.py` / `simulate.py` | `kaw` | Anti-Windup-Rate des Integrators |
| `main.py` | `umax` | Stellgrößenbegrenzung |

Die LMIs sind homogen in $(P,Y)$ bzw. $(Q,W)$; die Normierungen `P >> I` und `Q >> I` in `design_lmi` sind essenziell und dürfen nicht entfernt werden.

---

## Einschränkungen

- Nur quadratische Strecken ($m=p$); die Anti-Windup-Verstärkung verwendet $F_i^{-1}$.
- Skalarer Planungsparameter, $D=0$.
- Die LMI-Garantie setzt konstantes $B$ und $C$ sowie eine polytopische (konvexe) Parameterabhängigkeit voraus. Für $B(\theta)$ ist eine parameterabhängige Lyapunov-Formulierung nötig.
- Die LMI garantiert die Stabilität des **linearen** Kreises. Während der Sättigung folgt die Beschränktheit des Reglerzustands aus der Anti-Windup-Struktur. Globale Stabilität unter Sättigung (Kreiskriterium / Sektor-LMI) ist **nicht** implementiert und nur für offen-schleifig stabile Strecken möglich; bei instabilen Strecken gibt es nur regionale Garantien.
- Die Option `lqr` liefert keine Garantie für zeitvariantes θ(t).
- Die LQR-Verstärkung stabilisiert den *erweiterten* Kreis; $A+BF_x$ allein ist nicht zwingend stabil, was nur für die diagnostische Stabilität des Faktors $M$ relevant ist.

---

## Fehlersuche

| Symptom | Wahrscheinliche Ursache |
|---|---|
| `ImportError: cannot import name ... from 'plant'` | Funktion fehlt, falsch geschrieben oder Datei nicht gespeichert |
| `ValueError: inhomogeneous shape` in `plant.py` | Komma statt Dezimalpunkt in einer Matrix |
| `np.block`-Formfehler in `ctrl_mats` | falsches Aufteilen von `Fa`; richtig: `Fx=Fa[:, :n]`, `Fi=Fa[:, n:]` |
| Terminal mit Text überflutet | Debug-`print` in `ctrl_mats` vergessen (wird ~5 × pro Zeitschritt aufgerufen) |
| `FileNotFoundError: logs_<methode>.pkl` | zuerst `main.py <methode>` ausführen |
| cvxpy-Formfehler | Beobachterblock braucht `np.eye(n)`, Zustandsrückführungsblock `np.eye(na)` |
| LMI-Verstärkungen instabil / sinnlos | Normierung `P >> I` / `Q >> I` fehlt |
| Kein Plotfenster (VM / ohne Display) | die gespeicherte Datei `results_<methode>.png` öffnen |

---

## Literatur

> Die hier verwendete Struktur aus koprimer Faktorisierung und Anti-Windup ist nach der Standardliteratur umgesetzt; bitte prüfen Sie die genauen Formulierungen anhand der Quellen, bevor Sie diese zitieren.

- S. X. Ding, *Model-Based Fault Diagnosis Techniques*, Springer, 2013 (koprime Faktorisierung, beobachterbasierte Parametrierung stabilisierender Regler).
- K. Zhou, J. C. Doyle, K. Glover, *Robust and Optimal Control*, Prentice Hall, 1996 (doppelt koprime Faktorisierung, Youla-Parametrierung).
- H. H. Niemann, J. Stoustrup, Arbeiten zu Verstärkungsplanung und Umschaltung mittels Youla-Parametrierung (z. B. *Gain scheduling using the Youla parameterization*, IEEE CDC, 1999).
- M. V. Kothare, P. J. Campo, M. Morari, C. N. Nett, *A unified framework for the study of anti-windup designs*, Automatica, 1994 (Anti-Windup mit koprimer Faktorisierung).

---

## Lizenz

Wählen Sie eine Lizenz für Ihr Repository (z. B. MIT) und fügen Sie eine `LICENSE`-Datei hinzu.

**Autor:** *Ihr Name / Ihre Institution*
