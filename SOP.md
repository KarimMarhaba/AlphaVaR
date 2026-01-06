# 📋 AlphaVaR Dev-SOP (Positron)

**Wann nutzen?** Bei jedem neuen Code-Input (Feature, Fix, Test) durch die AI.

## Phase 1: Integration & Dependencies
*Vermeidet "Namespace missing" & "Global Variables" Fehler.*

- [ ] **Code einfügen:** Neue Datei in `R/` erstellen oder bestehende updaten.
- [ ] **Neue Pakete?** (z.B. neues `library(pkg)` oder `pkg::fun`)
    - 👉 Konsole: `usethis::use_package("paketname")`
- [ ] **Global Variables?** (Spaltennamen ohne Quotes, z.B. in `filter()`)
    - 👉 `R/globals.R` öffnen und Namen zur Liste hinzufügen.

## Phase 2: Dokumentation & Namespace
*Vermeidet "Function not found" Fehler.*

- [ ] **Roxygen Header:**
    - Steht `#' @export` über der Funktion? (Damit User sie sehen).
    - Steht `#' @importFrom` dabei? (Falls extern).
- [ ] **Generieren:**
    - 👉 Konsole: `devtools::document()`
    - *Check:* Hat sich die Datei `NAMESPACE` verändert?

## Phase 3: Testing (TDD)
*Vermeidet Regressionen.*

- [ ] **Test-Datei:** In `tests/testthat/` Datei erstellen/updaten.
- [ ] **Ausführen:**
    - 👉 Konsole: `devtools::test()`
    - *Ziel:* Alle Tests **PASS**.

## Phase 4: R CMD Check
*Publikationsreife prüfen.*

- [ ] **Check:**
    - 👉 Konsole: `devtools::check()`
- [ ] **Ergebnis:**
    - Errors: **0**
    - Warnings: **0** (Ausnahme: qpdf)
    - Notes: Minimieren (oft `globals.R` vergessen).

## Phase 5: Abschluss & Git
- [ ] **Lokal Installieren:** (Für eigene Nutzung)
    - 👉 Konsole: `devtools::install()`
- [ ] **Git Push:**
    - Source Control (links) -> Message ("Feat: ...") -> Commit -> Sync.
- [ ] **Context Update:**
    - Falls sich Struktur drastisch änderte: `DEVELOPER_CONTEXT.md` anpassen.

---

### 🚀 Quick-Command Referenz (Konsole)

```r
# 1. Neue Abhängigkeit (einmalig)
usethis::use_package("paketname")

# 2. Doku & Namespace update (Oft!)
devtools::document()

# 3. Alles laden (zum Rumprobieren)
devtools::load_all()

# 4. Tests
devtools::test()

# 5. Full Check
devtools::check()

# 6. Installieren
devtools::install()