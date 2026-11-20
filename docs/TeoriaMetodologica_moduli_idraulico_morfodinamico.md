# Ricostruzione batimetrica automatica – Approccio teorico dei moduli idraulico e morfodinamico

**Elisa Caccamo** · 8 ottobre 2026

> Versione Markdown del file `TeoriaMetodologica_moduli_idraulico_morfodinamico.docx`. Le equazioni sono state ricostruite in LaTeX a partire dal testo estratto; la figura del §4.3 non era presente nel testo ed è indicata come segnaposto.

---

## Stato del lavoro

Resoconto breve per giornata, dal più recente. Il dettaglio è nella §9; i punti aperti sono nella tabella della §9.4.

**20/11/2026, run della sera** (§9.14; resoconto in `Final/Output_folder/RESOCONTO.md`, run precedente in `_run_20-11/`). Run con la maschera del tratto 0 accorciata prima dell'ombra del ponte e la regola più stretta per i salti isolati. Nessun errore.

- **Punto 8 chiuso:** la regola corregge solo le due sezioni attese (ex 728 ed ex 1239, ora 726 e 1237), con gli stessi valori; la 492 (ora 490) resta al valore grezzo.
- **Punto 19 chiuso:** la risalita a fine tratto 0 è sparita (ultime 5 sezioni 260.10–260.02 m) e la pendenza del tratto 0 è quella della retta (2.0·10⁻⁴). Il tratto 0 accorciato finisce a 61–68 m dai tratti 1 e 2, oltre `distance_threshold` (51 m): la biforcazione 0 → 1, 2 è stata aggiunta a `forced_bifurcations` nel config.
- **Numerazione:** il tratto 0 ha 2 sezioni in meno, quindi tutti gli ID a valle scendono di 2 (nodo 15 → 13, 1070 → 1068, 1369 → 1367).
- **Validazione invariata** (nessuna stazione nei tratti 0–3): 48 GNSS RMSE 0.30 m, errore medio +0.22 m; 59 stazioni RMSE 0.40 m, +0.34 m; test del segno 5 / 7 / 36, p = 0.77.
- **Nuovo punto 21, $S_0$ al nodo 13:** il nodo è alla testa della rete e la finestra a monte copre solo il tratto 0 (114 m, quasi piatto), che dà una pendenza sotto la soglia minima; $S_0$ viene alzata alla soglia (1.22·10⁻³). Il nodo passa da asimmetrico ($\psi$ = 0.82) a bilanciato ($\psi$ = 0.73, una sola soluzione, $D_0$ = 2.53 m, $\beta$ = 10.9). Da decidere come trattare $S_0$ ai nodi di testa.
- Le ultime due sezioni di ogni tratto hanno la stessa quota per costruzione (in `calculate_widths` la prima e l'ultima sezione copiano la vicina): non è un'anomalia.
- Corretti i commenti "Chezy's eq" nel sistema di biforcazione del modulo idraulico (la formula è Manning).

**20/11/2026, decisioni dopo il run.**

- **Regola per i salti isolati più stretta:** una sezione è un picco solo se è un estremo locale e si stacca da entrambe le vicine di più di glitch_tol. Sui dati del 20/11 corregge solo i due picchi veri (728, 1239) e lascia i gradini (492, 14). Implementata in `Preprocessing_HydraulicModule.m` (§9.6); da verificare al prossimo run.
- **Fine del tratto 0 (punto 19):** il tratto 0 finisce subito a monte di un ponte, e i tratti 1 e 2 cominciano subito a valle. Elisa ha confermato che la maschera del tratto 0 finisce sull'ombra del ponte e l'ha accorciata (nuovi `I_000.tif` e `dtm_000.tif`, stessa griglia; quote del DTM ritagliato 259.92–260.41 m, la risalita a 261.36 m non c'è più). La prima esportazione della maschera era float32 con nodata −3.4·10³⁸ fuori dal poligono, che `centerline_from_mask` avrebbe scartato (somma negativa): `I_000.tif` è stata riesportata come int16 0/1 con nodata 0, come le altre (stessa griglia, stessi 16412 pixel). Requisito degli input: maschere int16 0/1. Il codice non cambia (la conversione aggiunta per poco tempo è stata tolta). Prossimo run con maschera nuova e regola più stretta.
- **s035:** resta sul tratto 20 (89% dei punti), assegnazione automatica confermata.
- **Introduzione e Discussione del paper:** la preparazione manuale degli input non va enfatizzata (Scaletta aggiornata).
- **s015-t30:** resta nella validazione, esclusa dal confronto tirante–profondità della §9.11.
- **Sezione 256:** segnalata (sponda unica specchiata, $W_2$ del nodo 244 da 78.3 a 88.1 m), come esempio della fragilità della stima della larghezza bankfull (punto 1).
- **Stazioni per la validazione principale:** in sospeso fino alla risposta di Carlo. **Sensitività:** rimandata a dopo la chiusura dei punti 8 e 19.
- **Pozze:** confronto con e senza pozze fatto (§9.11, "Con e senza pozze"): stazioni in pozza +0.34 m di errore medio, le altre +0.14 m.
- **Baseline semplice:** proposta registrata (punto 20).

**20/11/2026, run** (§9.13; resoconto in `Final/Output_folder/RESOCONTO.md`). Primo run con la regola unica per i salti isolati, s015 sul tratto 30, i nomi `sNNN-tCC` e $Q$ = 10.76 m³/s. Nessun errore.

- **Regola per i salti isolati:** corregge esattamente le 4 sezioni previste (14, 492, 728, 1239) e lascia intatti i gradini; nessuna stazione cambia errore medio di più di 4 cm.
- **Validazione identica al 19/11 entro 1 mm:** 48 GNSS RMSE 0.30 m, errore medio +0.22 m in entrambi i moduli; moduli non distinguibili (5 / 7 / 36, p = 0.77). Errore medio confermato dovuto soprattutto alle pozze (correlazione 0.90 con la profondità osservata).
- **Nomi e figure corretti** (punto 17 chiuso); nessun nodo cambia tipo o soluzione; portata 10.76 con effetto trascurabile.
- **Da decidere:**
  - sezione 492, piede di un gradino che la regola porta a metà (punto 8);
  - fine del tratto 0 al nodo 15, pelo libero che risale di 1.1–1.3 m nelle ultime sezioni, probabilmente per la barra tra i rami nel transetto (punto 19);
  - s015-t30, sezione di riferimento al nodo: va esclusa dalle analisi basate sulla sezione di riferimento;
  - insieme di stazioni e parametri per l'analisi di sensitività, di cui questo run è il riferimento.

**20/11/2026.** Corrette nel codice le anomalie minori segnalate dal run del 19/11 (punto 17 della §9.4) e valutato se serve ancora `segments_to_smooth` (§9.10, "Lisciamento").

- **Figure:** le sottostazioni si chiamano ora `sNNN-tCC` (trattino invece del trattino basso, che MATLAB rendeva come pedice); titoli ed etichette dicono "in situ" (GNSS + ADCP) invece di "GNSS". I nomi dei file di output restano invariati.
- **Tabella diagnostica dei tratti sospetti** (`Preprocessing_HydraulicModule.m`): ora è calcolata sulla $z_{DTM}$ grezza (prima del filtro e del lisciamento), su tutti i tratti di `segments_to_smooth` e su quelli già sospetti, con rugosità, salto massimo, effetto del lisciamento ed esito "smoothing not needed" / "smoothing removes it" / "to be inspected manually".
- **Lisciamento:** con il nodata corretto `segments_to_smooth` non serve più per le oscillazioni; toglie solo tre picchi isolati di 0.2–0.4 m (tratti 26, 32, 12), mentre picchi simili restano in tratti fuori dalla lista (34, 35, 37). Da decidere (punto 8): toglierlo, o sostituirlo con una regola uguale per tutti i tratti.
- **Lisciamento, scelta B (20/11, pomeriggio):** `segments_to_smooth` e il vecchio filtro sui salti (soglia 1 m) sono sostituiti da una sola regola per i salti isolati, uguale per tutti i tratti (`glitch_tol` = 0.20 m nel config; §9.6, "Regola per i salti isolati"). Con i dati del 19/11 corregge 4 sezioni (14, 492, 728, 1239) e non tocca i gradini verso valle che la mediana mobile appiattiva.
- **s015:** assegnata a mano al tratto madre 30 (sta al nodo 1070, tra la fine del 30 e l'inizio di 31 e 32); diventa `s015-t30`.
- **Errore medio positivo** (§9.11): viene soprattutto dalle pozze, dove il pelo libero è quasi piatto e il moto uniforme con la pendenza del tratto sottostima il tirante; residuo correlato 0.91 con la profondità osservata. Da chiarire il livello del volo usato per la portata (0.695 m nel config, 97 cm la lettura del 17/3 a San Benigno).
- **ADCP:** spiegato il problema della copertura e aggiunto alla domanda 6 per Carlo.
- **Prossimo run:** verificare i nomi delle sottostazioni, i titoli, la nuova regola per i salti isolati e la tabella diagnostica.
- **Decisioni del 20/11 (sera):**
  - lisciamento: questione chiusa con la regola unica per i salti isolati; resta solo la verifica al prossimo run (punto 8);
  - pozze e sezioni molto irregolari: limite del metodo, dichiarato e non correggibile con la rettangolarizzazione (punto 16, §9.11);
  - s015 al tratto 30: confermato.
- **Portata del volo (§9.12), chiusa:** 10.76 m³/s, portata giornaliera ARPA del 17/03/2025, coerente con il livello 0.97 m e la scala di deflusso valida dal 01/01 al 16/04/2025 (10.9 m³/s). La scala 2023, che dava 23.4 m³/s, non vale per il 2025. Nel config `Q_upstream` = 10.76 e `level_upstream` = 0.97 (solo valore iniziale del solutore). L'errore medio positivo non dipende dalla portata.
- **Analisi di sensitività:** piano approvato; la portata è fissata, si può impostare.

**19/11/2026, run con la correzione del nodata, gli ADCP e le sottostazioni** (§9.10; resoconto completo in `Final/Output_folder/RESOCONTO.md`). Pre-processing + modulo idraulico, modulo morfodinamico e `ModuleComparison` girano senza errori.

- **Correzione del nodata riuscita:** salti segnalati dal filtro da 149 a 2; i tratti di `segments_to_smooth` non oscillano più nemmeno senza lisciamento; a s006, s019, s020 $z_{DTM}$ coincide con il DTM completo entro 4 cm. Residui medi: s006 −1.47 → −0.21, s019 −0.77 → +0.27, s020 −0.45 → −0.04 m.
- **48 stazioni GNSS originali:** RMSE da 0.36 a 0.30 m in entrambi i moduli (senza s006 resta 0.30, come il 18/11); errore medio da +0.12 a +0.22 m, ora l'errore principale (punto 11).
- **Moduli non distinguibili:** RMSE MM − HM +0.4 cm [−0.2, +1.2]; MM 5 / HM 7 / 36 pari, p = 0.77. Il `compare_modules_gnss.m` corretto, sui CSV del 18/11, dà 3 / 8 / 37 (punto 14 chiuso).
- **Tutte le 59 stazioni** (con sottostazioni e ADCP): RMSE 0.40 m, errore medio +0.33 m; 6 / 10 / 43, p = 0.45. Gli 8 transetti ADCP sono il 14% delle stazioni ma il 51% dei punti, e hanno un residuo 0.3–0.6 m più alto delle GNSS vicine: **i transetti ADCP coprono soprattutto le parti profonde, dove non si poteva guadare**, quindi la loro media di stazione è più bassa della media di sezione (§9.9).
- **Rete:** portata conservata ovunque, 10.80 m³/s all'uscita; 1369 torna asimmetrico ($\beta$ = 14.9, soluzione quasi bilanciata); 15 cambia soluzione ($\psi$ da 0.41 a 0.82) solo perché $\psi_{hyd}$ passa da 0.49 a 0.67.
- **Sottostazioni:** nessuno dei due moduli riproduce il dislivello tra i rami dello stesso transetto (s003: ramo 12 0.4–0.6 m sotto gli altri nei dati; s007: i due canali alla stessa quota nei dati, gradino di 0.6 m nei modelli).
- **Anomalie minori nelle figure e nel log** (nessuna modifica fatta): pedici TeX nei nomi con "_", titoli "GNSS" anche con i dati ADCP, tabella diagnostica dei tratti sospetti che scrive "apply smoothing" per tutti.
- **Da decidere (scaletta in §5 del resoconto):** se servono ancora profilo monotono e `segments_to_smooth` (restano risalite di 10–55 cm e il salto di +1.1 m in fondo al tratto 0); origine dell'errore medio positivo (livello idrometrico 7/3 contro 17/3, scabrezza); come entrano gli ADCP nella validazione e con quale peso (domanda 6 a Carlo); regola di scelta della soluzione e larghezza bankfull ai nodi, da mettere nella sensitività.

**19/11/2026.** Analizzati i run del 18/11 sera (pre-processing + modulo idraulico, modulo morfodinamico, `ModuleComparison`) con la rete corretta: portata propagata solo dentro il tratto e tratti 15, 16, 29 senza portata (§9.8).

- **Portata conservata ovunque:** 10.80 m³/s all'uscita in entrambi i moduli, nessun nodo sbilanciato. Punti 6 e 7 della §9.4 chiusi.
- **Nodi:** 10 risolti dal sistema (366 e 817 esclusi); 1369 diventa bilanciato ($\beta$ = 11.7, sotto la soglia), gli altri 9 asimmetrici con $\psi$ entro ±0.02 dal run precedente.
- **Validazione GNSS (48 stazioni):** RMSE da 0.60 (DTM grezzo) a 0.36 m in entrambi i moduli; errore medio da +0.56 a +0.12 m. Senza s006: 0.30 m.
- **I due moduli non sono distinguibili:** differenza di RMSE +0.4 cm [−0.1, +0.9]; MM più vicino in 3 stazioni, HM in 8, 37 pari (p ≈ 0.23). Differiscono solo nelle 11 stazioni dei tratti a canali multipli.
- **$z_{DTM}$ resta l'errore dominante:** confermata alla sezione 170 la causa 1 della §9.6 (nodata −32767 dei DTM di segmento non convertito), con il filtro sui salti che sposta l'errore sulle sezioni sane. Correzione in attesa di decisione.
- **Corretto `compare_modules_gnss.m`:** il test del segno esclude le parità (|ΔRMSE| ≤ 1 mm, parametro `tie_tol`); la tabella delle stazioni riporta sia la differenza di RMSE sia quella dell'errore medio. Il resoconto del run contava come "MM più vicino" le stazioni con fondo MM più basso: i numeri corretti sono 3 contro 8 (punto 14 della §9.4). Va rilanciato `ModuleComparison`.
- **Domande per Carlo:** in un file a parte, `docs/Domande_per_Carlo.md` (classi delle stazioni, $\theta^*_{BF}$, portata bankfull, regola di scelta, ruolo del trasporto al volo, densità dei punti ADCP). La vecchia §8 è stata tolta da questo documento.
- **Dati ADCP nella validazione:** gli 8 transetti ADCP del 07/03/2025 entrano come stazioni a01–a08, senza distinzione dalle sezioni GNSS; quota ellissoidica verificata. Le sezioni GNSS s003 e s007, che attraversano più tratti, sono divise in sottostazioni assegnate a mano (s015 da definire). Nuova funzione `read_field_points.m` (§9.9). La densità dei punti ADCP resta in sospeso fino alla risposta di Carlo.
- **Correzione di $z_{DTM}$, primo passo (deciso il 19/11):** si tiene l'approccio attuale e si corregge solo il nodata. In `Preprocessing_HydraulicModule.m` il nodata −32767 dei DTM di segmento diventa NaN alla lettura; in `calculate_widths.m` la media della sezione esclude esplicitamente i NaN. Da rilanciare il pre-processing e guardare i risultati (salti segnalati dal filtro, oscillazioni dei tratti di `segments_to_smooth`, $z_{DTM}$ contro il DTM completo a s006, s019, s020) prima di decidere sul profilo monotono (§9.6).
- **Prossimi passi:** rilanciare i tre script (pre-processing, modulo morfodinamico, `ModuleComparison`) con la correzione del nodata, i dati ADCP e le sottostazioni; definire la divisione di s015; decidere come presentare il contributo del modulo morfodinamico nel paper.

---

## 1. Introduzione e notazione

Il framework corregge le quote del DTM nell'area bagnata dell'Orco ricostruendo la quota del fondo $\eta$ di ogni sezione trasversale rettangolarizzata. Due moduli condividono la stessa geometria e la stessa topologia di rete, ma usano principi diversi: il **modulo idraulico** impone solo vincoli di moto uniforme e di continuità della portata liquida; il **modulo morfodinamico** aggiunge il trasporto solido frazionario (due classi, grossolana C e fine F) e la continuità del sedimento, in condizioni di equilibrio.

In entrambi i moduli il calcolo procede da monte verso valle e ogni sezione viene classificata, tramite la matrice di adiacenza, in una di tre tipologie: canale singolo, nodo di biforcazione, confluenza. Le sezioni successive del documento seguono questa stessa griglia.

**Notazione.** Pedice 0 = canale alimentatore (sezione madre) o canale di valle della confluenza; pedici 1 e 2 = rami figli della biforcazione o rami affluenti della confluenza. Le grandezze sono:

| Simbolo | Significato | Unità |
| --- | --- | --- |
| $z_{DTM}$ | quota grezza del DTM sulla centerline | m s.l.m. |
| $\eta$ | quota del fondo ricostruita | m s.l.m. |
| $h$, $D$ | tirante ($h$ nel modulo idraulico, $D$ nel morfodinamico, come nel paper di Ragno) | m |
| $W$ | larghezza della sezione | m |
| $S$ | pendenza locale del tratto | – |
| $Q$, $q$ | portata liquida totale e per unità di larghezza ($q=Q/W$) | m³/s, m²/s |
| $q_{sx}$ | trasporto solido per unità di larghezza della frazione $x$ (C o F) | m²/s |
| $f_x$ | frazione della classe $x$ sulla superficie del letto ($f_C+f_F=1$) | – |
| $\theta_g$ | numero di Shields riferito al diametro medio geometrico $d_g$ | – |
| $\Delta$ | peso relativo sommerso del sedimento, $(\rho_s-\rho)/\rho$ | – |
| $W^m$ | larghezza osservata dalla maschera il giorno del volo LiDAR | m |
| $W^{bf}$ | larghezza bankfull teorica, larghezza rettangolare equivalente $A/D$ dal DTM corretto | m |
| $D_{i,s}$ | tirante di Shields del ramo $i$, usato solo per misurare la geometria bankfull del ramo | m |
| $\psi$ | quota di portata del ramo 1 in condizioni bankfull, $\psi=Q_1^{bf}/Q_0^{bf}$; al volo la portata si ripartisce con la stessa $\psi$ | – |
| $W_1/W_0$ | quota di larghezza del ramo 1: è la $\psi$ della soluzione bilanciata (portata ripartita in proporzione alle larghezze, $D_1=D_2=D_0$) | – |
| $\psi_{hyd}$ | quota di portata del ramo 1 nel modulo idraulico al volo, $Q_1^{hyd}/Q_0$, usata come riferimento per scegliere tra le soluzioni asimmetriche (§9.1) | – |
| $\psi_x$ | quota del flusso solido della classe $x$ (C o F) che va nel ramo 1, $Q_{s1,x}/Q_{s0,x}$ | – |
| $\beta$ | rapporto di forma della madre in condizioni bankfull, $\beta=W_0^{bf}/(2D_0)$ (semilarghezza su tirante). Sotto un valore critico (≈ 12–13 con rami uguali) il nodo ha solo la soluzione bilanciata; sopra compaiono due soluzioni asimmetriche e la bilanciata diventa instabile (Bolla Pittaluga et al., 2003) | – |
| $\theta^*_c$, $\theta^*_{BF}$ | numero di Shields critico (0.03) e bankfull ($1.62\,\theta^*_c$), riferiti a $d_{50}$ | – |
| $\theta_{bf}/\theta^*_c$ | Shields bankfull effettivo della madre diviso per quello critico, con $\theta_{bf}=D_0 S_0/(\Delta d_{50})$. Dice quale limite fissa $D_0$ (§9.1): 1.62 = tirante di Shields pieno $D_S$; 1.00 = tirante critico $D_c$ (sponde basse); valori intermedi = livello di sfioro del transetto | – |

---

## 2. Base comune ai due moduli

Entrambi i moduli partono dalla stessa rappresentazione della rete, costruita una sola volta.

- **Geometria.** Per ogni poligono di canale singolo si estrae la centerline (RivMAP, `centerline_from_mask`), la si ricampiona ogni $dx=10$ m e si traccia la sezione ortogonale. La larghezza $W$ è la distanza tra il primo e l'ultimo punto della sezione che cadono nella maschera binaria.

- **Pendenza.** $S$ è costante per tratto e si ricava nel pre-processing dalla quota del DTM lungo la centerline (quota media di sezione, dopo la rimozione degli outlier): si interpola $z(s)$ con una retta ai minimi quadrati, $z=a\,s+b$, dove $s$ è l'ascissa curvilinea lungo la centerline, e si pone $S=\max\left(1.6\cdot10^{-4},\,-a\right)$. Ai tratti con un solo punto si assegna $S=1.6\cdot10^{-4}$. Le pendenze sono calcolate una sola volta e usate da entrambi i moduli.

- **Topologia.** Due tratti sono connessi se la distanza tra la fine dell'uno e l'inizio dell'altro è inferiore a una soglia di distanza prefissata. La matrice antisimmetrica $a_{ij}$ vale $+1$ se il tratto $j$ alimenta il tratto $i$ e $-1$ nella posizione simmetrica.

- **Classificazione.** Righe e colonne della matrice sono gli ID di inizio e di fine di ogni tratto. Se la fine del tratto $j$ è vicina all'inizio del tratto $i$, $a_{\text{fine}_j,\,\text{inizio}_i}=+1$ e $a_{\text{inizio}_i,\,\text{fine}_j}=-1$. Una colonna di **inizio** con più di un $+1$ indica una sezione che riceve acqua da più tratti (le righe): è una **confluenza** (`conf_info`). Una colonna di **fine** con più di un $-1$ indica una sezione che alimenta più tratti: è una **biforcazione** (`bif_info`). Tutte le altre sezioni sono **canale singolo**.

La portata di monte proviene dalle stazioni idrometriche (Spineto, San Benigno) ed è l'unica condizione al contorno idraulica. Il risultato di ogni modulo è la quota $\eta=z_{DTM}-h$ su ogni sezione; le quote sono poi interpolate (scattered interpolant) sulla griglia del DTM e si produce la mappa delle differenze (DoD).

---

## 3. Modulo idraulico

Il modulo idraulico ricava il tirante $h$ di ogni sezione rettangolare dal moto uniforme di Manning (forma di Gauckler–Strickler, $k_s=1/n$, $n=0.03$ costante) con la portata liquida propagata da monte, e pone $\eta=z_{DTM}-h$. Il coefficiente di Chézy non compare: è legato a $n$ da $C=R^{1/6}/n$ (con $R\approx h$ per sezione larga), quindi con $n$ costante $C$ cresce con il tirante (con $n=0.03$: $C\approx30$ m$^{1/2}$/s per $h=0.5$ m, $\approx36$ m$^{1/2}$/s per $h=1.5$ m). Non c'è trasporto solido: l'unico vincolo di rete è la conservazione della portata liquida, più la continuità dell'energia nei nodi di biforcazione.

$$
\eta = z_{DTM} - h
$$

### 3.1 Canale singolo

In ogni sezione la portata $Q$ nota da monte viene inserita nella legge di moto uniforme per una sezione rettangolare larga (raggio idraulico $\approx h$). Con Gauckler–Strickler ($k_s=1/n$):

$$
Q = W\,h\,k_s h^{1/6}\sqrt{h\,S} = W\,k_s\sqrt{S}\,h^{5/3}
\quad\Longrightarrow\quad
h = \left(\frac{Q\,n}{W\sqrt{S}}\right)^{3/5}
$$

La portata resta invariata e passa alla sezione successiva. Poiché $S$ e $k_s=1/n$, con $n$ di Manning costante, sono costanti nel tratto, tra due sezioni consecutive dello stesso canale vale $W_1h_1^{5/3}=W_2h_2^{5/3}$: il tirante varia solo con la larghezza.

$$
h_2 = \left(\frac{W_1}{W_2}\right)^{3/5} h_1
$$

### 3.2 Nodo di biforcazione

Nella sezione madre la portata $Q_0$ è nota; vanno determinate la ripartizione $Q_1$, $Q_2$ e i tiranti $h_1$, $h_2$ dei rami figli. Il sistema ha 4 equazioni in 4 incognite ($h_1$, $h_2$, $Q_1$, $Q_2$): moto uniforme in ciascun ramo, conservazione della portata al nodo e uguaglianza dell'energia specifica all'imbocco dei due rami (Bernoulli, perdite trascurate).

$$
h_1 = \left(\frac{Q_1\,n}{W_1\sqrt{S_1}}\right)^{3/5} \tag{1}
$$

$$
h_2 = \left(\frac{Q_2\,n}{W_2\sqrt{S_2}}\right)^{3/5} \tag{2}
$$

$$
Q_0 = Q_1 + Q_2 \tag{3}
$$

$$
h_1 + \frac{Q_1^2}{2g\,W_1^2h_1^2} = h_2 + \frac{Q_2^2}{2g\,W_2^2h_2^2} \tag{4}
$$

Il sistema è non lineare e si risolve numericamente; ne seguono $\eta_1=z_{DTM,1}-h_1$ e $\eta_2=z_{DTM,2}-h_2$. La ripartizione è governata solo da geometria ($W$, $S$) e scabrezza: la quota del fondo non entra come incognita, perché è una conseguenza del tirante.

Se il sistema non ha soluzione fisica, la portata viene ripartita in proporzione alle larghezze dei rami.

### 3.3 Confluenza

Alla sezione di confluenza si somma la portata dei due affluenti, $Q_0=Q_1+Q_2$. La portata totale entra nella legge di moto uniforme della sezione di valle (eq. di §3.1 con $W_0$, $S_0$) e fornisce $h_0$ e $\eta_0$. Non si impone alcuna condizione energetica tra i rami affluenti.

---

## 4. Modulo morfodinamico – base teorica

Il modulo morfodinamico impone che la rete sia in equilibrio: in ogni ramo il trasporto solido di ciascuna classe granulometrica è quello compatibile con l'idraulica locale, e nei nodi i flussi di sedimento si conservano. Il riferimento è Ragno et al. (2023), che estende a una miscela bimodale il modello a due celle di Bolla Pittaluga, Repetto e Tubino (2003, "BRT").

### 4.1 Granulometria della miscela

La superficie del letto è una miscela di due classi, grossolana (C) e fine (F), con frazioni $f_C$ e $f_F=1-f_C$. Si usa la scala sedimentologica $\phi$ e il diametro medio geometrico $d_g$, con $d_{ref}=1$ mm:

$$
\phi = -\log_2 d\,[\text{mm}], \qquad
\phi_m = f_C\,\phi_C + f_F\,\phi_F, \qquad
d_g = d_{ref}\,2^{-\phi_m}
$$

**Dati di campo.** La miscela di riferimento è la media delle due curve granulometriche rilevate il 28/02/2019 (Pratoregio) e il 18/03/2019 (meandro di San Benigno). La soglia tra le classi è 2 mm: le frazioni $f_{0F}$ e $f_{0C}$ del canale di monte e i diametri rappresentativi $d_F$ e $d_C$ derivano da queste curve, e $d_{50}$ si ottiene per interpolazione della curva media. $d_{50}$ si usa solo per il tirante bankfull (§4.3); tutte le grandezze di trasporto usano $d_g$.

I diametri rappresentativi delle classi sono medie geometriche in scala $\phi$: ogni intervallo tra due setacci è rappresentato dal suo punto medio $\phi_k$, e la classe grossolana dalla media pesata sulle frazioni $f_k$. La classe fine, priva di setaccio inferiore, è rappresentata da $d_F=1$ mm.

$$
\phi_k = -\log_2\sqrt{d_{k-1}\,d_k}, \qquad
d_C = 2^{-\sum_k f_k\,\phi_k / f_C}, \qquad
d_F = 1\ \text{mm}
$$

### 4.2 Trasporto frazionario (Wilcock & Crowe, 2003)

Il trasporto per unità di larghezza della classe di diametro $d$ dipende dallo sforzo al fondo e da una funzione di trasporto che include l'effetto di nascondimento (hiding):

$$
q_s(\phi) = \frac{(\tau/\rho)^{3/2}}{g\,\Delta}\,G(\zeta), \qquad
\zeta = g_{hr}\,\frac{\theta_g}{\theta_r}, \qquad
\theta_g = \frac{\tau}{\rho\,g\,\Delta\,d_g}
$$

$$
G(\zeta) =
\begin{cases}
0.002\,\zeta^{7.5} & \zeta < 1.35 \\[4pt]
14.2\left(1 - \dfrac{0.894}{\zeta^{0.5}}\right)^{4.5} & \zeta \ge 1.35
\end{cases}
$$

$$
g_{hr} = \left(\frac{d}{d_g}\right)^{-b}, \qquad
b = \frac{0.67}{1 + \exp\left(1.5 - d/d_g\right)}
$$

$\theta_r=0.036$ (miscela ghiaiosa senza sabbia significativa) e $\Delta=1.65$. Lo sforzo al fondo in moto uniforme è $\tau=\rho g D S$. Il numero di Shields $\theta_g$ è uno solo per la miscela, calcolato su $d_g$ e comprensivo di $\Delta$; la dipendenza dalla singola classe entra soltanto tramite $g_{hr}$. Il trasporto della classe $x$ sulla superficie del letto è $f_x\,q_{sx}$.

Il diametro medio geometrico $d_g$ è calcolato con la composizione della superficie del tratto considerato ($f_F$, $f_C$) e quindi varia da ramo a ramo.

### 4.3 Condizioni bankfull

Il trasporto viene valutato in condizioni di piena a rive piene, perché a portate ordinarie il sedimento grossolano è quasi immobile e il sistema non avrebbe soluzione significativa. Il tirante bankfull del canale madre si ricava fissando il numero di Shields bankfull come multiplo di quello critico:

$$
\theta^*_{BF} = 1.62\,\theta^*_c, \quad \theta^*_c = 0.03
\quad\Longrightarrow\quad
D_0 = \frac{\theta^*_{BF}\,\Delta\,d_{50}}{S_0}
$$

Il coefficiente 1.62 è quello adottato per fiumi ghiaiosi a canale singolo, con $\theta^*_c=0.03$ e $d_{50}$ dalla curva granulometrica media. Lo stesso criterio fornisce il tirante bankfull di ciascun ramo figlio.

Con $D_0$ si ricava la larghezza bankfull $W_0^{bf}$ dal DTM corretto, nel quale l'area bagnata contiene il fondo ricostruito dal modulo idraulico. Lungo il transetto ortogonale alla centerline (direzione calcolata solo con i punti dello stesso tratto) si individuano le sponde come intersezioni del profilo con il livello $\eta_0+D_0$, si integra l'area $A_0$ tra le sponde e si definisce la larghezza rettangolare equivalente. La portata bankfull segue da Manning:

$$
W_0^{bf} = \frac{A_0}{D_0}, \qquad
Q_0^{bf} = \frac{W_0^{bf}\,D_0^{5/3}\sqrt{S_0}}{n}
$$

Se il livello bankfull $\eta_0+D_0$ non è contenuto entro il transetto (il terreno non lo raggiunge su nessuno dei due lati), le sponde e l'area si calcolano al livello massimo del profilo e la larghezza equivalente è $A_0$ divisa per la profondità a quel livello. È un'approssimazione, che tende a sottostimare la larghezza bankfull.

La larghezza bankfull non può essere minore di quella bagnata osservata il giorno del volo: se $W^{bf}<W^m$ si pone $W^{bf}=W^m$. In tal caso il rapporto $W_0^m/W_0^{bf}$ usato per tornare alle condizioni del volo vale 1. Lo stesso criterio si applica alle larghezze bankfull dei rami, prima della riscalatura su $W_0$.

Per ciascun ramo figlio si calcola con lo stesso criterio il tirante $D_{i,s}$, con la pendenza $S_i$ del ramo, e con esso la larghezza bankfull $W_i^{bf}$ dal DTM corretto. $D_{i,s}$ serve solo a definire la geometria: i tiranti del nodo $D_1$, $D_2$ sono incognite del sistema (§5.2). Per rispettare la chiusura del nodo $W_0=W_1+W_2$, le larghezze dei rami si riscalano conservando il loro rapporto:

$$
D_{i,s} = \frac{\theta^*_{BF}\,\Delta\,d_{50}}{S_i}, \qquad
W_i = W_0^{bf}\,\frac{W_i^{bf}}{W_1^{bf} + W_2^{bf}}
$$

**Ritorno alle condizioni del volo.** Il nodo fornisce la ripartizione bankfull $\psi$. Le portate reali dei rami sono $\psi\,Q_0$ e $(1-\psi)\,Q_0$, con $Q_0$ la portata della madre al momento del volo; la conservazione al nodo è esatta per costruzione. Per i sedimenti si usa la stessa regola per ciascuna classe $x$:

$$
\psi = \frac{Q_1^{bf}}{Q_0^{bf}} = \frac{q_1W_1}{q_0W_0}, \qquad
Q_1 = \psi\,Q_0, \qquad
Q_2 = (1-\psi)\,Q_0
$$

$$
\psi_x = \frac{f_{1x}\,q_{s1x}\,W_1}{f_{0x}\,q_{s0x}\,W_0}, \qquad
Q_{s1,x} = \psi_x\,Q_{s0,x}, \qquad
Q_{s2,x} = Q_{s0,x} - Q_{s1,x}
$$

Le larghezze reali dei rami si riducono tutte dello stesso fattore, così la configurazione del volo resta geometricamente simile a quella bankfull su cui sono calcolati $\psi$ e $\eta_i$:

$$
W_{i,\text{reale}} = W_i\,\frac{W_0^m}{W_0^{bf}}
$$

Le larghezze dei rami nella sezione di imbocco possono quindi differire da quelle delle maschere, che vicino al nodo separano male i due rami; più a valle i rami usano la larghezza osservata. Il fondo dei rami resta la quota $\eta_i$ del nodo (non dipende dalla portata) e il tirante reale all'imbocco è $D_i=z_{DTM,i}-\eta_i$. Più a valle nei rami vale la logica del canale singolo ($Q_i$ costante, $D$ dal moto uniforme, $\eta=z_{DTM}-D$).

> **[Figura]** *Sezione rettangolarizzata al nodo: bankfull equivalente, nodo, ritorno al volo.*
> *(immagine non inclusa nella versione testuale del documento originale)*

Il pannello a mostra come la sezione madre diventa un rettangolo di pari area $A_0$ al livello bankfull, mentre $W_0^m$ è misurata al livello del volo. Il pannello b è il nodo risolto in bankfull, con pelo libero comune e tiranti $D_1$, $D_2$ incogniti. Il pannello c riporta i rami alle condizioni del volo: larghezze ridotte dello stesso fattore, fondi invariati, portate ripartite con $\psi$.

### 4.4 Nodo BRT esteso alle miscele

A monte della biforcazione, su una lunghezza $\alpha W_0$ ($\alpha=4$), il canale madre è diviso in due celle. Il flusso laterale di sedimento tra le celle $q_{sy}$ ha due componenti: il trasporto per la deviazione della corrente ($q_y$) e il richiamo gravitativo lungo la pendenza trasversale, governato dal coefficiente di Ikeda $r$:

$$
q_{sy,x} = q_{s0,x}\left[\frac{q_y}{q_0} - \frac{r}{\sqrt{\theta_{g0}\,l_{hr,x}}}\,\frac{\eta_1 - \eta_2}{W_0/2}\right]
$$

Valori adottati: $\alpha=4$ (Ragno et al., 2021, 2023), $r=0.65$ (intervallo 0.3–1; Ikeda, 1982; Talmon et al., 1995). Per $l_{hr}$ si adotta la "grain independence", $l_{hr,x}=d_{g0}/d_x$, per cui $\theta_{g0}\,l_{hr,x}=\tau/(\rho g\Delta d_x)$ è il numero di Shields della singola classe: la pendenza trasversale agisce di più sulla classe grossolana, meno mobile. $\theta_{g0}$ è il numero di Shields bankfull del canale madre, con $d_{g0}$ calcolato con la composizione del tratto madre:

$$
\theta_{g0} = \frac{\tau_{bf0}}{\rho\,g\,\Delta\,d_{g0}}, \qquad
\tau_{bf0} = \rho\,g\,D_0\,S_0, \qquad
l_{hr,x} = \frac{d_{g0}}{d_x}
\quad\Longrightarrow\quad
\theta_{g0}\,l_{hr,x} = \frac{\tau_{bf0}}{\rho\,g\,\Delta\,d_x}
$$

Il trasporto in ingresso $q_{s0x}$ è la capacità di trasporto del canale madre in condizioni bankfull, con la composizione del tratto madre: dati di campo alla testata della rete, poi quella propagata lungo la rete (§5.1, §5.4). Nel paper la biforcazione è libera e simmetrica ($W_1=W_2=W_0/2$) e il sistema si chiude con: $q_1+q_2=2q_0$, due bilanci di cella per C e F, due continuità di sedimento per C e F e quota del pelo libero uguale nei tre rami ($H_0=H_1=H_2$).

---

## 5. Modulo morfodinamico – logica per tipologia di sezione

Il modulo morfodinamico riceve dal modulo idraulico le quote $\eta$ e i parametri di base, poi propaga verso valle sia la portata liquida sia i flussi solidi delle due classi.

### 5.1 Canale singolo

Il tirante deriva dal moto uniforme, come nel modulo idraulico. Con $D$ e $S$ si calcolano $\tau=\rho g D S$, $\theta_g$ e il trasporto frazionario $q_{sC}$, $q_{sF}$ di Wilcock & Crowe (§4.2). In equilibrio il flusso solido totale di ciascuna classe, $Q_{sx}=f_x\,q_{sx}\,W$, si conserva lungo il tratto: se la larghezza cambia, cambia il trasporto per unità di larghezza ma non il flusso totale. La portata liquida passa invariata alla sezione successiva.

Alla testata della rete $Q_{s,x}=f_x\,q_{sx}\,W$ nella prima sezione. Nei tratti successivi $Q_{s,x}$ e la composizione $f_x$ si ereditano dal tratto a monte e restano costanti lungo il tratto; il trasporto per unità di larghezza è $q_{sx}=Q_{s,x}/W$.

### 5.2 Nodo di biforcazione – il sistema completo (8 equazioni)

Incognite: $q_1$, $q_2$, $\eta_1$, $\eta_2$, $f_{1C}$, $f_{2C}$, $D_1$, $D_2$. Note: $q_0$, $W_0$, $W_1$, $W_2$, $\eta_0$, $D_0$, $f_{0C}$, $q_{s0x}$ e le leggi di trasporto. Per estendere il nodo BRT a rami di larghezza diversa si impone $W_0=W_1+W_2$, riscalando le larghezze bankfull dei rami (§4.3), e si usa il flusso laterale d'acqua:

$$
Q_y = \frac{1}{2}\left(Q_1 - Q_2 - Q_0\,\frac{W_1 - W_2}{W_0}\right)
$$

**Q o q.** Le larghezze $W_0$, $W_1$, $W_2$ sono note dal DTM e diverse tra loro, quindi i bilanci al nodo (continuità liquida e solida) si scrivono in portata totale, $Q=qW$ e $Q_{s,x}=f_x\,q_{sx}\,W$. La portata per unità di larghezza serve solo localmente, per $\tau$, $\theta_g$ e la legge di Wilcock & Crowe. In `fsolve` conviene usare $Q_1$ come incognita e calcolare $q_1=Q_1/W_1$ dentro il residuo; le equazioni sotto restano scritte in $q$ per confronto con gli appunti.

$$
q_1W_1 + q_2W_2 = q_0W_0 \tag{1 — continuità liquida}
$$

$$
\eta_1 + D_1 = \eta_2 + D_2 \tag{2 — pelo libero ramo 1 = ramo 2}
$$

$$
\eta_2 + D_2 = \eta_0 + D_0 \tag{3 — pelo libero ramo 2 = madre}
$$

$$
f_{0F}\,q_{syF} = \frac{1}{2\alpha}\left(f_{1F}\,q_{s1F}\,\frac{W_1}{W_0} - f_{2F}\,q_{s2F}\,\frac{W_2}{W_0} - f_{0F}\,q_{s0F}\,\frac{W_1 - W_2}{W_1 + W_2}\right) \tag{4}
$$

$$
f_{0C}\,q_{syC} = \frac{1}{2\alpha}\left(f_{1C}\,q_{s1C}\,\frac{W_1}{W_0} - f_{2C}\,q_{s2C}\,\frac{W_2}{W_0} - f_{0C}\,q_{s0C}\,\frac{W_1 - W_2}{W_1 + W_2}\right) \tag{5}
$$

$$
f_{1F}\,q_{s1F}\,W_1 + f_{2F}\,q_{s2F}\,W_2 = f_{0F}\,q_{s0F}\,W_0 \tag{6 — continuità solida, F}
$$

$$
f_{1C}\,q_{s1C}\,W_1 + f_{2C}\,q_{s2C}\,W_2 = f_{0C}\,q_{s0C}\,W_0 \tag{7 — continuità solida, C}
$$

$$
\frac{\eta_1 + \eta_2}{2} = \eta_0 \tag{8 — quota del nodo}
$$

Le (4) e (5) sono i bilanci di cella per le classi F e C. Per $W_1=W_2=W_0/2$ le (4)–(5), combinate con (6)–(7), ricadono esattamente nelle eq. (14b–14c) di Ragno et al. L'equazione (8) non è nel paper: sostituisce un'informazione che il paper ottiene dalla pendenza dei rami e fissa la quota del nodo sulla media dei due figli.

### 5.3 Riduzione del sistema a 4 equazioni in [q₁, η₁, f₁F, f₂F]

Quattro equazioni si risolvono in forma esplicita, lasciando $q_1$, $\eta_1$, $f_{1F}$ e $f_{2F}$ come incognite:

$$
q_2 = \frac{q_0W_0 - q_1W_1}{W_2}\ \ (1), \qquad
\eta_2 = 2\eta_0 - \eta_1\ \ (8)
$$

$$
D_1 = \eta_0 + D_0 - \eta_1\ \ (2{+}3), \qquad
D_2 = D_0 - \eta_0 + \eta_1\ \ (3)
$$

**Perché 4 incognite e non 2.** Si eliminano per sostituzione solo le equazioni che non contengono trasporti solidi: (1), (8), (2)–(3). Le continuità solide (6)–(7) sono lineari in $f_{1F}$ e $f_{2F}$ solo se $q_{s1F}$, $q_{s1C}$, $q_{s2F}$, $q_{s2C}$ sono numeri noti; la versione precedente le risolveva così, ottenendo una formula esplicita di $f_{1F}$ e riducendo il sistema a $[q_1,\eta_1]$. I trasporti però dipendono dalla composizione dei rami, per due vie:

$$
q_{six} = q_{six}(D_i,\,f_{iF}): \qquad
d_{g,i} = 2^{-(f_{iC}\,\phi_C + f_{iF}\,\phi_F)}
\ \Longrightarrow\
\theta_{g,i} = \frac{\tau_i}{\rho g\Delta\,d_{g,i}}, \quad
g_{hr,x} = \left(\frac{d_x}{d_{g,i}}\right)^{-b}
$$

La formula esplicita era quindi circolare: $f_{1F}$ compariva anche a secondo membro, dentro i $q_s$. Valutando i $q_s$ con una composizione fissata (ad esempio quella della madre) si otteneva un $f_{1F}$ che non soddisfa (6)–(7) una volta ricalcolati i trasporti con quel valore.

Il secondo motivo è la singolarità. Il determinante del sistema lineare (6)–(7) è proporzionale al denominatore della vecchia formula:

$$
\det \propto q_{s1F}\,q_{s2C} - q_{s1C}\,q_{s2F} = 0
\quad\Longleftrightarrow\quad
\frac{q_{s1F}}{q_{s1C}} = \frac{q_{s2F}}{q_{s2C}}
$$

Nel regime di basso trasporto ($\zeta<1.35$, $G=0.002\,\zeta^{7.5}$) lo sforzo si semplifica nel rapporto tra le classi, che dipende solo dalla composizione:

$$
\frac{q_{sF}}{q_{sC}} = \left(\frac{g_{hr,F}}{g_{hr,C}}\right)^{7.5}
$$

Con la stessa composizione nei due rami il rapporto è identico, il determinante è nullo e $f_{1F}$ è indeterminato; vicino a questa condizione il sistema è mal condizionato e $f$ può uscire da $[0,1]$. Sull'Orco la classe grossolana è spesso in questo regime. Il significato fisico è che la differenza di composizione tra i rami (il sorting) è un risultato dell'accoppiamento tra bilanci di cella (4)–(5) e continuità (6)–(7), non un dato da cui partire.

| Equazione | Eliminata per sostituzione | Motivo |
| --- | --- | --- |
| (1) | sì: $q_2$ | non contiene trasporti |
| (8) | sì: $\eta_2$ | non contiene trasporti |
| (2)–(3) | sì: $D_1$, $D_2$ | non contiene trasporti |
| (6)–(7) | no: $f_{1F}$, $f_{2F}$ restano incognite | i $q_s$ dipendono da $f$; sistema singolare per $\zeta<1.35$ |
| (4)–(5) | no: $q_1$, $\eta_1$ restano incognite | equazioni residuo del nodo |

Restano 4 equazioni, (4)–(7), in 4 incognite $[q_1,\eta_1,f_{1F},f_{2F}]$, risolte insieme con `fsolve`. A ogni iterazione il residuo calcola $D_1$, $D_2$ da $\eta_1$, poi $d_g$ e $g_{hr}$ da $f_{1F}$, $f_{2F}$, poi i $q_{six}$, e valuta (4)–(7): trasporti e composizioni restano sempre coerenti.

Con $\eta_2=2\eta_0-\eta_1$ si ha $\eta_1-\eta_2=2(\eta_1-\eta_0)$, e le (4)–(5) diventano le due equazioni residuo seguenti, risolte con `fsolve` insieme alle continuità solide (6)–(7):

$$
\frac{f_{0F}\,q_{s0F}}{W_0}\left[\frac{2q_1W_1 - q_0W_0 - q_0W_1 + q_0W_2}{2\alpha\,q_0} - \frac{4r\,(\eta_1 - \eta_0)}{\sqrt{\theta_{g0}\,l_{hr,F}}}\right]
= \frac{1}{2\alpha W_0}\Big[2f_{1F}\,q_{s1F}\,W_1 - f_{0F}\,q_{s0F}\,W_0 - f_{0F}\,q_{s0F}\,(W_1 - W_2)\Big]
$$

$$
\frac{f_{0C}\,q_{s0C}}{W_0}\left[\frac{2q_1W_1 - q_0W_0 - q_0W_1 + q_0W_2}{2\alpha\,q_0} - \frac{4r\,(\eta_1 - \eta_0)}{\sqrt{\theta_{g0}\,l_{hr,C}}}\right]
= \frac{1}{2\alpha W_0}\Big[2q_{s1C}\,W_1 - 2f_{1F}\,q_{s1C}\,W_1 - f_{0C}\,q_{s0C}\,W_0 - f_{0C}\,q_{s0C}\,(W_1 - W_2)\Big]
$$

I trasporti $q_{s1x}$ e $q_{s2x}$ non sono costanti: dipendono da $D_1(\eta_1)$, $D_2(\eta_1)$ tramite $\tau=\rho g D S$, e da $f_{1F}$, $f_{2F}$ tramite $d_g$ e $g_{hr}$. A ogni iterazione di `fsolve` si ricalcolano quindi $D_i$, poi $q_{six}$, poi i residui di (4)–(7). Ottenuti $q_1$, $\eta_1$, $f_{1F}$ e $f_{2F}$, si ricavano $q_2$, $\eta_2$, $D_1$, $D_2$ e i flussi solidi dei figli, e si torna alle condizioni del volo con le regole di §4.3 ($\psi$ per portate liquide e solide, fattore $W_0^m/W_0^{bf}$ per le larghezze).

### 5.4 Confluenza

La confluenza somma i flussi in entrata e non richiede un sistema non lineare:

- Portata liquida: $Q_0=Q_1+Q_2$; il tirante $D_0$ si ricava da Manning invertita e $\eta_0=z_{DTM}-D_0$.
- Flussi solidi totali per classe: $Q_{s,x}=f_{1x}\,q_{s1x}\,W_1+f_{2x}\,q_{s2x}\,W_2$.
- Trasporto per unità di larghezza a valle: $q_{s0x}=Q_{s,x}/W_0$.
- Composizione della superficie a valle: la frazione $f_F$ che, con la capacità di trasporto della sezione di valle, trasporta fine e grossolano nel rapporto dell'apporto (equazione sotto). Diventa la condizione di monte per il tratto successivo e per l'eventuale biforcazione seguente.

$$
\frac{f\,q_{sF}(f)}{(1-f)\,q_{sC}(f)} = \frac{Q_{s,F}}{Q_{s,C}}
$$

Il rapporto dei flussi $Q_{s,x}/\sum Q_s$ descrive la composizione del materiale trasportato, non quella della superficie del letto, che è la grandezza usata da Wilcock-Crowe.

---

## 6. Sintesi per tipologia di sezione

La differenza chiave è nel nodo di biforcazione: il modulo idraulico ripartisce la portata con un bilancio di energia, quello morfodinamico con un bilancio di sedimento che rende $\eta_1$ un'incognita.

| Sezione | Modulo idraulico | Modulo morfodinamico |
| --- | --- | --- |
| Canale singolo | Moto uniforme (Manning) → $h$; $Q$ propagata invariata; $\eta=z_{DTM}-h$ | Moto uniforme → $D$; $\tau=\rho g D S$ → $q_{sC}$, $q_{sF}$ (Wilcock & Crowe); $Q$ e $Q_{s,x}$ propagati invariati |
| Biforcazione | 4 eq. in ($h_1$, $h_2$, $Q_1$, $Q_2$): Manning × 2, $Q_0=Q_1+Q_2$, Bernoulli; condizioni normali | 8 eq. in ($q_1$, $q_2$, $\eta_1$, $\eta_2$, $f_{1C}$, $f_{2C}$, $D_1$, $D_2$) ridotte a 4 in ($q_1$, $\eta_1$, $f_{1F}$, $f_{2F}$): nodo BRT con hiding, continuità liquida e solida, pelo libero comune, $\eta_0=(\eta_1+\eta_2)/2$; condizioni bankfull, poi riscalate; ritorno alle condizioni del volo con $\psi$ |
| Confluenza | $Q_0=Q_1+Q_2$ → Manning → $h_0$, $\eta_0$ | $Q_0=Q_1+Q_2$ → $D_0$, $\eta_0$; $Q_{s,x}$ sommati; $q_{s0x}=Q_{s,x}/W_0$; $f_F$ di equilibrio con la capacità di trasporto di valle |
| Ingresso | $Q$ di stazione, $W$, $S$, $n$; pendenze $S$ calcolate nel pre-processing | Uscita del modulo idraulico, $d_{50}$, $\phi_C$, $\phi_F$, $f_{0C}$, $\Delta$, $r$, $\alpha$, $\theta_r$; pendenze $S$ calcolate nel pre-processing |
| Uscita | $\eta$ su ogni sezione → DTM corretto | $\eta$ su ogni sezione, ripartizione dei flussi solidi, composizione $f_x$ dei rami → DTM corretto |

---

## 7. Questioni aperte

*Questa sezione è il quadro dei problemi aperti a ottobre 2026, prima della revisione del modulo morfodinamico (§9). L'elenco aggiornato dei punti in sospeso è la §9.4.*

Questa sezione raccoglie i problemi emersi applicando il modulo morfodinamico (`2_ORCO_bifurcationsNEW_modifiche.mlx`) ai dati dell'Orco del volo LiDAR del 17/03/2025 ($Q=10.8$ m³/s a San Benigno). Il codice è coerente con le sezioni 4–5 e conserva portata liquida e flussi solidi a tutti i nodi, ma **nessuna delle 12 biforcazioni viene risolta con il sistema nodale**: in tutte interviene il ripiego (ripartizione proporzionale alle larghezze, §7.6). Le cause sono nei dati in ingresso al nodo, non nelle equazioni: per ciascuna si descrive cosa succede, dove, perché, e le scelte di metodo possibili.

### 7.1 Stato dei nodi

Valori stampati dalla tabella diagnostica del modulo morfodinamico (sezione 5 dello script). $D_{0}$, $D_{1,s}$, $D_{2,s}$ sono i tiranti di Shields $\theta^*_{BF}\Delta d_{50}/S$ della madre e dei rami; $W_0^m$ è la larghezza della maschera, $W_0^{bf}$ quella bankfull dal DTM; `exit` è il codice di uscita di `fsolve`.

| Bif | $S_0$ | $S_1$ | $S_2$ | $D_0$ (m) | $D_{1,s}$ (m) | $D_{2,s}$ (m) | $W_0^m$ (m) | $W_0^{bf}$ (m) | exit | Causa principale |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 15 | 1.6e-4 | 6.7e-4 | 2.8e-3 | 31.2 | 7.5 | 1.8 | 20 | 54.8 | −3 | $S_0$ al minimo; $W_2^{bf}$ non calcolabile |
| 106 | 3.9e-3 | 8.8e-3 | 2.7e-3 | 1.3 | 0.6 | 1.9 | 26 | 56.3 | 3 | converge, ma $q_2<0$: $S_1/S_2=3.3$ |
| 205 | 5.6e-4 | 4.7e-3 | 1.6e-4 | 8.9 | 1.1 | 31.2 | 36.5 | 83.4 | −3 | $S_0$ bassa, $S_2$ al minimo |
| 244 | 1.6e-4 | 1.1e-2 | 6.8e-3 | 31.2 | 0.4 | 0.7 | 44 | 127.8 | 0 | $S_0$ al minimo; $S_1/S_0=70$ |
| 366 | 2.8e-3 | 1.6e-4 | 5.4e-3 | 1.8 | 31.2 | 0.9 | 41.5 | 137.7 | 0 | $S_1$ al minimo; clamp con $D_0$ plausibile |
| 469 | 5.4e-3 | 8.8e-3 | 1.2e-2 | 0.9 | 0.6 | 0.4 | 13.5 | 29.6 | −3 | pendenze dei rami 1.6–2.1 volte $S_0$ |
| 584 | 3.2e-3 | 9.3e-3 | 8.8e-3 | 1.6 | 0.5 | 0.6 | 62 | 68.6 | 0 | rami 2.8–2.9 volte $S_0$; clamp |
| 645 | 1.5e-3 | 2.1e-3 | 8.7e-3 | 3.3 | 2.4 | 0.6 | 32.5 | 106.6 | −3 | $S_2/S_0=5.7$; clamp |
| 817 | 3.1e-3 | 1.6e-4 | 3.5e-3 | 1.6 | 31.2 | 1.4 | 31 | 92.2 | −3 | $S_1$ al minimo; $W_2^{bf}$ non calcolabile |
| 1070 | 3.5e-3 | 1.1e-2 | 5.0e-3 | 1.4 | 0.4 | 1.0 | 11.5 | 74.3 | 0 | $S_1/S_0=3.2$ |
| 1293 | 3.9e-4 | 1.8e-3 | 9.1e-3 | 12.7 | 2.8 | 0.5 | 41.5 | – | – | $S_0$ bassa; $W_0^{bf}$ non calcolabile |
| 1369 | 1.6e-4 | 2.3e-3 | 3.3e-3 | 31.2 | 2.1 | 1.5 | 36.5 | 74.6 | −3 | $S_0$ al minimo |

Quello che funziona: la portata è conservata a tutti i nodi e alle confluenze, nessuna sezione resta senza quota di fondo, e la differenza media rispetto al modulo idraulico è di 0.045 m (massimo 0.73 m). A portata del volo i flussi solidi sono dell'ordine di $10^{-16}$ m³/s, cioè praticamente nulli: è coerente con §4.3 (a portate ordinarie il grossolano è immobile) e non è un problema, perché il trasporto conta solo nel calcolo bankfull dei nodi.

### 7.2 Pendenze dei tratti non rappresentative ai nodi

**Cosa succede.** In 7 nodi su 12 almeno una delle tre pendenze (madre o rami) è al valore minimo $1.6\cdot10^{-4}$ (madre: 15, 244, 1369; rami: 205, 366, 817), e in 1293 la madre ha $S_0=3.9\cdot10^{-4}$. Dove le pendenze non sono al minimo, quelle di madre e rami differiscono comunque di un fattore 2–70 (244, 584, 645, 1070, 106).

**Dove.** Nel pre-processing (§2.2): la pendenza di ogni tratto è la retta ai minimi quadrati $z(s)$ sull'intero tratto, limitata a $1.6\cdot10^{-4}$. Il valore minimo compare quando la retta ha pendenza nulla o negativa (tratti corti, quote disturbate, rigurgiti a ponti o confluenze) o quando il tratto ha un solo punto.

**Perché crea il problema.** Le pendenze entrano nel nodo in due modi.

1. Il tirante bankfull è inversamente proporzionale alla pendenza:
$$
D = \frac{\theta^*_{BF}\,\Delta\,d_{50}}{S} \approx \frac{0.005\ \text{m}}{S}
\quad\Longrightarrow\quad
S=1.6\cdot10^{-4}\ \Rightarrow\ D\approx31\ \text{m}
$$
Con un tirante di 9–31 m il livello bankfull supera il terreno (§7.3) e la geometria del nodo non ha più significato.

2. Il nodo BRT è molto sensibile al rapporto tra le pendenze dei rami. Con il pelo libero comune i tiranti dei rami sono vicini, quindi lo sforzo $\tau=\rho g D S$ scala circa con $S$. Nel regime di basso trasporto di Wilcock-Crowe ($\zeta<1.35$):
$$
q_s \propto \tau^{3/2}\,\zeta^{7.5} \propto \tau^{9}
\quad\Longrightarrow\quad
\frac{q_{s1}}{q_{s2}} \sim \left(\frac{S_1}{S_2}\right)^{9}
$$
Una differenza del 10% tra $S_1$ e $S_2$ cambia il rapporto dei trasporti di un fattore 2.4; un rapporto di 3 lo cambia di un fattore $\sim 2\cdot10^4$. Con rapporti di pendenza di 3–70 le continuità solide (6)–(7) e i bilanci di cella (4)–(5) non hanno una soluzione con entrambi i rami attivi: `fsolve` non converge (exit 0 o −3, residuo $10^{0}$–$10^{5}$) oppure converge con un ramo chiuso ($q_2<0$ in 106).

Nella realtà, a poche decine di metri dal nodo, la pendenza del pelo libero di madre e rami non può differire di così tanto: le differenze vengono dalla regressione su tratti di lunghezza molto diversa.

**Scelte di metodo.**

- **A. Pendenza del nodo su una finestra di lunghezza fissa**. Per la madre la pendenza si calcola sul pelo libero $z_{DTM}$ di un tratto di lunghezza $L$ a monte del nodo, per ciascun ramo su un tratto di lunghezza $L$ a valle, attraversando se necessario i confini tra tratti (la topologia è nota). $L$ fissa per tutta la rete, ad esempio 300–500 m, oppure proporzionale alla larghezza, $L=k\,W_0^m$ con $k\approx10$. Le pendenze restano "dal DTM" (§2.2) ma sono stimate su una scala confrontabile per madre e rami. Il calcolo nel resto della rete non cambia.
- **B. Sostituire il valore minimo.** Quando la regressione dà $S\le1.6\cdot10^{-4}$ o il tratto è troppo corto, si usa la pendenza del tratto a monte (o a valle) o quella della finestra A, invece di $1.6\cdot10^{-4}$. Da applicare comunque, anche con A.
- **C. Modificare il valore di ripiego $1.6\cdot10^{-4}$.** Il minimo imposto nel pre-processing quando la regressione dà pendenza nulla o negativa si sostituisce con un valore scelto sui dati dell'Orco, ad esempio un percentile basso (5–10%) delle pendenze valide della rete o la pendenza media del corso d'acqua. Poiché $D\approx0.005\ \text{m}/S$, il valore di ripiego fissa direttamente il tirante massimo: con $S=10^{-3}$ si ottiene $D\approx5$ m invece di 31 m. È la modifica più semplice, ma resta un valore unico per tutta la rete e non corregge la differenza di pendenza tra madre e rami.

### 7.3 Tirante bankfull di Shields superiore all'altezza delle sponde

**Cosa succede.** In 13 sezioni il livello bankfull $\eta+D$ non è contenuto nel transetto di ±100 m (avviso *clamp used*): madri 15, 205, 244, 366, 584, 645, 1293, 1369 e rami 235, 377, 446, 646, 818. In quasi tutte una sola sponda viene trovata e l'altra è specchiata (*mirrored*), a 40–100 m dalla centerline. Le larghezze bankfull risultano 2–6 volte quelle della maschera (366: 137.7 m contro 41.5 m; 1070: 74.3 m contro 11.5 m; 645: 106.6 m contro 32.5 m).

**Dove.** In `W_bankfull_from_dtm` (§4.3): il livello $z_{target}=\eta+D$ supera il massimo del profilo del DTM almeno su un lato, quindi viene abbassato al massimo del profilo e le sponde si cercano a quella quota.

**Perché.** Due cause, che si sommano:

1. pendenze piccole (§7.2), che danno $D$ di 9–31 m;
2. anche dove $D$ è plausibile (1.4–1.8 m in 366, 584, 1070), le sponde reali sono più basse di $D$ rispetto al fondo ricostruito. Il criterio $\theta^*_{BF}=1.62\,\theta^*_c$ con $d_{50}=62$ mm è tarato su fiumi ghiaiosi a canale singolo con sponde ben definite; nei tratti a canali multipli dell'Orco le sponde dei singoli rami sono basse e a quel livello l'acqua allaga la piana.

**Scelte di metodo.**

- **A. Bankfull geomorfologico dal DTM** *(proposta principale)*. Sul transetto si cerca il livello di sfioro $z_{sfioro}$, cioè la quota più bassa a cui l'acqua esce dall'alveo: il minimo tra i massimi del profilo a sinistra e a destra della centerline (entro una distanza massima). Il tirante bankfull diventa
$$
D_{bf} = \min\left(D_{Shields},\; z_{sfioro}-\eta\right)
$$
Shields resta il riferimento teorico, la geometria lo limita. Con $D_{bf}$ si calcolano $W^{bf}$, $Q^{bf}$ e lo Shields bankfull effettivo $\theta_{bf}=D_{bf}S/(\Delta d_{50})$, da riportare per controllo.
- **B. Calibrare $\theta^*_{BF}$ sull'Orco.** Si stima $\theta^*_{BF}$ nelle sezioni a canale singolo con sponde ben definite (dove il livello di sfioro è chiaro) e si usa quel valore in tutta la rete al posto di $1.62\,\theta^*_c$. Mantiene un criterio unico, ma richiede di scegliere le sezioni di taratura.

### 7.4 Larghezze bankfull non calcolabili

**Cosa succede.** In alcune sezioni `W_bankfull_from_dtm` non restituisce una larghezza:

- 41 e 1293: la sponda specchiata cade fuori dai dati del DTM (*NaN values while interpolating*);
- 889: nessun incrocio oltre la semi-larghezza della maschera (*no crossing satisfies abs(s) ≥ halfW*);
- 1294: area negativa e nessuna sezione precedente nello stesso tratto.

**Effetto.** Nei nodi 15 e 817 manca la larghezza di un ramo e si usa $W_1=W_2=W_0/2$; nel nodo 1293 manca $W_0^{bf}$ e il sistema non viene impostato.

**Perché.** Sono conseguenze dirette di §7.3 (livello troppo alto, sponde specchiate lontane) e dei bordi del DTM.

**Scelte di metodo.**

- **A.** Se la larghezza bankfull non è calcolabile si usa la larghezza della maschera, $W^{bf}=W^m$: è coerente con la regola $W^{bf}\ge W^m$ di §4.3 e non introduce valori arbitrari.
- **B.** Quando una sponda è specchiata, l'area si calcola sul solo semi-transetto in cui la sponda è stata trovata e si raddoppia, invece di interpolare il profilo oltre i dati: è coerente con l'ipotesi di simmetria che sta dietro lo specchiamento.

Molti di questi casi dovrebbero sparire risolvendo §7.2 e §7.3.

### 7.5 Robustezza numerica del nodo

**Cosa succede.** Anche nei nodi con pendenze e tiranti plausibili (469, 584, 1070) `fsolve` termina con exit −3 (non riesce più a ridurre il residuo) o 0 (numero massimo di iterazioni), con residui di 10–80.

**Perché.** Con le pendenze attuali il sistema probabilmente non ha soluzione con entrambi i rami attivi (§7.2). In più il punto di partenza è sempre la biforcazione bilanciata ($q_1=q_0$, $\eta_1=\eta_0$, $f_{1F}=f_{2F}=f_{0F}$), che può essere lontano dalla soluzione quando le pendenze sono molto diverse.

**Scelte di metodo.**

- **A.** Rivalutare dopo aver risolto §7.2 e §7.3: con pendenze e tiranti coerenti il problema potrebbe sparire.

### 7.6 Ripiego attuale

Quando il nodo non viene risolto, il codice:

1. ripartisce portata liquida e flussi solidi in proporzione alle larghezze bankfull dei rami, o, se mancano, a quelle delle maschere ($\psi=W_1/W_0$);
2. assegna ai rami la composizione della madre;
3. calcola il fondo all'imbocco dei rami con il moto uniforme, come nel canale singolo.

In questi nodi il risultato è quindi idraulico (moto uniforme con ripartizione geometrica), non morfodinamico. Va dichiarato nei risultati finché le questioni sopra non sono risolte: allo stato attuale le quote ai nodi del modulo morfodinamico non contengono l'effetto del bilancio di sedimento.

### 7.7 Decisioni richieste

| # | Decisione | Proposta |
| --- | --- | --- |
| 1 | Pendenza ai nodi (§7.2) | Modificare il valore di ripiego $1.6\cdot10^{-4}$ (C): sostituirlo con un valore scelto sui dati dell'Orco (percentile basso, 5–10%, delle pendenze valide della rete, o pendenza media del corso d'acqua). Il valore fissa il tirante massimo ($S=10^{-3}\Rightarrow D\approx5$ m invece di 31 m). È la modifica più semplice, ma resta un valore unico per tutta la rete e non corregge la differenza di pendenza tra madre e rami |
| 2 | Tirante bankfull (§7.3) | Da decidere |
| 3 | Larghezze non calcolabili (§7.4) | Da decidere |
| 4 | Ramo che si chiude a bankfull (§7.1, nodo 106) | Da decidere |

---

## 9. Aggiornamento del modulo morfodinamico (novembre 2026)

Questa sezione raccoglie le modifiche introdotte nel codice il 12/10/2026 (`2_ORCO_bifurcationsNEW_modifiche_modifiche.mlx`) e il 18/11/2026 (`ORCO_bifurcations.m`) e i risultati dei run del 18/11/2026. Dove è in contrasto con le sezioni precedenti (in particolare §4.3, §5.2–5.3, §6 e §7), vale quanto scritto qui; le sezioni precedenti saranno aggiornate in seguito.

### 9.1 Modifiche al metodo

**Pendenza della madre al nodo.** $S_0$ si ricalcola con una retta ai minimi quadrati sul pelo libero $z_{DTM}$ entro $L=400$ m a monte del nodo, proseguendo nel tratto a monte se necessario (a una confluenza, nel ramo con la portata maggiore). Se la retta non è utilizzabile (meno di 5 punti, finestra più corta di 100 m, pendenza nulla o negativa) la finestra si raddoppia fino a 1600 m; se non basta si usa il 25° percentile delle pendenze dei tratti. Una pendenza valida si alza solo se è sotto il 5° percentile. Le pendenze dei rami non sono più dati del nodo.

**Tirante bankfull.** Il tirante di Shields si limita con il livello di sfioro del transetto e non scende mai sotto il tirante critico:

$$
D_0 = \max\Big(\min\big(D_S,\ z_{sfioro}-\eta_0\big),\ D_c\Big), \qquad
D_S = \frac{1.62\,\theta^*_c\,\Delta\,d_{50}}{S_0}, \qquad
D_c = \frac{\theta^*_c\,\Delta\,d_{50}}{S_0}
$$

$z_{sfioro}$ è il più basso dei due punti più alti del terreno, uno per lato, tra $W^m/2$ e $3W^m$ dalla centerline. Dove il tirante si ferma a $D_c$ (sponde basse) la larghezza resta quella misurata al livello di sfioro. Le larghezze dei rami si misurano con il tirante di Shields calcolato con $S_0$; se non sono calcolabili si usa il rapporto delle larghezze delle maschere.

**Calcolo della larghezza bankfull al nodo (verificato sul codice il 19/11/2026, `W_bankfull_from_dtm.m` e `MorphodynamicModule.m`).** La larghezza non è legata al tirante critico: si misura al livello bankfull di Shields, limitato dalle sponde.

1. Livello di partenza: $\eta_0+D_S$, con $\eta_0$ fondo del modulo idraulico e $D_S=1.62\,\theta^*_c\Delta d_{50}/S_0$.
2. Si traccia il transetto perpendicolare alla centerline sul DTM corretto dal modulo idraulico e si cerca il livello di sfioro $z_{sfioro}$ tra $W^m/2$ e $3W^m$ dalla centerline. Se il livello di Shields supera $z_{sfioro}$, il livello diventa $z_{sfioro}$ − 1 cm.
3. Le sponde sono le prime intersezioni del livello con il terreno oltre $W^m/2$, una per lato; se si trova un solo lato, l'altro si ottiene per simmetria.
4. L'area $A$ tra le sponde sotto il livello dà la larghezza rettangolare equivalente $W^{bf}=A/D$, con $D$ profondità del livello sul fondo; si impone $W^{bf}\ge W^m$.
5. Se la profondità al livello di sfioro è minore del tirante critico $D_c$, il tirante del nodo diventa $D_0=D_c$, ma larghezza e area restano quelle misurate al livello di sfioro (sponde basse, $\theta_{bf}/\theta^*_c=1$).
6. Per i rami si ripete lo stesso calcolo alla loro prima sezione, con lo stesso $D_S$ (pendenza della madre); le due larghezze si riscalano poi in modo che $W_1+W_2=W_0$.

**Sistema del nodo in forma BRT (6 equazioni).** Come in Bolla Pittaluga et al. (2003) e Ragno et al. (2023), il pelo libero è comune al nodo e i due rami, di uguale lunghezza, terminano allo stesso livello: hanno quindi la stessa pendenza di equilibrio $S_b$, che è un'incognita. Incognite $[\psi, D_1, D_2, S_b, f_{1F}, f_{2F}]$ con $\psi=Q_1/Q_0$; equazioni:

- moto uniforme in ciascun ramo, $\psi Q_0 = W_1 D_1^{5/3}\sqrt{S_b}/n$ e $(1-\psi) Q_0 = W_2 D_2^{5/3}\sqrt{S_b}/n$;
- bilanci di cella per F e C (§4.4, eq. 4–5 di §5.2), con $\eta_1-\eta_2=D_2-D_1$;
- continuità solida per F e C (eq. 6–7 di §5.2).

**Cosa si intende per nodo in forma BRT.** È il modello a due celle di Bolla Pittaluga, Repetto e Tubino (2003, "BRT"), esteso alle miscele da Ragno et al. (2023) e qui a rami di larghezza diversa. Ipotesi:

1. **equilibrio:** in ciascun ramo il trasporto di ogni classe è pari alla capacità di trasporto, e i flussi solidi si conservano al nodo;
2. **condizioni bankfull:** la configurazione del nodo si forma a rive piene (§4.3);
3. **moto uniforme in ciascun ramo**, sezione rettangolare larga, Manning con $n$ costante;
4. **pelo libero comune al nodo:** $\eta_0+D_0=\eta_1+D_1=\eta_2+D_2$, quindi $\eta_1-\eta_2=D_2-D_1$ e $\eta_i=\eta_0+D_0-D_i$;
5. **rami di uguale lunghezza che terminano allo stesso livello:** i due rami hanno la stessa pendenza di equilibrio $S_b$, che è un'incognita;
6. **due celle di lunghezza $\alpha W_0$ a monte del nodo** ($\alpha=4$): il flusso laterale di sedimento tra le celle è dovuto alla deviazione della corrente e al richiamo gravitativo lungo la pendenza trasversale (coefficiente di Ikeda $r=0.65$, "grain independence", §4.4);
7. **larghezze note**, con $W_1+W_2=W_0$; la soluzione bilanciata di riferimento è la ripartizione proporzionale alle larghezze, $\psi=W_1/W_0$.

Equazioni, nelle 6 incognite $\psi$, $D_1$, $D_2$, $S_b$, $f_{1F}$, $f_{2F}$ (con $Q_1=\psi Q_0$, $Q_2=(1-\psi)Q_0$ e $Q_{si,x}=f_{ix}\,q_{six}(D_i,S_b,f_{iF})\,W_i$):

$$
\psi\,Q_0 = \frac{W_1\,D_1^{5/3}\sqrt{S_b}}{n}, \qquad
(1-\psi)\,Q_0 = \frac{W_2\,D_2^{5/3}\sqrt{S_b}}{n} \tag{1–2}
$$

$$
f_{0x}\,q_{sy,x} = \frac{1}{2\alpha W_0}\left(Q_{s1,x} - Q_{s2,x} - Q_{s0,x}\,\frac{W_1 - W_2}{W_1 + W_2}\right), \qquad x = F, C \tag{3–4}
$$

$$
q_{sy,x} = q_{s0,x}\left[\frac{Q_y}{\alpha W_0\,q_0} - \frac{r}{\sqrt{\theta_{g0}\,l_{hr,x}}}\,\frac{2\,(D_2 - D_1)}{W_0}\right], \qquad
Q_y = \frac{1}{2}\left(Q_1 - Q_2 - Q_0\,\frac{W_1 - W_2}{W_0}\right)
$$

$$
Q_{s1,x} + Q_{s2,x} = Q_{s0,x}, \qquad x = F, C \tag{5–6}
$$

Le (1)–(2) sono il moto uniforme nei rami, le (3)–(4) i bilanci di cella per le classi fine e grossolana, le (5)–(6) la continuità solida. La continuità liquida è implicita nell'uso di $\psi$.

L'equazione (8), $\eta_0=(\eta_1+\eta_2)/2$, e le pendenze imposte ai rami non si usano più. Motivo: il sistema di §5.2 non conteneva il moto uniforme nei rami. Con le pendenze misurate (rami 1.5–5 volte più ripidi della madre) non aveva soluzione; con la stessa pendenza nei tre rami dava una seconda soluzione in cui un ramo portava fino a 2.5 volte la portata compatibile con il suo tirante, e con rami di uguale larghezza non dava mai le soluzioni asimmetriche di BRT. La forma a 6 equazioni le ritrova: con rami uguali, sotto un rapporto di forma critico $\beta=W_0/(2D_0)\approx13$ esiste solo la soluzione bilanciata, sopra compaiono due soluzioni asimmetriche speculari.

**Scelta della soluzione.** `fsolve` parte da 55 configurazioni e trova tutte le soluzioni fisiche. La soluzione bilanciata ($\psi=W_1/W_0$, $D_1=D_2=D_0$, $S_b=S_0$) esiste sempre: se è l'unica si adotta; se esistono soluzioni asimmetriche la bilanciata è instabile (BRT) e si adotta l'asimmetrica più vicina a $\psi_{hyd}$, la quota del ramo 1 nel modulo idraulico al volo ($Q_i^{hyd}\propto W_i^m h_i^{5/3}\sqrt{S_i}$ alla prima sezione dei rami).

**Ritorno al volo.** $\psi$ e $\psi_x$ ripartiscono portata e flussi solidi del volo come in §4.3; la prima sezione di ogni ramo si calcola con il moto uniforme alla sua portata, come ogni altra sezione. Le quote bankfull $\eta_i=\eta_0+D_0-D_i$ restano una diagnostica: usate direttamente al volo davano fondi sopra il pelo libero in 5–7 nodi.

### 9.2 Risultati del run di riferimento (larghezze bankfull)

Volo del 17/03/2025, `node_width_source = 'dtm'`: al nodo tutte le larghezze sono quelle bankfull calcolate sul DTM corretto (§4.3), cioè l'area $A$ tra le sponde al livello bankfull divisa per il tirante bankfull, $W^{bf}=A/D_0$. È il run di riferimento: tirante, larghezza, portata e trasporto del nodo sono tutti in condizioni bankfull. Tutti i 12 nodi sono risolti dal sistema (residuo $\le10^{-15}$) e in tutti esistono la soluzione bilanciata e due asimmetriche, quindi tutti risultano asimmetrici. La pendenza di equilibrio dei rami è 0.61–0.99 volte $S_0$; in tutti i nodi il ramo con il fondo più basso ha la superficie più grossolana ($f_F$ minore), come atteso dal richiamo gravitativo lungo la pendenza trasversale. I dislivelli bankfull tra i rami arrivano a 1.4 m (584, 645). In 9 madri su 12 il tirante bankfull è quello critico ($\theta_{bf}=\theta_c$). Portata e flussi solidi si conservano a tutti i nodi; nessuna sezione resta senza quota.

### 9.3 Quanto il risultato dipende dalla larghezza (test di sensitività)

Per misurare la sensibilità del nodo alla larghezza, il run è stato ripetuto sostituendo le larghezze bankfull con quelle delle maschere del volo (`node_width_source = 'mask'`), lasciando invariati pendenze, tiranti bankfull, equazioni e regola di scelta. Questo secondo run **non è un'alternativa al metodo**: combina un tirante bankfull con una larghezza misurata a portata ordinaria (10.8 m³/s), quindi la madre non è in condizioni bankfull e il trasporto solido al nodo non è giustificato. Serve solo a capire quanto il risultato dipende dalla larghezza; il risultato del metodo è quello del run con le larghezze bankfull (§9.2).

| Bif | $\beta$ bankfull | $\beta$ maschera | soluzione bankfull | soluzione maschera | $\psi$ bankfull | $\psi$ maschera | $\psi_{hyd}$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 15 | 31.1 | 8.0 | asimmetrica | bilanciata | 0.41 | 0.64 | 0.49 |
| 106 | 23.0 | 12.3 | asimmetrica | asimmetrica | 0.35 | 0.56 | 0.54 |
| 205 | 13.9 | 7.1 | asimmetrica | bilanciata | 0.47 | 0.40 | 0.74 |
| 244 | 31.2 | 11.3 | asimmetrica | asimmetrica | 0.80 | 0.72 | 0.76 |
| 366 | 30.4 | 10.8 | asimmetrica | bilanciata | 0.25 | 0.43 | 0.15 |
| 469 | 15.4 | 6.5 | asimmetrica | bilanciata | 0.72 | 0.68 | 0.67 |
| 584 | 33.6 | 24.0 | asimmetrica | asimmetrica | 0.51 | 0.75 | 0.66 |
| 645 | 28.3 | 8.0 | asimmetrica | bilanciata | 0.46 | 0.45 | 0.35 |
| 817 | 35.1 | 12.1 | asimmetrica | bilanciata | 0.17 | 0.26 | 0.08 |
| 1070 | 14.4 | 9.0 | asimmetrica | bilanciata | 0.76 | 0.81 | 0.84 |
| 1293 | 18.6 | 11.1 | asimmetrica | bilanciata | 0.69 | 0.76 | 0.67 |
| 1369 | 12.0 | 4.6 | asimmetrica | bilanciata | 0.26 | 0.30 | 0.27 |

| Indicatore | larghezze bankfull (riferimento) | larghezze della maschera (sensitività) |
| --- | --- | --- |
| $\beta$ | 12–35 | 4.6–24 |
| nodi asimmetrici | 12 su 12 | 3 su 12 (106, 244, 584) |
| scostamento massimo dalla soluzione bilanciata, $\lvert\psi-W_1/W_0\rvert$ | 0.39 | 0.10 |
| dislivello massimo tra i rami, $\lvert\eta_i-\eta_0\rvert$ | 1.40 m | 0.21 m |
| fondo rispetto al modulo idraulico (media / massimo) | −0.006 / 0.30 m | 0.023 / 0.58 m |

**Conclusioni.**

1. **Il nodo è molto sensibile alla larghezza.** Il rapporto di forma critico ($\approx12$–13) cade tra i valori di $\beta$ dei due run: con le larghezze bankfull tutti i nodi sono sopra la soglia e asimmetrici, con larghezze più piccole (quelle della maschera, 1.5–4 volte minori) 9 nodi su 12 scenderebbero sotto la soglia e diventerebbero bilanciati. Il tipo di soluzione dipende quindi in modo critico da come si misura la larghezza bankfull.
2. **L'effetto sulla ripartizione è rilevante.** La quota di portata del ramo 1 differisce tra i due run in media di 0.11 e fino a 0.24 (584), 0.23 (15) e 0.21 (106); l'effetto sul fondo al volo è piccolo in media (≤ 0.02 m) ma arriva a 0.3–0.6 m localmente.
3. **Conseguenza per il metodo.** Il risultato da usare è quello con le larghezze bankfull, l'unico coerente con l'ipotesi che la morfologia dei nodi si formi in condizioni bankfull. La sensibilità mostra però che la stima della larghezza bankfull va curata: al livello di sfioro la larghezza dal DTM comprende anche barre e canali adiacenti, tipici dei tratti a canali multipli, e può sovrastimare la larghezza del canale che alimenta il nodo. Una via per renderla più robusta è cercare le sponde solo nel canale attivo che alimenta il nodo, senza estendere la ricerca fino a $3W^m$. Nel paper il test va riportato come analisi di sensitività sulla larghezza bankfull.

### 9.4 Punti in sospeso

Questa tabella è l'elenco aggiornato di tutti i punti ancora aperti dei due moduli; la §7 resta come quadro della versione di ottobre. Quando un punto si chiude, il suo stato diventa "chiuso" con il riferimento alla sezione che lo documenta.

| # | Punto | Stato | Prossimo passo |
| --- | --- | --- | --- |
| 1 | **Stima della larghezza bankfull**, che controlla il rapporto di forma $\beta$ e quindi il tipo di soluzione al nodo (§9.3) | aperto; confermato dal run del 18/11 sera: a 1369 $\beta$ passa da 12.0 a 11.7 e il nodo diventa bilanciato (§9.8). Nel run del 19/11, con $z_{DTM}$ corretta, 1369 torna asimmetrico ($\beta$ = 14.9) e $\beta$ cambia di 2–5 in altri quattro nodi (§9.10). Run del 20/11: alla sezione 256 (ramo del nodo 244) una sola sponda viene trovata e specchiata, e $W_2^{dtm}$ passa da 78.3 a 88.1 m senza cambiare la soluzione: esempio della fragilità della stima | sponde cercate solo nel canale attivo che alimenta il nodo; da mettere nell'analisi di sensitività |
| 2 | **Sponde basse:** in 9 madri su 12 il tirante bankfull è quello critico (5 su 10 nel run del 18/11 sera); $\theta^*_{BF}=1.62\,\theta^*_c$ (fonte da documentare) | aperto | eventuale calibrazione sull'Orco |
| 3 | **Portata bankfull calcolata localmente** a ogni nodo, non coerente lungo la rete (per esempio 584: 216 m³/s, il nodo 645 su un suo ramo: 481 m³/s) | aperto | — |
| 4 | **Regola di scelta della soluzione:** nei nodi 205, 584 e 645 nessuna soluzione asimmetrica è vicina a $\psi_{hyd}$ | aperto (invariato nel run del 18/11 sera). Run del 19/11: il nodo 15 salta da $\psi$ = 0.41 a 0.82 solo perché $\psi_{hyd}$ passa da 0.49 a 0.67 (§9.10) | domanda 4 a Carlo; regola di scelta e larghezza bankfull nell'analisi di sensitività |
| 5 | **Dati ADCP** (8 transetti del 07/03/2025) | integrati nel codice il 19/11/2026 come stazioni a01–a08, senza distinzione dalle sezioni GNSS (§9.9); run del 19/11: residuo 0.3–0.6 m più alto delle GNSS vicine perché i transetti coprono soprattutto le parti profonde, non guadabili; errore medio su tutte le stazioni da +0.22 a +0.33 m (§9.9, §9.10) | decidere come entrano nella validazione (con la media di stazione insieme alle GNSS, solo separati, o solo dove completano una sezione GNSS) e con quale peso (domanda 6 a Carlo, `adcp_min_spacing` = 0) |
| 6 | **Portata non conservata nel modulo idraulico:** il tratto 6 riceveva la portata del tratto 5 e da lì la rete portava 9.99 invece di 10.80 m³/s (§9.5) | **chiuso** (§9.8): bilancio chiuso a tutti i nodi in entrambi i moduli (sbilancio massimo 4.6·10⁻¹² m³/s), 10.80 m³/s all'uscita | — |
| 7 | **Rami 16–15 e 29** (§9.7): nel modello tolgono portata alla rete; al volo sono quasi fermi (29 staccato dal tratto 22, 16–15 parte di un meandro bagnato solo in parte), mentre a portate più alte sono attivi | **chiuso per le condizioni del volo** (§9.8): 366 e 817 "branch without discharge" in entrambi i moduli; 127 559 pixel bagnati (6.4%) restano al DTM grezzo; nessuna stazione GNSS coinvolta | resta aperta, come estensione, la rete alla portata bankfull (limiti del paper) |
| 8 | **Quota di partenza $z_{DTM}$** (§9.6): oscilla e in alcuni tratti è fino a 1.6 m più bassa del DTM completo; causa i forti residui negativi di s006, s019, s020 | causa confermata alla sezione 170 (§9.8): nodata dei DTM di segmento, propagato dal filtro sui salti. **Corretto nel codice il 19/11/2026** (§9.6, "Correzione adottata"): nodata −32767 → NaN alla lettura dei DTM di segmento, media della sezione sui soli valori non NaN. **Run del 19/11** (§9.10): salti segnalati da 149 a 2, tratti di `segments_to_smooth` senza oscillazioni, $z_{DTM}$ entro 4 cm dal DTM completo a s006, s019, s020 **20/11:** `segments_to_smooth` e filtro a 1 m sostituiti da una regola unica per i salti isolati (`glitch_tol` = 0.20 m, §9.6) **Deciso il 20/11:** questione chiusa con la regola unica per i salti isolati **Run del 20/11 (§9.13):** la regola corregge esattamente le 4 sezioni previste; caso limite alla 492, piede di un gradino portato a metà **20/11 (decisione):** regola resa più stretta (picco solo se si stacca da entrambe le vicine); sui dati del 20/11 corregge solo 728 e 1239 **Run del 20/11 sera (§9.14):** corrette solo le ex 728 e 1239; la 492 resta al valore grezzo. **Chiuso il 20/11** | — |
| 9 | **Stazioni GNSS su più tratti** s003, s007, s015, s035 | s003 (tratti 9, 11, 12) e s007 (tratti 6, 7) divise in sottostazioni assegnate a mano il 19/11/2026 (§9.9); nel run del 19/11 assegnate correttamente; nessuno dei due moduli riproduce il dislivello tra i rami (§9.10); s015 assegnata a mano al tratto madre 30 (sta al nodo 1070) il 20/11. Run del 20/11: s015-t30 ha la sezione di riferimento al nodo 1070, con il 63% dei punti sui tratti 31 e 32; resta ambigua solo s035. Deciso il 20/11: s015 resta nella validazione ed è esclusa dal confronto tirante–profondità (§9.11) | s035 resta sul tratto 20 (deciso il 20/11): punto chiuso |
| 10 | **Classi delle stazioni:** con la nuova rete 37 stazioni su 48 sono nel canale singolo; restano 5 a valle di una biforcazione, 1 a valle di una confluenza, 5 negli altri tratti a canali multipli, con intervalli molto larghi (§9.8) | in attesa della risposta di Carlo (domanda 1). Run del 19/11 con ADCP e sottostazioni: 44 / 8 / 1 / 6 stazioni; ambigue s015 e s035 | due classi, distanza dal nodo, o quattro classi solo descrittive; se il dislivello tra i rami delle sottostazioni diventa una metrica |
| 11 | **Errore sistematico positivo** del DTM corretto: sceso a +0.12 m (per punto) dopo i punti 6 e 7; +0.15–0.16 m senza s006, +0.30 m negli altri tratti a canali multipli (§9.8) | aperto; **dopo il punto 8 sale a +0.22 m** sulle 48 stazioni GNSS (+0.33 m con gli ADCP) ed è ora l'errore principale (§9.10). **Analisi del 20/11** (§9.11): viene soprattutto dalle pozze (pelo libero quasi piatto, moto uniforme con la pendenza del tratto); residuo correlato 0.91 con la profondità osservata, −0.50 con il rapporto tra pendenza locale e del tratto. La portata del volo è confermata (§9.12), quindi l'errore non dipende da lei; nelle pozze è un limite del metodo (punto 16) | rilanciare con la nuova regola per i salti isolati e rivedere l'errore residuo; scabrezza e larghezze nell'analisi di sensitività |
| 12 | **Intervalli di confidenza:** `validation_gnss_mean.m` ricampiona i punti, e gli intervalli che ne escono sono troppo stretti; quelli corretti ricampionano le stazioni | implementato in `ModuleComparison.m`: RMSE, errore medio, mediana e NMAD con intervallo per stazione, riduzione dell'RMSE rispetto al DTM grezzo e differenza tra i moduli con intervallo appaiato | nel paper riportare questi intervalli al posto di quelli di `validation_gnss_mean.m` (vedi PaperBAT, 5.5). Run del 19/11: anche in `ModuleComparison.m` RMSE ed errore medio di gruppo sono calcolati sui punti, quindi un transetto ADCP pesa 5–6 stazioni GNSS dentro ogni campione bootstrap; decidere se pesare le stazioni |
| 13 | **Confronto dei due moduli da ripetere** (§9.5) | ripetuto dopo i punti 6 e 7 (§9.8) e dopo il punto 8 (§9.10): differenza non significativa in tutte le combinazioni di stazioni (+0.2/+0.5 cm, p 0.45–0.77) | decidere come presentarlo nel paper (classi, domanda 1) |
| 14 | **Test del segno e parità:** `compare_modules_gnss.m` (righe 111–112 e 170–171) conta come vittorie le stazioni con RMSE uguale nei due moduli a meno di 10⁻¹¹ (37 su 48), che finiscono a HM o a MM per rumore di arrotondamento; nel log risulta "MM migliore in 9/48, p = 4·10⁻⁵". La colonna `MM_minus_HM_cm` della tabella delle stazioni è la differenza dell'errore medio, non dell'RMSE | **corretto nel codice** (19/11/2026): `opts.tie_tol` = 1 mm (impostato in `ModuleComparison.m`); `closer` vale MM, HM o `tie`; colonne `dRMSE_MM_HM_cm` (differenza di RMSE) e `dME_MM_HM_cm` (differenza dell'errore medio, ex `MM_minus_HM_cm`); `n_ties` nella tabella dei gruppi; barre grigie per le parità nel pannello B. Con i CSV del 18/11 sera: MM più vicino in 3, HM in 8, 37 pari, p ≈ 0.23 | **chiuso** (§9.10): sui CSV del 18/11 la funzione corretta dà 3 / 8 / 37, p = 0.23; sui dati del 19/11, 5 / 7 / 36 (48 GNSS) e 6 / 10 / 43 (tutte) |
| 15 | **Trasporto solido al volo praticamente nullo:** i flussi solidi di MM alla portata del volo sono 10⁻¹⁶–10⁻²¹ m³/s, quindi le ripartizioni $\psi_x$ al volo agiscono su quantità nulle | da verificare | confermare che è atteso (10.8 m³/s è molto sotto la soglia di moto) e dichiararlo: al volo il modulo morfodinamico agisce solo tramite $\psi$ |
| 16 | **Pozze e sezioni molto irregolari:** il fondo rettangolarizzato con il moto uniforme e la pendenza del tratto non riproduce le pozze a pelo libero quasi piatto né i dislivelli locali tra rami vicini (s003-t12, s007-t06; §9.10–9.11) | **limite del metodo, deciso il 20/11/2026**: si dichiara, non si corregge | riportarlo tra i limiti nella Discussion |
| 17 | **Anomalie minori di figure e log** (run del 19/11): nomi delle sottostazioni resi come pedici TeX, titoli "GNSS" anche con gli ADCP, tabella diagnostica dei tratti che scriveva "apply smoothing" per tutti | **chiuso il 20/11/2026**: verificato nel run del 20/11 (trattino, "in situ", niente pedici, tabella diagnostica coerente) | — |
| 18 | **Portata del volo** (§9.12): livello a San Benigno 0.97 m il 17/03/2025; con la scala 2023 corrisponde a 23.4 m³/s, mentre il config usa 10.8 m³/s (fonte non annotata) e un livello di 0.695 m ricavato con una regressione lineare sulla scala | **chiuso il 20/11/2026** (§9.12): 10.76 m³/s, portata giornaliera ARPA del 17/03/2025, coerente con 0.97 m e la scala valida dal 01/01 al 16/04/2025 (10.9 m³/s); la scala 2023 non vale per il 2025 | — |
| 19 | **Fine del tratto 0 al nodo 15** (§9.13): il pelo libero è quasi piatto (260.03–260.13 m) e risale a 260.26, 261.36, 261.36 m nelle ultime tre sezioni; la regola corregge solo la 14. Il tratto 0 ha pendenza minima in tutto il modulo idraulico e il nodo 15 ricava $S_0$ allargando la finestra | in chiusura. **20/11:** la maschera del tratto 0 finisce sull'ombra del ponte stradale (confermato da Elisa); i tratti 1 e 2 cominciano a valle del ponte. **Run del 20/11 sera (§9.14):** maschera accorciata, risalita sparita, pendenza del tratto 0 = 2.0·10⁻⁴; biforcazione 0 → 1, 2 forzata nel config. **Chiuso il 20/11**; resta il punto 21 | requisito della segmentazione: niente ponti o strutture dentro le maschere |
| 20 | **Baseline semplice** (proposta del 20/11/2026): oggi l'unico confronto è con il DTM grezzo | proposta | DTM grezzo meno un tirante uniforme da Manning con la portata totale, senza rete e senza nodi; validata con le stesse stazioni e gli stessi intervalli. Mostra quanto contano la rete e la ripartizione ai nodi (PaperBAT, attività 2) |
| 21 | **$S_0$ al nodo 13 (ex 15), testa della rete** (§9.14): la finestra a monte copre solo il tratto 0 accorciato (114 m, pendenza 2.0·10⁻⁴), e $S_0$ è alzata alla soglia minima 1.22·10⁻³. Il nodo diventa bilanciato ($\psi$ = 0.73 contro 0.82), $D_0$ = 2.53 m. Nessuna stazione nei tratti 0–3 | aperto | decidere come trattare i nodi di testa: tenere la soglia e dichiararlo; usare la pendenza attraverso il ponte (fine del tratto 0 → inizio dei rami); oppure ripartire la portata con $\psi_{hyd}$ senza risolvere il nodo. Da mettere nella sensitività |

### 9.5 Validazione GNSS: modulo idraulico e modulo morfodinamico

> **Nota (19/11/2026).** I numeri di questa sezione vengono dal run del 18/11 pomeriggio, prima della correzione della portata e dei rami 15, 16, 29. I numeri aggiornati sono nella §9.8; questa sezione resta come documentazione del metodo e del run precedente.

Il confronto tra i due moduli si ripete con `Final/ModuleComparison.m` (dopo i due moduli; figure e tabelle in `Output_folder/04_Validation/ModuleComparison`).

**Dati e metodo.** I DTM corretti dal modulo idraulico e dal modulo morfodinamico (run di riferimento con le larghezze bankfull) sono confrontati con i punti GNSS del 07/03/2025. Una **stazione** è una sezione trasversale misurata con il GNSS: tutti i punti `sNNN_pMMM` con lo stesso `sNNN`. Le stazioni sono 48, con 956 punti dentro l'area corretta. La numerazione delle stazioni non segue la rete da monte a valle: per esempio s001 e s002 sono nel tratto 14, a valle di s004–s006 nel tratto 8. Per ogni punto il residuo è la quota del DTM nella posizione del punto meno la quota media della sua stazione (quota ellissoidica), perché il modello assegna una sola quota a tutta la sezione. Un residuo positivo indica un DTM più alto del fondo misurato.

Ogni punto è assegnato alla sezione della centerline più vicina, e la stazione al tratto della maggioranza dei suoi punti. In quattro stazioni i punti cadono su due tratti (s003: 46% dei punti sul tratto assegnato, s007: 64%, s015: 39%, s035: 89%), quindi la loro classificazione è incerta. Il tipo di tratto viene dalla topologia della rete: la frazione di filoni vale 1 alla testa, si dimezza in ogni ramo di una biforcazione e si somma alle confluenze.

- **Canale singolo** (frazione 1, tratti 8 e 14): 8 stazioni.
- **Tratto che inizia a una biforcazione:** 18 stazioni. Dodici sono nel tratto 30, lungo 1.8 km, per cui "a valle della biforcazione" vale davvero solo per le prime.
- **Tratto a canali multipli che inizia a una confluenza:** 10 stazioni.
- **Altri tratti a canali multipli:** 12 stazioni.

Gli intervalli al 95% sono calcolati con il bootstrap delle **stazioni**, non dei punti, perché i punti di una stessa sezione condividono la quota di riferimento e non sono indipendenti. Ricampionando i punti si ottengono intervalli troppo stretti (0.35–0.41 m per l'RMSE del modulo idraulico, contro 0.27–0.50 m).

**Risultati** (metri; RMSE con intervallo al 95% tra parentesi quadre):

| Tipo di tratto | Stazioni / punti | DTM grezzo ME / RMSE | Modulo idraulico ME / RMSE | Modulo morfodinamico ME / RMSE | RMSE morfodinamico − idraulico | Stazioni in cui il morfodinamico è più vicino |
| --- | --- | --- | --- | --- | --- | --- |
| Canale singolo | 8 / 147 | +0.44 / 0.47 | −0.05 / 0.54 [0.10–0.93] | −0.07 / 0.55 [0.10–0.95] | +0.007 [−0.012, +0.018] | 5 / 8 |
| Canale singolo senza s006 | 7 / 130 | +0.47 / 0.49 | +0.13 / 0.21 [0.07–0.31] | +0.11 / 0.20 [0.08–0.30] | −0.007 [−0.012, +0.008] | 5 / 7 |
| A valle di una biforcazione | 18 / 384 | +0.54 / 0.57 | +0.24 / 0.32 [0.21–0.42] | +0.26 / 0.33 [0.23–0.44] | +0.015 [+0.009, +0.021] | 1 / 18 |
| A valle di una confluenza | 10 / 213 | +0.54 / 0.54 | +0.14 / 0.35 [0.25–0.46] | +0.16 / 0.35 [0.25–0.45] | +0.000 [−0.009, +0.011] | 2 / 10 |
| Altri tratti a canali multipli | 12 / 212 | +0.69 / 0.76 | +0.24 / 0.37 [0.15–0.51] | +0.26 / 0.40 [0.17–0.54] | +0.027 [+0.016, +0.035] | 2 / 12 |
| Tutte | 48 / 956 | +0.56 / 0.60 | +0.17 / 0.38 [0.27–0.50] | +0.19 / 0.39 [0.28–0.51] | +0.013 [+0.007, +0.018] | 10 / 48 |
| Tutte senza s006 | 47 / 939 | +0.56 / 0.60 | +0.20 / 0.32 [0.26–0.39] | +0.22 / 0.34 [0.27–0.40] | +0.013 [+0.005, +0.018] | 10 / 47 |

**Confronto stazione per stazione.** Il modulo idraulico ha l'RMSE più basso in 38 stazioni su 48 (test del segno a due code, p ≈ 6·10⁻⁵). La differenza di RMSE complessiva, +1.3 cm (intervallo +0.7…+1.8 cm), esclude lo zero: è piccola ma sistematica, non casuale. È più marcata nei tratti a valle di una biforcazione, dove il modulo morfodinamico è più vicino ai GNSS in una sola stazione su 18, e negli altri tratti a canali multipli (+2.7 cm). Nel canale singolo il modulo morfodinamico è più vicino in 5 stazioni su 8, per una ragione che non dipende dal nodo (punto 1 qui sotto).

**Dove agisce il modulo morfodinamico.** Le quote di fondo dei due moduli lungo la centerline coincidono nei tratti 0 e 3, a monte della prima biforcazione. Più a valle differiscono in tutti i tratti, non solo nel primo tratto dei rami: il modulo morfodinamico cambia la ripartizione della portata a ogni nodo, e con essa la portata e il tirante di tutta la rete a valle. Le portate di ogni tratto, ricavate dalle quote di fondo dei due moduli con $Q = W\,h^{5/3}S^{1/2}/n$ e $h = z_{DTM} - z_{fondo}$, mostrano però due problemi che pesano su questo confronto.

1. **Il modulo idraulico non conserva la portata.** Il tratto 6 riceve 5.00 m³/s, cioè la portata del tratto 5, invece dei 5.80 m³/s del tratto 4 da cui nasce. Da lì in poi la rete del modulo idraulico porta 9.99 m³/s invece di 10.80 m³/s (−7.5%). La causa più probabile, coerente con questi valori, è nel ciclo delle sezioni di `Preprocessing_HydraulicModule.m`: nei blocchi del canale singolo e della confluenza la portata viene scritta anche nella riga successiva (`Q_map(i+1)`), anche quando quella riga è la prima sezione di un altro tratto, e quel tratto conserva poi il valore sbagliato. Il modulo morfodinamico propaga la portata per tratti e la conserva a tutti i nodi. Per questo, nei tratti a canale singolo 8 e 14, i due moduli differiscono di −2.5 e −1.6 cm: il morfodinamico, con la portata giusta, ha un fondo più profondo. Sono queste le differenze delle stazioni s001, s002, s004–s006 e s028–s030, e una parte delle differenze a valle. La propagazione è stata corretta nel codice dopo questa validazione (§9.4, punto 6): i numeri di questa sezione sono quelli del run precedente.
2. **Due rami non hanno un tratto a valle.** Il ramo 16–15 e il ramo 29 terminano senza ricollegarsi alla rete: la fine del tratto 15 dista circa 100 m dal tratto più vicino, quella del 29 circa 140 m, mentre la soglia di collegamento è 51 m. La loro portata esce dalla rete: 1.47 + 0.70 m³/s nel modulo idraulico, 2.65 + 1.39 m³/s nel morfodinamico. A valle del tratto 30 la rete porta così 7.83 m³/s (idraulico) e 6.75 m³/s (morfodinamico) invece di 10.80 m³/s. Meno portata significa tiranti minori e fondo più alto: è una possibile causa dell'errore positivo, più forte nella rete a canali multipli. L'analisi dei due rami (§9.7) mostra che al momento del volo sono quasi fermi e non portano una portata apprezzabile.

Questi valori sono ricavati a ritroso dalle quote di fondo. Vanno confermati con le portate salvate dai due moduli (`Q_map` nei workspace), che `Final/check_discharge_conservation.m` controlla nodo per nodo.

**I forti residui negativi** delle stazioni s006, s019 e s020 (in media −1.44, −0.71 e −0.36 m) non dipendono dal tirante calcolato né dal rilievo GNSS. Il fondo corretto è $z_{DTM} - h$, dove $z_{DTM}$ è la quota della sezione della centerline (media del DTM di segmento lungo la sezione, dentro la maschera). Vicino a queste stazioni, e a s032, $z_{DTM}$ è 0.6–1.6 m più bassa del DTM grezzo misurato lungo la sezione GNSS e oscilla di 0.7–1.5 m su ±3 sezioni (60 m). Sulle 48 stazioni, la differenza tra DTM grezzo e $z_{DTM}$ della sezione ha una correlazione di −0.78 con il residuo del modulo idraulico. A s006, per esempio, il DTM grezzo lungo la sezione è piatto (±4 cm) e sta 0.20 m sopra la media GNSS; $z_{DTM}$ vale invece 253.90 m, 1.56 m più in basso, e il tirante del modulo (0.58 m) viene sottratto da lì. Il filtro sulle oscillazioni isolate di $z_{DTM}$ interviene solo sopra 1 m di curvatura, e i tratti 8, 18 e 33 non sono tra quelli lisciati (`segments_to_smooth`). Cause possibili e proposta di correzione: §9.6.

**Conclusioni.**

1. **Entrambi i moduli migliorano nettamente il DTM grezzo:** l'RMSE scende da 0.60 a 0.38–0.39 m e l'errore sistematico da +0.56 a +0.17–0.19 m.
2. **Il modulo idraulico è leggermente ma sistematicamente più vicino ai GNSS:** +1.3 cm di RMSE (+0.7…+1.8 cm), più vicino in 38 stazioni su 48. La differenza è piccola rispetto all'errore complessivo e all'incertezza dell'RMSE stesso (0.27–0.50 m).
3. **Questo confronto non isola l'effetto del nodo:** il modulo idraulico perde portata dal tratto 6 in giù, e in entrambi i moduli la portata dei rami 16–15 e 29 esce dalla rete. Il confronto va ripetuto dopo aver corretto la propagazione della portata e verificato i due rami.
4. **L'errore dominante resta sistematico e positivo** (DTM corretto troppo alto), più forte nella rete a canali multipli. Le cause da verificare sono la portata persa nei rami senza uscita, la scabrezza (n = 0.03 costante) e le larghezze della maschera nei rami stretti.
5. **I forti residui negativi** (s006, s019, s020) vengono dalle oscillazioni di $z_{DTM}$ lungo la centerline, non dal tirante né dal rilievo, e vanno corretti nel pre-processing.
6. **Conseguenza.** Allo stato attuale il contributo del modulo morfodinamico non è un miglioramento dell'accuratezza del fondo, ma la coerenza fisica della ripartizione di acqua e sedimento ai nodi (equilibrio, selezione granulometrica). Per misurarne l'effetto sull'accuratezza servono le correzioni dei punti 3 e 5 e dati che distinguano la ripartizione tra i rami, come le misure ADCP di portata nei rami.

### 9.6 Quota di partenza del fondo ($z_{DTM}$)

**Che cos'è.** In entrambi i moduli il fondo di una sezione è $\eta = z_{DTM} - D$: $z_{DTM}$ è la quota del pelo libero al momento del volo LiDAR (17/03/2025) e $D$ è il tirante calcolato. $z_{DTM}$ viene calcolata nel pre-processing da `calculate_widths.m`: è la media del DTM di segmento (`DTM_segments/dtm_*.tif`), campionato con interpolazione bilineare lungo la sezione dal primo all'ultimo pixel della maschera, esclusi i valori negativi. Poi un filtro corregge i salti isolati (curvatura oltre 1 m) e i tratti elencati in `segments_to_smooth` vengono lisciati con una mediana mobile.

**Il problema.** Un errore su $z_{DTM}$ passa tale e quale nel fondo. Vicino alle stazioni s006, s019, s020 e s032, $z_{DTM}$ è 0.6–1.6 m più bassa del DTM LiDAR completo letto lungo la sezione GNSS, e oscilla di 0.7–1.5 m tra sezioni distanti 10 m. A s006, per esempio, le sezioni 168–172 valgono 254.23, 254.85, 253.90, 255.40 e 254.42 m, mentre il DTM completo lungo la sezione GNSS è piatto a 255.45 m (±4 cm). Sulle 48 stazioni, la differenza tra DTM completo e $z_{DTM}$ ha una correlazione di −0.78 con il residuo del modulo idraulico: è la causa principale dei forti residui negativi (§9.5).

**Non dipende dalla data del rilievo.** I punti GNSS sono del 07/03/2025, il volo LiDAR del 17/03/2025, con un livello idrometrico più alto al momento del volo (dal grafico dell'idrometro, circa 86 il 7 marzo e circa 97 il 17 marzo). Questo può spiegare una parte della differenza tra DTM grezzo e fondo GNSS, ma non la differenza tra $z_{DTM}$ e il DTM completo, che vengono dallo stesso volo: il problema è nel modo in cui si ricava la quota della sezione dal dato LiDAR. La data diversa conta invece se il fondo è cambiato tra il 7 e il 17 marzo (il colmo del 15–16 marzo è modesto).

**Cause possibili, da verificare.**

1. **Nodata dei DTM di segmento.** Il pre-processing non converte in NaN il valore di nodata dei DTM di segmento (lo fa solo per il DTM completo). Con l'interpolazione bilineare, un punto vicino a un pixel senza dato mescola la quota vera con il nodata (per esempio −32767): l'esclusione dei valori negativi toglie i casi peggiori, ma non quelli in cui il nodata ha un peso piccolo, che restano abbassati di decimetri o metri. Spiegherebbe perché gli errori grandi sono tutti verso il basso.
2. **Riempimento dei vuoti.** I DTM di segmento sono riempiti (i vuoti sull'acqua sono interpolati) e il riempimento può dare quote più basse del pelo libero.
3. **Media lungo tutta la sezione.** La media comprende tutto ciò che sta tra il primo e l'ultimo pixel della maschera, comprese isole o barre in mezzo.

**Verifica alla sezione 170 (18/11/2026).** Confermata la causa 1, con un effetto aggiuntivo del filtro sui salti; le cause 2 e 3 non sono confermate in questo punto. Dettagli nella §9.8.

**Correzione adottata (19/11/2026).** Si tiene l'approccio attuale (media del DTM di segmento lungo la sezione) e si corregge solo la causa 1:

1. in `Preprocessing_HydraulicModule.m` il nodata −32767 di ogni DTM di segmento diventa NaN subito dopo la lettura (`dtm(dtm == -32767) = NaN`), come già si faceva per il DTM completo;
2. in `calculate_widths.m` la media della sezione esclude esplicitamente i NaN. L'interpolazione bilineare restituisce NaN appena uno dei quattro pixel vicini è NaN, quindi i punti contaminati dal nodata escono dalla media invece di abbassarla.

Il filtro sui salti isolati e `segments_to_smooth` restano invariati. Dopo il rilancio del pre-processing si guarda se segnalano ancora salti e oscillazioni; i punti 1 e 2 della proposta di correzione qui sopra (DTM completo, mediana sui pixel bagnati, profilo monotono) si valutano solo dopo.

**Esito (run del 19/11/2026, §9.10).** Salti segnalati dal filtro da 149 a 2; i tratti di `segments_to_smooth` non oscillano più; a s006, s019 e s020 $z_{DTM}$ coincide con il DTM completo entro 4 cm e i residui medi diventano −0.21, +0.27 e −0.04 m. La correlazione del punto 3 della proposta qui sotto non è stata ricalcolata. Il profilo monotono (punto 2 della proposta) servirebbe ora solo per poche risalite di 10–55 cm e per il salto di +1.1 m in fondo al tratto 0: da decidere (punto 8 della §9.4).

**Regola per i salti isolati (20/11/2026).** Il filtro sui salti (curvatura oltre 1 m) e la lista `segments_to_smooth` (mediana mobile su tratti scelti a mano) sono sostituiti da una sola regola, uguale per tutti i tratti. Per ogni sezione interna $i$ di un tratto, con $d_1=z_i-z_{i-1}$ e $d_2=z_i-z_{i+1}$, la sezione è un salto isolato se è un estremo locale stretto ($d_1 d_2>0$) e se si stacca da **entrambe** le vicine di più della soglia: $\min(|d_1|,|d_2|)>\max(\text{glitch\_tol},\,3\sigma_e)$, con $\sigma_e=1.4826\,\mathrm{MAD}(e)$ del tratto, $e_i=z_i-(z_{i-1}+z_{i+1})/2$ e glitch_tol = 0.20 m. Un salto isolato è sostituito dalla media delle due vicine, sui valori grezzi. La condizione su entrambi i lati lascia intatti i gradini. Con i dati del 20/11 la regola corregge 2 sezioni: 728 (tratto 26, +0.36 e +0.39 m dalle vicine) e 1239 (tratto 35, +0.55 e +0.51 m). La prima versione della regola, provata nel run del 20/11 (§9.13), chiedeva solo $|e_i|>$ soglia e correggeva anche la 492 (piede di un gradino di −0.71 m, che portava a metà gradino) e la 14 (salto alla fine del tratto 0): per questo è stata resa più stretta.

**Limite che resta.** La media comprende tutti i punti del transetto tra il primo e l'ultimo pixel della maschera del segmento (passo 0.5 m). Se in mezzo alla sezione c'è un'isola o una barra fuori dalla maschera, la sua quota entra nella media e alza $z_{DTM}$. Alla sezione 170 l'effetto è trascurabile (media e mediana differiscono di 1 cm), ma non è stato verificato su tutta la rete; va dichiarato tra i limiti del metodo.

**Proposta di correzione.**

1. Calcolare $z_{DTM}$ dal DTM LiDAR completo, con il nodata convertito in NaN, solo sui pixel bagnati della sezione (maschera dell'area bagnata), con la mediana invece della media.
2. Imporre che il pelo libero non risalga verso valle lungo ogni tratto (regressione monotona di $z_{DTM}$ lungo il tratto, con lo stesso livello ai nodi). Sostituirebbe la lista manuale `segments_to_smooth` e il filtro a soglia 1 m.
3. Verifica: la correlazione tra (DTM completo − $z_{DTM}$) e il residuo deve sparire, e s006, s019 e s020 non devono più avere residui tra −0.4 e −1.4 m. La tabella dei profili di `ModuleComparison.m` riporta già, per le sezioni vicine a queste stazioni, $z_{DTM}$ e il DTM completo letto sull'asse.

### 9.7 Rami 16–15 e 29

La maschera e la quota del pelo libero lungo la centerline mostrano che i due rami non sono biforcazioni che portano via acqua.

- **Ramo 29.** Il pelo libero sale da 241.60 m, alla giunzione con i tratti 28 e 30, a 243.80 m all'altro estremo (pendenza 3.6·10⁻³, simile a quella dell'asta principale). Se c'è un deflusso, va quindi verso la giunzione: il ramo 29 non toglie acqua alla rete. L'estremo di monte è a circa 140 m dal tratto 22, dove il pelo libero è a 244.9–245.0 m, circa 1.1 m più in alto; al momento del volo i due non sono collegati (vedi sotto). Nel modello invece il tratto è girato (è in `segments_to_flip`) e trattato come un ramo della biforcazione 817, che porta via 0.70 m³/s (modulo idraulico) o 1.39 m³/s (modulo morfodinamico).
- **Ramo 16–15.** Il pelo libero è praticamente orizzontale, da 251.15 a 251.24 m su circa 720 m, alla stessa quota della giunzione (fine del tratto 14: 251.12 m; inizio del tratto 17: 251.11 m), e l'estremo del tratto 15 non è vicino a nessun altro tratto. È un ramo quasi fermo, riempito dal rigurgito della giunzione, che non porta portata. Nel modello invece riceve 1.47 m³/s (idraulico) o 2.65 m³/s (morfodinamico) alla biforcazione 366.

**Osservazioni sul campo e sulle immagini.** Al momento del volo entrambi i rami contengono acqua ma sono quasi fermi, senza una portata apprezzabile. Il ramo 29 non è collegato al tratto 22, ma tra i due c'è una zona attiva: oltre una certa portata i due si collegano. Il ramo 16–15 fa parte di un ampio meandro, bagnato solo in parte al momento del volo. La pendenza del pelo libero lungo il ramo 29 (3.6·10⁻³) può indicare un deflusso molto piccolo, per esempio per filtrazione dal tratto 22, ma non una portata paragonabile a quella che il modello gli assegna.

Ne segue che la rete cambia con la portata: alle condizioni del volo i due rami non portano portata, mentre alla portata formativa (bankfull), su cui lavora il nodo del modulo morfodinamico, sono attivi e collegati in modo diverso (il ramo 29 alimentato dal tratto 22, il ramo 16–15 come parte del meandro). I nodi 366 e 817, come sono costruiti ora, non rappresentano né l'una né l'altra condizione.

**Conseguenze e scelte da fare.**

- **Ramo 16–15:** portata nulla, quindi tutta la portata del tratto 14 va al tratto 17 e il nodo 366 non è una biforcazione. Il fondo del ramo non si può ricavare dal moto uniforme (con $Q = 0$ il tirante è nullo): va deciso come stimarlo, per esempio con il tirante della giunzione, oppure lasciandolo non corretto. Nel ramo non ci sono stazioni GNSS.
- **Ramo 29:** alle condizioni del volo la portata è nulla o trascurabile, quindi tutta la portata del tratto 28 va al tratto 30 e il nodo 817 non è una biforcazione; anche qui il fondo non si ricava dal moto uniforme.
- **Modulo morfodinamico:** va deciso come trattare i nodi 366 e 817 alla portata bankfull, quando i due rami sono attivi con una topologia diversa da quella del volo (escluderli dal sistema del nodo, oppure ricostruire la rete alla portata bankfull).
- In tutti i casi la rete a valle dei tratti 17 e 30 torna a portare tutta la portata del volo.

**Scelte adottate (18/11/2026).** Non ci sono altri rami fermi nella rete. Nel codice di `Final` i tratti 15, 16 e 29 sono elencati in `no_flow_segments` nel config, e i due moduli li trattano così:

1. **portata nulla al volo:** alle biforcazioni 366 e 817 tutta la portata (e, nel modulo morfodinamico, tutto il trasporto solido con la composizione della madre) va all'altro ramo; il sistema del nodo non viene risolto (tipo di nodo "branch without discharge");
2. **fondo non ricostruito:** tirante e fondo dei tre tratti non vengono calcolati e i loro pixel bagnati restano con il DTM grezzo (`exclude_no_flow_pixels.m`); nel paper vanno dichiarati come rami fermi al momento del volo, non ricostruiti (non ci sono stazioni GNSS);
3. **portata bankfull:** la ricostruzione della rete alla portata formativa, in cui i rami sono attivi, resta un'estensione da discutere tra i limiti.

Con la nuova rete i tratti 17, 20, 28, 30, 33 e 35 tornano a canale singolo (portano tutta la portata), per cui la classificazione delle stazioni per tipo di tratto della §9.5 cambia: per esempio le 12 stazioni del tratto 30 passano da "a valle di una biforcazione" a "canale singolo".

### 9.8 Run del 18/11/2026 (sera): rete corretta

Run con la portata propagata solo dentro il tratto (punto 6 della §9.4) e i tratti 15, 16, 29 senza portata (`no_flow_segments`, §9.7); config con `node_width_source = 'dtm'`. I tre script (`Preprocessing_HydraulicModule`, `MorphodynamicModule`, `ModuleComparison`) terminano senza errori. Warning non bloccanti nel modulo morfodinamico: tirante di sfioro sotto quello critico (quindi $D_{bf}=D_c$) in 17 sezioni; sezione 484 con NaN nel profilo; larghezza bankfull del ramo 2 della Bif 469 non calcolabile (si usa il rapporto delle maschere). Figure e tabelle in `Final/Output_folder` (resoconto completo: `Output_folder/RESOCONTO.md`).

**Bilancio di portata.** Nessun nodo sbilanciato in nessuno dei due moduli: lo sbilancio massimo è 4.6·10⁻¹² m³/s (HM, nodo 205), in MM è nullo; la portata è costante dentro ogni tratto. I tratti senza valle sono esattamente quelli attesi: 15 e 29 (Q = 0) e 42 (uscita, 10.80 m³/s in entrambi i moduli; prima HM 9.99 → 7.83, MM 6.75). Le portate dei rami differiscono tra i moduli soprattutto nei tratti 10 (HM 2.86, MM 5.75 m³/s), 9 (7.94 contro 5.05), 4/6 (5.81 contro 3.81) e 22 (7.08 contro 5.49).

**Nodi.** 366 e 817 sono "branch without discharge" in entrambi i moduli. Gli altri 10 nodi sono risolti dal sistema (residuo ≤ 1.9·10⁻¹⁵).

Come leggere la tabella (simboli nella §1): $\beta=W_0^{bf}/(2D_0)$ è il rapporto di forma bankfull della madre e decide quante soluzioni ha il nodo; $\psi$ è la quota di portata del ramo 1 nella soluzione adottata; $W_1/W_0$ è la $\psi$ della soluzione bilanciata, quindi $\psi\neq W_1/W_0$ indica una soluzione asimmetrica; $\psi_{hyd}$ è la quota del ramo 1 nel modulo idraulico, usata per scegliere tra le asimmetriche; $\theta_{bf}/\theta^*_c$ dice quale limite fissa il tirante bankfull della madre (1.62: Shields pieno; 1.00: tirante critico, sponde basse; intermedio: livello di sfioro).

| Bif | $\beta$ | $\psi$ | $W_1/W_0$ | $\psi_{hyd}$ | soluzioni | tipo | $\theta/\theta_c$ madre | run prec. (§9.3): $\beta$ / $\psi$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 15 | 31.1 | 0.41 | 0.69 | 0.49 | 3 | asimmetrica | 1.00 | 31.1 / 0.41 |
| 106 | 23.0 | 0.35 | 0.59 | 0.54 | 3 | asimmetrica | 1.62 | 23.0 / 0.35 |
| 205 | 13.7 | 0.47 | 0.35 | 0.74 | 3 | asimmetrica | 1.10 | 13.9 / 0.47 |
| 244 | 30.9 | 0.80 | 0.66 | 0.76 | 3 | asimmetrica | 1.25 | 31.2 / 0.80 |
| 469 | 14.2 | 0.71 | 0.68 | 0.67 | 3 | asimmetrica | 1.20 | 15.4 / 0.72 |
| 584 | 33.3 | 0.51 | 0.79 | 0.66 | 3 | asimmetrica | 1.00 | 33.6 / 0.51 |
| 645 | 28.0 | 0.46 | 0.72 | 0.35 | 3 | asimmetrica | 1.00 | 28.3 / 0.46 |
| 1070 | 16.3 | 0.78 | 0.80 | 0.83 | 3 | asimmetrica | 1.00 | 14.4 / 0.76 |
| 1293 | 18.4 | 0.69 | 0.60 | 0.67 | 3 | asimmetrica | 1.00 | 18.6 / 0.69 |
| 1369 | **11.7** | 0.23 | 0.23 | 0.27 | **1** | **bilanciata** | 1.62 | 12.0 / 0.26 |

1369 diventa bilanciato: con $\beta$ = 11.7, sotto la soglia critica (≈ 12–13), esiste solo la soluzione $\psi=W_1/W_0$. È un'ulteriore conferma della sensibilità a $\beta$ (§9.3, punto 1 della §9.4). Il tirante bankfull è quello critico in 5 madri su 10 (prima 9 su 12); 106 e 1369 sono al tirante di Shields, 205, 244 e 469 al livello di sfioro. In 205, 584 e 645 $\psi$ resta lontano da $\psi_{hyd}$ (punto 4).

**Trasporto solido al volo.** I flussi solidi alla portata del volo sono dell'ordine di 10⁻¹⁶–10⁻²¹ m³/s (per esempio Bif 15: $Q_{s,F}$ della madre 1.4·10⁻¹⁶ m³/s): a 10.8 m³/s si è molto sotto la soglia di moto. Al volo, quindi, il modulo morfodinamico agisce sul fondo solo attraverso la ripartizione della portata $\psi$; le ripartizioni $\psi_x$ dei flussi solidi agiscono su quantità nulle (punto 15 della §9.4).

**Area non corretta.** 127 559 pixel bagnati (6.4% dell'area bagnata) restano al DTM grezzo, gli stessi nei due moduli: sono i tratti 15, 16, 29 (sezioni 367–445 e 818–888). Non contengono stazioni GNSS; i punti GNSS dentro l'area corretta restano 956.

**Validazione GNSS** (48 stazioni, 956 punti; m; intervalli al 95% dal bootstrap delle stazioni). Classi delle stazioni con la nuova rete: canale singolo 37 stazioni / 752 punti; a valle di una biforcazione 5 / 110 (s003, s015, s016, s032, s046); a valle di una confluenza 1 / 20 (s048); altri tratti a canali multipli 5 / 74 (s007, s008, s044, s045, s047).

| Gruppo | st / pt | RMSE grezzo | RMSE HM | RMSE MM | ME grezzo / HM / MM | RMSE MM − HM |
| --- | --- | --- | --- | --- | --- | --- |
| Canale singolo | 37 / 752 | 0.61 [0.53–0.68] | 0.38 [0.24–0.52] | 0.38 [0.24–0.52] | +0.57 / +0.12 / +0.12 | +0.000 |
| Canale singolo senza s006 | 36 / 735 | 0.61 [0.54–0.69] | 0.30 [0.23–0.38] | 0.30 [0.23–0.38] | +0.58 / +0.16 / +0.16 | 0.000 |
| A valle di una biforcazione | 5 / 110 | 0.47 [0.37–0.61] | 0.15 [0.12–0.18] | 0.15 [0.11–0.21] | +0.43 / +0.04 / +0.05 | 0.000 [−0.030, +0.034] |
| A valle di una confluenza | 1 / 20 | 0.35 | 0.04 | 0.04 | +0.35 / +0.03 / +0.03 | 0.000 |
| Altri tratti a canali multipli | 5 / 74 | 0.74 [0.36–1.03] | 0.44 [0.09–0.60] | 0.47 [0.10–0.64] | +0.63 / +0.30 / +0.29 | +0.035 [−0.005, +0.047] |
| **Tutte** | 48 / 956 | 0.60 [0.53–0.67] | **0.36 [0.25–0.48]** | **0.36 [0.25–0.49]** | +0.56 / +0.12 / +0.12 | **+0.004 [−0.001, +0.009]** |
| Tutte senza s006 | 47 / 939 | 0.60 [0.53–0.68] | 0.30 [0.23–0.36] | 0.30 [0.23–0.37] | +0.56 / +0.15 / +0.15 | +0.004 [−0.001, +0.011] |

La riduzione dell'RMSE rispetto al DTM grezzo (intervallo appaiato per stazione) è +0.24 m su tutte le stazioni (HM [+0.11, +0.34], MM [+0.10, +0.34]) e +0.31 / +0.30 m senza s006 (HM [+0.26, +0.35], MM [+0.26, +0.35]). Rispetto al run precedente (§9.5) l'errore medio per punto scende da +0.17 a +0.12 m (HM) e da +0.19 a +0.12 m (MM); l'RMSE per stazione da 0.38/0.39 a 0.36/0.36 m.

**Confronto tra i moduli.** La differenza di RMSE MM − HM su tutte le stazioni è +0.4 cm [−0.1, +0.9]: l'intervallo comprende lo zero, e i due moduli non sono distinguibili (prima +1.3 cm [+0.7, +1.8]). Nelle 37 stazioni del canale singolo i due moduli hanno la stessa portata e lo stesso fondo, quindi lo stesso RMSE. Escludendo queste parità (|ΔRMSE| ≤ 1 mm), MM è più vicino ai GNSS in 3 stazioni (s003 −6.0, s047 −6.5, s015 −0.9 cm) e HM in 8 (s044 +9.4, s046 +8.4, s045 +5.1, s008 +4.6, s007 +4.5, s032 +1.1, s016 +0.9, s006 +0.3 cm): test del segno p ≈ 0.23. Il "9/48, p = 4·10⁻⁵" del log conta le parità numeriche (punto 14 della §9.4) e non va usato. Anche il resoconto del run sbagliava il conteggio ("MM più vicino in 7, HM in 4"): usava la colonna `MM_minus_HM_cm`, che è la differenza dell'errore medio (fondo MM più basso o più alto), non dell'RMSE. Dal 19/11 `compare_modules_gnss.m` esclude le parità (`tie_tol` = 1 mm) e riporta le due differenze in colonne separate.

La mappa delle differenze MM − HM (1 740 047 pixel corretti) ha media −0.001 m e range [−0.24, +0.16] m; il 72.6% dei pixel ha |dZ| ≤ 5 mm e il P95 di |dZ| è 0.06 m. Le differenze si concentrano nei tratti a canali multipli (tratti 1–2, 4–13, 21–27; rami 31–32, 34/36, 38–39); la differenza media più grande è nel tratto 10 (−18.6 cm), dove la portata di MM è doppia di quella di HM.

**Residui negativi forti e $z_{DTM}$.** I residui di s006 (−1.47 m), s019 (−0.77) e s020 (−0.45) non cambiano: non dipendevano dalla portata. Alle sezioni assegnate a queste stazioni $z_{DTM}$ è più bassa del DTM completo sull'asse di 1.56 m (170), 1.18 m (1160) e 0.61 m (1176), a sezioni alterne: dove $z_{DTM}$ è giusta coincide con il DTM completo entro 2–3 cm.

Il controllo alla sezione 170 (`check_zDTM_section170.png`, `log_zDTM_check_170*.txt`) mostra il meccanismo:

1. Il DTM di segmento `dtm_008` contiene 141 870 pixel con il nodata −32767, non convertito in NaN. Alle sezioni 169 e 171 l'interpolazione bilineare vicino a questi pixel dà valori positivi ma falsi (211.91 e 107.50 m), che il filtro sui valori negativi non toglie e che abbassano la media di 1.32 e 1.84 m.
2. La quota grezza della 170 è invece giusta (255.47 m).
3. Il filtro sulle oscillazioni isolate segnala 169–172 e sostituisce ogni sezione con la media delle due vicine grezze: la 170 diventa (254.23 + 253.57)/2 = 253.90 m. Il filtro sposta così l'errore sulle sezioni sane.

La causa 1 della §9.6 è confermata; le cause 2 (riempimento: DTM di segmento uguale al DTM completo entro ±2 cm) e 3 (media lungo la sezione: media e mediana differiscono di 1 cm) non lo sono alla 170. Le sezioni 1160, 1162, 1174 e 1176 mostrano lo stesso schema, ma non sono state verificate pixel per pixel. Una nota minore: la griglia di `interp2` è costruita con `linspace` sui limiti della maschera (passo 0.5016 m invece di 0.5 m), con uno sfasamento fino a circa 0.25 m ai bordi, trascurabile su sezioni piatte.

**Conclusioni.**

1. Con la rete corretta la portata si conserva ovunque (punti 6 e 7 chiusi) e l'errore medio scende a +0.12 m in entrambi i moduli.
2. I due moduli danno lo stesso RMSE (0.36 m, 0.30 m senza s006) e la loro differenza non è significativa. Differiscono solo nelle 11 stazioni fuori dal canale singolo, troppo poche per un confronto per classe.
3. L'errore dominante resta $z_{DTM}$ (nodata dei DTM di segmento più filtro sui salti): va corretto prima di ripetere il confronto.
4. Il contributo del modulo morfodinamico resta la coerenza fisica della ripartizione ai nodi, non l'accuratezza del fondo (§9.5, conclusione 6); al volo agisce solo tramite $\psi$, perché il trasporto solido è nullo.

### 9.9 Dati di validazione: GNSS e ADCP (19/11/2026)

**Due modalità di rilievo, una sola validazione.** I dati di campo del 07/03/2025 vengono da due modalità: sezioni trasversali rilevate a guado con antenna GNSS (`punti_bat.csv`, 48 sezioni) e transetti ADCP (`punti_bat_ADCP.csv`, 8 transetti, 1024 punti). Nella validazione non si distinguono: ogni sezione GNSS e ogni transetto ADCP è una stazione, confrontata con la sua quota media come descritto nella §9.5. Nel paper si dichiara solo che i dati provengono dalle due modalità. Le stazioni passano da 48 a 56 (le sottostazioni qui sotto ne aggiungono altre).

**Lettura dei dati ADCP** (`read_field_points.m`). Il file ha le colonne `N_UTM32N, E_UTM32N, Z_BOTTOM` e un quarto campo compilato solo sulla prima riga di ogni transetto (per esempio `20250307151356_-1sec`): da lì inizia un nuovo transetto. I transetti diventano le stazioni a01–a08, nell'ordine del file. Sono lunghi 14–41 m, nella parte di valle della rete, vicino alle stazioni GNSS s009, s010, s015, s016, s022 e s023.

**Quota di riferimento.** $Z_{BOTTOM}$ è una quota ellissoidica, lo stesso riferimento del DTM. Verifica del 19/11/2026: nei transetti 3, 4 e 7 la quota ADCP differisce di 4–9 cm da quella ellissoidica dei punti GNSS più vicini (28–43 m), mentre la quota ortometrica è più bassa di circa 48 m.

**Copertura dei transetti (19/11/2026).** I transetti ADCP coprono soprattutto le parti profonde della sezione, dove non si poteva guadare; le sezioni GNSS vanno invece da sponda a sponda. Nel run del 19/11 (§9.10) le stazioni ADCP hanno un errore medio 0.3–0.6 m più alto delle stazioni GNSS vicine (tratto 30: +0.35/+0.61 contro +0.08/+0.20; tratto 35: +0.44/+0.74 contro +0.05/+0.18), anche se la quota dei singoli punti concorda entro 4–9 cm. La validazione confronta il fondo corretto, quasi uniforme sulla sezione, con la media della stazione: un transetto che misura solo la parte profonda ha una media più bassa della media di sezione, quindi un residuo più alto per effetto della copertura, non della quota. Con gli ADCP l'errore medio su tutte le stazioni sale da +0.22 a +0.33 m. Va deciso come usarli (punto 5 della §9.4).

**Densità dei punti (in sospeso).** L'ADCP registra un punto ogni 0.1–0.3 m e resta fermo qualche secondo all'inizio dei transetti (49–271 punti per transetto); il GNSS ha un punto ogni 1.5 m circa (mediana 20 punti per stazione). Nelle statistiche calcolate su tutti i punti insieme un transetto ADCP pesa quindi molto più di una stazione GNSS. Il parametro `adcp_min_spacing` del config dirada i punti ADCP; per ora vale 0 (tutti i punti), in attesa della risposta di Carlo (`docs/Domande_per_Carlo.md`, domanda 6).

**Sezioni GNSS su più tratti.** Alcune sezioni GNSS attraversano più rami. Come stazione unica, sono confrontate con una media che mescola fondi diversi e assegnate al tratto della maggioranza dei punti. Vengono quindi divise in sottostazioni, una per tratto, ognuna confrontata con la media dei suoi punti. Il tratto è assegnato a mano (`gnss_station_split` nel config) e non dalla sezione della centerline più vicina:

| Stazione | Punti | Tratto | Sottostazione |
| --- | --- | --- | --- |
| s003 | p001–p013 | 9 | s003-t09 |
| s003 | p014–p021 | 11 | s003-t11 |
| s003 | p022–p029 | 12 | s003-t12 |
| s007 | p001–p016 | 6 | s007-t06 |
| s007 | p017–p025 | 7 | s007-t07 |

s015 è ancora da definire. s007 attraversa i due canali (tratti 6 e 7), separati da una barra centrale, che confluiscono nel tratto 8 dove si trova s006. Le sottostazioni dello stesso transetto misurano il fondo di rami diversi nello stesso punto della rete: confrontano direttamente il dislivello tra i rami, che dipende dalla ripartizione della portata. Un caso simile, ma con stazioni già separate, è quello di s046 (tratto 24) e s047 (tratto 26), che confluiscono nel tratto 27 dove si trova s048.

### 9.10 Run del 19/11/2026: correzione del nodata, ADCP e sottostazioni

Run di pre-processing + modulo idraulico, modulo morfodinamico e `ModuleComparison` con la correzione del nodata (§9.6), i dati ADCP e le sottostazioni (§9.9). Tutti senza errori. Riferimento: run del 18/11 sera (§9.8). Resoconto completo, con le tabelle, in `Final/Output_folder/RESOCONTO.md`.

**Effetto della correzione del nodata su $z_{DTM}$.**
- Salti segnalati dal filtro: da 149 a 2 (sezione 14, ultimo salto di +1.11 m del tratto 0; sezione 1239, tratto 35, stazione a08). Salti grezzi oltre 0.5 m tra sezioni consecutive: da 106 a 4.
- La quota grezza cambia solo verso l'alto; $z_{DTM}$ finale cambia di più di 5 cm in 122 sezioni (media +0.08 m).
- Tratti di `segments_to_smooth` (12, 19, 26, 30, 32, 38, 39): range prima del lisciamento da 4.7–11.3 m a 0.8–1.4 m (30: 5.7 m, ma è la pendenza di 182 sezioni). Rugosità grezza mediana 3.6 cm, come gli altri tratti (2.9 cm): non oscillano più. Il lisciamento cambia il range di pochi centimetri, tranne nel tratto 26 (−0.39 m).
- Risalite verso valle oltre 10 cm rimaste nella $z$ grezza: tratti 29 (senza portata), 35, 34, 37, 0 e 22, di 10–55 cm.
- Sezioni 168–172, 1157–1163, 1173–1178: le sezioni alterne abbassate di 0.6–1.6 m sono sparite; $z_{DTM}$ coincide con il DTM completo sull'asse entro 4 cm (7 cm alla 172).

| Stazione | Errore medio HM, 18/11 | 19/11 |
| --- | --- | --- |
| s006 | −1.47 | −0.21 |
| s019 | −0.77 | +0.27 |
| s020 | −0.45 | −0.04 |

**Rete e nodi.**
- Bilancio di portata chiuso a tutti i nodi in entrambi i moduli; 10.80 m³/s all'uscita 42; 366 e 817 senza portata; 127 559 pixel (6.4%) al DTM grezzo, come prima.
- **1369** torna asimmetrico ($\beta$ da 11.7 a 14.9), con la soluzione scelta ($\psi$ = 0.230) praticamente uguale alla bilanciata (0.235); portata bankfull da 1118 a 638 m³/s.
- **15** cambia soluzione: $\psi$ da 0.41 a 0.82, perché $\psi_{hyd}$ passa da 0.49 a 0.67 con la nuova $z_{DTM}$. Ramo 1 nel modulo morfodinamico da 4.47 a 8.82 m³/s, ramo 2 da 6.33 a 1.98.
- $\beta$ cambia di 2–5 a 205 (13.7 → 15.3), 584 (→ 38.3), 1070 (→ 19.2) e 1293 (→ 15.8): la correzione cambia le pendenze e quindi le grandezze bankfull (punto 1 della §9.4). Tirante critico ancora in 5 madri su 10.

**Validazione** (intervalli al 95% dal bootstrap delle stazioni).

| Insieme | Stazioni | RMSE HM | RMSE MM | Errore medio HM / MM | RMSE MM − HM | MM / HM / pari, p |
| --- | --- | --- | --- | --- | --- | --- |
| 48 GNSS originali (18/11 sera, §9.8) | 48 | 0.36 | 0.36 | +0.12 / +0.12 | +0.4 cm [−0.1, +0.9] | 3 / 8 / 37, p = 0.23 |
| 48 GNSS originali | 48 | 0.30 [0.23–0.36] | 0.30 [0.23–0.37] | +0.22 / +0.22 | +0.4 cm [−0.2, +1.2] | 5 / 7 / 36, p = 0.77 |
| 51 GNSS (con sottostazioni) | 51 | 0.30 [0.23–0.36] | 0.30 [0.24–0.37] | +0.22 / +0.22 | +0.5 cm [−0.1, +1.2] | 6 / 9 / 36, p = 0.61 |
| Tutte (GNSS + ADCP) | 59 | 0.40 [0.32–0.47] | 0.40 [0.32–0.47] | +0.33 / +0.33 | +0.2 cm [−0.0, +0.6] | 6 / 10 / 43, p = 0.45 |

- **Precisione:** sulle 48 stazioni l'RMSE scende da 0.36 a 0.30 m; senza s006 resta 0.30, come il 18/11. Il guadagno viene tutto dalle stazioni con $z_{DTM}$ sbagliata. NMAD da 0.19 a 0.15 m; riduzione dell'RMSE rispetto al DTM grezzo +0.30 [+0.28, +0.33] m (era +0.24 [+0.11, +0.34]).
- **Errore sistematico:** l'errore medio sale da +0.12 a +0.22 m (stazioni con errore medio negativo da 10 a 3). Prima era in parte compensato dalle sezioni abbassate dal nodata; ora è l'errore principale (punto 11 della §9.4).
- **Moduli:** non distinguibili con nessuna combinazione di stazioni. Negli "altri tratti a canali multipli" con le sottostazioni (6 stazioni) la differenza è +4.1 cm [+0.1, +5.8], a favore di HM, ma su 6 stazioni.
- **ADCP:** residuo 0.3–0.6 m più alto delle GNSS vicine, per la copertura dei transetti (§9.9). Un transetto ADCP pesa 5–6 stazioni GNSS anche nelle statistiche di gruppo di `ModuleComparison`, che calcolano RMSE ed errore medio sui punti di ogni campione bootstrap.
- Mappa MM − HM: range [−0.19, +0.12] m, 69.6% dei pixel entro 5 mm, P95 0.071 m.

**Sottostazioni: dislivello tra i rami dello stesso transetto** (m).

| Coppia | GNSS | DTM grezzo | HM | MM |
| --- | --- | --- | --- | --- |
| s003: t09 − t12 | +0.40 | +0.11 | −0.14 | −0.05 |
| s003: t11 − t12 | +0.63 | +0.26 | +0.14 | +0.08 |
| s007: t06 − t07 | +0.02 | +0.60 | +0.63 | +0.75 |

Nessuno dei due moduli riproduce i dislivelli misurati: il ramo 12 di s003 è il più basso nei dati ma non nei modelli (residuo +0.46 / +0.43 m), e a s007 i modelli tengono il gradino di 0.6 m del DTM grezzo (s007-t06 +0.52 / +0.57 m). Il dislivello modellato è dominato da $z_{DTM}$ dei rami, non dalla ripartizione della portata.

**Anomalie minori:** nelle figure i nomi con "_" diventavano pedici TeX; i titoli dicevano "GNSS" anche con i dati ADCP; la tabella diagnostica dei tratti sospetti nel log del pre-processing scriveva "OK - oscillation, apply smoothing" per tutti i tratti, anche quelli che non oscillano più e il 35, che non è in `segments_to_smooth`. Corrette nel codice il 20/11/2026 (punto 17 della §9.4).

**Lisciamento (valutazione del 20/11/2026, su `z_matrix` di `PRE_workspace_2025.mat`).** Rugosità = deviazione standard della $z$ grezza attorno alla sua mediana mobile su 7 sezioni.

| Tratto | In `segments_to_smooth` | Rugosità (m) | Salto massimo grezzo (m) | Cambio massimo prodotto dal lisciamento (m) |
| --- | --- | --- | --- | --- |
| 26 | sì | 0.070 | 0.39 | 0.38 (un picco isolato alla terza sezione: 244.95 contro 244.59 e 244.56) |
| 32 | sì | 0.043 | 0.26 | 0.26 |
| 12 | sì | 0.087 | 0.42 | 0.21 |
| 39 | sì | 0.036 | 0.18 | 0.13 |
| 30 | sì | 0.007 | 0.17 | 0.07 |
| 19 | sì | 0.023 | 0.71 (un solo gradino) | 0.05 |
| 38 | sì | 0.008 | 0.11 | 0.04 |
| 34 | no | 0.083 | 0.41 | — |
| 35 | no | 0.085 | 0.55 | — (0.53 dal filtro sui salti, sezione 1239) |
| 37 | no | 0.062 | 0.25 | — |
| 0 | no | 0.301 | 1.11 | — (0.55 dal filtro sui salti, sezione 14) |
| altri tratti | no | 0.005–0.058 (mediana 0.03) | 0.01–0.30 | — |

Con il nodata corretto nessun tratto oscilla più: il lisciamento della lista toglie al massimo tre picchi isolati di 0.2–0.4 m (tratti 26, 32, 12) e sposta gli altri tratti di pochi centimetri. Picchi della stessa entità restano nei tratti 34, 35 e 37, che non sono nella lista. La lista manuale non segue quindi più un criterio. Due possibilità: toglierla (si accettano picchi isolati fino a circa 0.4 m, con effetto solo locale sul fondo), oppure sostituirla con una regola uguale per tutti i tratti, per esempio il filtro sui salti isolati con una soglia più bassa di 1 m. Il tratto 0 resta un caso a parte (salto di +1.1 m in fondo al tratto di testa).

**Scelta (20/11/2026): regola uguale per tutti i tratti** (§9.6, "Regola per i salti isolati"). I salti massimi dei tratti 12, 32, 34 e 37 sono gradini verso valle e non vengono toccati.

**Conclusioni.**
1. La causa 1 della §9.6 era quella dominante: corretto il nodata, $z_{DTM}$ non ha più sezioni abbassate né oscillazioni.
2. Il profilo monotono e `segments_to_smooth` probabilmente non servono più per le oscillazioni; restano da valutare le poche risalite di 10–55 cm e il salto di +1.1 m in fondo al tratto 0 (punto 8).
3. L'errore principale è ora l'errore medio positivo (+0.22 m sulle GNSS): il fondo corretto è più alto di quello misurato (punto 11).
4. I due moduli restano indistinguibili; il nodo 15 mostra quanto la soluzione dipenda da $\psi_{hyd}$ (punto 4).

### 9.11 Errore medio positivo: analisi del 20/11/2026

Analisi sui dati del run del 19/11 (51 stazioni GNSS, sottostazioni comprese; modulo idraulico). Per ogni stazione si confronta il tirante del modello $h$ alla sezione di riferimento con la profondità osservata $d_{obs}=z_{DTM}-\bar z_{GNSS}$, cioè il pelo libero del volo meno la media GNSS della stazione. Il residuo della stazione è $h$ troppo piccolo: residuo $=d_{obs}-h$.

**Il tirante del modello varia poco, la profondità osservata molto.** $h$ va da 0.10 a 0.60 m (mediana 0.33), $d_{obs}$ da 0.16 a 1.16 m (mediana 0.53). Il residuo ha correlazione 0.91 con $d_{obs}$ e −0.20 con $h$: l'errore positivo è concentrato nelle stazioni più profonde, non distribuito in modo uniforme.

**Non basta una scabrezza o una portata diversa.** Per avere $h=d_{obs}$, con il moto uniforme ($h\propto (nQ)^{3/5}$), $n$ o $Q$ andrebbero moltiplicati per 1.72 in mediana, ma con un intervallo interquartile 1.25–3.23. Un unico fattore toglierebbe l'errore medio, ma renderebbe troppo profonde le stazioni oggi già corrette o negative (s006, s007-t07, s003-t09, s044).

**Le stazioni più profonde sono in tratti a pelo libero quasi piatto.** La pendenza locale del pelo libero (retta su ±3 sezioni, circa 60 m) ha correlazione −0.41 con il residuo (−0.50 con il rapporto tra pendenza locale e pendenza del tratto usata dal modello).

| Terzo delle stazioni per residuo | Residuo (m) | $d_{obs}$ (m) | $h$ (m) | Pendenza locale | Pendenza del tratto | Larghezza (m) |
| --- | --- | --- | --- | --- | --- | --- |
| basso | +0.02 | 0.38 | 0.36 | 3.0·10⁻³ | 2.7·10⁻³ | 26 |
| medio | +0.15 | 0.49 | 0.33 | 2.8·10⁻³ | 3.4·10⁻³ | 35 |
| alto | +0.38 | 0.66 | 0.30 | 1.3·10⁻³ | 3.4·10⁻³ | 39 |

Nel terzo con il residuo più alto la pendenza locale è meno della metà di quella del tratto in 9 stazioni su 17 (3 su 17 nel terzo basso), e le sezioni sono più larghe. Sono le pozze a monte delle soglie (sequenza pool–riffle): il pelo libero è quasi orizzontale e l'acqua è profonda e lenta. Il moto uniforme con la pendenza media del tratto dà lì un tirante troppo piccolo, e il fondo risulta troppo alto. Esempio: nel tratto 6 (s007-t06) il pelo libero è piatto tra le sezioni 127 e 135 (256.45–256.55 m) e scende di 0.6 m entro la sezione 140; il modello dà 0.21 m di tirante, i GNSS 0.74 m.

**Livello idrometrico.** A San Benigno (lettura giornaliera delle 12:00) il livello è 86 cm il 7/3 (rilievo GNSS) e 97 cm il 17/3 (volo), con un colmo di 104 cm il 16/3. Il confronto non dipende direttamente dalla data del rilievo: il modello usa il pelo libero e la portata del volo, e $d_{obs}$ è riferita al pelo libero del volo; conta solo se il fondo è cambiato tra le due date. Resta però da verificare la portata del volo: §9.12.

**Con e senza pozze (20/11/2026).** Criterio fissato prima del confronto: una stazione è "in pozza" se la pendenza locale del pelo libero (±3 sezioni) è meno della metà della pendenza del tratto usata dal modello. Errore medio per stazione del modulo idraulico, 50 stazioni GNSS con sottostazioni (s015 esclusa, §9.13), dati del 19/11 (uguali al 20/11 entro 1 mm):

| Gruppo | Stazioni | Errore medio (media per stazione) | Mediana | RMSE sui punti |
| --- | --- | --- | --- | --- |
| in pozza | 16 | +0.34 m | +0.28 m | 0.46 m |
| non in pozza | 34 | +0.14 m | +0.13 m | 0.22 m |
| tutte | 50 | +0.20 m | — | — |

La differenza tra i due gruppi è significativa (test di Mann–Whitney, p = 0.008). Le pozze spiegano quindi gran parte dell'errore medio, ma non tutto: anche fuori dalle pozze il fondo resta 0.14 m più alto di quello misurato. Questa parte residua va guardata nell'analisi di sensitività (scabrezza, larghezze della maschera). Nel paper il confronto serve come diagnostica del limite delle pozze; la validazione principale resta su tutte le stazioni.

**Conclusioni.**
1. L'errore medio positivo non è un errore uniforme: viene soprattutto dalle pozze, dove il moto uniforme con la pendenza del tratto sottostima il tirante. È un limite del metodo (moto uniforme, nessun rigurgito) da dichiarare.
2. La portata del volo è confermata (§9.12): l'errore medio positivo non dipende da lei. Scabrezza e larghezze entrano nell'analisi di sensitività.
3. Decisione del 20/11/2026: l'errore nelle pozze e nelle sezioni molto irregolari è un limite del metodo, che ricostruisce un fondo rettangolarizzato con il moto uniforme; si dichiara nel paper e non si corregge (punto 16 della §9.4).

### 9.12 Portata del volo e scala di deflusso (20/11/2026)

**Conclusione.** La portata del volo è confermata: 10.76 m³/s, la portata giornaliera pubblicata da ARPA Piemonte per San Benigno il 17/03/2025. Corrisponde al livello del volo (0.97 m) con la scala di deflusso valida nel 2025. Nel config si passa da 10.8 a 10.76 m³/s, una differenza dello 0.4%, trascurabile sui tiranti. L'errore medio positivo quindi non dipende dalla portata.

**Scala valida al momento del volo.** ARPA pubblica due scale per il 2025 (`scala di deflusso_SAN BENIGNO ORCO_2025.csv`): una valida dal 01/01 al 16/04/2025 e una dal 17/04/2025. Per il volo del 17/03 vale la prima. Interpolata, dà:

| Livello a San Benigno | Data | Portata, scala 01/01–16/04/2025 (m³/s) | Portata, scala 2023 (m³/s) |
| --- | --- | --- | --- |
| 0.86 m | 07/03/2025 (rilievo GNSS) | 5.3 | 15.7 |
| 0.97 m | 17/03/2025 (volo) | **10.9** | 23.4 |
| 1.04 m | 16/03/2025 (colmo) | 15.5 | 29.1 |

La scala 2023 non vale per il 2025: a parità di portata il livello è circa 0.2 m più alto (per esempio 9.5 m³/s: 0.95 m nel 2025, circa 0.75 m nel 2023), probabilmente per il cambiamento del fondo alla stazione dopo le piene del 2024. La stima di 23.4 m³/s ricavata in un primo momento con la scala 2023 è quindi da scartare.

**Portata letta o portata dalla scala.** La portata pubblicata da ARPA non è misurata in continuo: è il livello misurato trasformato con la scala di deflusso valida in quel periodo. Leggere la portata pubblicata o applicare la scala al livello dà quindi lo stesso risultato, a meno di piccole differenze: la portata giornaliera è la media delle portate durante il giorno, mentre la scala applicata al livello giornaliero dà la portata di quel livello, e la scala non è lineare. Qui le due strade danno 10.76 e 10.9 m³/s. Si usa il dato giornaliero pubblicato, che è il valore validato da ARPA. L'ora del volo non è nota, quindi non si possono usare i dati orari.

**`level_upstream`.** Il valore 0.695 m del config era stato ricavato da 10.8 m³/s con una retta di regressione sulla scala 2023 (`FORECAST.LINEAR`), che non rappresenta una scala non lineare. Nel codice serve solo come valore iniziale del tirante per `fsolve` alle biforcazioni del modulo idraulico, quindi non ha influito sui risultati. Si sostituisce con il livello letto, 0.97 m, documentandone l'uso.

**Nota sul file `SAN BENIGNO ORCO_march2025.csv`.** Nella colonna delle portate i valori non corrispondono ai livelli dello stesso giorno: per esempio il 17/03 riporta 6.48 m³/s, mentre 10.76 m³/s compare sul 15/03; con la scala 2025 solo il 16/03 torna (15.3 contro 15.5 m³/s). La colonna delle portate va riscaricata da ARPA prima di usarla.

**Conseguenze.** La portata resta uno dei parametri dell'analisi di sensitività, con l'incertezza tipica della scala di deflusso (circa ±10–20%) e non con l'intervallo 10.8–23.4 m³/s. L'errore medio positivo resta spiegato soprattutto dalle pozze (§9.11, limite del metodo); scabrezza e larghezze entrano nella sensitività.

### 9.13 Run del 20/11/2026: regola unica per i salti isolati, s015, portata 10.76 m³/s

Run dei tre script con la regola unica per i salti isolati (§9.6), s015 assegnata al tratto 30, i nomi `sNNN-tCC`, i titoli "in situ" e $Q$ = 10.76 m³/s (§9.12). Nessun errore. Riferimento: run del 19/11 (§9.10), salvato in `Output_folder/_run_19-11/`. Resoconto completo: `Final/Output_folder/RESOCONTO.md` (in locale).

**Regola per i salti isolati.** Corregge esattamente le 4 sezioni previste: 14 (+0.55 m), 492 (−0.38), 728 (+0.38), 1239 (+0.53). In tutte le altre 1560 sezioni $z_{DTM}$ coincide con la quota grezza. I gradini verso valle dei tratti 12, 32, 34 e 37 non sono toccati; i picchi dei tratti 26 e 35 sono tolti. Rispetto al 19/11 $z_{DTM}$ cambia di più di 5 cm in 22 sezioni, tutte nei tratti prima lisciati (12, 19, 30, 32, 39); il fondo del modulo idraulico segue $z_{DTM}$, i tiranti cambiano meno di 1 cm.

**Due casi locali, senza stazioni.**

- **Sezione 492 (tratto 19).** Non è un picco ma il piede di un gradino di −0.71 m: le sezioni 489–494 valgono 250.16, 250.14, 249.86, 249.15, 249.20, 249.13. La 492 è appena 5 cm sotto la 493, quindi è un minimo locale e la regola la porta a metà gradino (249.53 m, fondo +0.33 m).
- **Fine del tratto 0 (sezioni 13–15).** Il tratto di testa è quasi piatto (260.03–260.13 m nelle sezioni 1–12), poi risale a 260.26, 261.36 e 261.36 m nelle ultime tre; il tratto 1, a valle del nodo 15, parte da 259.88 m. La regola corregge la 14 (a 260.81 m) ma non la 15, che è l'ultima sezione del tratto e la copia della 14 (bordo di `calculate_widths`). Il tratto 0 ha pendenza minima ($1.6\cdot10^{-4}$) in tutto il modulo idraulico, perché la retta sul pelo libero risale, e il nodo 15 ricava $S_0$ allargando la finestra. La causa più probabile è che i transetti al nodo attraversino la barra tra i due rami (limite della media su tutta la sezione, §9.6). Da verificare.

**Tabella diagnostica.** Esiti coerenti con rugosità e salti: nessuna oscillazione, tranne il tratto 0, segnato "TO BE INSPECTED MANUALLY" per il salto a fine tratto, non per un'oscillazione.

**Nomi, figure, nodi, portata.** Sottostazioni con il trattino, figure senza pedici e con "in situ"; s015-t30 sul tratto 30 (`manual`); resta ambigua solo s035. Bilancio chiuso, uscita 10.76 m³/s. Nessun nodo cambia tipo o soluzione; 244 cambia di poco ($W_1/W_0$ da 0.66 a 0.63) per la sezione 256, dove il profilo ora ha una sola sponda e viene specchiato (nuovo warning; $W_2^{dtm}$ da 78.3 a 88.1 m). La portata 10.76 riduce i tiranti dello 0.2%.

**Validazione: identica al 19/11 entro 1 mm.**

| Insieme | Stazioni | RMSE HM / MM (m) | Errore medio HM / MM (m) | RMSE MM − HM | MM / HM / pari, p |
| --- | --- | --- | --- | --- | --- |
| A: GNSS originali | 48 | 0.30 / 0.30 | +0.22 / +0.22 | +0.5 cm [−0.1, +1.2] | 5 / 7 / 36, p = 0.77 |
| B: GNSS con sottostazioni | 51 | 0.30 / 0.31 | +0.22 / +0.22 | +0.5 cm [−0.0, +1.2] | 6 / 9 / 36, p = 0.61 |
| C: tutte, ADCP compresi | 59 | 0.40 / 0.40 | +0.34 / +0.34 | +0.2 cm [+0.0, +0.6] | 6 / 10 / 43, p = 0.45 |

Riduzione dell'RMSE rispetto al DTM grezzo su A: +0.30 [+0.28, +0.33] m. Con s015 nel canale singolo le classi di A diventano 38 / 4 / 1 / 5 stazioni.

**Errore medio positivo.** L'analisi della §9.11 è confermata: correlazione del residuo 0.90 con la profondità osservata, −0.42 con il rapporto tra pendenza locale e pendenza del tratto; nel terzo con il residuo più alto la pendenza locale è meno della metà di quella del tratto in 9 stazioni su 17.

**s015-t30.** La sezione di riferimento è la 1069, alla fine del tratto 30 (nodo 1070), dove il tirante è 0.64 m; il 63% dei punti cade sui tratti 31 e 32. L'errore medio della stazione (+0.03 m, calcolato con il DTM corretto nei punti) è affidabile, mentre il confronto tirante–profondità alla sezione di riferimento (−0.16 m) non lo è: s015 va esclusa dalle analisi che usano la sezione di riferimento (§9.11).

**Conclusioni.**

1. La regola per i salti isolati funziona: il punto 8 si può chiudere dopo aver deciso i due casi locali (492 e fine del tratto 0).
2. Il run è il riferimento per l'analisi di sensitività: tutte le correzioni sono dentro e le statistiche sono stabili.

### 9.14 Run del 20/11/2026 (sera): maschera del tratto 0 e regola più stretta

Input: maschera e DTM del tratto 0 accorciati prima dell'ombra del ponte (`I_000.tif` int16 0/1); regola per i salti isolati con la condizione su entrambe le vicine (§9.6); biforcazione 0 → 1, 2 aggiunta a `forced_bifurcations`, perché il tratto 0 accorciato finisce a 61.5 e 67.6 m dai tratti 2 e 1, oltre `distance_threshold` = 51 m. Il codice è quello del 20/11, a parte i commenti. Il tratto 0 ha 13 sezioni invece di 15, quindi tutti gli ID a valle scendono di 2.

**Tratto 0.** Asse da 134.6 a 113.6 m. Ultime 5 sezioni: z grezzo 260.10, 260.06, 260.02, 260.02, 260.02 m (prima: 260.10, 260.08, 260.26, 261.36, 261.36). La pendenza del tratto è quella della retta, 2.0·10⁻⁴ (prima era portata a `slope_min` = 1.6·10⁻⁴ perché la retta saliva verso valle). Il fondo idraulico non risale più a fine tratto.

**Salti isolati.** Corrette solo 726 (ex 728, tratto 26: +0.36 e +0.39 m dalle vicine) e 1237 (ex 1239, tratto 35: +0.55 e +0.51 m). La 490 (ex 492) resta a 249.15 m: la pendenza del tratto 19 passa da 0.01132 a 0.01146, il fondo idraulico alla 490 scende di 0.38 m, nessun tirante cambia di più di 3 mm.

**Nodo 13 (ex 15).**

| | 20/11 (nodo 15) | 20/11 sera (nodo 13) |
| --- | --- | --- |
| $S_0$ al nodo | 2.60·10⁻³ (valore di ripiego) | 1.22·10⁻³ (finestra di 400 m alzata alla soglia) |
| $D_0$ | 1.18 m | 2.53 m |
| $W_0$ bankfull | 74.8 m | 55.1 m |
| $\beta$ | 31.6 | 10.9 |
| $\psi_{hyd}$ / $\psi$ | 0.67 / 0.82 | 0.67 / 0.73 |
| Soluzioni / tipo | 3 / asimmetrica | 1 / bilanciata |
| $Q_1$ / $Q_2$ nel modulo morfodinamico | 8.79 / 1.97 m³/s | 7.91 / 2.85 m³/s |

Il nodo è la testa della rete: la finestra a monte si ferma all'inizio del tratto 0, che è quasi piatto (rigurgito del ponte), e la pendenza è sotto la soglia minima. Con la $S_0$ più bassa il tirante al nodo raddoppia, $\beta$ scende sotto il valore critico e resta solo la soluzione bilanciata. Il fondo del modulo morfodinamico cambia di +1.6 cm sul tratto 1 e −4.6 cm sul tratto 2 (massimo 5.3 cm); il modulo idraulico non cambia. Gli altri nodi sono identici al 20/11. Punto 21 della §9.4.

**Validazione.** Identica al 20/11 entro 1 mm, perché nessuna stazione cade nei tratti 0–3 né sulla sezione 490: 48 GNSS RMSE 0.30 m, errore medio +0.22 m, riduzione +0.300 [+0.278, +0.325] m; 59 stazioni RMSE 0.40 m, +0.34 m; MM − HM +0.5 cm [−0.1, +1.2], test del segno 5 / 7 / 36, p = 0.77. Mappa MM − HM: P95 da 0.071 a 0.066 m.

**Sezioni di estremità.** In ogni tratto la prima e l'ultima sezione copiano quota e larghezza della vicina (`calculate_widths`, "Edge case handling"): le ultime due quote uguali del tratto 0 vengono da qui.

**Conclusioni.** Punti 8 e 19 chiusi. Resta da decidere il trattamento di $S_0$ al nodo di testa (punto 21). Il run della sera sostituisce quello del 20/11 come riferimento per la sensitività.

---

## Registro modifiche

| Data | Modifica |
| --- | --- |
| 09/10/2026 | Conversione da `.docx` a Markdown; equazioni ricostruite in LaTeX. |
| 09/10/2026 | Spostato nel repo GitHub (`docs/`): da qui in poi la versione di riferimento è questa. |
| 09/10/2026 | Aggiunta la sezione 7 "Questioni aperte" e la sezione 8 "Domande per Carlo". |
| 18/11/2026 | Aggiunta la sezione 9 "Aggiornamento del modulo morfodinamico": modifiche introdotte nel codice il 12/10/2026 e il 18/11/2026, risultati del run del 18/11/2026 e test di sensitività alla larghezza. |
| 18/11/2026 | Aggiunta la §9.5: validazione GNSS del modulo idraulico e del modulo morfodinamico, per tipo di tratto e profilo per profilo. |
| 18/11/2026 | Corretta la §9.5: intervalli di confidenza per stazione invece che per punto; confronto stazione per stazione (il modulo idraulico è leggermente ma sistematicamente più vicino ai GNSS, non indistinguibile); quattro tipi di tratto; portata non conservata nel modulo idraulico e rami senza uscita; origine dei forti residui negativi in $z_{DTM}$. |
| 18/11/2026 | La §9.4 diventa la tabella di tutti i punti in sospeso (voci 6–13 nuove). Corretta nel codice la propagazione della portata del modulo idraulico (punto 6). |
| 18/11/2026 | Aggiunte la §9.6 (quota di partenza $z_{DTM}$: problema, cause possibili, proposta) e la §9.7 (rami 16–15 e 29); aggiornati i punti 7 e 8 della §9.4; nota in testa alla §7 che rimanda alla §9.4. |
| 18/11/2026 | §9.7 e punto 7 della §9.4: aggiunte le osservazioni sul campo sui rami 16–15 e 29 (quasi fermi al volo, attivi a portate più alte). |
| 18/11/2026 | §9.7: scelte adottate per i rami 15, 16 e 29 (portata nulla, fondo non ricostruito, nessun sistema del nodo), implementate in `Final`. §9.4: aggiornati i punti 7, 10 e 12. |
| 19/11/2026 | Aggiunta in testa la sezione "Stato del lavoro" (resoconto per giornata). Aggiunta la §9.8: risultati dei run del 18/11 sera con la rete corretta (bilancio di portata, nodi, validazione GNSS, confronto dei moduli, verifica di $z_{DTM}$ alla sezione 170). §9.4: chiusi i punti 6 e 7, aggiornati 1, 2, 4, 8, 10, 11, 13; nuovi punti 14 (parità nel test del segno) e 15 (trasporto al volo nullo). Nota in testa alla §9.5 e verifica in §9.6. |
| 19/11/2026 | Tolta la §8 "Domande per Carlo": le domande aperte sono in `docs/Domande_per_Carlo.md`. Aggiunte la §9.9 (dati ADCP nella validazione, sottostazioni di s003 e s007) e in §9.1 il calcolo della larghezza bankfull al nodo verificato sul codice. §9.4: aggiornati i punti 5, 9 e 10. Codice: `read_field_points.m` (nuova), `validation_gnss_mean.m`, `station_reach_class.m`, `gnss_cross_section_profiles.m`, i tre script principali e il config (`adcp_file`, `adcp_min_spacing`, `gnss_station_split`). |
| 19/11/2026 | Correzione di $z_{DTM}$, primo passo: nodata dei DTM di segmento convertito in NaN (`Preprocessing_HydraulicModule.m`), media della sezione sui soli valori non NaN (`calculate_widths.m`). §9.6 ("Correzione adottata" e limite delle isole nella media), §9.4 punto 8, stato del 19/11. |
| 19/11/2026 | Corretto `compare_modules_gnss.m` (parità escluse dal test del segno con `tie_tol` = 1 mm, colonne `dRMSE_MM_HM_cm` e `dME_MM_HM_cm`): §9.4 punto 14, §9.8, stato del 19/11. Corretto in §9.8 il conteggio delle stazioni (MM più vicino in 3, HM in 8). Aggiunti alla notazione (§1) $\beta$, $W_1/W_0$, $\psi_{hyd}$, $\psi_x$, $\theta^*_c$, $\theta^*_{BF}$, $\theta_{bf}/\theta^*_c$, con la guida alla lettura della tabella dei nodi in §9.8. La §8 diventa 8.1 (ottobre) e 8.2 (domande del 19/11/2026). |
| 20/11/2026 | Run della sera (§9.14): punti 8 e 19 chiusi, nuovo punto 21 ($S_0$ al nodo di testa 13); biforcazione 0 → 1, 2 forzata nel config; commenti "Chezy's eq" corretti nel modulo idraulico. |
| 20/11/2026 | Notazione della scabrezza allineata al codice: solo $n$ di Manning (Gauckler–Strickler $k_s=1/n$), Chézy tolto dal testo; legame $C=R^{1/6}/n$ nel §3. |
| 20/11/2026 | Maschera del tratto 0 accorciata da Elisa e riesportata come int16 0/1 (codice invariato). |
| 20/11/2026 | Punto 19: la maschera del tratto 0 finisce sull'ombra di un ponte, Elisa la accorcia; s035 resta sul tratto 20; punto 9 chiuso. |
| 20/11/2026 | Decisioni dopo il run: regola per i salti isolati più stretta (§9.6, codice); confronto con e senza pozze (§9.11); §9.4: punti 1, 8, 9 e 19 aggiornati (ponte a valle del tratto 0), nuovo punto 20 (baseline semplice). |
| 20/11/2026 | Run del 20/11: nuova §9.13 e stato del 20/11 (run); §9.4: punto 17 chiuso, punti 8 e 9 aggiornati, nuovo punto 19 (fine del tratto 0 al nodo 15). |
| 20/11/2026 | §9.12 riscritta: portata del volo confermata (10.76 m³/s, ARPA, coerente con la scala 01/01–16/04/2025); la stima di 23.4 m³/s con la scala 2023 è scartata. Config: `Q_upstream` = 10.76, `level_upstream` = 0.97. Punto 18 chiuso, punto 11 aggiornato. |
| 20/11/2026 | Nuova §9.12 (portata del volo e scala di deflusso 2023); §9.11 rimanda alla §9.12 e registra la decisione sulle pozze; §9.4: punto 8 deciso, punto 11 aggiornato, nuovi punti 16 (limite: pozze e sezioni irregolari) e 18 (portata del volo). |
| 20/11/2026 | Regola unica per i salti isolati al posto di `segments_to_smooth` e del filtro a 1 m (§9.6, §9.10); s015 assegnata al tratto 30; nuova §9.11 (analisi dell'errore medio positivo, livello idrometrico); §9.4 punti 8, 9, 11; domanda 6 per Carlo riscritta (copertura ADCP). |
| 20/11/2026 | Stato del 20/11. Codice: sottostazioni `sNNN-tCC`, titoli ed etichette "in situ", tabella diagnostica dei tratti sospetti su $z_{DTM}$ grezza. §9.10: anomalie minori corrette e valutazione del lisciamento (tabella dei tratti). §9.4: punto 8 (prossimo passo), nuovo punto 17. |
| 19/11/2026 | Run del 19/11 con la correzione del nodata, gli ADCP e le sottostazioni: nuova §9.10; "Stato del lavoro" del 19/11 (run); §9.9, copertura dei transetti ADCP (soprattutto le parti profonde, non guadabili); §9.6, esito della correzione. §9.4: chiuso il punto 14, aggiornati 1, 4, 5, 8, 9, 10, 11, 12, 13. Resoconto completo in `Final/Output_folder/RESOCONTO.md`. |
