# Ricostruzione batimetrica automatica – Approccio teorico dei moduli idraulico e morfodinamico

**Elisa Caccamo** · 8 ottobre 2026

> Versione Markdown del file `TeoriaMetodologica_moduli_idraulico_morfodinamico.docx`. Le equazioni sono state ricostruite in LaTeX a partire dal testo estratto; la figura del §4.3 non era presente nel testo ed è indicata come segnaposto.

---

## Stato del lavoro

Resoconto breve per giornata, dal più recente. Il dettaglio è nella §9; i punti aperti sono nella tabella della §9.4.

**19/11/2026.** Analizzati i run del 18/11 sera (pre-processing + modulo idraulico, modulo morfodinamico, `ModuleComparison`) con la rete corretta: portata propagata solo dentro il tratto e tratti 15, 16, 29 senza portata (§9.8).

- **Portata conservata ovunque:** 10.80 m³/s all'uscita in entrambi i moduli, nessun nodo sbilanciato. Punti 6 e 7 della §9.4 chiusi.
- **Nodi:** 10 risolti dal sistema (366 e 817 esclusi); 1369 diventa bilanciato ($\beta$ = 11.7, sotto la soglia), gli altri 9 asimmetrici con $\psi$ entro ±0.02 dal run precedente.
- **Validazione GNSS (48 stazioni):** RMSE da 0.60 (DTM grezzo) a 0.36 m in entrambi i moduli; errore medio da +0.56 a +0.12 m. Senza s006: 0.30 m.
- **I due moduli non sono distinguibili:** differenza di RMSE +0.4 cm [−0.1, +0.9]; MM più vicino in 3 stazioni, HM in 8, 37 pari (p ≈ 0.23). Differiscono solo nelle 11 stazioni dei tratti a canali multipli.
- **$z_{DTM}$ resta l'errore dominante:** confermata alla sezione 170 la causa 1 della §9.6 (nodata −32767 dei DTM di segmento non convertito), con il filtro sui salti che sposta l'errore sulle sezioni sane. Correzione in attesa di decisione.
- **Da decidere ora:** correzione di $z_{DTM}$ (§9.6), gestione delle parità nel test del segno, classi delle stazioni (troppo poche fuori dal canale singolo).

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
| $\psi$ | quota di portata del ramo 1 in condizioni bankfull, $\psi=Q_1^{bf}/Q_0^{bf}$ | – |

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

Il modulo idraulico ricava il tirante $h$ di ogni sezione rettangolare dal moto uniforme (Chezy/Manning) con la portata liquida propagata da monte, e pone $\eta=z_{DTM}-h$. Non c'è trasporto solido: l'unico vincolo di rete è la conservazione della portata liquida, più la continuità dell'energia nei nodi di biforcazione.

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

- Portata liquida: $Q_0=Q_1+Q_2$; il tirante $D_0$ si ricava da Chezy/Manning invertita e $\eta_0=z_{DTM}-D_0$.
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
| Canale singolo | Moto uniforme (Chezy/Manning) → $h$; $Q$ propagata invariata; $\eta=z_{DTM}-h$ | Moto uniforme → $D$; $\tau=\rho g D S$ → $q_{sC}$, $q_{sF}$ (Wilcock & Crowe); $Q$ e $Q_{s,x}$ propagati invariati |
| Biforcazione | 4 eq. in ($h_1$, $h_2$, $Q_1$, $Q_2$): Chezy × 2, $Q_0=Q_1+Q_2$, Bernoulli; condizioni normali | 8 eq. in ($q_1$, $q_2$, $\eta_1$, $\eta_2$, $f_{1C}$, $f_{2C}$, $D_1$, $D_2$) ridotte a 4 in ($q_1$, $\eta_1$, $f_{1F}$, $f_{2F}$): nodo BRT con hiding, continuità liquida e solida, pelo libero comune, $\eta_0=(\eta_1+\eta_2)/2$; condizioni bankfull, poi riscalate; ritorno alle condizioni del volo con $\psi$ |
| Confluenza | $Q_0=Q_1+Q_2$ → Chezy → $h_0$, $\eta_0$ | $Q_0=Q_1+Q_2$ → $D_0$, $\eta_0$; $Q_{s,x}$ sommati; $q_{s0x}=Q_{s,x}/W_0$; $f_F$ di equilibrio con la capacità di trasporto di valle |
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

## 8. Domande per Carlo

1. Al nodo, la pendenza di madre e rami si può calcolare su una finestra di lunghezza fissa $L$ invece che con la regressione sull'intero tratto (§7.2 A)? Se sì, con quale $L$?
2. Quando la regressione dà una pendenza nulla o negativa, è meglio usare la pendenza del tratto adiacente (§7.2 B) o modificare il valore di ripiego $1.6\cdot10^{-4}$ (§7.2 C)? In questo secondo caso, con quale criterio va scelto il nuovo valore?
3. Il tirante bankfull di Shields ($\theta^*_{BF}=1.62\,\theta^*_c$) supera l'altezza delle sponde in molti tratti a canali multipli. È meglio limitarlo con il livello di sfioro del transetto (§7.3 A) o calibrare $\theta^*_{BF}$ sull'Orco (§7.3 B)?
4. Quando la larghezza bankfull non è calcolabile, si può usare quella della maschera e, se si trova una sola sponda, calcolare l'area sul semi-transetto e raddoppiarla (§7.4)?
5. Se a bankfull un ramo si chiude ($q_2\le0$, nodo 106) ma è bagnato il giorno del volo, va considerato un risultato fisico o un segnale di dati in ingresso non affidabili?
6. Per i nodi che restano non risolti, è accettabile la ripartizione proporzionale alle larghezze (§7.6), dichiarando nei risultati che lì le quote di fondo sono idrauliche e non morfodinamiche?

---

## 9. Aggiornamento del modulo morfodinamico (novembre 2026)

Questa sezione raccoglie le modifiche introdotte nel codice il 12/10/2026 (`2_ORCO_bifurcationsNEW_modifiche_modifiche.mlx`) e il 18/11/2026 (`ORCO_bifurcations.m`) e i risultati dei run del 18/11/2026. Dove è in contrasto con le sezioni precedenti (in particolare §4.3, §5.2–5.3, §6, §7 e §8), vale quanto scritto qui; le sezioni precedenti saranno aggiornate in seguito.

### 9.1 Modifiche al metodo

**Pendenza della madre al nodo.** $S_0$ si ricalcola con una retta ai minimi quadrati sul pelo libero $z_{DTM}$ entro $L=400$ m a monte del nodo, proseguendo nel tratto a monte se necessario (a una confluenza, nel ramo con la portata maggiore). Se la retta non è utilizzabile (meno di 5 punti, finestra più corta di 100 m, pendenza nulla o negativa) la finestra si raddoppia fino a 1600 m; se non basta si usa il 25° percentile delle pendenze dei tratti. Una pendenza valida si alza solo se è sotto il 5° percentile. Le pendenze dei rami non sono più dati del nodo.

**Tirante bankfull.** Il tirante di Shields si limita con il livello di sfioro del transetto e non scende mai sotto il tirante critico:

$$
D_0 = \max\Big(\min\big(D_S,\ z_{sfioro}-\eta_0\big),\ D_c\Big), \qquad
D_S = \frac{1.62\,\theta^*_c\,\Delta\,d_{50}}{S_0}, \qquad
D_c = \frac{\theta^*_c\,\Delta\,d_{50}}{S_0}
$$

$z_{sfioro}$ è il più basso dei due punti più alti del terreno, uno per lato, tra $W^m/2$ e $3W^m$ dalla centerline. Dove il tirante si ferma a $D_c$ (sponde basse) la larghezza resta quella misurata al livello di sfioro. Le larghezze dei rami si misurano con il tirante di Shields calcolato con $S_0$; se non sono calcolabili si usa il rapporto delle larghezze delle maschere.

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
| 1 | **Stima della larghezza bankfull**, che controlla il rapporto di forma $\beta$ e quindi il tipo di soluzione al nodo (§9.3) | aperto; confermato dal run del 18/11 sera: a 1369 $\beta$ passa da 12.0 a 11.7 e il nodo diventa bilanciato (§9.8) | sponde cercate solo nel canale attivo che alimenta il nodo |
| 2 | **Sponde basse:** in 9 madri su 12 il tirante bankfull è quello critico (5 su 10 nel run del 18/11 sera); $\theta^*_{BF}=1.62\,\theta^*_c$ (fonte da documentare) | aperto | eventuale calibrazione sull'Orco |
| 3 | **Portata bankfull calcolata localmente** a ogni nodo, non coerente lungo la rete (per esempio 584: 216 m³/s, il nodo 645 su un suo ramo: 481 m³/s) | aperto | — |
| 4 | **Regola di scelta della soluzione:** nei nodi 205, 584 e 645 nessuna soluzione asimmetrica è vicina a $\psi_{hyd}$ | aperto (invariato nel run del 18/11 sera) | — |
| 5 | **Validazione con i dati ADCP** (quella GNSS è nella §9.5) | da fare | — |
| 6 | **Portata non conservata nel modulo idraulico:** il tratto 6 riceveva la portata del tratto 5 e da lì la rete portava 9.99 invece di 10.80 m³/s (§9.5) | **chiuso** (§9.8): bilancio chiuso a tutti i nodi in entrambi i moduli (sbilancio massimo 4.6·10⁻¹² m³/s), 10.80 m³/s all'uscita | — |
| 7 | **Rami 16–15 e 29** (§9.7): nel modello tolgono portata alla rete; al volo sono quasi fermi (29 staccato dal tratto 22, 16–15 parte di un meandro bagnato solo in parte), mentre a portate più alte sono attivi | **chiuso per le condizioni del volo** (§9.8): 366 e 817 "branch without discharge" in entrambi i moduli; 127 559 pixel bagnati (6.4%) restano al DTM grezzo; nessuna stazione GNSS coinvolta | resta aperta, come estensione, la rete alla portata bankfull (limiti del paper) |
| 8 | **Quota di partenza $z_{DTM}$** (§9.6): oscilla e in alcuni tratti è fino a 1.6 m più bassa del DTM completo; causa i forti residui negativi di s006, s019, s020 | causa confermata alla sezione 170 (§9.8): nodata dei DTM di segmento, propagato dal filtro sui salti; stesso schema a sezioni alterne a s019 (1160, 1162) e s020 (1174, 1176) | decidere e implementare la proposta della §9.6 (DTM completo con nodata NaN, mediana sui pixel bagnati, profilo monotono; cosa fare del filtro sui salti e di `segments_to_smooth`); poi rilanciare i tre script |
| 9 | **Stazioni GNSS ambigue** s003, s007, s015, s035 (punti su due tratti) | da verificare | assegnare il tratto a mano sulla mappa `GNSS_station_assignment.png` |
| 10 | **Classi delle stazioni:** con la nuova rete 37 stazioni su 48 sono nel canale singolo; restano 5 a valle di una biforcazione, 1 a valle di una confluenza, 5 negli altri tratti a canali multipli, con intervalli molto larghi (§9.8) | da decidere | due classi (canale singolo, 37, contro canali multipli, 11) oppure differenza tra i moduli in funzione della distanza dal nodo |
| 11 | **Errore sistematico positivo** del DTM corretto: sceso a +0.12 m (per punto) dopo i punti 6 e 7; +0.15–0.16 m senza s006, +0.30 m negli altri tratti a canali multipli (§9.8) | aperto | rivedere dopo il punto 8; poi scabrezza (n = 0.03 costante) e larghezze della maschera nei rami stretti |
| 12 | **Intervalli di confidenza:** `validation_gnss_mean.m` ricampiona i punti, e gli intervalli che ne escono sono troppo stretti; quelli corretti ricampionano le stazioni | implementato in `ModuleComparison.m`: RMSE, errore medio, mediana e NMAD con intervallo per stazione, riduzione dell'RMSE rispetto al DTM grezzo e differenza tra i moduli con intervallo appaiato | nel paper riportare questi intervalli al posto di quelli di `validation_gnss_mean.m` (vedi PaperBAT, 5.5) |
| 13 | **Confronto dei due moduli da ripetere** (§9.5) | ripetuto dopo i punti 6 e 7 (§9.8): differenza non significativa | ripetere dopo il punto 8 |
| 14 | **Test del segno e parità:** `compare_modules_gnss.m` (righe 111–112 e 170–171) conta come vittorie le stazioni con RMSE uguale nei due moduli a meno di 10⁻¹¹ (37 su 48), che finiscono a HM o a MM per rumore di arrotondamento; nel log risulta "MM migliore in 9/48, p = 4·10⁻⁵". La colonna `MM_minus_HM_cm` della tabella delle stazioni è la differenza dell'errore medio, non dell'RMSE | da correggere | escludere le parità con una tolleranza (1 mm: MM più vicino in 3, HM in 8, 37 pari, p ≈ 0.23); rinominare la colonna o aggiungere la differenza di RMSE |
| 15 | **Trasporto solido al volo praticamente nullo:** i flussi solidi di MM alla portata del volo sono 10⁻¹⁶–10⁻²¹ m³/s, quindi le ripartizioni $\psi_x$ al volo agiscono su quantità nulle | da verificare | confermare che è atteso (10.8 m³/s è molto sotto la soglia di moto) e dichiararlo: al volo il modulo morfodinamico agisce solo tramite $\psi$ |

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

**Nodi.** 366 e 817 sono "branch without discharge" in entrambi i moduli. Gli altri 10 nodi sono risolti dal sistema (residuo ≤ 1.9·10⁻¹⁵):

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

**Confronto tra i moduli.** La differenza di RMSE MM − HM su tutte le stazioni è +0.4 cm [−0.1, +0.9]: l'intervallo comprende lo zero, e i due moduli non sono distinguibili (prima +1.3 cm [+0.7, +1.8]). Nelle 37 stazioni del canale singolo i due moduli hanno la stessa portata e lo stesso fondo, quindi lo stesso RMSE. Escludendo queste parità (|ΔRMSE| ≤ 1 mm), MM è più vicino ai GNSS in 3 stazioni (s003 −6.0, s047 −6.5, s015 −0.9 cm) e HM in 8 (s044 +9.4, s046 +8.4, s045 +5.1, s008 +4.6, s007 +4.5, s032 +1.1, s016 +0.9, s006 +0.3 cm): test del segno p ≈ 0.23. Il "9/48, p = 4·10⁻⁵" del log conta le parità numeriche (punto 14 della §9.4) e non va usato.

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
