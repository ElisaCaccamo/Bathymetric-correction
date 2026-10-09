# Ricostruzione batimetrica automatica – Approccio teorico dei moduli idraulico e morfodinamico

**Elisa Caccamo** · 8 ottobre 2026

> Versione Markdown del file `TeoriaMetodologica_moduli_idraulico_morfodinamico.docx`. Le equazioni sono state ricostruite in LaTeX a partire dal testo estratto; la figura del §4.3 non era presente nel testo ed è indicata come segnaposto.

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
>
> **TODO (Elisa):** l'immagine è solo nel `.docx` originale. Caricarla come PNG in `docs/fig/` (oppure passarla a Claude), poi sostituire questo segnaposto con `![Sezione rettangolarizzata al nodo](fig/nodo_bankfull.png)`.

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

- **A. Pendenza del nodo su una finestra di lunghezza fissa** *(proposta principale)*. Per la madre la pendenza si calcola sul pelo libero $z_{DTM}$ di un tratto di lunghezza $L$ a monte del nodo, per ciascun ramo su un tratto di lunghezza $L$ a valle, attraversando se necessario i confini tra tratti (la topologia è nota). $L$ fissa per tutta la rete, ad esempio 300–500 m, oppure proporzionale alla larghezza, $L=k\,W_0^m$ con $k\approx10$. Le pendenze restano "dal DTM" (§2.2) ma sono stimate su una scala confrontabile per madre e rami. Il calcolo nel resto della rete non cambia.
- **B. Sostituire il valore minimo.** Quando la regressione dà $S\le1.6\cdot10^{-4}$ o il tratto è troppo corto, si usa la pendenza del tratto a monte (o a valle) o quella della finestra A, invece di $1.6\cdot10^{-4}$. Da applicare comunque, anche con A.
- **C. Pendenza unica del nodo** *(test di sensibilità)*. $S_0=S_1=S_2=S_{nodo}$, ad esempio la pendenza della finestra centrata sul nodo. La ripartizione dipende allora solo da larghezze e composizione: è un limite inferiore dell'effetto delle pendenze, utile per capire quanto la soluzione dipende da esse, ma rinuncia all'informazione sulla diversa pendenza dei rami.

Impatto sul documento: §2.2 (definizione della pendenza ai nodi).

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
- **C. Lasciare il criterio di Shields e segnalare i nodi con clamp** come non risolvibili. È la situazione attuale: i nodi vanno nel ripiego.

Impatto sul documento: §4.3 (definizione del tirante bankfull).

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
- **B. Continuazione sulle pendenze.** Si risolve prima il nodo con pendenza unica (§7.2 C), poi si portano gradualmente $S_1$, $S_2$ ai valori reali usando ogni soluzione come punto di partenza della successiva. Se a un certo passo la soluzione si perde (un ramo si chiude), quello è il limite fisico del modello per quel nodo.
- **C. Ramo che si chiude** ($q_2\le0$, come in 106). Il modello dice che, a quelle pendenze, l'equilibrio bankfull ha un solo ramo attivo, ma il ramo è bagnato il giorno del volo. Va deciso se trattarlo come esito fisico (ramo inattivo a bankfull, ripartizione da un altro criterio alle condizioni del volo) o come segnale di pendenze non affidabili.

### 7.6 Ripiego attuale

Quando il nodo non viene risolto, il codice:

1. ripartisce portata liquida e flussi solidi in proporzione alle larghezze bankfull dei rami, o, se mancano, a quelle delle maschere ($\psi=W_1/W_0$);
2. assegna ai rami la composizione della madre;
3. calcola il fondo all'imbocco dei rami con il moto uniforme, come nel canale singolo.

In questi nodi il risultato è quindi idraulico (moto uniforme con ripartizione geometrica), non morfodinamico. Va dichiarato nei risultati finché le questioni sopra non sono risolte: allo stato attuale le quote ai nodi del modulo morfodinamico non contengono l'effetto del bilancio di sedimento.

### 7.7 Decisioni richieste

| # | Decisione | Proposta |
| --- | --- | --- |
| 1 | Pendenza ai nodi (§7.2) | Finestra di lunghezza fissa $L$ per madre e rami (A) + sostituzione del valore minimo (B); pendenza unica (C) come test di sensibilità |
| 2 | Valore di $L$ | 300–500 m, oppure $L=10\,W_0^m$ |
| 3 | Tirante bankfull (§7.3) | $D_{bf}=\min(D_{Shields},\,z_{sfioro}-\eta)$ (A) |
| 4 | Larghezze non calcolabili (§7.4) | $W^{bf}=W^m$ (A); area sul semi-transetto ×2 con sponde specchiate (B) |
| 5 | Ramo che si chiude (§7.5 C) | Da decidere dopo aver applicato 1–3 |

L'ordine consigliato è 1 → 3 → 4, rieseguendo dopo ogni passo la tabella diagnostica di §7.1 per vedere quanti nodi passano da ripiego a soluzione.

---

## Registro modifiche

| Data | Modifica |
| --- | --- |
| 09/10/2026 | Conversione da `.docx` a Markdown; equazioni ricostruite in LaTeX. |
| 09/10/2026 | Spostato nel repo GitHub (`docs/`): da qui in poi la versione di riferimento è questa. |
| 09/10/2026 | Aggiunta la sezione 7 "Questioni aperte": nodi di biforcazione non risolti sui dati 2025, cause e scelte di metodo. |
