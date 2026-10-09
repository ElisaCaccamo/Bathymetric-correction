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

## Registro modifiche

| Data | Modifica |
| --- | --- |
| 09/10/2026 | Conversione da `.docx` a Markdown; equazioni ricostruite in LaTeX. |
| 09/10/2026 | Spostato nel repo GitHub (`docs/`): da qui in poi la versione di riferimento è questa. |
